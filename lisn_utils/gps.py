'''
Concatenation and decimation of GPS binary files
Suported formats:
 - Novatel
 - Leica
 - Ashtech

@author: Juan C. Espinoza
@contact: jucar.espinoza@gmail.com

'''
import os, re, gzip, math, traceback
import tarfile, zipfile, itertools
import numpy as np
from cStringIO import StringIO
from subprocess import Popen, PIPE
from math import sin, cos, asin, acos
from warnings import warn
from struct import pack, unpack
from urllib2 import urlopen
from datetime import datetime, timedelta
from ftplib import FTP
from utility import Data, EpochDict, open_file, PI, deg2rad, rad2deg, rglob
from gpsdatetime import GPSDateTime
from plotter import plot_data_vars, plot_rnx, plot_s4
from headers import *

try:
    from _crc import CRC32
except ImportError:
    warn('Error importing _crc.so extension NVD conversion will not be available')


#globals
jday_gps0   = 2444245
w_l1        = 0.190293672983469
w_l2        = 0.2442102134245683

kdays     = [31,28,31,30,31,30,31,31,30,31,30,31]

localpath = './'

rnx_codes = ('C1','C2','P1','P2','L1','L2','S1','S2')
tec_vars = ['eqTEC', 'ltime', 'ele', 'lat', 'lon', 'sTEC', 'azi']

mystrip   = lambda s  : s.strip()
truth     = lambda x  : True
to3float  = lambda s  : tuple([tofloat(s[k*14:(k+1)*14]) for k in (0,1,2)])

def limit(x):
    '''
    Limit the value of x to 10^10 (RINEX limitation)
    '''
    if x<0: return -(abs(x)%1000000000)
    return x%1000000000

def gdate(jd):
    '''
    Computes the Gregorian Calendar date (year, month, day)
    given the Julian date (jd)
    Translate to python from write_daily_nvd_file.c of Robert Sheehan
    <sheehanp@bc.edu>.
    '''
    l= jd+68569
    n= 4*l/146097
    l= l-(146097*n+3)/4
    i= 4000*(l+1)/1461001
    l= l-1461*i/4+31
    j= 80*l/2447
    k= l-2447*j/80
    l= j/11
    j= j+2-12*l
    i= 100*(n-49)+i+l

    return (int(i), int(j), int(k))

def jd(yy, mm, dd):
    '''
    Computes the Julian date (jd) given a Gregorian Calendar date
    (year, month, day).
    Translate to python from write_daily_nvd_file.c of Robert Sheehan
    <sheehanp@bc.edu>.
    '''
    return dd-32075+1461*(yy+4800+(mm-14)/12)/4+ \
                367*(mm-2-(mm-14)/12*12)/12-3*((yy+4900+(mm-14)/12)/100)/4

def fmod(x):
    fmodx = x
    while fmodx < 0 or fmodx > 2*PI:
        if fmodx >= 0:
            fmodx = fmodx - 2*PI
        else:
            fmodx = fmodx + 2*PI
    return fmodx

def min2date(mins):
    '''
    Converts time in minutes from january 1th into day of year, hour and minutes
    '''
    mm  = mins%60
    dif = (mins - mm)/60
    hh = dif%24
    doy = (dif - hh)/24+1
    return doy, hh, mm

def lla2xyz(lat, lon, alt):
    '''
    Convert lat, lon, height in WGS84 to ECEF (X,Y,Z)
    lat and long given in decimal degrees,altitude should be given in meters
    '''
    lat = lat*deg2rad
    lon = lon*deg2rad
    a   = 6378137.0
    f   = 1/298.257223563
    e2  = 2*f-f**2
    chi = math.sqrt(1-e2*(math.sin(lat))**2)
    X   = (a/chi +alt)*math.cos(lat)*math.cos(lon)
    Y   = (a/chi +alt)*math.cos(lat)*math.sin(lon)
    Z   = (a*(1-e2)/chi + alt)*math.sin(lat)
    return X,Y,Z

def xyz2lla(x, y, z):
    '''
    Calculate geodetic location from x,y and z
    lat & lon in degrees and alt in same units as .
    '''
    a   = 6378137.0
    e   = 8.1819190842622e-2
    b   = math.sqrt(a**2*(1-e**2))
    ep  = math.sqrt((a**2-b**2)/b**2)
    p   = math.sqrt(x**2+y**2)
    th  = math.atan2(a*z, b*p)
    lon = math.atan2(y, x)
    lat = math.atan2((z+ep**2*b*math.sin(th)**3), (p-e**2*a*math.cos(th)**3))
    N   = a/math.sqrt(1-e**2*math.sin(lat)**2)
    alt = p/math.cos(lat)-N
    
    return lat*rad2deg, lon*rad2deg, alt

    '''
    r = sqrt(x**2+y**2+z**2)
    if r < 6370000:
        r = 6370000
    xlat = atan(z/(sqrt(x**2+y**2)))
    xlat = xlat*rad2deg
    xlon = atan(y/x)
    if x >= 0.:
        xlon = xlon*rad2deg
    if x < 0.:
        if y >= 0.:
            xlon = (xlon*rad2deg) + 180.
        if y < 0.:
            xlon = (xlon*rad2deg) - 180.
    RE = 6378390*(0.99832 + 0.00168*cos(2.0*xlat/rad2deg))
    r -= RE
    return xlat, xlon, r
    '''

def TRNN(glat, glon, R):
    '''
    GEOCENTRIC TO GEODETIC CONVERSION
    GEOCENTRIC glat (RAD), glon (RAD), RADIUS (KM)
    SGD IS GEODETIC LAT (RAD), LONG (RAD), RADIUS (KM)
    SGR IS GEODETIC LAT (DEG), LONG (DEG), ALTITUDE (KM)
    translate to python from gpselem.f
    '''

    XC = [0,0,0]
    SGD = [0,0,0]
    SGR = [0,0,0]
    A = 6378.135
    B = 6356.75
    F     = (A - B) / A
    EP2   = 2.0 * F - F * F
    EP4   = EP2 * EP2
    C8    = 1.0 - EP2
    C5    = 2.0 * A**2 * EP2 * C8
    C6    = A**2 * C8 * C8
    C1    = -EP4 / C8
    C2    = -2.0 * EP2
    C3    = A**2 * EP4
    RO    = R * math.sin(glat)
    XC[0] = RO * math.cos(glon)
    XC[1] = RO * math.sin(glon)
    XC[2] = R * math.cos(glat)
    RO2   = RO * RO
    Z2    = XC[2] * XC[2]
    Z0    = A * XC[2] * math.sqrt(C8/(C8*RO2+Z2))
    E1    = 4.0 * C1
    D2    = C2 * XC[2]
    E2    = 3.0 * D2
    D3    = C3 - C8 * Z2 - RO2
    E3    = 2.0 * D3
    D4    = C5 * XC[2]
    D5    = C6 * Z2
    TEST = 1
    while (TEST-0.00001>0):
        FZ0   = (((C1*Z0+D2)*Z0+D3)*Z0+D4)*Z0+D5
        FPZ   = ((E1*Z0+E2)*Z0+E3)*Z0+D4
        Z1 = Z0 - FZ0/FPZ
        TEST=abs(Z1-Z0)
        Z0 = Z1
    SGD[1]=glon
    SGD[0]=math.atan(math.sqrt(C8*(A**2*C8-Z1*Z1))/Z1)
    if(SGD[0]<0):
        SGD[0]=SGD[0]+PI
    SPLT=math.cos(SGD[0])
    SGD[2]=XC[2]/SPLT-A*C8/math.sqrt(1.-EP2*SPLT*SPLT)
    GLT=PI/2-SGD[0]
    ALR=A*(.998320047+.001683494*math.cos(2.*GLT)-.000003549*math.cos(4.*GLT)+\
            0.000000008*math.cos(6.*GLT))
    SGR[0]=GLT*rad2deg
    SGR[1]=SGD[1]*rad2deg
    SGR[2]=SGD[2]
    SGD[2]=SGD[2]+ALR
    return SGD,SGR

def SILL(EL,AZ,XH,XTLON,XTLAT):
    '''
    THIS ROUTINE COMPUTES SUB-IONOSPHERIC LATITUDE AND
    LONGITUDE, FROM STATION COORDINATES, AND IONOSPHERIC
    HEIGHT.
    translate to python from gpselem.f
    '''

    ELR   = deg2rad*EL
    AZR   = deg2rad*AZ
    STLAT = deg2rad*XTLAT
    STLON = deg2rad*XTLON
    RE    = 6378.39*(0.99832 + 0.00168*math.cos(2.0*STLAT ))
    CHI   = math.asin(math.cos(ELR)*RE/(RE + XH))
    A     = 1.5707963 - CHI - ELR
    ARG = math.sin(STLAT)*math.cos(A)+math.cos(STLAT)*math.sin(A)*math.cos(AZR)    
    SILAT = math.asin(ARG)
    ARG = math.sin(A)*math.sin(AZR)/math.cos(SILAT)
    DLONG = math.asin(ARG)
    SILON = STLON+DLONG
    RANGE = (RE+XH)*math.sin(A)/math.cos(ELR)    
    SILAT = rad2deg*SILAT
    SILON = rad2deg*SILON
    return SILAT,SILON,RANGE

def GTODDBL(XLAT, XLON, SLAT, SLON, XALT=200.0):
    '''
    '''
    
    STLAT = SLAT*deg2rad
    RE = 6378.39*(.99832 + .00168*cos(2.*STLAT))
    C = abs(SLON - XLON)
    AS = deg2rad*(90.0 - SLAT)
    BS = deg2rad*(90.0 - XLAT)
    COSC  = cos(AS)*cos(BS) + sin(AS)*sin(BS)*cos(C*deg2rad)
    if COSC >= 0.999999999: COSC = 0.999999999
    if COSC <= -.999999999: COSC = -0.999999999
    HDIST = acos(COSC)*(RE + XALT)
    return abs(HDIST)

def cat_files(file_out, source, run_date=None, interval=10, extension='',
              old_extension='', label='', header='', gz=False, daily=True,
              fix_novatel=False, recursive=False):
    '''
    Concatenate files from "source" (also can deal with gziped files) and create
    a new decimated file depending of extension.
    
    Inputs:
        file_out  : file name for the new file
        source    : path where the files are located can be also a list of files
        run_date  : date to filter source files [datetime]
        interval  : interval in seconds
        extension : type of file: nvd, lb2, obs, txt
        old_extension : extension of the files in source
        label     : label used to filter files
        header    : TODO
        gz        : Gzip the new file
        daily     : Used for scintillation files        
    Output:
        Boolean   : True if the new file has valid records

    '''
    if isinstance(source, list):
        all_files = [f for f in source]
        file_list = [f for f in source]
    elif isinstance(source, basestring):
        if isinstance(run_date, (datetime, GPSDateTime)):
            run_date0 = run_date-timedelta(1)
            run_date1 = run_date+timedelta(1)
            filter  = '*%s*%s*.%s' % (label, run_date.strftime('%y%m%d'),
                                      old_extension)
            filter0 = '*%s*%s*.%s' % (label, run_date0.strftime('%y%m%d'),
                                      old_extension)
            filter1 = '*%s*%s*.%s' % (label, run_date1.strftime('%y%m%d'),
                                      old_extension)
            file_list0 = rglob(source, filter0, recursive=recursive)
            file_list1 = rglob(source, filter1, recursive=recursive)
            file_list0.sort()
            file_list1.sort()
        else:
            filter = '*.%s' % old_extension
        file_list = rglob(source, filter, recursive=recursive)
        file_list.sort()
        if isinstance(run_date, datetime):
            all_files = file_list0[-20:]+file_list+file_list1[:5]
        else:
            all_files = [f for f in file_list]
        if len(file_list) == 0:
            print 'No files found at:', source
            return False
    else:
        print 'Invalid source'
        return False

    print 'Creating %s from *%s*.%s files' % (file_out, label, old_extension)
    dec_file     = DecimateFile(file_out, extension, 
                                '%s.%s' % (label, old_extension), run_date)
    last_gps_sec = 0
    last_sec     = -1
    line         = ''
    for file in all_files:
        #print ' Processing: %s' % file
        try:
            fileo = open_file(file)
        except IOError:
            print ' Error reading file:', file, 'skipping'
            continue
        if extension=='nvd':
            last_gps_sec = dec_file.add_nvd(fileo, last_gps_sec, interval, 
                                            fix_novatel=fix_novatel)
        elif extension=='lb2':
            last_gps_sec = dec_file.add_lb2(fileo, last_gps_sec, interval)
        elif extension=='obs':
            last_gps_sec = dec_file.add_obs(fileo, last_gps_sec, interval)
        else:
            if 'txt' in old_extension:
                line = open_file(file_list[0]).readline()
            last_sec = dec_file.add_file(fileo, last_sec, interval, line, daily)
        fileo.close()
    dec_file.close(gz)
    if dec_file.nrec<>0:
        return True
    return False

def nvd_to_rnx(input_file, marker='abcd', rnx_path=localpath, date=None,
                obscodes=rnx_codes, interval=1, alt_name=None, 
                fix_novatel=False, **kwargs):
    '''
    Read a Novatel binary file (nvd) and create a Rinex observation file.

    Inputs:
        input_file: binary novatel file
        marker: Four char station name for rinex header.
        rnx_path: path where the rinex file will be created.
        obscodes: list of rinex codes
        interval: interval for decimated data.
        kwargs: aditional information for rinex header.
    Outputs: Rinex filename if created and number of records
    '''

    def decode_id43(buffer):
        '''
        Unpack id43 records from buffer, return  a list of records objects
        '''
        num_obs = unpack('=L', buffer[:4])[0]
        ind = 4
        rec_list = []
        last_rec = ReadRecord()
        for i in range(num_obs):
            record1 = ReadRecord('BB', buffer[ind:ind+2])
            record2 = ReadRecord('=dfdffffL', buffer[ind+4:ind+44])
            record1.vals += record2.vals
            record1.id43()
            #print record1.prn, record1
            if record1.prn <> last_rec.prn:
                last_rec = record1
            else:
                last_rec.update(record1)
                rec_list.append(last_rec)
            ind += 44
        return rec_list

    if input_file.lower().endswith('.nvd.gz') or \
        input_file.lower().endswith('.nvd'):
        try:
            rawfile = open_file(input_file)
            print 'Opening file', input_file
        except IOError:
            print 'Error reading file:', input_file
            return None, 0
    else:
        print input_file, 'Invalid nvd file'
        return None, 0

    nrec = 0
    intervals = set()
    sync = [0, 0, 0]
    last_obs = datetime(1,1,1)

    while True:
        try:
            ch = unpack('B', rawfile.read(1))[0]
        except:
            break
        sync = [sync[1], sync[2], ch]
        #print old_pos, sync
        if sync == [0xAA, 0x44, 0x12]:
            file_pos = rawfile.tell()
            buff0 = rawfile.read(1)
            hdr_len = unpack('B', buff0)[0]
            buff1 = rawfile.read(hdr_len-4)
            header = ReadRecord('=HBBHHBBHlLHH', buff1)
            header.id43_header()
            if fix_novatel and header.week < 1024:
                header.week += 1024
            week_sec = header.milli_sec/1000
            gps_sec = header.week*SECONDS_IN_WEEK + week_sec
            dow = week_sec/86400
            daysec = week_sec%86400
            hh = daysec/3600
            mm = (daysec%3600)/60
            ss = daysec%60
            rec_date = gdate(header.week*7 + dow + jday_gps0) + (hh, mm, ss)
            buff2 = rawfile.read(header.msg_len)
            buffer = '\xaa\x44\x12'+buff0+buff1+buff2
            crc_1 = unpack('=L', rawfile.read(4))[0]
            crc_2 = CRC32(buffer) & 0xFFFFFFFF

            if crc_1 <> crc_2:
                rawfile.seek(file_pos, 0)
                print 'CRC_error', crc_1, crc_2, rec_date
                continue

            is_dec = gps_sec%interval == 0

            if header.msg_id == 43 and is_dec:
                obs_list = decode_id43(buff2)
                if nrec == 0 and len(obs_list) <> 0:
                    first_obs = GPSDateTime(*rec_date)
                    date_obs  = GPSDateTime(*rec_date)
                    if alt_name:
                        rnx_name = '%s.%02do' % (alt_name, date_obs.year%100)
                    else:
                        rnx_name = '%s%03d0.%02do' % (marker, date_obs.doy,
                                                      date_obs.year%100)
                    print 'Creating rinex file:', os.path.join(rnx_path, rnx_name)
                    obscodes = [x for x in obscodes if x in obs_list[0]]
                    Rinex = MkRinex(marker, interval, obscodes)
                    file_tmp = os.tmpfile()
                elif nrec <> 0 and len(obs_list) <> 0:
                    intervals.add(GPSDateTime(*rec_date)-date_obs)
                    date_obs = GPSDateTime(*rec_date)
                    if date_obs<=last_obs:
                        continue
                    last_obs = date_obs
                else:
                    continue
                prn_list = ['G%02d' % rec.prn for rec in obs_list]
                sat_line = Rinex.mk_satline(date_obs, prn_list)
                dat_line = ''
                for rec in obs_list:
                    obs_values = [rec[x] for x in obscodes]
                    lli_values = [0 for x in obs_values]
                    snr_values = [Rinex.snr_flag(rec, x) for x in obscodes]
                    dat_line += Rinex.mk_recordline(obs_values, lli_values,
                                snr_values)
                file_tmp.write(sat_line+dat_line)
                nrec += 1
    if nrec<>0:
        kwargs['endtime'] = last_obs
        if intervals:
            if min(intervals) <> interval:
                Rinex.interval = min(intervals)
        header = Rinex.mk_header(first_obs, **kwargs)
        if not os.path.exists(rnx_path):
            os.makedirs(rnx_path)
        file_tmp.seek(0,0)
        file_rnx = open(os.path.join(rnx_path, rnx_name), 'w')
        file_rnx.write(header)
        file_rnx.write(file_tmp.read())
        rawfile.close()
        file_tmp.close()
        file_rnx.close()
        print 'Rinex records:', nrec
        return os.path.join(rnx_path, rnx_name), nrec
    else:
        print 'No records found for given date'
        return None, 0

def lb2_to_rnx(input_file, marker='site', rnx_path=localpath, date=None,
                obscodes=rnx_codes, interval=1, alt_name=None, wpos=False,
                **kwargs):
    '''
    Read a Leica file lb2 and create a Rinex observation file.

    Inputs:
        input_file: binary leica file (id37 or lb2)
        marker: Four char station name for rinex header.
        rnx_path: path where the rinex file will be created.
        obscodes: list of rinex codes
        interval: interval for decimated data
        kwargs: aditional information for rinex header.
    Outputs: Rinex filename if created and number of records
    '''

    def decode_id3(buffer, ndata):
        '''
        Unpack id3 records from buffer, return  a list of records objects
        '''
        ind = 6
        obs_list = []
        while ind<ndata:
            if len(buffer[ind:ind+3])<>3:
                ind += len(buffer[ind:ind+3])
                continue
            record = ReadRecord('=Bh', buffer[ind:ind+3])
            record.id3(buffer, ind+3)
            if sum(record.values())<>0:
                obs_list.append(record)
            ind = record.ind
            #print record.ind
        if ind==ndata:
            return obs_list
        else:
            return None

    def decode_id4(buffer, ndata):
        '''
        Unpack id4 records from buffer, return  a list of records objects
        '''
        nsat = unpack('B', buffer[0])[0]
        ind = 1
        pos_list = []
        while ind<ndata:
            record = ReadRecord('=Bhh', buffer[ind:ind+5])
            record.id4()
            pos_list.append(record)
            ind += 5
        if ind==ndata:
            return pos_list
        else:
            return None

    def decode_id85(buffer, ndata):
        '''
        Unpack id85 records from buffer, return  a list of records objects
        '''
        record = ReadRecord('B', buffer[0])
        record.id85(buffer)
        #if record.ind == ndata:
        return record

    def decode_id37(buffer, header, ind):
        '''
        Unpack id37 records from buffer, return  a list of records objects
        '''
        obs_list = []
        for i in range(header.nprn):
            record = ReadRecord('=Bdhdh', buffer[ind:ind+21])
            record.id37()
            if sum(record.values())<>0:
                obs_list.append(record)
            ind += 24
        return obs_list

    if input_file.lower().endswith('.lb2') or \
        input_file.lower().endswith('.lb2.gz'):
        try:
            rawfile = open_file(input_file)
            print 'Opening file', input_file
        except IOError:
            print 'Error reading file:', input_file
            return None, 0
    else:
        print input_file, 'Invalid lb2 file'
        return None, 0

    file_name = input_file.split('/')[-1]
    if date:
        file_date = GPSDateTime(date)
    else:
        try:
            file_date = GPSDateTime(file_name.split('.lb2')[0], filename=True)
        except:
            print 'Filename does not have a valid date stamp: yymmdd'
            return None, 0

    week = file_date.GPSweek
    nrec = 0
    prec = 0
    sync = [0, 0]
    intervals = set()
    first_jd = 0
    first_obs = 0
    last_obs = datetime(1,1,1)
    
    while True:
        try:
            ch = unpack('B', rawfile.read(1))[0]
        except:
            break
        sync = [sync[1], ch]
        if sync == [0x9C, 0xAE]:
            file_pos = rawfile.tell()
            buff0 = rawfile.read(2)
            rec_len = unpack('H', buff0)[0]
            if rec_len >8192:
                rawfile.seek(file_pos, 0)
                continue
            buff1 = rawfile.read(1)
            rec_id = unpack('B', buff1)[0]
            buff2 = rawfile.read(rec_len-5)
            buff3 = rawfile.read(2)
            checksum1 = unpack('h', buff3)[0]
            num_bytes = 3+rec_len-5
            checksum2 = sum(unpack('%iB' % num_bytes, buff0+buff1+buff2))*-1

            '''
            if checksum1 <> checksum2:
                print 'Checsum error'
                rawfile.seek(file_pos, 0)
                continue
            '''

            save_rnx = False
            save_pos = False
            ndata = rec_len-5

            if rec_id == 0x02:
                #obs_list = decode_id2(buff2, ndata)
                continue
            elif rec_id == 0x03:
                if ndata <=6:
                    rawfile.seek(file_pos, 0)
                    continue
                header = ReadRecord('=3chB', buff2[:6])
                header.id3_header()
                is_dec = header.gpstow%(interval*10) == 0
                week_sec = header.gpstow/10
                dow = week_sec/86400
                day_sec = week_sec % 86400
                hms = (day_sec/3600, (day_sec % 3600)/60, day_sec % 60)
                rec_date = gdate(week*7 + dow + jday_gps0) + hms

                obs_list = decode_id3(buff2, ndata)
                save_rnx = True
                save_pos = False
            elif rec_id == 0x04:
                pos_list = decode_id4(buff2,ndata)
                save_pos = True
                save_rnx = False
            elif rec_id == 0x85:
                record = decode_id85(buff2, ndata)
                if record.date[2]>80:
                    yr = record.date[2]+1900
                else:
                    yr = record.date[2]+2000
                date_id85 = (yr, record.date[1], record.date[0])
                jd_id85 = jd(*date_id85)
                continue
            elif rec_id == 0x37:
                station_prefix = buff2[:2]
                i0 = 2
                if buff2[3:6] == 'ver':
                    print 'ver True'
                    i0 += 10
                header = ReadRecord('11B', buff2[i0:i0+11])
                header.id37_header()
                is_dec = header.gps_sec%interval == 0
                rec_date = header.gps_date+header.dhms[1:]

                obs_list = decode_id37(buff2, header, i0+13)
                save_rnx = True
                save_pos = False
            else:
                #print 'Lb2 id %d unknown, skipping' % rec_id
                continue

            if save_rnx and is_dec:
                if nrec == 0 and len(obs_list)<>0:
                    first_obs = GPSDateTime(*rec_date)
                    date_obs  = GPSDateTime(*rec_date)
                    if date_obs<file_date:
                        continue
                    elif date_obs>(file_date+timedelta(2)):
                        continue
                    if alt_name:
                        rnx_name = '%s.%02do' % (alt_name, date_obs.year%100)
                    else:
                        rnx_name = '%s%03d0.%02do' % (marker, date_obs.doy,
                                                      date_obs.year%100)
                    print 'Creating rinex file:', os.path.join(rnx_path,
                                                                rnx_name)
                    obscodes = [x for x in obscodes if x in obs_list[0]]
                    Rinex = MkRinex(marker, interval, obscodes)
                    file_tmp = os.tmpfile()
                elif nrec <> 0 and len(obs_list)<>0:
                    intervals.add(GPSDateTime(*rec_date)-date_obs)
                    date_obs = GPSDateTime(*rec_date)
                    if date_obs<=last_obs:
                        continue
                    last_obs = date_obs
                else:
                    continue

                prn_list = ['G%02d' % rec.prn for rec in obs_list]
                sat_line = Rinex.mk_satline(date_obs, prn_list)
                dat_line = ''
                for rec in obs_list:
                    obs_values = [rec[x] for x in obscodes]
                    lli_values = [0 for x in obs_values]
                    snr_values = [Rinex.snr_flag(rec, x) for x in obscodes]
                    dat_line += Rinex.mk_recordline(obs_values, lli_values,
                                snr_values)
                file_tmp.write(sat_line+dat_line)
                nrec += 1

            if save_pos and wpos:
                if prec==0:
                    pos_name = '%s_%02d%02d%02d.pos' % (marker, rec_date[0]%100,
                        rec_date[1], rec_date[2])
                    pos_path = os.path.join(os.path.dirname(rnx_path), 'posit')
                    print 'Creating position file:', os.path.join(pos_path,
                                                                    pos_name)
                    if not os.path.exists(pos_path):
                        os.makedirs(pos_path)
                    pos_file = open(os.path.join(pos_path, pos_name), 'w')
                    pos_file.write('HH MM SS PRN  az el')
                pos_line = ''
                for rec in pos_list:
                    pos_line += '%02d %02d %02d %3d %3d %2d\n' % (rec_date[3],
                                rec_date[4], rec_date[5],rec.prn, rec['az'],
                                rec['el'])
                pos_file.write(pos_line)
                prec +=1
    if prec<>0:
        print 'pos records:', prec
        pos_file.close()

    if nrec<>0:
        kwargs['endtime'] = last_obs
        if intervals:
            if min(intervals) <> interval:
                Rinex.interval = min(intervals)
        header = Rinex.mk_header(first_obs, **kwargs)
        if not os.path.exists(rnx_path):
            os.makedirs(rnx_path)
        file_tmp.seek(0,0)
        file_rnx = open(os.path.join(rnx_path, rnx_name), 'w')
        file_rnx.write(header)
        file_rnx.write(file_tmp.read())
        rawfile.close()
        file_tmp.close()
        file_rnx.close()
        print 'Rinex records:', nrec
        return os.path.join(rnx_path, rnx_name), nrec
    else:
        print 'No records found for given date'
        return None, 0

def obs_to_rnx(input_file, marker='site', rnx_path=localpath, date=None,
                obscodes=rnx_codes, interval=1, alt_name=None, pcode='P1', **kwargs):
    '''
    Read a observable ascii file (obs) and create a Rinex observation file.

    Inputs:
        input_file: observables file (GPS-Scinda)
        marker: Four char station name for rinex header.
        rnx_path: path where the rinex file will be created.
        obscodes: list of rinex codes
        interval: interval for decimated data.
        kwargs: aditional information for rinex header.
    Outputs: Rinex filename if created and number of records
    '''

    def decode_obs(rawfile, last, pcode):
        '''
        Unpack obs records from rawfile, return  a list of records objects
        '''
        first = False
        obs_list = last[-1:]
        if len(obs_list)<>0:
            millisec = obs_list[-1].millisec
        while True:
            buff = rawfile.readline()
            if not buff:
                last = []
                break
            dataline = ReadRecord(None, buff)
            if len(dataline.vals)<>11: continue
            dataline.obs_line(pcode)
            if len(obs_list)<>0:
                if millisec < dataline.millisec:
                    last = [dataline]
                    break
            millisec = dataline.millisec
            if dataline.prn not in [dl.prn for dl in obs_list]:
                obs_list.append(dataline)
        return obs_list, last

    if input_file.lower().endswith('.obs') or \
        input_file.lower().endswith('.obs.gz'):
        try:
            rawfile = open_file(input_file)
            print 'Opening file', input_file
        except IOError:
            print 'Error reading file:', input_file
            return None, 0
    else:
        print input_file, 'Invalid obs file'
        return None, 0

    file_name = input_file.split('/')[-1]
    
    if date:
        file_date = date
    else:
        try:
            file_date = GPSDateTime(file_name, filename=True)
        except:
            print 'Filename does not have a valid date stamp: yymmdd'
            return None, 0

    week = file_date.GPSweek
    nrec = 0
    millisec = 0
    intervals = set()
    last = []
    #last_obs = datetime(1,1,1)

    while True:
        obs_list, last = decode_obs(rawfile, last, pcode)
        #print obs_list, last
        if len(last)==0:
            break
        week_sec = obs_list[0].millisec/1000
        gps_sec = week*SECONDS_IN_WEEK + week_sec
        dow = week_sec/86400
        daysec = week_sec%86400
        hh = daysec/3600
        mm = (daysec%3600)/60
        ss = daysec%60
        rec_date = gdate(week*7 + dow + jday_gps0) + (hh, mm, ss)
        if not gps_sec%interval == 0:
            continue
        if nrec == 0 and len(obs_list) <> 0:
            first_obs = GPSDateTime(*rec_date)
            date_obs  = GPSDateTime(*rec_date)
            last_obs  = GPSDateTime(*rec_date)
            if date_obs<file_date:
                continue
            elif date_obs>(file_date+timedelta(1)):
                continue
            if alt_name:
                rnx_name = '%s.%02do' % (alt_name, date_obs.year%100)
            else:
                rnx_name = '%s%03d0.%02do' % (marker, date_obs.doy,
                                              date_obs.year%100)
            print 'Creating rinex file:', os.path.join(rnx_path, rnx_name)
            obscodes = [x for x in obscodes if x in obs_list[0]]
            Rinex = MkRinex(marker, interval, obscodes)
            file_tmp = os.tmpfile()
        elif nrec<>0 and len(obs_list)<>0:            
            date_obs = GPSDateTime(*rec_date)
            if date_obs<=last_obs:
                continue
            intervals.add(date_obs-last_obs)
            last_obs = date_obs
        else:
            continue

        prn_list = ['G%02d' % rec.prn for rec in obs_list]
        sat_line = Rinex.mk_satline(date_obs, prn_list)
        dat_line = ''
        for rec in obs_list:
            obs_values = [rec[x] for x in obscodes]
            lli_values = [0 for x in obs_values]
            snr_values = [Rinex.snr_flag(rec, x) for x in obscodes]
            dat_line += Rinex.mk_recordline(obs_values, lli_values, snr_values)
        file_tmp.write(sat_line+dat_line)
        nrec += 1
    if nrec<>0:
        kwargs['endtime'] = last_obs
        if intervals:
            if min(intervals) <> interval:
                Rinex.interval = min(intervals)
        header = Rinex.mk_header(first_obs, **kwargs)
        if not os.path.exists(rnx_path):
            os.makedirs(rnx_path)
        file_tmp.seek(0,0)
        file_rnx = open(os.path.join(rnx_path, rnx_name), 'w')
        file_rnx.write(header)
        file_rnx.write(file_tmp.read())
        rawfile.close()
        file_tmp.close()
        file_rnx.close()
        print 'Rinex records:', nrec
        return os.path.join(rnx_path, rnx_name), nrec
    else:
        print 'No records found for given date'
        return None, 0

def btog(c):
    if c in (None, '', ' '):
        return 'G'
    return c.upper()

def toint(x):
    if x is None or x.strip() == '':
        return 0
    return int(x)

def tofloat(x):
    if x is None or x.strip() == '':
        return 0.
    return float(x)

def choose(a, b):
    if a is not None and b in (' ', None):
        return a
    return b.replace('&', ' ')

def value(thing):
    '''
    Ensure that arbitrary attributes can be set on `thing'.
    E.g. foo = value(foo); foo.bar = 'qux'
    '''

    if type(thing)in (float, int):
        return my_float(thing)
    elif type(thing)==str:
        thing = my_str(thing)
    return thing

def fitlin(X, Y):
    '''
    '''
    n = float(len(X))
    xbar = sum(X)/n
    ybar = sum(Y)/n
    sumxx = sum([i**2 for i in X])
    sumxy = sum([X[i]*Y[i] for i in range(len(X))])
    b = (sumxy - n*xbar*ybar)/(sumxx-n*xbar**2)
    a = ybar - b*xbar
    return a, b

def fix_jumps(Y, X, vmax, j=None, id=None):
    '''
    Fix jumps of a list Y
    '''

    if j:
        start = j
        end = j+1
    else:
        start = 1
        end = len(Y)

    for i in range(start, end):
        diff = abs(Y[i] - Y[i-1])
        if diff>vmax:
            if i>=start+4:
                i1 = i-10
                if i1<1: i1 = 1
                i2 = i-1
                x1 = [X[x]-X[i1] for x in range(i1, i2)]
                y1 = [Y[x] for x in range(i1,i2)]
                a, b = fitlin(x1, y1)
                ynew  = a + b*(X[i]-X[i1])
                ydiff = ynew - Y[i]
            else:
                psum = sum([Y[x] for x in range(i)])
                ynew  = psum/i
                ydiff = ynew - Y[i]
            #if id:
            #    print 'Jump found: id=%s, x=%.2f, y1=%.3f, y0=%.3f, diff =%.3f, ydiff =%.3f'% \
            #    (id, X[i], Y[i], Y[i-1], diff, ydiff)
            for k in range(i,len(Y)):
                Y[k] += ydiff
    return Y

def sign(a, b):
    if b>=0:
        return abs(a)
    else:
        return -abs(a)

def mnbrak(ax, bx, F, args):
    '''
    Find a range where a zero of the function F(x, *args) is located
    ax, bx are just guest (initial values)
    '''

    golden = 1.618
    glimit = 100.
    tiny = 1e-18
    dum = 0

    fa = F(ax, *args)
    fb = F(bx, *args)

    if fb>fa:
        dum, ax, bx = ax, bx, dum
        dum, fb, fa = fb, fa, dum

    cx = bx+golden*(bx-ax)
    fc = F(cx, *args)

    while fb>fc:
        r = (bx-ax)*(fb-fc)
        q = (bx-cx)*(fb-fa)
        u = bx-((bx-cx)*q-(bx-ax)*r)/(2.0*sign(max(abs(q-r), tiny), q-r))
        ulim = bx+glimit*(cx-bx)
        if (bx-u)*(u-cx) > 0.0:
            fu = F(u, *args)
            if fu < fc:
                ax, bx = bx, u
                fa, fb = fb, fu
                break
            elif fu > fb:
                cx=u
                fc=fu
                break
            u=cx+golden*(cx-bx)
            fu = F(u, *args)
        elif (cx-u)*(u-ulim) > 0.0:
            fu = F(u, *args)
            if fu < fc:
                bx, cx = cx,u
                u=cx+golden*(cx-bx)
                fb, fc = fc, fu
                fu=F(u, *args)
        elif (u-ulim)*(ulim-cx) >= 0.0:
            u=ulim
            fu = F(u, *args)
        else:
            u=cx+golden*(cx-bx)
            fu = F(u, *args)
        ax, bx, cx = bx, cx, u
        fa, fb, fc = fb, fc, fu
    return (ax, bx, cx), (fa, fb, fc)

def brent(val, fval, F, args, itmax=20, tol=1e-2):
    '''
    Find a zero of a function F(x, *args) using brent method.
    val = (ax, bx, cx)
    fval = (F(ax), F(bx), F(cx))

    val: is a range where is located the Zero.
    '''

    cgold = 0.381966
    zeps = 1e-10
    e = 0

    ax, bx, cx = val
    fa, fb, fc = fval

    a=min(ax, cx)
    b=max(ax, cx)

    x, w, v = bx, bx, bx
    fx, fw, fv = fb, fb, fb

    for i in range(itmax):
        xm = .5*(a+b)
        tol1 = tol*abs(x)+zeps
        tol2 = 2*tol1
        if abs(x-xm)<=(tol2-0.5*(b-a)):
            return x, fx, i
        if abs(e) > tol1:
            r=(x-w)*(fx-fv)
            q=(x-v)*(fx-fw)
            p=(x-v)*q-(x-w)*r
            q=2.0*(q-r)
            if q > 0.0: p = -p
            q=abs(q)
            etemp=e
            e=d
            if (abs(p) >= abs(0.5*q*etemp)) or (p<=q*(a-x)) or (p >=q*(b-x)):
                if x>=xm: e = a-x
                else: e = b-x
                d=cgold*e
            else:
                d=p/q
                u=x+d
                if (u-a < tol2) or (b-u < tol2):
                    d=sign(tol1, xm-x)
        else:
            if x>=xm: e = a-x
            else: e = b-x
            d=cgold*e
        if abs(d)>=tol1: u = x+d
        else: u = x+sign(tol1, d)
        fu = F(u, *args)

        if (fu <= fx):
            if (u >= x): a=x
            else: b=x
            v, w, x = w, x, u
            fv, fw, fx = fw, fx, fu
        else:
            if (u < x): a=u
            else: b=u
            if (fu <= fw) or (w == x):
                v, w = w, u
                fv, fw = fw, fu
            elif (fu <= fv) or (v == x) or (v == w):
                v=u
                fv=fu
    print 'Max iteration reach'
    return x, fx, i

class ReadRecord(dict):
    '''
    Class to unpack binary string or buffer into data values.
    '''
    def __init__(self, format=None, buffer=None):
        if format is not None and buffer is not None:
            self.vals = unpack(format, buffer)
        elif format is None and buffer is not None:
            self.vals = [x.strip() for x in buffer.split(' ') if len(x)<>0]
        self.prn = 0

    def id43_header(self):
        '''
        unpack id43 header (nvd file)
        '''
        self.msg_id        = self.vals[0]    #unsigned short
        self.msg_type      = self.vals[1]    #char
        self.port_add      = self.vals[2]    #char
        self.msg_len       = self.vals[3]    #unsigned short
        self.sequence      = self.vals[4]    #unsigned short
        self.idle_time     = self.vals[5]    #char
        self.time_status   = self.vals[6]    #char
        self.week          = self.vals[7]    #unsigned short
        self.milli_sec     = self.vals[8]    #long
        self.receiver_stat = self.vals[9]    #unsigned long
        self.reserved      = self.vals[10]   #unsigned short
        self.sw_version    = self.vals[11]   #unsigned short

    def id43(self):
        '''
        unpack id43 record of nvd file
        '''

        class track_status(object):
            '''
            Class to hold track status of a nvd record
            '''
            def __init__(self, track_stat):
                self.phase_lock = (track_stat >> 10) & 0x01
                self.parity     = (track_stat >> 11) & 0x01
                self.code_lock  = (track_stat >> 12) & 0x01
                self.group      = (track_stat >> 20) & 0x01
                self.freq       = (track_stat >> 21) & 0x03
                self.code       = (track_stat >> 23) & 0x07

        track = track_status(self.vals[9])

        self.prn        = self.vals[0]    #unsigned short
        self.reserved   = self.vals[1]    #unsigned short
        if track.freq == 0:
            self['C1']  = self.vals[2]    #double
            self['L1']  = limit(self.vals[4]*-1) #double
            self['S1']  = self.vals[7]    #float
        if track.freq == 1:
            self['P2']  = self.vals[2]    #double
            self['L2']  = limit(self.vals[4]*-1) #double
            self['S2']  = self.vals[7]    #float
        self.psr_std    = self.vals[3]    #float
        self.adr_std    = self.vals[5]    #float
        self.dopp       = self.vals[6]    #float
        self.locktime   = self.vals[8]    #float

    def id3_header(self):
        '''
        '''
        temp = self.vals[0:3]+('\x00',)
        self.gpstow        = unpack('=I', ''.join(temp))[0]
        self.accum         = self.vals[3]
        self.iflag         = self.vals[4]

    def id3(self, buffer, ind):
        '''
        '''
        self.dopp = self.vals[0]
        L1denb = (self.dopp >> 3) & 0x01
        L2denb = (self.dopp >> 4) & 0x01
        stat = self.vals[1]
        chan = (stat & 0x07);
        L1trk = ((stat >> 4) & 0x03)
        L2trk = ((stat >> 6) & 0x03)
        L1lck = ((stat >> 9) & 0x01)
        L2lck = ((stat >> 10) & 0x01)
        self.prn = ((stat >> 11) & 0x1f)+1
        self.L1qual = 0
        self['L1'] = 0.
        self.L1code = 0
        self.L1contrk = 0
        self.L1dopp = 0
        self.L1costi = 0
        self.L1cn0 = 0
        if L1trk == True:
            data = unpack('=BdhB', buffer[ind:ind+12])
            ind += 12
            self.L1qual = data[0]
            self['L1'] = limit(data[1])
            self.L1code = limit(data[2])
            self.L1contrk = data[3]
            if L1denb == True:
                data = unpack('=L', buffer[ind:ind+4])
                ind += 4
                self.L1dopp = data[0]
            self.L1costi = self.L1qual & 0x7
            self.L1cn0 = self.L1qual >> 3
        self.L2qual = 0
        self['L2'] = 0.
        self.L2code = 0
        self.L2contrk = 0
        self.L2dopp = 0
        if L2trk == True:
            data = unpack('=BdhB', buffer[ind:ind+12])
            ind += 12
            self.L2qual = data[0]
            self['L2'] = limit(data[1])
            self.L2code = limit(data[2])
            self.L2contrk = data[3]
            if L2denb == True:
                data = unpack('=L', buffer[ind:ind+4])
                ind += 4
                self.L2dopp = data[0]
            self.L2costi = self.L2qual & 0x7
            self.L2cn0 = self.L2qual >> 3
        self['C1'] = (float(self['L1']) + self.L1code/32.) * w_l1
        self['P2'] = (float(self['L2']) + self.L2code/32.) * w_l2
        self.ind = ind

    def id4(self):
        self.prn = self.vals[0]
        self['az'] = self.vals[1]
        self['el'] = self.vals[2]

    def id37_header(self):
        self.pc_date       = self.vals[0:3]    #char
        self.pc_time       = self.vals[3:6]    #char
        self.dhms          = self.vals[6:10]   #char
        self.nprn          = self.vals[10]     #char
        if self.pc_date[0] < 90:
            juld_pc = jd(self.pc_date[0]+2000, self.pc_date[1], self.pc_date[2])
        else:
            juld_pc = jd(self.pc_date[0]+1900, self.pc_date[1], self.pc_date[2])
        idow_pc = (juld_pc % 7) + 1
        idow_gps = self.dhms[0]
        del_dow = idow_pc - idow_gps
        if del_dow == 0:
            juld_gps = juld_pc
        else:
            if del_dow > 1:
                del_dow = del_dow - 7
            if del_dow < -1:
                del_dow = del_dow + 7
            juld_gps = juld_pc - del_dow
        idsec = 3600*self.dhms[1] + 60*self.dhms[2] + self.dhms[3]
        idsec_pc = 3600*self.pc_time[0] + 60*self.pc_time[1] + self.pc_time[2]
        gps_day = juld_gps - jday_gps0
        self.gps_sec = gps_day * 86400 + idsec
        self.gps_date = gdate(juld_gps)

    def id37(self):
        '''
        unpack id37 record of lb2 file
        '''
        self.prn     = self.vals[0] + 1   #unsigned short
        self['L1']   = limit(self.vals[1])       #double
        self.code1   = self.vals[2]       #short
        self['L2']   = limit(self.vals[3])       #double
        self.code2   = self.vals[4]       #short
        #self.el      = self.vals[4]      #unsigned short
        #self.az      = self.vals[5]      #unsigned short
        self['C1']   = (self['L1'] + self.code1/32.0) * w_l1
        self['P2']   = (self['L2'] + self.code2/32.0) * w_l2

    def id85(self, buffer):
        '''
        unpack id85 record of lb2 file
        '''
        ind = 1
        iflag = self.vals[0]
        if iflag and 0x1:
            self.date = unpack('3B', buffer[1:4])
            self.time = unpack('3B', buffer[4:7])
            self.offset = unpack('h', buffer[7:9])
            ind+=8
        if (iflag>>1) and 0x1:
            ind+=1
            llh = unpack('=3d', buffer[ind:ind+24])
            self.lat = llh[0]
            self.lon = llh[1]
            self.hgt = llh[2]
            ind += 24
        self.ind = ind

    def obs_line(self, pcode):
        '''
        unpack obs line record of obs file
        '''
        self.millisec = int(self.vals[0])
        #self.warning_CA = int(self.vals[1])
        #self.warning_P1 = int(self.vals[2])
        #self.warning_P2 = int(self.vals[3])
        self['S1'] = float(self.vals[4])-30
        self['S2'] = float(self.vals[5])-30
        self[pcode] = float(self.vals[6])
        self['P2'] = float(self.vals[7]) + float(self.vals[6])
        self['L1'] = limit(float(self.vals[8]))
        self['L2'] = limit(float(self.vals[9]))+ limit(float(self.vals[8]))*60/77
        self.prn = int(self.vals[10])

class DecimateFile(object):
    '''
    Class to create a new file with decimated values.
    '''

    def __init__(self, file_name=None, extension='dat', old_extension='', date=None):
        if file_name is None:
            file_name = os.path.join(localpath, 'abcd.'+extension)
        self.file_name  = file_name
        self.old_extension   = old_extension
        self.fileo      = os.tmpfile()
        self.data       = ''
        self.nrec       = 0
        self.first_line = None
        self.fix_line   = ''
        self.date       = date
        self.bad_doy    = None
        self.doy        = 0
        self.last_prns  = []

    def add_nvd(self, filei, gps_sec_last=0, interval=10, fix_novatel=False):
        '''
        Add new raw nvd data files to the decimated new file.
        Based in write_daily_nvd_file.c of Robert Sheehan <sheehanp@bc.edu>

        Inputs:
            filei: input file object
            gps_sec_last: seconds (gps time) to start append new data
            interval: decimate value in seconds

        Outputs:
            gps_sec_last: seconds (gps time) for the last data block appened.
        '''
        nrec = 0
        sync = [0, 0, 0]
        while True:
            try:
                ch = unpack('B', filei.read(1))[0]
            except:
                break
            sync = [sync[1], sync[2], ch]
            #print old_pos, sync
            if sync == [0xAA, 0x44, 0x12]:
                file_pos = filei.tell()
                buff0 = filei.read(1)
                hdr_len = unpack('B', buff0)[0]
                buff1 = filei.read(hdr_len-4)
                try:
                    header = ReadRecord('=HBBHHBBHlLHH', buff1)
                except:
                    filei.seek(file_pos)
                    continue
                header.id43_header()
                if fix_novatel and header.week<1024:
                    header.week += 1024
                week_sec = header.milli_sec/1000
                gps_sec = header.week*SECONDS_IN_WEEK + week_sec
                dow = week_sec/86400
                daysec = week_sec%86400
                hh = daysec/3600
                mm = (daysec%3600)/60
                ss = daysec%60
                rec_date = gdate(header.week*7 + dow + jday_gps0) + (hh, mm, ss)
                fdate = GPSDateTime(*rec_date[:3])

                if self.date:
                    if fdate < self.date:
                        #4 for crc 3 for sync
                        filei.seek(file_pos+hdr_len+header.msg_len+(4-3))
                        continue
                    elif fdate > self.date:
                        print 'Truncate file:', filei.name
                        break

                if header.msg_len < 8000:
                    buff2 = filei.read(header.msg_len)
                    buff = '\xAA\x44\x12'+buff0+buff1+buff2
                    try: crc_1 = unpack('=L', filei.read(4))[0]
                    except: continue
                    crc_2 = CRC32(buff) & 0xFFFFFFFF
                    if crc_1 <> crc_2:
                        filei.seek(file_pos, 0)
                        print 'CRC_error:', crc_1, crc_2, rec_date
                        continue
                    else:
                        is_dec = gps_sec%interval == 0
                        if ((header.msg_id == 43) and is_dec):
                            if gps_sec > gps_sec_last:
                                gps_sec_last = gps_sec
                                self.fileo.write(buff)
                                self.fileo.write(pack('=L', crc_1))
                                nrec += 1
                else:
                    filei.seek(file_pos, 0)
        self.nrec += nrec
        return gps_sec_last

    def add_lb2(self, fileo, gps_sec_last=0, interval=10):
        '''
        Add new binary id37 data files to the decimated new file,
        based in write_daily_id37_file.c of Robert Sheehan <sheehanp@bc.edu>.

        Inputs:
            fileo: binary id37 object file
            gps_sec_last: seconds (gps time) to start append new data
            interval: id37 sample interval in seconds
            run_date: datetime to restrict data for the given date.

        Outputs:
            gps_sec_last: seconds (gps time) for the last data block appened.
        '''

        nrec = 0
        sync = [0, 0]
        while True:
            try:
                ch = unpack('B', fileo.read(1))[0]
            except:
                break
            sync = [sync[1], ch]
            #print old_pos, sync
            if sync == [0x9C, 0xAE]:
                file_pos = fileo.tell()
                buff0 = fileo.read(2)
                if len(buff0)<2: break
                rec_len = unpack('H', buff0)[0]
                if rec_len >8192:
                    fileo.seek(file_pos, 0)
                    continue
                buff1 = fileo.read(1)
                if len(buff1)<1: break
                rec_id = unpack('B', buff1)[0]
                buff2 = fileo.read(rec_len-5)
                if len(buff2)<(rec_len-5): break
                buff3 = fileo.read(2)
                if len(buff3)<2: break
                checksum1 = unpack('h', buff3)[0]
                num_bytes = 3+rec_len-5
                checksum2 = -1*sum(unpack('%iB' % num_bytes, buff0+buff1+buff2))
                #print checksum1, checksum2
                if checksum1 <> checksum2:
                    print 'Checsum error'
                    fileo.seek(file_pos, 0)
                    continue
                if rec_id == 0x37:
                    station_prefix = buff2[:2]
                    i0 = 2
                    if buff2[3:6] == 'ver':
                        print 'ver True'
                        i0 += 10
                    #format = '%iB' % len(buff2[i0:])
                    format = '11B'
                    header = ReadRecord(format, buff2[i0:i0+11])
                    header.id37_header()
                    fdate = datetime(*header.gps_date)
                    if self.date:
                        if fdate < self.date:
                            fileo.seek(file_pos+rec_len, 0)
                            continue
                        elif fdate > self.date:
                            print 'Truncate file:', fileo.name
                            break
                    is_dec = header.gps_sec%interval == 0
                    if is_dec and (header.gps_sec>gps_sec_last):
                        gps_sec_last = header.gps_sec
                        self.fileo.write('\x9C\xAE'+buff0+buff1+buff2)
                        self.fileo.write(pack('h', checksum1))
                        nrec += 1
        self.nrec += nrec
        return gps_sec_last

    def add_obs(self, filei, gps_sec_last=0, interval=10):
        '''
        Add new obs data files to the decimated new file.

        Inputs:
            filei: input file object
            gps_sec_last: seconds (gps time) to start append new data
            interval: decimate value in seconds

        Outputs:
            gps_sec_last: seconds (gps time) for the last data block appened.
        '''

        nrec = 0
        if not self.date:
            try:
                week = GPSDateTime(filei.name.split('/')[-1], filename=True).GPSweek
            except:
                print '%s does not have a valid datestamp yymmdd, skipping...' % filei.name
                return gps_sec_last
        else:
            week = self.date.GPSweek
            gps_sec_min = week*SECONDS_IN_WEEK+int(self.date.strftime('%w'))*86400
            gps_sec_max = gps_sec_min+86400

        while True:
            buff = filei.readline()
            if not buff:
                break
            values = [x.strip() for x in buff.split(' ') if len(x)<>0]
            if len(values)<>11:
                continue
            try:
                week_sec = int(values[0])/1000
                dum = int(values[1])
                dum = int(values[2])
                dum = int(values[3])
            except:
                continue
            gps_sec = week*SECONDS_IN_WEEK + week_sec

            if self.date and gps_sec>=gps_sec_max:
                continue
            if self.date and gps_sec<gps_sec_min:
                continue            
            '''
            dow = week_sec/86400
            daysec = week_sec%86400
            hh = daysec/3600
            mm = (daysec%3600)/60
            ss = daysec%60
            rec_date = gdate(week*7 + dow + jday_gps0) + (hh, mm, ss)
            fdate = GPSDateTime(*rec_date[:6])
            #print fdate
            '''
            
            if gps_sec>gps_sec_last:
                self.last_prns = []
            is_dec = gps_sec%interval == 0            
            if gps_sec >= gps_sec_last and is_dec:                
                gps_sec_last = gps_sec
                if values[-1] not in self.last_prns:
                    self.last_prns.append(values[-1])
                    self.fileo.write(buff)
                    nrec += 1
        self.nrec += nrec
        return gps_sec_last

    def add_file(self, filei, last_sec, interval, line, daily=True):
        '''
        Add new text data files to the decimated new file.

        Inputs:
            filei: input file object
            gps_sec_last: seconds (gps time) to start append new data
            interval: decimate value in seconds

        Outputs:
            gps_sec_last: seconds (gps time) for the last data block appened.
        '''

        if self.first_line is None:
            if not self.date:
                self.date = GPSDateTime(filei.name.split('/')[-1], filename=True)
            if 'txt' in self.old_extension:
                self.doy        = self.date.doy                
                self.first_line = line
                if len(line.split())==7:
                    self.fix_line   = line[3:15].split(' ')
                    if not self.bad_doy:
                        self.bad_doy = int(line.split()[2]+'%s' % int(line.split()[3]))
                    self.fix_line = '%s%s%s' % (self.fix_line[0], self.fix_line[1],
                                                int(self.fix_line[2]))
                else:
                    print line
                    if not self.bad_doy:
                        self.bad_doy = int(line.split()[0][4:])
                    self.fix_line = line.split()[0]
            elif 's4' in self.old_extension or 'pos' in self.old_extension:
                self.first_line = '%4d %02d %02d\n' % (self.date.year,
                                    self.date.month, self.date.day)
                self.fix_line = '%02d %003d ' % (self.date.year%100,
                                                 self.date.doy)
            elif 'scn' in self.old_extension:
                self.first_line = 'T'
                self.fix_line = 'T %02d %02d %02d' % (self.date.year%100,
                                    self.date.month, self.date.day)

        if self.fix_line == '000000' and (self.date.strftime('%y%m%d') \
                not in filei.name):
            print 'skipping,', filei.name
            return last_sec
        nrec = 0
        if 'txt' in self.old_extension and len(line.split())==7:
            buff = filei.readline()
        while True:
            buff = filei.readline()
            if not buff:
                break
            #if re.findall('[^\x0a\x20-\x7f]', buff):
            #    print 'Error chars'
            #    continue
            try:
                values = [x.strip() for x in buff.split(' ') if x]
                if 's4.gz' in self.old_extension:
                    sec = int('%02d%05d' % (int(values[1]), int(values[2])))
                elif 'txt' in self.old_extension:
                    str_sec = values[1]
                    sec = int(values[0][4:])*100000 + (int(str_sec[:2])*60*60 + int(str_sec[2:4])*60 + int(str_sec[4:6]))
                    
                    if int(values[0][4:])>self.bad_doy:
                        self.bad_doy += 1
                        self.doy += 1
                        #sec += 100000
                    nsats = int(values[2])
                elif 'pos' in self.old_extension:
                    sec = int(values[2]+values[3])
                elif 'scn' in self.old_extension:
                    if values[0]<>'T':
                        continue
                    sec = int(''.join(values[1:]))
                else:
                    self.fileo.write(buff)
            except:
                traceback.print_exc(4)
                continue
            if self.fix_line in buff and self.nrec<>0 and sec==0:
                print 'first obs found'
                self.data  = ''
                last_sec   = -1
                nrec = 0
            if (self.fix_line in buff or not daily) and sec>last_sec:            
                if 'txt' in self.old_extension:
                    #if nrec>1 and sec-last_sec>3000: continue
                    fmt0 = '%02d %03d% 6d% 3d'
                    fmt1 = ['% 3d %3.2f %4.1f %3.1f' for i in range(nsats)]
                    fmt  = fmt0+''.join(fmt1)+'\n'
                    str_sec = values[1]
                    vals = (GPSDateTime(filei.name.split('/')[-1], filename=True).year%100, self.doy, sec-int(values[0][4:])*100000, nsats)
                    vals += tuple([float(s.replace(',', '.')) for s in values[3:] if s])
                    if len(fmt.split())<>len(vals):
                        continue
                    self.data += fmt % vals
                elif 'scn' in self.old_extension:
                    self.data += buff
                    while True:
                        nline = filei.readline()
                        if 'T' in nline:
                            filei.backline()
                            break
                        elif nline:
                            self.data += nline
                        else:
                            break
                else:
                    self.data += buff.strip()+'\n'
                last_sec = sec
                nrec += 1
        self.nrec += nrec
        return last_sec

    def close(self, gz=False, header=None):
        '''
        Save and close files, create path if not exist
        '''
        self.fileo.seek(0)
        if self.nrec <> 0:
            path = os.path.dirname(self.file_name)
            if not os.path.exists(path):
                if path<>'':
                    os.makedirs(os.path.dirname(self.file_name))
            if gz:
                new_file = gzip.open(self.file_name+'.gz', 'wb')
            else:
                new_file = open(self.file_name, 'wb')
            if self.data:
                new_file.write(self.data)
            else:
                new_file.write(self.fileo.read())
            new_file.close()
        self.fileo.close()
        print 'Saved %d records' % self.nrec

class Yuma(dict):
    '''
    Class to read an hold data from almanacs files (download if neccesary)
    in yuma format.
    Search almanac file for the given date (datetime) and path.
    '''
    def __init__(self, date, path=None, N=31, local=False):
        '''
        '''
        
        self.date = GPSDateTime(date)
        if not path:
            path = os.path.join(localpath, 'almanac')
        self.path = path
        if not os.path.exists(path):
            os.makedirs(path)
        for i in range(N):
            idate = self.date - timedelta(i)
            idoy = idate.doy
            fname = '%d%03d.ALM' % (idate.year, idoy)
            if os.path.exists(os.path.join(path, fname)):
                break
            fname = '%d%03d.alm' % (idate.year, idoy)
            if os.path.exists(os.path.join(path, fname)):
                break
            if not local:
                fname = self.download(idate)
                if fname:
                    break
            fname = False
        if not fname:
            raise RuntimeError('No almanac file found for %s in %s' \
                               % (self.date.date(), path))
        else:
            fid = open_file(os.path.join(path, fname))
            while True:
                try:
                    line = fid.next()
                except StopIteration:
                    break
                if not line: continue
                if line[:2]=='**':
                    ID = 'G%02d' % int(fid.next().split(':')[1])
                    self[ID]={}
                    self[ID]['HEALTH'] = float(fid.next().split(':')[1])
                    self[ID]['ECC']    = float(fid.next().split(':')[1])
                    self[ID]['TOA']    = float(fid.next().split(':')[1])
                    self[ID]['INC']    = float(fid.next().split(':')[1])
                    self[ID]['RRAND']  = float(fid.next().split(':')[1])
                    self[ID]['SMA']    = float(fid.next().split(':')[1])
                    self[ID]['RAND']   = float(fid.next().split(':')[1])
                    self[ID]['AOP']    = float(fid.next().split(':')[1])
                    self[ID]['AM']     = float(fid.next().split(':')[1])
                    self[ID]['AF0']    = float(fid.next().split(':')[1])
                    self[ID]['AF1']    = float(fid.next().split(':')[1])
                    week = int(fid.next().split(':')[1])
                    realweek = idate.GPSweek
                    if week<realweek:
                        week+=1024*(realweek/1024)
                    self[ID]['WK'] = week

    def download(self, date):
        '''
        '''
        toa_1 = ['', '380928', '466944', '552960', '032768', '118784', '']
        toa_2 = ['319488', '405504', '503808', '589824', '061440',
                 '147456', '233472']

        if date<datetime(1996,3,31): toa = toa_1
        else: toa = toa_2

        URL1 = 'http://www.navcen.uscg.gov/?Do=getAlmanac&almanac=%03d&year=%d&format=yuma' \
                % (date.doy, date.year)
        URL2 = 'http://celestrak.com/GPS/almanac/Yuma/%d/almanac.yuma.week%04d.%s.txt' \
                % (date.year, (date.GPSweek)%1024, toa[date.GPSdow])
        URL3 = 'https://gps.afspc.af.mil/gps/archive/%d/almanacs/yuma/wk%d%s.alm' \
                % (date.year, date.GPSweek, date.strftime('%a').lower())
        urls = [URL1, URL2, URL3]

        for url in urls:
            err = True
            try:
                print 'Downloading almanac file from %s' % url
                data = urlopen(url)
                line = data.readline()
                if '****' in line and 'Week' in line:
                    err = False
                    break
            except:
                err = True

        if err: return False
        fname = '%d%03d.ALM' % (date.year, date.doy)
        new_alm = open(os.path.join(self.path, fname), 'w')
        new_alm.write(line+data.read())
        new_alm.close()
        data.close()
        return fname

    def sat_pos(self, prns, dt, user_pos=(0, 0, 0)):
        '''
        Translate from gps.c of "The Essential GNSS Project"
        Copyright (c) 2007 author: doxygen
        '''
        
        # iterate to include the Sagnac correction
        # since the range is unknown, an approximate of 70 ms is good enough
        # to start the iterations so that 2 iterations are enough

        ret   = []
        rango = 0.07*LIGHT_SPEED
        userX, userY, userZ = user_pos
        for prn in prns:
            for i in xrange(2):
                Xs, Ys, Zs = self.calc_sat_pos(prn, dt, rango)
                dx         = Xs - userX
                dy         = Ys - userY
                dz         = Zs - userZ
                rango      = math.sqrt( dx*dx + dy*dy + dz*dz )
            ret.append([Xs, Ys, Zs])

        return ret

    def calc_sat_pos(self, prn, dt, estimate_range):
        '''
        '''
        
        toe = self[prn]['TOA']
        m0  = self[prn]['AM']
        delta_n = 0
        tot = dt.GPSweek*SECONDS_IN_WEEK+dt.GPStow
        tk  = tot - (self[prn]['WK']*SECONDS_IN_WEEK+toe)
        # compute the corrected mean motion
        a   = self[prn]['SMA']**2
        n   = math.sqrt(GPS_GRAVITY_CONSTANT/(a*a*a))
        n  += delta_n
        # kepler's equation for eccentric anomaly
        M = m0 + n*tk
        E = M
        for i in range(8):
            E = M + self[prn]['ECC']*math.sin(E)

        cosE = math.cos(E)
        sinE = math.sin(E)
        # true anomaly
        v = math.atan2(math.sqrt(1.0-self[prn]['ECC']**2)*sinE, (cosE - self[prn]['ECC']))
        #argument of latitude
        u = v + self[prn]['AOP']
        # radius in orbital plane
        r = a * (1.0 - self[prn]['ECC']*math.cos(E))
        # orbital inclination
        i = self[prn]['INC']

        #argument of latitude correction
        #! no correction
        cosu = math.cos(u)
        sinu = math.sin(u)
        x_op = r*cosu
        y_op = r*sinu

        omegak = self[prn]['RAND'] + (self[prn]['RRAND'] - GPS_WGS84_EARTH_ROTATION_RATE)*tk - GPS_WGS84_EARTH_ROTATION_RATE*(toe + estimate_range/LIGHT_SPEED )

        # compute the WGS84 ECEF coordinates,
        # vector r with components x & y is now rotated using, R3(-omegak)*R1(-i)
        cos_omegak = math.cos(omegak)
        sin_omegak = math.sin(omegak)
        cosi = math.cos(i)
        sini = math.sin(i)
        x = x_op * cos_omegak - y_op * sin_omegak * cosi
        y = x_op * sin_omegak + y_op * cos_omegak * cosi
        z = y_op * sini

        return x,y,z

    def get_sat_clock_correction(self, prns, dt):
        '''
        '''
        ret = []
        for prn in prns:
            toe = self[prn]['TOA']
            toc = self[prn]['TOA']
            m0  = self[prn]['AM']
            delta_n = 0
            tot = dt.GPSweek*SECONDS_IN_WEEK+dt.GPStow
            tk  = tot - (self[prn]['WK']*SECONDS_IN_WEEK+toe)
            tc  = tot - (self[prn]['WK']*SECONDS_IN_WEEK + toc)
            # compute the corrected mean motion
            a   = self[prn]['SMA']**2
            n   = math.sqrt(GPS_GRAVITY_CONSTANT/(a*a*a))
            n  += delta_n
            # kepler's equation for eccentric anomaly
            M = m0 + n*tk
            E = M
            for i in xrange(8):
                E = M + self[prn]['ECC']*math.sin(E)

            #relativistic correction
            d_tr  = GPS_CLOCK_CORRECTION_RELATIVISTIC_CONSTANT*self[prn]['ECC']*self[prn]['SMA']*math.sin(E)
            d_tr *= LIGHT_SPEED

            #clock correction
            d_tsv  = self[prn]['AF0']+self[prn]['AF1']*tc#+af2*tc*tc
            ret.append(d_tsv*LIGHT_SPEED+d_tr)

            #clock drift
            #clock_drift = self[prn]['AF1']*LIGHT_SPEED

        return ret

class gpselaz(dict):
    '''
    Class to calculate el, az, lat and lon of satellites, uses yuma almanac info
    '''
    def __init__(self, date, path=None):
        yuma = Yuma(date, path)
        U = 398603.2
        self.date = GPSDateTime(date)
        self.elemnts = {}
        for id in yuma:
            IDate = GPSDateTime(GPSweek=yuma[id]['WK'])
            ragren = self.ra_gren(IDate)
            SMA = yuma[id]['SMA']*yuma[id]['SMA']/1000.0
            #***********************************************************
            #*    CALCULATE PERIOD IN MINUTES FROM SMA                 *
            #*    FORMULA   SMA=(U**1/3)*((86400.0/TWO*PI*N)**2/3)     *
            #*              SMA=SEMI-MAJOR AXIS(KM)                    *
            #*              U=398603.2                                 *
            #*              N=PERIOD(REVS/DAY)                         *
            #***********************************************************
            N = 86400.0/(2*PI)*math.sqrt(U/SMA**3)
            N = 1440.0/N
            #calculate argument of perigee
            AOP = fmod(yuma[id]['AOP'])*rad2deg
            #calculate inclination
            INC = fmod(yuma[id]['INC'])*rad2deg
            RAAN = (yuma[id]['RAND']*rad2deg + ragren)*deg2rad
            RAAN = fmod(RAAN)*rad2deg
            #calculate mean anomaly at 00 hrs
            XMO = 1440.0/N*360.0/86400.0/rad2deg*yuma[id]['TOA']
            AM = fmod(yuma[id]['AM'] - XMO)*rad2deg
            self.elemnts[id] = [IDate, SMA, N, yuma[id]['ECC'], INC, AM, AOP,
                                RAAN]

    def update(self, flat, flon, falt, step=600):
        '''
        Set de location step in seconds
        '''
        #self = {}
        self.step = step
        jyear     = self.date.year
        jmonth    = self.date.month
        jday      = self.date.day
        ndoy      = self.date.doy

        #NAVSTAR PASS PLANNER
        ST2 = [0, 0, 0]
        ST2[0] = (90-flat)*deg2rad
        ST2[1] = flon*deg2rad
        ST2[2] = falt/1000+6370         #earth radius = 6370 km
        #Change IRPT setting according to time step
        IRPT = 86400/step

        if jyear%4 == 0: YRDAYS = 366
        else: YRDAYS = 365

        TRSS  = (ndoy-1)*24.0*60.0
        SILAT = 0
        SILON = 0

        for prn, elem in self.elemnts.items():
            #prn = 'G%02d' % sat
            self[prn] = {}
            self[prn]['hour'] = []
            self[prn]['az'] = []
            self[prn]['el'] = []
            self[prn]['lat'] = []
            self[prn]['lon'] = []
            [TEP,IEY,AO,EO,ONO,OPO,XIO,XLO,XNO,THGR] = self.KBNV(*elem)
            TW  = TRSS - TEP
            TWS = TW*60.0
            x365 = 365
            if IEY%4 == 0: x365 = 366
            for I in range(int(IRPT)):
                RNGP = 0
                TW = TWS/60.0
                TP = TW + TEP
                [JDY,IHR,AMN] = min2date(TP)
                JDYP = JDY
                if (JDY > YRDAYS): JDYP = JDY - YRDAYS
                IMN = int(AMN)
                FMN = AMN - float(IMN)
                SC = FMN*60.0
                TW2W = TW
                if (TW2W < -28800.0):
                    TW2W = TW2W + x365*1440
                [BGR,GCR] = self.SGPN(TW2W,AO,EO,ONO,OPO,XIO,XLO,XNO,THGR)
                [SL2,D] = self.SPGM(ST2,BGR)
                EL2    = 90.0 - SL2[0]*rad2deg
                A1     = EL2
                RNG    = SL2[2]
                if (A1 < -3.0):
                    RNGP = RNG
                    TWS  = TWS + step
                    continue
                ELS    = A1
                SATLAT = GCR[0]
                SATLON = GCR[1]
                SATHT  = GCR[2]
                ELR    = EL2*deg2rad
                AZR    = (SL2[1])
                AZS    = AZR*rad2deg
                if (AZS < 0.0):
                    AZS = AZS + 360.0
                RDOT   = (RNG-RNGP)/step
                XLAM   = 300.0/1575.42
                DOP    = -(RDOT/XLAM)
                if (ELS > 0):
                    [SILAT,SILON,RANGE] = SILL(ELS,AZS,350.0,flon,flat)
                self[prn]['hour'].append(IHR+IMN/60.)
                self[prn]['az'].append(AZS)
                self[prn]['el'].append(ELS)
                self[prn]['lat'].append(SILAT)
                self[prn]['lon'].append(SILON)
                RNGP = RNG
                TWS  = TWS + step

    def find(self, prn, itime):
        '''
        Find azimuth, elevation, latitud and longitud, given PRN and time.
        '''
        factor = 3600/self.step
        xaz, xel, xlat, xlon = 0, 0, 0, 0
        eatimes = self[prn]['hour']
        kk = len(eatimes)-1
        dif = [abs(itime-t) for t in eatimes]
        idx = dif.index(min(dif))
        if idx==0:
            i0, i1, i2 = 0, 0, 1
        elif idx==kk:
            i0, i1, i2 = kk, kk-1, kk
        else:
            idif1 = abs(itime - eatimes[idx-1])
            idif2 = abs(itime - eatimes[idx+1])
            i0, i1, i2 = idx - 1, idx -1, idx
            if idif2 < idif1:
                i0, i1, i2 = idx, idx, idx +1
        az1 = self[prn]['az'][i1]
        az2 = self[prn]['az'][i2]
        if (az1-az2) > 330.0:
            if az1 < az2:
                az1 = az1 + 360.0
            elif az1 > az2:
                az1 = az1 - 360.0
        xaz = az1 + (az2 - az1)*(itime-eatimes[i1])*factor
        xel = self[prn]['el'][i0] + (self[prn]['el'][i2] - \
                self[prn]['el'][i1])*(itime-eatimes[i0])*factor
        xlat = self[prn]['lat'][i0] + (self[prn]['lat'][i2] - \
                self[prn]['lat'][i1])*(itime-eatimes[i0])*factor
        xlon = self[prn]['lon'][i0] + (self[prn]['lon'][i2] - \
                self[prn]['lon'][i1])*(itime-eatimes[i0])*factor
        return xaz, xel, xlat, xlon

    def ra_gren(self, idate):
        '''
        Returns the right ascension at Greenwich on idate at 00 hours.
        '''
        dayjul0 = 2433282.5
        sec1 = 8640184
        sec2 = 0.812866
        sec3 = 0.093104
        sec4 = 6.2e-6
        consecs = 24110.54841

        yr = idate.year - 1900
        mn = idate.month
        dy = idate.day

        if mn <= 2:
            mn += 9
            yr -= 1
        else:
            mn -= 3

        #days elapsed from FEB. 29, 1900 TO JAN. 1, 1950 = 18204
        dayjul = ((1461*yr)/4) + (((153*mn)+2)/5) + dy - 18204 + dayjul0
        T = (dayjul - 2451545.)/36525.
        AA1 = sec1*T + T*(sec2 + T*(sec3 - sec4*T))
        ragren = consecs/240. + AA1/240.
        ragren = ragren%360
        if ragren < 0:
            ragren = ragren + 360
        return ragren

    def KBNV(self, idate, AO, PERM, EO, XI, G, OPO, ONO):
        '''
        Disc file read of GPS Keplerian elements
        elem = [IDate,SMA,N,ECC,INC,AM,AOP,RAAN]
        '''

        THGS = 99.24377
        SIDR = 0.985647346
        SDAY = 7304
        IED = idate.doy
        IEY = idate.year-1900
        if EO == 0:
            return [None for x in range(12)]
        else:
            IEPT = IED - 1 + (IEY - 50)*365 + (IEY - 49)/4
            FEPT = float(IEPT)
            THGR = ((FEPT - SDAY)*SIDR + THGS)*deg2rad
            THGR = fmod(THGR)
            TE   = (IED - 1)*1440.0
            XNO  = 2*PI/PERM
            ONO  = ONO*deg2rad
            XIO  = XI*deg2rad
            OPO  = OPO*deg2rad
            G    = G*deg2rad
            XLO  = fmod(ONO + OPO + G)
            AO   = AO/6378.166
            return TE,IEY,AO,EO,ONO,OPO,XIO,XLO,XNO,THGR

    def SGPN(self, T,AO,EO,ONO,OPO,XIO,XLO,XNO,THGR):
        '''
        THIS IS A VERSION OF THE NORAD
        SIMPLIFIED GENERAL PERTURBATIONS (SGP)
        CODE FOR THE PROPAGATION OF NORAD MEAN
        ORBITAL ELEMENTS TO TIME T...
        ...INITIALIZATION BY SUBROUTINE KBNO
        VERSION FOR NAG ELEMENTS...
        '''

        GCN=[0,0,0]
        XJ2 =0.0010826158
        SIDD = 0.4375269512707e-2
        XL=XLO+T*XNO
        XL=fmod(XL)
        ON=fmod(ONO)
        OP=fmod(OPO)
        E=EO
        SOP=math.sin(OP)
        COP=math.cos(OP)
        AXN=E*COP
        AYN=E*SOP
        SI=math.sin(XIO)
        CI=(XIO)
        ON=fmod(ON)
        OP=fmod(OP)
        XL=fmod(XL)
        U=fmod(XL-ON)
        TP1=XJ2*(1.0-1.5*SI*SI)/(2.0*AO*AO*(1.0-E*E))
        A=AO*(1.0+TP1)
        IT=50
        TP1=U
        TP2=1.0
        E6A=1.0e-8*TP1+1.0e-6
        while (abs(TP2)>=E6A):
            SINEO=math.sin(TP1)
            COSEO=math.cos(TP1)
            IT=IT-1
            if (IT<=0):
                err=True
            TP5=1.0-COSEO*AXN-SINEO*AYN
            TP5=(U-(-SINEO*AXN+COSEO*AYN+TP1))/TP5
            if (1.0<abs(TP5)):
                TP5=abs(TP5)/TP5
            TP2=TP5
            TP1=TP2+TP1
        ESE=-AYN*COSEO+AXN*SINEO
        ECE=AXN*COSEO+AYN*SINEO
        ESQ=AXN*AXN+AYN*AYN
        E=math.sqrt(ESQ)
        R=A*(1.0-ECE)
        P=A*(1.0-ESQ)
        TP1=1.0+math.sqrt(1.0-ESQ)
        SSU=ESE*AXN/TP1
        CSU=(ESE/TP1*AYN-AXN+COSEO)/R*A
        SSU=A*(-SSU-AYN+SINEO)/R
        SU=self.ACTN(SSU,CSU)
        C2SU=2.0*SU
        S2SU=math.sin(C2SU)
        C2SU=math.cos(C2SU)
        TP7=XIO
        SI=math.sin(TP7)
        CI=math.cos(TP7)
        S2I=math.sin(2.0*TP7)
        FA=3.0*XJ2/(2.0*P*P)
        DR=math.atan(SSU/CSU)
        RD=math.atan(SSU/CSU)
        RD=abs(math.cos(RD))
        DLR=RD*FA/3.0*SI*SI*P
        DLSU=-(-7.0*SI*SI+6.0)/12.0*S2SU*FA
        DLON=FA*S2SU/2.0*CI
        DLXI=S2I*C2SU/4.0*FA
        R=R+DLR
        SU=fmod(SU+DLSU)
        ON=fmod(ON+DLON)
        SOP=AYN/E
        COP=AXN/E
        OP=self.ACTN(SOP,COP)
        XI=XIO+DLXI
        SON=math.sin(ON)
        CON=math.cos(ON)
        SI=math.sin(XI)
        CI=math.cos(XI)
        SSU=math.sin(SU)
        CSU=math.cos(SU)
        X=CSU
        Y=CI*SSU
        Z=SI*SSU
        HRA=THGR+T*SIDD
        HRA=fmod(HRA)
        GCN[1]=self.ACTN(Y,X)+ON-HRA
        GCN[1]=fmod(GCN[1])
        XY=math.sqrt(X*X+Y*Y)
        GCN[0]=self.ACTN(XY,Z)
        GCN[2]=R*6378.166
        [GCD,GCR] = TRNN(GCN[0], GCN[1], GCN[2])
        return GCD,GCR

    def SPGM(self, RO, RT):

        RB=[0,0,0]
        CTO = math.cos(RO[0])
        STO = math.sin(RO[0])
        CTT = math.cos(RT[0])
        STT = math.sin(RT[0])
        CHD = CTO*CTT+STO*STT*math.cos(RT[1]-RO[1])
        SHD = math.sqrt(1.0-CHD*CHD)
        GAM = math.atan2(SHD,CHD)
        SBE = STT*math.sin(RT[1]-RO[1])/SHD
        CBE = (CTT-CTO*CHD)/(STO*SHD)
        RB[1] = math.atan2(SBE,CBE)
        RB[2] = math.sqrt(RO[2]*RO[2]+RT[2]*RT[2]-2.0*RO[2]*RT[2]*CHD)
        SEL   = RT[2]*SHD/RB[2]
        CEL   = (RT[2]*CHD-RO[2])/RB[2]
        RB[0] = math.atan2(SEL,CEL)
        #print 'spgm ',RO,RT
        return RB,GAM

    def ACTN(self, SX, CX):
        #ACTN = 0.0
        if (CX == 0.0):
            if (SX == 0.0):
                return 0.0
            elif (SX > 0.0):
                return PI/2
            else:
                return 3*PI/2
        elif (CX > 0.0):
            if (SX == 0.0):
                return 0.0
            elif (SX > 0.0):
                return math.atan(SX/CX)
            else:
                return 2.0*PI + math.atan(SX/CX)
        else:
            return PI + math.atan(SX/CX)

class S4Data(Data):
    '''
    '''

    def __init__(self, filename, station={}, vars='all'):
        '''
        '''
        self.data = {}
        self.epoch = GPSDateTime(2100,1,1)
        self.prns = set()
        self.bad_prn = []
        self.rec_bias = 0
        self.station = station
        if 'code' not in station:
            self.station['code'] = filename.split('/')[-1][:4]
        self.arcs = {}
        fo = open_file(filename)
        self.header = fo.next()
        if len(self.header.split())>6:
            fo.reset()
        while True:
            try:
                line = fo.next()
            except StopIteration:
                break
            if line[0] == '#' or len(line)==0: continue
            rec = S4Record(line)
            if rec.epoch is None:
                continue
            self.data[rec.epoch] = rec
            self.rec = rec           
            self.checkbreak(rec.epoch-self.epoch)
        
        self.dates = self.data.keys()
        self.dates.sort()
        self.intervals = set([self.dates[i]-self.dates[i-1] for i in range(1, len(self.dates))])
        self.interval = min(self.intervals)
        self.endphase(self.prns, True)
        self.check()
        fo.close()
        self.date = self[-1].epoch
        self.vars = set(['epoch', 's4', 'ele', 'azi'])

    def checkbreak(self, delta):
        '''
        '''
        self.prns.update(self.rec)
        
        for prn in set(self.rec).union(self.arcs):           
            if prn not in self.rec:
                self.endphase(prn)
                continue
            if delta>300:
                self.breakphase(prn)
                continue
            if prn not in self.arcs:
                self.arcs[prn] = [[len(self.data)-1, None]]
                continue
            if self.arcs[prn][-1][1] is not None:
                self.arcs[prn] += [[len(self.data)-1, None]]
                continue
            if prn not in self.arcs or self.arcs[prn][-1][1] is not None:
                continue

    def check(self):
        '''
        '''
        dum = {}
        for prn in self.arcs:
            dum[prn] = [arc for arc in self.arcs[prn] \
                                   if arc[1]-arc[0]>=10]
        self.arcs = dum

    def iterlist(self, prn, keys=['s4']):
        '''
        '''
        def chooser(rec, prn, keys):
            tup = (rec.epoch,)
            for key in keys:
                tup += (rec[prn][key],)
            return tup

        rec_idx = []
        for arc in self.arcs[prn]:
            rec_idx += [x for x in range(arc[0], arc[1])]
        
        for x in rec_idx:
            yield chooser(self[x], prn, keys)
        return
    
    def calcll(self):
        '''
        '''
        if 'latitude'not in self.station or 'longitude' not in self.station:
            raise RuntimeError('Missing station latitude and longitude')
        SLAT, SLON = self.station['latitude'], self.station['longitude']
        
        for prn in self.arcs:
            for rec in self.get_records(prn):
                xlat,xlon,rang = SILL(rec[prn]['ele'], rec[prn]['azi'], 350.0, SLON, SLAT)
                rec[prn]['lat'] = xlat
                rec[prn]['lon'] = xlon
        
        self.vars.add('lat')
        self.vars.add('lon')
    
    def plot(self, X, Y, **kwargs):
        '''
        '''
        return plot_data_vars(self, X, Y, **kwargs)

class S4Record(dict):
    '''
    '''
    def __init__(self, line):
        '''
        '''
        cols = [x.strip() for x in line.split(' ') if len(x)<>0]
        self.epoch = None
        try:
            tot_sec = int(cols[2])
            fyear   = int(cols[0])
            if fyear<100:
                fyear = fyear+2000
            else:
                fyear = fyear+1900
            fdoy    = int(cols[1])
            fhour   = tot_sec/3600
            fmin    = (tot_sec-(fhour*3600))/60
            fsec    = tot_sec-(fhour*3600)-(fmin*60)
            if len(cols)<>int(cols[3])*4+4:
                return
            for nsat in range(int(cols[3])):
                prn = 'G%02d' % int(cols[nsat*4+4])
                if prn not in self:
                    self[prn] = EpochDict()
                self[prn]['s4'] = float(cols[nsat*4+5].replace(',','.'))
                self[prn]['azi'] = float(cols[nsat*4+6].replace(',','.'))
                self[prn]['ele'] = float(cols[nsat*4+7].replace(',','.'))
        except (IndexError, ValueError):
            return
        self.epoch = GPSDateTime(year=fyear, doy=fdoy, hour=fhour, minute=fmin,
                                 second=fsec)
                

    def __getitem__(self, index):
        '''
        Allow you to access GPS satellites, eg record['G13'], as
        simply record[13].  For GLONASS or Galileo, you must use the full code.
        '''
        if index == 'epoch':
            return self.epoch
        elif isinstance(index, (int, long, float)):
            return dict.__getitem__(self, 'G%02d' % index)
        else:
            return dict.__getitem__(self, index)

class my_float(float):
    pass

class my_str(str):
    pass

class Stats(object):
    '''
    Class to calculate statistical values from a list
    '''
    def __init__(self, L):
        self.n = float(len(L))
        self.sum = sum(L)
        if len(L)>1:
            self.average = self.sum/self.n
            self.variance = sum([(x-self.average)**2 for x in L])/(self.n-1)
            self.std_dev = self.variance**0.5
        elif len(L)==1:
            self.average = L[0]
            self.variance = 0
            self.std_dev = 0
        else:
            raise RuntimeError('List is empty')


class MkRinex(object):
    '''
    Class useful to create a rinex file
    '''

    def __init__(self, marker, interval, obscodes):
        self.marker   = marker
        self.interval = interval
        self.obscodes = obscodes

    def mk_header(self, firsttime, **kwargs):
        '''
        Return a string with a rinex header

        first_obs = date for the first observation [datetime]
        '''

        program   = '%s.py v%4.2f' % ('lisn_utils', 2)
        firsttime = firsttime.timetuple()[:6]
        obscodes  = '%6i' % len(self.obscodes)        
        for n,obs in enumerate(self.obscodes):
            if n>1 and n%9==0:
                obscodes += '# / TYPES OF OBSERV\n'
                obscodes += '%6s' % ''
            obscodes += '%6s' % obs        
        
        fmt = '%%%is' % ((9-n%9-1)*6+19) 
        obscodes +=  fmt % '# / TYPES OF OBSERV'
        run_date   = GPSDateTime.now().strftime('%Y/%m/%d %H:%M:%S')

        run_by       = kwargs.get('run_by', lisn_header['run_by'])
        observer     = kwargs.get('observer', lisn_header['observer'])
        agency       = kwargs.get('agency', lisn_header['agency'])
        receivernum  = kwargs.get('receivernum', '')
        receivertype = kwargs.get('receivertype', '')
        receiverver  = kwargs.get('receiverver', '')
        antennanum   = kwargs.get('antennanum', '')
        antennatype  = kwargs.get('antennatype', '')
        xyz          = kwargs.get('xyz', (0, 0, 0))+('APPROX POSITION XYZ', )
        delta        = kwargs.get('delta', (0, 0, 0))+('ANTENNA: DELTA H/E/N', )
        interval     = kwargs.get('interval', self.interval)
        endtime      = kwargs.get('endtime')
        comments     = kwargs.get('comment', '')
        comments     = lisn_comment+comments
        comment      = ''
        for line in StringIO(comments):
            if line:
                comment += '%-60sCOMMENT\n' % line[:60].strip()

        header  = '%9s%38s%33s\n' % ('2.10', 'OBSERVATION DATA    G (GPS)',
                                     'RINEX VERSION / TYPE')
        header += '%-20s%-20s%-20sPGM / RUN BY / DATE\n' % (program[:20],
                                                            run_by[:20],
                                                            run_date)
        header += '%s' % comment
        header += '%-60sMARKER NAME\n' % self.marker
        header += '%-20s%-40sOBSERVER / AGENCY\n' % (observer[:20], agency[:40])
        header += '%-20s%-20s%-20sREC # / TYPE / VERS\n' % (receivernum,
                                                            receivertype,
                                                            receiverver)
        header += '%-20s%-20s%-20sANT # / TYPE\n' % (antennanum, antennatype, '')
        header += '%14.4f%14.4f%14.4f%37s\n' % xyz
        header += '%14.4f%14.4f%14.4f%38s\n' % delta
        header += '%6i%6i%68s\n' % (1, 1, 'WAVELENGTH FACT L1/2')
        header += '%s\n' % obscodes
        header += '%10.3f%58s\n' % (interval, 'INTERVAL')
        header += '%6d%6d%6d%6d%6d%13.7f     GPS         TIME OF FIRST OBS\n' % firsttime
        if endtime:
            endtime = endtime.timetuple()[:6]
            header += '%6d%6d%6d%6d%6d%13.7f     GPS         TIME OF LAST OBS\n' % endtime
        header += '%73s\n' % 'END OF HEADER'
        return header

    def mk_satline(self, epoch, prn_list):
        '''
        Create rinex observation record
        '''

        sat_line = ' %02d %2d %2d %2d %2d%11.7f  0%3d' % (epoch.year%100,
                        epoch.month, epoch.day, epoch.hour, epoch.minute,
                        epoch.second, len(prn_list))
        for i, prn in enumerate(prn_list):
            if i>0 and i%12==0:
                sat_line += '\n % 34s' % prn
            else:
                sat_line += '%s' % prn
        return sat_line + '\n'

    def mk_recordline(self, obsvalues, LLI, SNR):
        '''
        Create a rinex observations given obs value, LLI (lost of lock
        indicator) and SNR flag

        output: rinex record in text format
        '''
        format = ''
        nrec   = 0

        for code in self.obscodes:
            #if code not in rnx_codes:
            #    continue
            if nrec%5==0 and nrec>1:
                format+='\n'
            format += '%14.3f%s%s'
            nrec+=1
        format+='\n'
        values = [x[i] for x in zip(obsvalues, LLI, SNR) for i in range(3)]
        return format % tuple(values)

    def snr_flag(self, rec, code=None):
        '''
        '''

        if code=='L1' and 'S1' in self.obscodes and 'S1' in rec:
            snr = rec['S1']
        elif code=='L2' and 'S2' in self.obscodes and 'S2' in rec:
            snr = rec['S2']
        else:
            return 0
        if snr<0:
            return 0
        flags = [0, 26, 28, 32, 36, 39, 42, 45, 49]
        value = [x for x in flags if snr>=x][-1]
        return flags.index(value)+1

class Hatanaka(object):
    '''
    Class to handled hatanaka compress and uncompress files
    '''
    def __init__(self):
        '''
        '''

        ###check for rinex to crinex command
        if os.path.exists('/usr/local/bin/rnx2crx'):
            self.rnx2crx = '/usr/local/bin/rnx2crx'
        else:
            self.rnx2crx = None

        if os.path.exists('/usr/local/bin/crx2rnx'):
            self.crx2rnx = '/usr/local/bin/crx2rnx'
        else:
            self.crx2rnx = None

        if os.path.exists('/usr/bin/compress'):
            self.zip_cmd = '/usr/bin/compress'
        elif os.path.exists('/usr/local/bin/compress'):
            self.zip_cmd = '/usr/local/bin/compress'
        elif os.path.exists('/bin/compress'):
            self.zip_cmd = '/bin/compress'
        else:
            self.zip_cmd = None
            print 'compress command not found, please install ncompress '

        if os.path.exists('/usr/bin/uncompress'):
            self.unzip_cmd = '/usr/bin/uncompress'
        elif os.path.exists('/usr/local/bin/uncompress'):
            self.unzip_cmd = '/usr/local/bin/uncompress'
        elif os.path.exists('/bin/uncompress'):
            self.unzip_cmd = '/bin/uncompress'
        else:
            self.unzip_cmd = None
            print 'uncompress command not found, please install ncompress '

    def compress(self, file, Z=False, rem=False):
        '''
        Create a CRINEX file from RINEX
        '''
        if self.rnx2crx:
            print 'Creating crinex file: %s' % file[:-1]+'d'
            os.system('cat %s | %s - > %s' % (file, self.rnx2crx,
                                                file[:-1]+'d'))
            if os.path.exists(file[:-1]+'d'):
                new_file = file[:-1]+'d'
                if Z and self.zip_cmd:
                    print 'Compressing file : %s.Z' % new_file
                    os.system('%s -f %s' % (self.zip_cmd, new_file))
                    new_file += '.Z'
                if rem:
                    os.remove(file)
                return new_file
            else:
                return file
        else:
            print 'RNX2CRX not found'
            return file

    def uncompress(self, file):
        '''
        Create RINEX file from CRINEX
        '''
        
        if file.endswith('.Z') and self.unzip_cmd:
            os.system('%s -f %s' % (self.unzip_cmd, file))
            if os.path.exists(file[:-2]):
                file = file[:-2]
                
        if self.crx2rnx and os.path.exists(file) and file.endswith('d'):
            print 'Creating rinex file: %s' % file[:-1]+'o'
            os.system('cat %s | %s - > %s' % (file, self.crx2rnx,
                                                file[:-1]+'o'))
            if os.path.exists(file[:-1]+'o'):
                return file
        else:
            print 'Invalid file: %s' % file
            return False

def versioncheck(ver):
    '''
    Given RINEX format version ver, verify that this program can handle it.
    '''
    nums = ver.split('.')
    if not 0 < len(nums) < 3:
        raise ValueError('RINEX Version not parsable')
    elif int(nums[0]) != 2:
        raise IOError('RINEX File not version 2; unsupported')
    elif len(nums) > 1 and int(nums[1]) > 11:
        warn('RINEX minor version more recent than program.')
    return ver.strip()


def crxcheck(ver):
    '''
    Check whether Compact RINEX version is known to this program.
    '''
    if ver != '1.0':
        raise ValueError('CRINEX version ' + ver + ' not supported.')
    return ver.strip()


def iso(c):
    '''
    Ensure that the character c is `O' (for RINEX observation data.)
    '''
    if c.upper() != 'O':
        raise IOError('RINEX File is not observation data')
    return c.upper()

def headertime(s):
    '''
    Parse RINEX header times into GPSDateTime objects
    '''
    if not s.strip():
        return None
    year   = int(s[0:6])
    if year<1900: year += 1900
    month  = toint(s[6:12])
    day    = toint(s[12:18])
    hour   = toint(s[18:24])
    minute = toint(s[24:30])
    second = tofloat(s[30:43])
    usec   = (second - int(second)) * 1000000
    return GPSDateTime(year, month, day, hour, minute, int(second), int(usec))

def recordtime(s, baseyear=None):
    '''
    Parse RINEX time epoch into GPSDateTime objects
    the latter has two digit years which can be disambiguation with `baseyear'.
    '''
    
    if not s.strip():
        return None
    year = int(s[0:3])
    if baseyear is not None:
        year += (int(baseyear)/100)*100
    elif year < 80:
        year += 2000
    else:
        year += 1900
    month  = toint(s[3:6])
    day    = toint(s[6:9])
    hour   = toint(s[9:12])
    minute = toint(s[12:15])
    second = tofloat(s[15:26])
    usec   = (second - int(second)) * 1000000
    return GPSDateTime(year, month, day, hour, minute, int(second), int(usec))

class wavelength(object):
    '''
    Parse RINEX WAVELENGTH FACT L1/2 headers

    These headers specify 1: Full cycle ambiguities (default),
    2: half cycle ambiguities (squaring), or 0: does not apply,
    either globally or for particular satellites.
    This is only valid for GPS satellites on frequencies L1 or L2.
    '''
    # If prn list is empty (numsats = 0), L1/2 ambiguity applies to all
    # satellites.  Otherwise, it applies to satellites given in the prnlist;
    # continuation lines are allowed.
    # Ambiguity information is valid until the next 'global' WAVELENGTH FACT,
    # or until that prn is reset.
    def __init__(self):
        '''Set all satellites to default wavelength ambiguity, 1.'''
        self.waveinfo = dict([('G%02d' % prn, (1, 1)) for prn in range(1, 33)])

    def __call__(self, s):
        '''Update wavelength ambiguities with information from a new header.'''
        l1amb   = toint(s[0:6])
        l2amb   = toint(s[6:12])
        numsats = toint(s[12:18])
        if not numsats:  # This is a `global' line
            self.waveinfo = dict([('G%02d' % prn, (l1amb, l2amb))
                for prn in range(1, 33)])
        else:
            for p in range(numsats):
                prn = btog(s[21 + 6 * p]) + '%02d' % \
                                              toint(s[22 + 6 * p : 24 + 6 * p])
                self.waveinfo[prn] = (l1amb, l2amb)
        return self.waveinfo.copy()


class obscode(object):
    '''
    Parse RINEX # / TYPES OF OBSERV headers, specifying observation types.

    These header list observation codes which will be listed in this file.
    Continuation lines are necessary for more than 9 observation types.
    It is possible to redefine this list in the course of a file.
    '''
    # There must be `numtypes' many observation codes, possibly over two lines.
    # Continuation lines have blank `numtypes'.
    def __init__(self):
        self.numtypes = None

    def __call__(self, s):
        nt = toint(s[0:6])
        if self.numtypes is not None and not nt:  # continuation line
            if len(self.obstypes) >= self.numtypes:
                raise RuntimeError('Observation code headers seem broken.')
            for ot in range(min(self.numtypes - len(self.obstypes), 9)):
                self.obstypes += [s[6 * ot + 10 : 6 * ot + 12]]
        elif nt:
            self.numtypes = nt
            self.obstypes = []
            for ot in range(min(nt, 9)):
                self.obstypes += [s[6 * ot + 10 : 6 * ot + 12]]
        else:
            raise RuntimeError('Observation type code continuation header '
                               'without beginning!')
        return self.obstypes[:]


class satnumobs(object):
    '''
    Parse RINEX PRN / # OF OBS headers.

    These headers list how many of each observation type were recorded for
    each satellite included in the file.  If present, there will be one for
    each satellite in the file (as reported in the # OF SATELLITES header.)
    If there are more than 9 observation types, a continuation line will be
    necessary for each satellite.
    This program will determine this information anyway, and check against the
    header if it is supplied.
    '''
    def __init__(self):
        self.sno = {}
        self.prn = None

    def __call__(self, s):
        '''Return a dictionary, by satellite PRN code, of observation counts.

        The counts are a list in the same order as obscode().
        '''
        prn = s[0:3]
        if prn.strip() == '' and self.prn is not None:  # continuation line
            pass
        elif prn.strip() != '':
            prn = btog(prn[0]) + '%02d' % toint(prn[1:])
            self.prn = prn
            if prn in self.sno:
                warn('Repeated # OF OBS for PRN ' + prn + ', why?')
            else:
                self.sno[prn] = []
        else:
            raise RuntimeError('PRN / # OF OBS continuation without beginning!')
        for no in range(9):
            obs = s[no * 6 + 3: no * 6 + 9]
            if obs.strip() == '':
                break
            else:
                self.sno[self.prn] += [toint(obs)]
        return self.sno

class SatBias(dict):
    '''
    Read (download if necesary) satellite biases from UNIBE and create a
    dictionary by prn, bias values.
    '''

    ftp_add = 'ftp.unibe.ch'

    def __init__(self, date, obscodes=['C1'], path=None, factor=-TECUns):

        self.date = GPSDateTime(date)
        if not path:
            path = os.path.join(localpath, 'biases')
        if not os.path.exists(path):
            os.makedirs(path)
        self.p1c1 = os.path.join(path, self.date.strftime('p1c1%y%m.dcb'))
        self.p1p2 = os.path.join(path, self.date.strftime('p1p2%y%m.dcb'))

        try:
            f_p1p2 = open(self.p1p2)
            f_p1c1 = open(self.p1c1)
        except IOError:
            print 'Attempt to download DCB files...'
            if self.download():
                f_p1p2 = open(self.p1p2)
                f_p1c1 = open(self.p1c1)
            else:
                #raise RuntimeError('Satellite bias files not available.')
                print 'Satellite bias files not available.'
                self.default()
                return

        for line in f_p1p2:
            if line[0] == 'G':
                self[line[:3]] = float([x for x in line.split(' ') \
                                        if len(x)<>0][1])*factor

        if 'C1' in obscodes:
            for line in f_p1c1:
                if line[0] == 'G':
                    self[line[:3]] -= float([x for x in line.split(' ') \
                                             if len(x)<>0][1])*factor

    def download(self):
        '''
        Try to donwload bias files from unibe
        '''
        ftp_path = '/aiub/CODE/'
        actual = GPSDateTime().year == self.date.year and \
                 GPSDateTime().month == self.date.month
        if not actual:
            ftp_path = os.path.join(ftp_path, `self.date.year`)
            label = '%02d%02d' % (self.date.year%100, self.date.month)
            ext = '.Z'
        else:
            label = ''
            ext = ''

        try:
            ftp = FTP(self.ftp_add)
            ftp.login()
            ftp.cwd(ftp_path)
            ftp.retrbinary('RETR P1C1%s.DCB%s' % (label, ext),
                               open('%s%s' % (self.p1c1, ext), 'wb').write)
            ftp.retrbinary('RETR P1P2%s.DCB%s' % (label, ext),
                               open('%s%s' % (self.p1p2, ext), 'wb').write)
            ftp.close()
        except:            
            print 'Error downloading files'
            return False
        if not actual:
            os.system('gunzip %s%s' % (self.p1p2, ext))
            os.system('gunzip %s%s' % (self.p1c1, ext))
        return True

    def default(self):
        '''
        '''
        for i in range(33):
            self['G%02d' % i] = 0

class RNXData(Data):
    '''
    A RNXData object is primarily an iterable of records, one for each epoch.
    in chronological order; each record is a dictionary by satellite id
    of dictionaries by observation code of values.

    RNXData.header['name'] also gives access to RINEX header values.
    '''
    TYPE = 'RNX'
    
    def __init__(self, filename, check=True, pcode=None, headeronly=False,
                 station=False, verbose=False):
        self.header = {}
        self.satsystem = None
        self.date = None
        self.data = {}
        self.rec_bias = 0
        self.tzoffset = 0
        self.prns = set()
        self.inmotion = False
        self.arcs = {}
        self.station = {}
        self.bad_prn = []
        self.vars = set()
        self.verbose = verbose
        self.epoch = GPSDateTime(2100,1,1)
    
        if verbose: 
            print 'Parsing file: %s' % filename

        fid = open_rinex(filename, verbose)

        procheader(fid, RINEX, self.header, 0)
        baseyear = self.timesetup(filename)
        if headeronly:
            fid.close()
            return self
        if 'is_crx' in self.header:
            record = records_crx(baseyear)
        else:
            record = records_rnx(baseyear)
    
        while True:
            try:
                record.update(fid)
            except StopIteration:
                break
    
            if record.flag == 5:
                procheader(fid, RINEX, self.header, len(self),
                           xrange(record.numrec), record.epoch)
            elif record.flag == 4:
                procheader(fid, RINEX, self.header, len(self),
                           xrange(record.numrec))
            elif 2 <= record.flag <= 3:
                self.inmotion = record.flag == 2
                procheader(fid, RINEX, self.header, len(self),
                           xrange(record.numrec))
            elif 0 <= record.flag <= 1 or record.flag==6:
                self.data[record.epoch] = RNXRecord(record, record.offset(fid))
                for prn in record.prnlist(fid):
                    dataline = record.get_obs(fid, prn, self.obscodes())
                    self.update(record.epoch, prn, dataline)
                self.rec = self.data[record.epoch]               
                self.checkbreak(record.epoch-self.epoch, record.flag==6)
                self.epoch = record.epoch
        fid.close()
        del self.epoch
        del self.rec
        self.dates = self.data.keys()
        self.dates.sort()
        self.intervals = set([self.dates[i]-self.dates[i-1] for i in range(1, len(self.dates))])
        self.interval = min(self.intervals)
        if check:
            self.endphase(self.prns, last=True)
            self.check()
        if pcode:
            for x in self.header['obscodes']:
                i = self.header['obscodes'][x].index('C1')
                self.header['obscodes'][x][i] = 'P1'
        self.vars = set(self.obscodes())
        self.vars.add('epoch')
        if station:
            self.station = station
        else:
            self.station['code'] = self.header['marker'][0]
            self.station['location'] = xyz2lla(*self.header['markerpos'][0])                
    
    def info(self):
        '''
        #TODO
        '''
        types = {'O': 'Observation', 'N':'Navigation'}
        s  = 'Rinex %s Data\n' % types[self.header['filetype']]
        s += ' marker    : %s\n' % self.header['marker'][0]
        s += ' obscodes  : %s\n' % self.obscodes()
        s += ' # records : %d\n\n' % len(self)


        prns = self.arcs.keys()
        prns.sort()
        for prn in prns:
            s += ' PRN %s:\n' % prn
            for arc in self.arcs[prn]:
                s += '  %s - %s\n' % (self[arc[0]].epoch, self[arc[1]].epoch)

        return s
    
    def update(self, dt, prn, obs):
        '''
        Add new observables info of a prn to last record.
        '''
        self.data[dt][prn] = obs  

    def obscodes(self, index=-1):
        '''
        Return (current) list of observation codes stored in this RNXData.
        '''
        if 'obscodes' not in self.header:
            raise RuntimeError('RINEX file did not define data records')

        obs = self.header['obscodes']

        if index == 0:
            index = min(obs)
        elif index == -1:
            index = max(obs)
        else:
            index = max([k for k in obs if k <= index])
        return obs[index]    

    def checkbreak(self, delta, slip=False):
        '''
        Check times slip and bad records slips.

        Checks at last record added.  This should be called for each record
        inserted, after all its values have been added.

        Also update attributes dates, prns, intervals
        '''
        
        self.prns.update(self.rec) 
        
        #check for breaks
        for prn in set(self.rec).union(self.arcs):
            if slip:
                self.breakphase(prn)
                continue
            
            bad = self.rec.badness(prn)
            if bad:
                self.endphase(prn)
                continue
            
            if delta>480:
                self.breakphase(prn)
                continue

            if not bad and prn not in self.arcs:
                self.arcs[prn] = [[len(self.data)-1, None]]
                continue
            if not bad and self.arcs[prn][-1][1] is not None:
                self.arcs[prn] += [[len(self.data)-1, None]]
                continue            

    def timesetup(self, filename=None):
        '''
        Add significant attributes to object, also set up year disambiguation.
        '''
        if self.satsystem is None and 'satsystem' in self.header:
            self.satsystem = self.header['satsystem']
        baseyear = None
        if 'firsttime' in self.header:
            baseyear = self.header['firsttime'].year
        if 'endtime' in self.header:
            if baseyear is None:
                baseyear = self.header['endtime'].year
        '''
        if filename:
            filename = filename.split('/')[-1]
            if re.search('[0-9]{6}', filename):
                fdate = time.strptime(re.search('[0-9]{6}', filename).group(),
                    '%y%m%d')
                self.date = GPSDateTime(*fdate[:3])
        '''
        return baseyear

    def iterlist(self, prn, obscode=None, filter=None):
        '''
        Returns an iterator over the list of records given, prn and obscodes.
        '''
        def chooser(obs, rec, prn):
            if obs == 'epoch':
                return rec['epoch'].hours
            elif obs in rec[prn]:
                return rec[prn][obs]
            else:
                return None

        rec_idx = []
        for arc in self.arcs[prn]:
            if filter:
                rec_idx += [x for x in range(arc[0], arc[1]) \
                            if self[x][prn][filter[0]]>filter[1]]
            else:
                rec_idx += [x for x in range(arc[0], arc[1])]

        if isinstance(obscode, (list, tuple, set, dict)):
            if not obscode:
                obscode = None
            elif isinstance(obscode, (tuple, set, dict)):
                obscode = [o for o in obscode]
        elif isinstance(obscode, str):
            obscode = [obscode]
        if obscode is None:
            obscode = list(self.allobs)

        for x in rec_idx:
            yield [chooser(obs, self[x], prn) for obs in obscode]

    def check(self):
        '''
        Check header values of file with found and correct them.
        Also save useful attributes.
        '''        
        #checking headers values
        if 'interval' in self.header:
            if self.header['interval'] != self.interval:
                warn('INTERVAL ' + str(self.header['interval']) + \
                        ' does not match minimum observation interval ' +\
                        str(min(self.intervals)))
                self.header['interval'] = self.interval
        else:
            self.header['interval'] = self.interval

        if self.header['interval']<=0:
            warn('Zero or negative interval found...check')
            long = 1800/30
        else:
            long = 1800/self.interval
        for prn in self.arcs:
            self.arcs[prn] = [arc for arc in self.arcs[prn] \
                                   if arc[1]-arc[0]>long]

        #save date attribute
        self.date = self.dates[0].date()
        #save tzoffset attribute
        if 'longitude' in self.station:
            self.tzoffset = self.station['longitude']/15.
        else:
            lat, lon, height = xyz2lla(*self.header['markerpos'][0])
            self.tzoffset = lon/15.

    def merge(self, other, check=True):
        '''
        Merge current RNXData object with another.
        '''
        old_size = len(self)
        self.rec_bias  = 0
        self.epoch = GPSDateTime(2100,1,1)
        self.arcs = {}
        tmp = self.data
        self.data = {}
        
        for obs in other.header['obscodes'][0]:
            if obs not in self.header['obscodes'][0]:
                self.header['obscodes'][0].append(obs)

        for rec in other:
            if rec.epoch in tmp:
                tmp[rec.epoch].update(rec)
            else:
                tmp[rec.epoch] = rec
            self.prns.update(rec)

        self.dates = tmp.keys()
        self.dates.sort()
        if check:
            self.intervals = set([self.dates[i]-self.dates[i-1] for i in range(1, len(self.dates))])
            self.interval = min(self.intervals)
            for dt in self.dates:
                self.data[dt] = tmp[dt]
                self.rec = tmp[dt]
                self.checkbreak(dt-self.epoch, self.rec.flag==6)
                self.epoch = dt
            self.endphase(self.prns, True)
            self.check()
            del self.epoch
            del self.rec
        else:
            self.data = tmp
        del tmp
        print 'Merging old:%d + new:%d -> %d records' % (old_size, len(other), len(self))

    def getgaps(self):
        '''
        Return a list of gaps presents
        '''
        interval = min(self.intervals)
        gaps  = []
        for i in range(1,len(self)):
            delta = self[i].epoch-self[i-1].epoch
            if delta>interval:
                gaps.append([self[i-1].epoch, self[i].epoch, delta])
        return gaps

    def calcpos(self, obscode=None, local=False, samples=20):
        '''
        Calculate receiver position using least squares method,
        -- NO CORRECTIONS ARE APPLIED TO SOLUTION --
        '''

        if obscode:
            obscode = obscode.upper()
        else:
            if 'C1' in self.header['obscodes'][0]:
                obscode = 'C1'
            elif 'P2' in self.header['obscodes'][0]:
                obscode = 'P2'
            else:
                obscode = 'P1'
        if obscode not in self.header['obscodes'][0]:
            raise RuntimeError('Obscode "%s" missing: receiver position can not be calculate' % obscode)
        
        ephe = Yuma(self.date, local=local)
        
        ret  = []
        step = len(self)/samples

        if 'location' in self.station:
            xyz = lla2xyz(*self.station['location'])
        else:
            xyz = self.header['markerpos'][0]

        for i in range(0, len(self), step):
            rec = self[i]
            prns = [prn for prn in rec.keys() if prn in ephe]
            if len(prns)<4: continue
            XYZ_rec = np.tile(lla2xyz(*xyz), (len(prns), 1))
            for j in xrange(3):
                #calculate satellite position
                XYZ_sat  = np.array(ephe.sat_pos(prns, rec.epoch, XYZ_rec[0]))
                R_sat    = ((XYZ_sat[:,0]-XYZ_rec[:,0])**2+(XYZ_sat[:,1]-XYZ_rec[:,1])**2+(XYZ_sat[:,2]-XYZ_rec[:,2])**2)**0.5
                R_sat_i  = np.tile(R_sat, (3, 1)).T
                #form matrix system
                A = np.column_stack([(XYZ_sat-XYZ_rec)/R_sat_i, np.ones(len(prns))])
                Y = R_sat-np.array(rec.iterlist(prns, obscode))-np.array(ephe.get_sat_clock_correction(prns, rec.epoch))
                #solve linear system
                DXYZ = np.linalg.solve(np.dot(A.T, A), np.dot(A.T, Y))
                #recalculate reciver position
                XYZ_rec += np.tile(DXYZ[:3], (len(prns), 1))
            ret.append(XYZ_rec[0])
            #break
        ret = np.array(ret)
        #print ret[:,0], ret[:,1], ret[:,2]
        return ret[:,0].mean(), ret[:,1].mean(), ret[:,2].mean()

    def calctec(self, level=5, rec_bias=None, bad_prn=[], path=localpath, iterations=25):
        '''
        Calculate TEC, STEC and eqTEC and append as observation to records.
        the calculation is only for those records in arcs.
        level:
            0 = only calculate atec rtec ele azi lat lon
            1 = correct bad points
            2 = fix jumps
            3 = correct cycle slips #TODO
            4 = calc RX bias and eqTEC
            5 = leveling zero TEC
        '''
        if len(self)<50 or not self.arcs:
            raise RuntimeError('Not enough data to calculate TEC')
        self.bad_prn = bad_prn
        elem = gpselaz(self.date, os.path.join(path, 'almanac'))
        if 'location' in self.station:
            elem.update(*self.station['location'])
        else:
            elem.update(*xyz2lla(*self.header['markerpos'][0]))
        nprns = []
        for prn in self.arcs:
            if prn not in elem or not self.arcs[prn]:
                nprns.append(prn)
        for prn in nprns:
            del self.arcs[prn]

        for prn in self.arcs:
            for x in self.get_index(prn):
                [xazi, xele, xlat, xlon] = elem.find(prn, self[x].epoch.hours)
                self[x][prn]['ele'] = xele
                self[x][prn]['azi'] = xazi
                self[x][prn]['lat'] = xlat
                self[x][prn]['lon'] = xlon
                self[x][prn]['ltime'] = self[x].epoch.hours + xlon/15.
                self[x][prn]['rTEC'] = self[x].ptec(prn)
                self[x][prn]['aTEC'] = self[x].ctec(prn)
        
        self.vars = self.vars.union(['ele', 'azi', 'lat', 'lon', 'ltime', 'rTEC', 'aTEC'])
        
        if level>=1:
            if self.verbose: print 'Correcting bad points'
            for prn, arclist in self.arcs.items():
                for arc in arclist:
                    [self.correct_points(prn, x) for x in range(*arc)[1:-1]]
        if level>=2:
            if self.verbose: print 'Fixing jumps'
            for prn, arclist in self.arcs.items():
                for arc in arclist:
                    rtec = [self[x][prn]['rTEC'] for x in range(arc[0],arc[1]+1)]
                    times = [self[x].epoch.doys for x in range(arc[0],arc[1]+1)]
                    rtec_index = [x for x in range(arc[0],arc[1]+1)]
                    new_rtec = fix_jumps(rtec, times, 1.173, id=prn)
                    for x, val in zip(rtec_index, new_rtec):
                        self[x][prn]['rTEC'] = val
        if level>=3:
            if self.verbose: print 'Processing data...'
            sat_bias = SatBias(self.date, self.obscodes(0), os.path.join(path, 'biases'))
            for prn in self.arcs:
                for arc in self.arcs[prn]:
                    suma = [(self[x][prn]['aTEC'] - \
                             self[x][prn]['rTEC'])*sin(self[x][prn]['ele']*deg2rad) \
                                for x in range(arc[0],arc[1]+1)]
                    xcount = [sin(self[x][prn]['ele']*deg2rad) \
                              for x in range(arc[0],arc[1]+1)]
                    avediff = sum(suma)/sum(xcount)
                    
                    for x in range(arc[0],arc[1]+1):
                        self[x][prn]['TEC'] = self[x][prn]['rTEC'] + avediff
            self.calctec2(self.arcs.keys(), sat_bias, 0)
            for prn, arclist in self.arcs.items():
                for arc in arclist:
                    self.correct_slips(prn, arc)
            self.vars = self.vars.union(['TEC', 'eqTEC', 'sTEC'])
        if level>=4:
            bias_file = '%s.bias' % self.header['marker'][0]
            bias_file = os.path.join(path, 'biases', bias_file)
            bias_list = {}
            old_rec_bias = None
            write_bias = True
            leveling = True            
            if os.path.exists(bias_file):                
                for line in open(bias_file):
                    bias_list[GPSDateTime(line.split(',')[0], fmt='%Y %m %d')] = float(line.split(',')[1].strip())            
            if rec_bias:
                if self.verbose: print 'Calculating eqTEC using RX bias =', rec_bias
                self.rec_bias = rec_bias
                self.calctec2(self.arcs.keys(), sat_bias, self.rec_bias)
                if level>=5:
                    if self.verbose: print 'Leveling bias for negative TEC'
                    self.zerolevel(sat_bias, iterations=iterations)
                    if abs(rec_bias-self.rec_bias)>3:
                        if self.verbose: 
                            print 'zero leveling bias difference to high using=%s' % rec_bias
                        self.rec_bias = rec_bias
                        self.calctec2(self.arcs.keys(), sat_bias, self.rec_bias)
            else:
                if self.verbose: print 'Calculating RX bias'
                self.rec_bias = self.calc_rec_bias(sat_bias)
                
                if bias_list:
                    dt_list = [dt for dt in bias_list]
                    dt_list.sort()
                    for dt_bias in dt_list:
                        if dt_bias<self.date:
                            old_rec_bias = bias_list[dt_bias]
                            old_dt = dt_bias
                        else:
                            break
                if self.rec_bias and old_rec_bias:
                    if abs(old_rec_bias-self.rec_bias)>10:
                        warn('Using RX bias=%.1f, date: %s' % (old_rec_bias,
                                                               dt))
                        self.rec_bias = old_rec_bias
                        #write_bias = False
                elif not self.rec_bias and old_rec_bias:
                    warn('Using RX bias=%.1f, date: %s' % (old_rec_bias, dt))
                    self.rec_bias = old_rec_bias
                    write_bias = False
                    leveling = False
                elif not self.rec_bias and not old_rec_bias:
                    raise RuntimeError('RX bias can not be estimated')
                self.calctec2(self.arcs.keys(), sat_bias, self.rec_bias)
                if level>=5 and leveling:
                    self.zerolevel(sat_bias, iterations=iterations)
                    #if old_rec_bias and abs(old_rec_bias-self.rec_bias)>8:
                    #    self.rec_bias = old_rec_bias

                if write_bias:
                    bias_list[self.date] = self.rec_bias
                    dt_list = [dt for dt in bias_list]
                    dt_list.sort()
                    fbias = open(bias_file, 'w')
                    [fbias.write('%d %02d %02d, %.1f\n' % (dt.year,
                        dt.month, dt.day, bias_list[dt])) for dt in dt_list]
                    fbias.close()

    def correct_points(self, prn, x):
        '''
        Correct points too distant
        '''

        if abs(self[x][prn]['rTEC']-self[x-1][prn]['rTEC'])>1.173:
            if abs(self[x+1][prn]['rTEC']-self[x-1][prn]['rTEC'])<1.173:
                self[x][prn]['rTEC'] = 0.5*(self[x+1][prn]['rTEC'] + \
                                        self[x-1][prn]['rTEC'])

    def badrecords(self, prn, ttype, tmax, tdiff):
        '''
        Return a index list with bad values for tec calculation
        '''
        bad = set()
        for arc in self.arcs[prn]:
            rec_idx = range(arc[0],arc[1]+1)
            bad.union(set([x for x in rec_idx if self[x][prn][ttype]>tmax]))
            tec_tmp = [self[x][prn][ttype] for x in rec_idx if x not in bad]

            if len(tec_tmp)<>0:
                ave = Stats(tec_tmp).average
            else:
                continue

            bad.union(set([x for x in rec_idx if abs(self[x][prn][ttype]-ave)>tdiff]))
            tmp = [self[x][prn][ttype] for x in rec_idx \
                   if x not in bad]

            if len(tmp)<>0:
                stats = Stats(tmp)
                ave   = stats.average
                stdev = stats.std_dev
                s1    = ave+stdev*4.
                s2    = ave-stdev*4.
                bad.union(set([x for x in rec_idx if (s1<self[x][prn][ttype]) \
                           or (self[x][prn][ttype]<s2)]))
        return bad

    def correct_slips(self, prn, arc):
        '''
        Detect and correct cycle slips taken 10 continues samples and ...
        '''
        return

        nsamp = 6
        rec_idx = range(arc[0]+1, arc[1]+1)
        for x in rec_idx[1:]:
            s1 = x-nsamp/2
            if s1<arc[0]+1:
                s1=arc[0]+1
            s2 = nsamp+s1
            if s2>arc[1]:
                s2 = arc[1]
                s1 = s2-nsamp
            pdiff = [self[x][prn]['rTEC']-self[x-1][prn]['rTEC'] \
                     for x in range(s1,x)+range(x,s2)]
            cdiff = [self[x][prn]['aTEC'] for x in range(s1,x)+range(x,s2)]
            pstat = Stats(pdiff)
            cstat = Stats(cdiff)
            tdiff = self[x][prn]['rTEC']-self[x-1][prn]['rTEC']
            if abs(tdiff)>1.173:
                #print prn, self[x][prn]['rTEC'], tdiff, pdiff, cdiff
                plim_dw = pstat.average-4.*pstat.std_dev
                plim_up = pstat.average+4.*pstat.std_dev
                if plim_up<tdiff or tdiff<plim_dw:
                    #potential cycle slip occurr
                    gdiff = abs(self[x][prn]['TEC'] - self[x][prn]['aTEC'])
                    clim_up = cstat.average + cstat.std_dev
                    clim_dw = cstat.average - cstat.std_dev
                    if clim_dw<self[x][prn]['TEC']<clim_up and abs(tdiff)>1.9:
                        continue
                    if self[x][prn]['ele']<30:
                        #fix as a jump
                        rtec = [self[k][prn]['rTEC'] for k in range(arc[0],arc[1]+1)]
                        secs = [self[k].epoch.seconds for k in range(arc[0],arc[1]+1)]
                        rtec = fix_jumps(rtec, secs, 1.173, x, id=prn)
                        for k in range(x, len(rtec)):
                            self[x+arc[0]][prn]['rTEC'] = rtec[x]
                        continue
                    islips = abs(tdiff/1.173)
                    xadd = 1.173*abs(tdiff)/tdiff*float(islips)
                    if self.verbose:
                        print 'Correcting cycleslip prn=%s, idx=%d, tecdif=%.3f' % \
                                (prn, x, tdiff)
                    for y in range(x, arc[1]):
                        self[y][prn]['rTEC'] = self[y][prn]['rTEC'] - xadd

    def calc_rec_bias(self, sat_bias):
        '''
        Receiver Bias estimation, using mnbrak and brent methods,
        minimizing sum[var('eqTEC')] betwen 03:00 and 06:00 LT.
        '''

        t1 = 3 - self.tzoffset
        t2 = 6 - self.tzoffset

        idx_list = [x for x in range(len(self)) \
                    if t1<self[x].epoch.hours<t2]

        if len(idx_list)==0:
            return False
        else:
            prn_list1 = []
            for x in idx_list:
                prn_list1 += self[x].keys()
            prn_list1 = set(prn_list1)

            prn_list = []
            for prn in prn_list1:
                bad = [self[x].badness(prn) for x in idx_list]
                if True in bad or prn in self.bad_prn or prn not in self.arcs:
                    continue
                prn_list += [prn]
            #find a range where there is a minimun
            val, fval = mnbrak(100, 50, self.fvarsum,
                                   (idx_list, prn_list, sat_bias))
            #find a minimum of fvarsum function
            bias = brent(val, fval, self.fvarsum,
                             (idx_list, prn_list, sat_bias))
            if self.verbose:
                print 'RX bias = %.2f, sum_var(eqTEC) = %.2f after %d iterations' \
                        % (int(bias[0]), bias[1], bias[2])
        if int(bias[0])==0:
            return 1
        return int(bias[0])

    def correct_sat_bias(self, qprn):
        '''
        '''
        if qprn not in self.arcs:
            print 'No correction can be done for prn=%s' % qprn
            return None
        intersc = []
        arcs = self.arcs[qprn]
        for arc in arcs:
            #check for an arc bigger than 30 minutes
            if self[arc[1]].epoch-self[arc[0]].epoch>30*60:
                xstr = self[arc[0]].epoch.hours + self.tzoffset
                xend = self[arc[1]].epoch.hours + self.tzoffset + 1
                #self.cls2app()
                xdmin = 5.0
                hdmin = 500.0

                qprn_idx = [x for x in range(arc[0],arc[1]+1) if xstr<self[x][qprn]['ltime']<xend \
                            and self[x][qprn]['ele']>20]

                arc_prn = set([prn for x in qprn_idx for prn in self[x].keys() \
                               if prn<>qprn])

                for aprn in arc_prn:
                    distmn = 10000
                    aprn_idx = [x for x in self.get_index(aprn) if \
                                xstr<self[x][aprn]['ltime']<xend and self[x][aprn]['ele']>20]

                    idxs = [x for x in aprn_idx if x in qprn_idx]

                    for x in idxs:
                        dist = abs(self[x][qprn]['lat']-self[x][aprn]['lat'])
                        if dist < distmn:
                            hdist = GTODDBL(self[x][qprn]['lat'], self[x][qprn]['lon'],
                                            self[x][aprn]['lat'], self[x][aprn]['lon'])
                            distmn = dist
                            idx = x
                    #print qprn, aprn, idx, distmn, hdist
                    if distmn<xdmin and hdist<hdmin:
                        intersc += [(aprn, idx)]
                #####
        #print intersc
        if intersc:
            sum1, sum2, sum3 = 0, 0, 0
            for bub in intersc:
                iprn, x = bub
                sf   = 1.0/cos(asin(0.94092*cos(self[x][qprn]['ele']*deg2rad)))
                sum1 = sum1 + self[x][qprn]['eqTEC']/sf
                sum2 = sum2 + self[x][iprn]['eqTEC']/sf
                sum3 = sum3 + 1.0/sf**2
            bias_crrc = (sum1-sum2)/sum3
            print 'Correcting bias for %s: ' % qprn , bias_crrc
            return bias_crrc
        else:
            print 'No correction can be done for', qprn
            return None

    def get_index(self, prn):
        '''
        Return a list of index given a prn
        '''
        if prn in self.arcs:
            return [x for arc in self.arcs[prn] for x in range(arc[0], arc[1]+1)]
        else:
            return []

    def fvarsum(self, rec_bias, idx_list, prn_list, sat_bias):
        '''
        calculate sum(std_dev(eqTEC)) per sat for a range of epochs 'idx_list'
        '''
        self.calctec2(prn_list, sat_bias, rec_bias)
        sd = 0

        for x in idx_list:
            prns = self[x].keys()
            eqTEC_list = [self[x][prn]['eqTEC'] for prn in prns \
                         if prn in prn_list and 'eqTEC' in self[x][prn]]
            if eqTEC_list:
                sd += Stats(eqTEC_list).std_dev
        return sd

    def calctec2(self, prns, sat_bias, rec_bias):
        '''
        TEC from carrier phase (rTEC) is smooth but ambiguous.
        TEC from code (aTEC), (C1, C2) or encrypted code (P1, P2) is
        absolute but noisy.  Phase data needs to be fitted to pseudorange data.
        This function calculates the TEC for each PRN in each record, where
        possible.
        Also calculate equivalent TEC using satellite bias and receiver bias
        corrections.
        '''
        
        for prn in prns:
            for arc in self.arcs[prn]:
                for x in range(arc[0],arc[1]+1):
                    self[x][prn]['sTEC'] = self[x][prn]['TEC'] - \
                                                sat_bias[prn]-rec_bias
                    sf = 1.0/cos(asin(0.94092*cos(self[x][prn]['ele']*deg2rad)))
                    self[x][prn]['eqTEC'] = self[x][prn]['sTEC']/sf

    def zerolevel(self, sat_bias, iterations=20):
        '''
        Correct receiver bias by leveling zero TEC values.
        '''
        cnt_zero = 0
        cnt_high = 0
        while True:
            prnszero = set()
            for prn in self.arcs.keys():
                if prn in self.bad_prn: continue
                try:
                    min_tec = self.get_min('eqTEC', prns=(prn,))
                    if min_tec<0: prnszero.add(prn)
                except ValueError:
                    continue
            if len(prnszero)==0 and cnt_zero==0:
                if cnt_high>iterations:
                    print 'Max iteration reach'
                    break
                self.rec_bias += 1
                self.calctec2(self.arcs.keys(), sat_bias, self.rec_bias)
                cnt_zero = 0
                cnt_high += 2
            elif len(prnszero)==0:
                break
            elif len(prnszero)>10:
                self.rec_bias -= 5
                self.calctec2(self.arcs.keys(), sat_bias, self.rec_bias)
                cnt_zero = 0
            elif cnt_zero > 6:
                break
            else:
                self.rec_bias -= 1
                cnt_zero += 1
                self.calctec2(self.arcs.keys(), sat_bias, self.rec_bias)

    def vtec(self, max_ele=25):
        '''
        Experimental total vertical tec calculation
        '''
        def calcvTEC(rec, max_ele):
            tec = 0
            count = 0
            for prn in rec.keys():
                if 'eqTEC' in rec[prn] and rec[prn]['ele']>=max_ele and \
                                                prn not in self.bad_prn:
                    tec += rec[prn]['eqTEC']
                    count += 1
            if count:
                return tec/count
            else:
                return 0
        tim = [rec.epoch for rec in self]
        VTEC = [calcvTEC(rec, max_ele) for rec in self]
        #eqTEC = fix_jumps(eqTEC, tim, 3, id='eqTEC')
        return [(tim[x], VTEC[x]) for x in range(len(VTEC)) \
                    if not tim[x].minute%5 and tim[x].second==0]

    def save(self, path='./', marker='code', alt_name=False, hours=None,
             seconds=None, **kwargs):
        '''
        Save data into a new Rinex file use kwargs for header values.
        '''

        def attr(value, atr):
            if hasattr(value, atr): 
                return getattr(value, atr)
            else: 
                return 0
        
        if hours:
            date = self[-1].epoch - timedelta(1./24*hours)
            if date<self[0].epoch:
                date = self[0].epoch
        else:
            date = self.date

        if alt_name:
            rnx_name = '%s.%02do' % (marker, date.year%100)
        else:
            rnx_name = '%s%003d0.%02do' % (marker, date.doy, date.year%100)

        obscodes = self.obscodes(0)
        rnx = MkRinex(marker, min(self.intervals), obscodes)
        file_tmp = os.tmpfile()
        nrec = 0
        for rec in self:
            if rec.epoch<date:
                continue
            if not hours and rec.epoch.day<>date.day:
                break
            if seconds and rec.epoch.second not in seconds:
                continue
            if nrec==0:
                first_rec = rec
            prns = rec.keys()
            prns.sort()
            sat_line = rnx.mk_satline(rec.epoch, prns)
            dat_line = ''            
            for prn in prns:
                obs_values = [rec[prn][code] for code in obscodes]
                lli_values = [attr(rec[prn][code], 'lli') for code in obscodes]
                snr_values = [attr(rec[prn][code], 'str') for code in obscodes]
                dat_line += rnx.mk_recordline(obs_values, lli_values,
                                snr_values)
            file_tmp.write(sat_line+dat_line)
            endtime = rec.epoch
            nrec += 1

        if nrec<>0:
            kwargs['endtime'] = endtime
            kwargs['interval'] = self.interval
            header = rnx.mk_header(first_rec.epoch, **kwargs)
            if not os.path.exists(path):
                os.makedirs(path)
            file_tmp.seek(0,0)
            file_rnx = open(os.path.join(path, rnx_name), 'w')
            file_rnx.write(header)
            file_rnx.write(file_tmp.read())
            file_tmp.close()
            file_rnx.close()
            print 'Rinex records saved %d of %d' % (nrec, len(self))
            return os.path.join(path, rnx_name), nrec
        else:
            file_tmp.close()
            print 'No records found for given date'
            return None, 0

    def save_data(self, vars=tec_vars, path=localpath, prns=None, gz=False,
                  alt_name=False):
        '''
        Save observables data into a new column text file.
        '''
        if not os.path.exists(path):
            os.makedirs(path)

        if 'network' in self.station:
            code = self.station['code']
            fullname = self.station['fullname']
            location = self.station['location']
            receiver = self.station['info']['receiver']
        else:
            code = self.header['marker'][0]
            fullname = '%s,' % self.header['marker'][0]
            location = xyz2lla(*self.header['markerpos'][0])
            receiver = self.header['receivertype']

        if alt_name:
            fname = '%s.dat' % alt_name
        else:
            fname = '%s_%s.dat' % (code, self.date.strftime('%y%m%d'))

        fname = os.path.join(path, fname)

        if gz:
            fname += '.gz'
            fd = gzip.open(fname, 'wb')
            print 'Creating TEC file: %s' % fname
        else:
            fd = open(fname, 'w')
            print 'Creating TEC file: %s' % fname
            
        header = '%9s%s\n' % (' ', fullname)
        header += '%13.4f%13.4f%13.4f\n' % tuple(location)
        header += '%s\n' % receiver
        header += 'RINEX\n'
        header += '%5d%5d%5d' % (self.date.year,
                                       self.date.month, self.date.day)
        fd.write(header)
        if prns:
            prns = [prn for prn in prns if prn in self.arcs]
        else:
            prns = [prn for prn in self.arcs.keys() if prn not in self.bad_prn]
        prns.sort()
        fmt = '\n%5d    0%4d%4d%4d' + ''.join(['%9.3f' for x in vars])
        for prn in prns:
            nrec = 0
            data = ''
            arcs = self.arcs[prn]
            for arc in arcs:
                for x in range(arc[0], arc[1]+1):
                    #values = [self[x][prn][var] for var in vars]
                    values = []
                    for var in vars:
                        if var in self[x][prn]:
                            values.append(self[x][prn][var])
                        elif var == 'time':
                            values.append(self[x].epoch.hours)
                        else:
                            values.append(-1)
                    nrec += 1
                    data += fmt % ((nrec, self[x].epoch.hour,
                           self[x].epoch.minute, self[x].epoch.second)\
                           +tuple(values))
            prn_line = '\nPRN %7d%7d' % (int(prn[1:]), nrec)
            fd.write(prn_line+data)
        fd.close()
        return fname

    def plot(self, X, Y, **kwargs):
        '''
        Plot vars
        '''
        
        return plot_data_vars(self, X, Y, **kwargs)

    def plot_epochs(self, figname=None, localtime=False):
        '''
        PLOT data availavility
        '''

        plot_rnx(self, figname, localtime)

class RNXRecord(dict):
    '''
    A record of observations (many satellites, many channels) at a given epoch.

    Has atributtes epoch, flag and clock offset in addition to a dictionary
    (by PRN code) of dictionaries (by RINEX observation code, e.g. C1, L2)
    of values.
    Can access as record.epoch, record[13], record['G17'], or iteration.
    '''

    def __init__(self, rec, clkoffset):
        self.epoch = rec.epoch
        self.flag = rec.flag
        self.clockoffset = clkoffset

    def __getitem__(self, index):
        '''
        Allow you to access GPS satellites, eg record['G13'], as
        simply record[13].  For GLONASS or Galileo, you must use the full code.
        '''
        if index == 'epoch':
            return self.epoch
        elif isinstance(index, (int, long, float)):
            return dict.__getitem__(self, 'G%02d' % index)
        else:
            return dict.__getitem__(self, index)

    def __contains__(self, index):
        '''
        Allow containment tests (eg if 13 in record:) for abbreviated GPS PRNs.
        '''
        if isinstance(index, (int, long, float)):
            return dict.__contains__(self, 'G%02d' % index)
        return dict.__contains__(self, index)

    def __eq__(self, other):
        '''
        Equal comparisson reference to the epoch of the record
        '''
        
        return self.epoch==other.epoch

    def __ne__(self, other):
        '''
        Not equal comparisson reference to the epoch of the record
        '''

        return self.epoch!=other.epoch

    def __lt__(self, other):
        '''
        Comparisson reference to the epoch of the record
        '''

        return self.epoch<other.epoch

    def __le__(self, other):
        '''
        Comparisson reference to the epoch of the record
        '''

        return self.epoch<=other.epoch

    def __gt__(self, other):
        '''
        Comparisson reference to the epoch of the record
        '''

        return self.epoch>other.epoch

    def __ge__(self, other):
        '''
        Comparisson reference to the epoch of the record
        '''

        return self.epoch>=other.epoch

    def get_max(self, key):
        '''
        '''
        return max([self[prn][key] for prn in self.keys() if key in self[prn]])
    
    def get_min(self, key):
        '''
        '''
        return min([self[prn][key] for prn in self.keys() if key in self[prn]])

    def iterlist(self, prns, obscode):
        '''
        '''
        return [self[prn][obscode] for prn in prns]

    def ptec(self, prn):
        '''
        Phase (carrier) TEC, if observations are available.

        Convert cycles of L1, L2 to ns (divide by frequency)
        Then TEC_P (TECU) = (L1(ns) - L2(ns)) * TECU/ns
        1 TECU = 10*16 el / m**2
        Suffers from integer ambiguity (`constant' offset which changes after
        each cycle slip.)  Should be dealt with in `arcs' between cycle slips.
        '''
        L1ns = self[prn]['L1']/F1
        L2ns = self[prn]['L2']/F2
        return (L1ns - L2ns) * TECUns

    def ctec(self, prn):
        '''
        Code (pseudorange) TEC, if observations are available.

        p(ns) = p(m) / .3; TEC_C (TECU) = (p2(ns) - p1(ns) - HWCAL) * TECU/ns
        HWCAL is the hardware calibration in ns
        Suffers from satellite bias, receiver bias, multipath and noise.
        '''
        if 'C1' in self[prn] and 'P2' in self[prn]:
            return (self[prn]['P2'] - self[prn]['C1']) * TECUns/(LIGHT_SPEED*1e-9)
        if 'P1' in self[prn] and 'P2' in self[prn]:
            return (self[prn]['P2'] - self[prn]['P1']) * TECUns/(LIGHT_SPEED*1e-9)
        if 'C1' in self[prn] and 'C2' in self[prn]:
            return (self[prn]['C2'] - self[prn]['C1']) * TECUns/(LIGHT_SPEED*1e-9)        

    def badness(self, prn):
        '''
        Indicate how `bad' this record is for TEC calculation
        '''

        if prn not in self:
            return True  # Very bad!
        elif 'L1' not in self[prn] or 'L2' not in self[prn]:
            return True
        elif 'C1' not in self[prn] and 'P1' not in self[prn]:
            return True
        elif 'C2' not in self[prn] and 'P2' not in self[prn]:
            return True
        elif abs(self.ctec(prn))>990 or abs(self.ptec(prn))>99999990:
            return True
        elif self[prn]['L2']==0 or self[prn]['P2']==0:
            return True
        return False    

class PRNrecord(EpochDict):
    '''
    '''
    
    def __getitem__(self, index):
        '''
        '''
        if index in self:
            return dict.__getitem__(self, index)
        else:
            val = value(0)
            val.lli = ' '
            val.str = ' '
            return val

class records_cnt(object):
    '''
    Class to parse rinex records only for count purpose used with count_rinex
    function.
    '''
    def __init__(self, iscrx):
        self.line = ''
        self.iscrx = iscrx

    def update(self, fid):
        self.line = self.getline(fid)
        self.numrec = toint(self.line[29:32])

    def getline(self, fid):
        if not self.iscrx:
            return fid.next()
        else:
            line = fid.next()
            ln = fid.next()
            try: xyz=line[0]
            except IndexError:
                warn('Duplicate record found, epoch=' + \
                     self.epoch.strftime('%y/%m/%d %H:%M:%S'))
                for n in range(self.numrec+2):
                    line = fid.next()
            if line[0] == '&':
                return line.replace('&', ' ')
            else:
                return ''.join(map(choose, self.line, line))

    def prnlist(self, fid):

        prnlist = []
        if not self.iscrx:
            line = self.line
            for z in range(self.numrec):
                s = z % 12
                if z and not s:
                    line = fid.next()
                prn = btog(line[32 + s * 3]) + '%02d' % \
                                            toint(line[33 + s * 3 : 35 + s * 3])
                prnlist += [prn]
            return prnlist
        else:
            for s in range(self.numrec):
                prn = btog(self.line[32 + s * 3]) + '%02d' % \
                                    toint(self.line[33 + s * 3 : 35 + s * 3])
                prnlist += [prn]
            return prnlist

    def dataline(self, fid, numobs):

        line = fid.next()
        if not self.iscrx and numobs>5:
            line = fid.next()
        return

class records_rnx(object):
    '''
    Parse record headers (epoch lines) in standard RINEX:
    Combine continuation lines if necessary.
    '''
    def __init__(self, baseyear):
        self.line = ''
        self.baseyear = baseyear
        self.epoch = None
        #self.intervals = set()

    def update(self, fid):
        '''
        Process a new epoch line.
        '''
        self.line     = self.getline(fid)
        self.oldepoch = self.epoch
        self.epoch    = recordtime(self.line[:26], self.baseyear)
        self.numrec   = toint(self.line[29:32])
        
        while self.oldepoch and self.epoch and self.epoch<=self.oldepoch:
            warn('Duplicate or earlier record found skiping, epoch= %s' %
                 self.epoch)
            n = 0
            if hasattr(self, 'is_crx'):
                n = 1
            for i in range(self.numrec+n):
                lin = fid.next()
            self.line = self.getline(fid)
            self.epoch = recordtime(self.line[:26], self.baseyear)
            self.numrec   = toint(self.line[29:32])

        #if self.epoch is not None and self.oldepoch is not None:
        #    self.intervals.add(self.epoch - self.oldepoch)
        self.numrec = toint(self.line[29:32])
        self.flag = toint(self.line[28])

    def getline(self, fid):
        return fid.next()

    def prnlist(self, fid):
        '''
        Return the list of PRNs (satellite IDs) included in this epoch line.
        May consume extra lines if there are more than 12 PRNs.
        '''
        prnlist = []
        line = self.line
        for z in range(self.numrec):
            s = z % 12
            if z and not s:
                line = fid.next()
            prn = btog(line[32 + s * 3]) + '%02d' % \
                                        toint(line[33 + s * 3 : 35 + s * 3])
            prnlist += [prn]
        return prnlist

    def get_obs(self, fid, prn, obscodes):
        obsline = PRNrecord()
        for x in range(len(obscodes)):
            ind = x%5
            if not ind:
                line = fid.next()
            val = value(tofloat(line[ind*16 : ind*16+14]))
            val.lli = toint(line[ind*16+14 : ind*16+15])
            val.str = toint(line[ind*16+15 : ind*16+16])
            obsline[obscodes[x]] = val
        return obsline

    def offset(self, fid):
        '''
        Return receiver clock offset (optionally) included at end of epoch line.
        '''
        return tofloat(self.line[68:])

class records_crx(records_rnx):
    '''
    Parse record headers in Compact RINEX:
    each line only contains differences from the previous.
    '''
    def __init__(self, baseyear):
        self.is_crx = True
        self.data = {}
        self.offsetval = None
        records_rnx.__init__(self, baseyear)

    def getline(self, fid):
        self.offsetval = None
        line = fid.next()
        if not line:
            return ''.join(map(choose, self.line, line))
        elif line[0] == '&':
            return line.replace('&', ' ')
        else:
            return ''.join(map(choose, self.line, line))

    def prnlist(self, fid):
        prnlist = []
        for s in range(self.numrec):
            prn = btog(self.line[32 + s * 3]) + '%02d' % \
                                       toint(self.line[33 + s * 3 : 35 + s * 3])
            prnlist += [prn]
        return prnlist

    def get_obs(self, fid, prn, obscodes):
        obsline = EpochDict()
        numobs = len(obscodes)

        if prn not in self.data:
            self.data[prn] = {}
            self.data[prn]['arcs'] = [dataArc() for n in range(numobs)]
            self.data[prn]['llis'] = [charArc() for n in range(numobs)]
            self.data[prn]['strs'] = [charArc() for n in range(numobs)]

        line = fid.next()
        vals = line.split(' ', numobs)

        for c, v in enumerate(vals[:numobs]):
            if len(v) >= 2 and v[1] == '&':
                self.data[prn]['arcs'][c] = dataArc(toint(v[0]))
                self.data[prn]['arcs'][c].update(toint(v[2:]))
            elif v.rstrip():
                self.data[prn]['arcs'][c].update(toint(v))
            elif v.rstrip():
                raise ValueError('Uninitialized data arc.')

            obsline[obscodes[c]] = value(self.data[prn]['arcs'][c].get())

        if len(vals) > numobs:
            for c, l in enumerate(vals[numobs][0:numobs*2:2]):
                self.data[prn]['llis'][c].update(l)
                obsline[obscodes[c]].lli = self.data[prn]['llis'][c].get()
            for c, s in enumerate(vals[numobs][1:numobs*2:2]):
                self.data[prn]['strs'][c].update(s)
                obsline[obscodes[c]].str = self.data[prn]['strs'][c].get()
        '''
        for x in range(numobs):
            val = value(self.data[prn]['arcs'][x].get())
            val.lli = self.data[prn]['llis'][x].get()
            val.str = self.data[prn]['strs'][x].get()
            obsline[obscodes[x]] = val
        '''
        return obsline

    def offset(self, fid):
        if self.offsetval is not None:
            return self.offsetval
        line = fid.next()
        if len(line) >= 2 and line[1] == '&':
            self.offsetArc = dataArc(toint(line[0]))
            self.offsetArc.update(toint(line[2:]))
        elif line.rstrip() and 'offsetArc' in self.__dict__:
            self.offsetArc.update(toint(line))
        elif line.rstrip():
            raise ValueError('Uninitialized clock offset data arc.')
        else:
            return 0.
        return self.offsetArc.get()/1000000000

class dataArc(object):
    '''
    Numeric records in Compact RINEX are Nth-order differences
    from previous records.
    Difference order is usually 3. Fields are separated by space.
    LLI and STR are kept separately at the end of the line in one character
    each.
    '''
    def __init__(self, order=3):
        self.order = order
        self.data = []
        self.index = 0

    def update(self, value):
        if self.index < self.order:
            self.data.append(value)
            self.index += 1
        else:
            self.data[self.order - 1] += value

        for diff in range(self.index - 2, -1, -1):
            self.data[diff] += self.data[diff + 1]
        return self.data[0]

    def get(self):
        if len(self.data):
            return self.data[0]/1000.
        else:
            return 0

class charArc(object):
    '''
    LLI and STR records in Compact RINEX only record changes from the previous
    record; space indicated no change.
    '''
    def __init__(self):
        self.data = '0'
    def update(self, char):
        self.data = ''.join(map(choose, self.data, char))
    def get(self):
        return toint(self.data)

class TECData(RNXData):
    '''
    '''
    
    TYPE = 'TEC'
    BREAK_TIME = 600
    
    def __init__(self, filename, satsystem='G', station={}):
        header = {}
        site = {}
        site['info'] = {}
        self.arcs = {}
        self.prns = set()
        self.data = {}
        self.satsystem = None
        self.rec_bias = 0
        self.bad_prn = []    
        self.fo = open_file(filename)
        site['fullname'] = self.fo.next().strip()
        site['location'] = tuple(self.fo.values(fmt=float)[:3])
        site['info']['receiver'] = self.fo.next().strip()
        header['receivertype'] = site['info']['receiver']
        header['marker'] = {0: filename.split('/')[-1][:4]}
        site['lisn_code'] = header['marker'][0]
        site['code'] = header['marker'][0]
        self.fo.next()
        self.date = GPSDateTime(*self.fo.values(fmt=int))
        
        while True:
            try:
                vals = self.fo.values()
            except StopIteration:
                break
            if 'PRN' in vals:
                prn = '%s%02d' % (satsystem, int(vals[1]))
                self.prns.add(prn)
                continue
            self.add(prn, vals)
        
        self.header = header
        if not station:
            self.station = site
        else:
            self.station = station
        self.dates = self.data.keys()
        self.dates.sort()
        self.intervals = set([self.dates[i]-self.dates[i-1] for i in range(1, len(self.dates))])
        self.interval = min(self.intervals)
        [self.checkbreak(i) for i in range(len(self))]
        self.endphase(len(self)-1, self.prns)
        if len(self)>0:
            self.vars = set(self[0].KEYS)
            self.vars.add('epoch')
        
    def add(self, prn, vals, secs=None):
        '''
        '''
        dt = self.date.replace(hour=int(vals[2]), minute=int(vals[3]), second=int(vals[4]))
        
        if dt not in self.data:
            self.data[dt] = TECrecord(dt, prn, vals)
        else:
            self.data[dt].update(TECrecord(dt, prn, vals))
        
    def checkbreak(self, n):
        '''
        '''
        
        for prn in set(self[n]).union(self.arcs):
            if n>0:
                if self[n].epoch-self[n-1].epoch>1800:
                    self.breakphase(n, prn)
                    continue
            
            bad = prn not in self[n] 
                                    
            if not bad and prn not in self.arcs:
                self.arcs[prn] = [[n, None]]
                continue
            if not bad and self.arcs[prn][-1][1] is not None:
                self.arcs[prn] += [[n, None]]
                continue
            if prn not in self.arcs or self.arcs[prn][-1][1] is not \
                    None:
                continue
            if bad:
                self.endphase(n, prn)
                continue
            
    def breakphase(self, n, prn):
        '''
        Begin new phase-connected-arc for satellite prn.
        '''
        if isinstance(prn, (list, tuple, set, dict)):
            [self.breakphase(p) for p in prn]
        elif prn not in self.arcs:
            self.arcs[prn] = [[n, None]]
        elif prn not in self[n]:
            self.endphase(n, prn)
        else:
            self.endphase(n, prn)
            self.arcs[prn] += [[n, None]] 
            
    def endphase(self, n, prn):
        '''
        End current phase-connected-arc, if any, for satellite prn.
        Ends arc just before the current record.
        '''

        if isinstance(prn, (list, tuple, set, dict)):
            [self.endphase(n, p) for p in prn]
        elif prn in self.arcs and self.arcs[prn][-1][1] is None:
            self.arcs[prn][-1][1] = n-1          
                
class TECrecord(RNXRecord):
    '''
    '''
    
    KEYS = ['eqTEC', 'ltime', 'ele', 'lat', 'lon', 'sTEC', 'azi']
    
    def __init__(self, epoch, prn, vals):
        self.epoch = epoch
        self[prn] = EpochDict(zip(self.KEYS, [float(x) for x in vals[5:]]))

class BUBData(list):
    '''
    '''
    
    TYPE = 'BUB'
        
    def __init__(self, filename):
        self.codes = set()
        self.prns = set()
        fo = open_file(filename)
        cnt = 0
        while True:
            try:
                values = fo.values()
            except StopIteration:
                break
            if ':' in values:
                dt = GPSDateTime(int(values[2]), int(values[3]), int(values[4]))
            elif 'bubble' in values:
                continue
            elif cnt!=int(values[0]):
                self.append(BUBRecord(dt, values))
                self.codes.add(self[-1].code)
                self.prns.add(self[-1].prn)
                cnt += 1
            else:
                self[-1].add(values)

    def __str__(self):
        '''
        '''
        return '%sData[records=%s]' % (self.TYPE, len(self))

class BUBRecord(dict):
    '''
    '''
    
    def __init__(self, date, values):
        self.epoch = date + int(values[3])
        self.prn = 'G%02d' % int(values[2])
        self.code = values[1]
        self['seconds'] = [int(values[3])]
        self['lat'] = [float(values[4])]
        self['lon'] = [float(values[5])]
        self['eqTEC'] = [float(values[6])]
        self['dTEC'] = [float(values[7])]
        self['bTEC'] = [float(values[8])]
        
    def add(self, values):
        self['seconds'].append(int(values[3]))
        self['lat'].append(float(values[4]))
        self['lon'].append(float(values[5]))
        self['eqTEC'].append(float(values[6]))
        self['dTEC'].append(float(values[7]))
        self['bTEC'].append(float(values[8]))

class header(object):
    '''
    For each RINEX header type, this holds a list of field objects
    which are defined in the associated line.
    This is for header values which should only occurr once.
    '''

    class field(object):
        '''
        Describes a value in a RINEX header: variable name, position in the
        line, and how to interpret it.
        '''
        def __init__(self, name, start, stop, convert=mystrip):
            self.name = name
            self.start = start
            self.stop = stop
            self.convert = convert

        def read(self, line):
            return self.convert(line[self.start:self.stop])
        '''
        def __deepcopy__(self, memo={}):
            #Fix deepcopying in Python 2.4.
            newfield = self.__class__(deepcopy(self.name, memo),
                              deepcopy(self.start, memo),
                              deepcopy(self.stop, memo), self.convert)
            memo[id(self)] = newfield
            return newfield
        '''

    def __init__(self, field_args):
        self.values = [header.field(*fargs) for fargs in field_args]

    def read(self, line, header, recordnum, epoch=None):
        for field in self.values:
            header[field.name] = field.read(line)

class listheader(header):
    '''
    For multiple header values

    The value may change for different observation records.
    If multiple instances are at the same record number, the last is used.
    They are accessed by record number; whichever value is valid for that
    record is returned.
    '''
    def read(self, line, header, recordnum, epoch=None):
        for field in self.values:
            if field.name not in header:
                header[field.name] = {}
            if field.name == 'comment':
                if recordnum not in header[field.name]:
                    header[field.name][recordnum] = []
                header[field.name][recordnum] += [field.read(line)]
            else:
                header[field.name][recordnum] = field.read(line)
                if epoch is not None:
                    header[field.name][recordnum].epoch = epoch

RINEX = {
    'CRINEX VERS   / TYPE' : header([('crnxver', 0, 3, crxcheck),
                              ('is_crx', 0, 0, truth)]),
    'CRINEX PROG / DATE  ' : header([('crnxprog', 0, 20),
                              ('crxdate', 40, 60),
                              ('is_crx', 0, 0, truth)]),
    'RINEX VERSION / TYPE' : header([('rnxver', 0, 9, versioncheck),
                              ('filetype', 20, 21, iso),
                              ('satsystem', 40, 41, btog)]),
    'PGM / RUN BY/ DATE  ' : header([('rnxprog', 0, 20),
                              ('run_by', 20, 40),
                              ('filedate', 40, 60)]),
    'PGM / RUN BY / DATE ' : header([('rnxprog', 0, 20),
                              ('run_by', 20, 40),
                              ('filedate', 40, 60)]),
    'COMMENT             ' : listheader([('comment', 0, 60)]),
    'MARKER NAME         ' : listheader([('marker', 0, 60)]),
    'MARKER NUMBER       ' : listheader([('markernum', 0, 20)]),
    'APPROX POSITION XYZ ' : listheader([('markerpos', 0, 42, to3float)]),
    'OBSERVER / AGENCY   ' : header([('observer', 0, 20),
                              ('agency', 20, 60)]),
    'REC # / TYPE / VERS ' : header([('receivernum', 0, 20),
                              ('receivertype', 20, 40),
                              ('receiverver', 40, 60)]),
    'ANT # / TYPE        ' : listheader([('antennanum', 0, 20),
                              ('antennatype', 20, 40)]),
    'ANTENNA: DELTA H/E/N' : listheader([('antennashift', 0, 42, to3float)]),
    'WAVELENGTH FACT L1/2' : listheader([('ambiguity', 0, 53, wavelength())]),
    '# / TYPES OF OBSERV ' : listheader([('obscodes', 0, 60, obscode())]),
    'INTERVAL            ' : header([('interval', 0, 10, tofloat)]),
    'TIME OF FIRST OBS   ' : header([('firsttime', 0, 43, headertime),
                              ('firsttimesys', 51, 54)]),
    'TIME OF LAST OBS    ' : header([('endtime', 0, 43, headertime),
                              ('endtimesys', 51, 54)]),
    'RCV CLOCK OFFS APPL ' :
                      listheader([('receiverclockcorrection', 0, 6, toint)]),
    'LEAP SECONDS        ' : listheader([('leapseconds', 0, 6, toint)]),
    '# OF SATELLITES     ' : header([('numsatellites', 0, 6, toint)]),
    'PRN / # OF OBS      ' :
                      header([('obsnumpersatellite', 3, 60, satnumobs())])
}

def procheader(fid, RINEX, header, numrec, numlines=itertools.repeat(0),
                        epoch=None):
    '''
    Parse Rinex headers into a dictionary.
    '''

    for n in numlines:
        try:
            line = '%-80s' % fid.next()  # pad spaces to 80
        except StopIteration:
            break
        label = line[60:]
        if label.strip() == 'END OF HEADER':
            break
        if label in RINEX:
            RINEX[label].read(line, header, numrec, epoch)
        else:
            warn('Header line ' + label + ' unrecognized; ignoring')

def open_rinex(filename, verbose=False):
    '''
    Open a rinex file (rinex, crinex, gziped y/o tarred), return an open_file
    object
    '''

    if tarfile.is_tarfile(filename):
        if verbose: print 'Unpacking tarfile.'
        zfile = tarfile.open(filename)  # Automatically handles tar.gz,bz2
        zfile = zfile.extractfile(zfile.next())
    elif zipfile.is_zipfile(filename):
        zfile = zipfile.ZipFile(filename)
        obs = [z.filename for z in zfile.filelist if re.search('\.[0-9]{2}[OoDd]$', z.filename)]
        if obs:
            zfile = StringIO(zfile.open(obs[0]).read())
        else:
            raise IOError('Invalid rinex file in ZIP')        
    elif filename.lower().endswith('.gz'):
        if verbose: print 'Gunzipping file.'
        zfile = gzip.open(filename)
        zfile.name = filename[:filename.rfind('.gz')]
    elif filename.endswith('.Z'):
        zfile = StringIO(Popen(["zcat", filename],
                         stdout=PIPE).communicate()[0])
    else:
        zfile = open(filename)
    if not re.search('\.[0-9]{2}[OoDd]', filename.split('/')[-1]):
        if verbose: warn(filename + ' Invalid Filename format')
    return open_file(zfile)

def count_rinex(filename):
    '''
    Simplified version of read_rinex that return the number of records of the
    file.
    '''
    fid = open_rinex(filename)
    header = {}
    procheader(fid, RINEX, header, 0)
    record = records_cnt('is_crx' in header)
    count = 0
    while True:
        try:
            record.update(fid)
        except StopIteration:
            break
        count += 1
        for prn in record.prnlist(fid):
            record.dataline(fid, header['obscodes'][0])
    fid.close()
    return count

def read_rinex(filename, check=True, pcode=None, headeronly=False,
               station=False, verbose=False):
    '''
    keep for old compatibility.
    '''
    
    return RNXData(filename, check, pcode, headeronly,
                   station, verbose)

if __name__=='__main__':
    pass