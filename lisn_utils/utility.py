# -*- coding: utf-8 -*-
'''
Miscellaneous classes and functions.
These are not very specific in usage, however could be useful anywhere.

@author: Juan C. Espinoza
@contact: jucar.espinoza@gmail.com

'''
import os, time, math, gzip, json, shutil, glob
import traceback, fnmatch
import numpy as np
import smtplib
from email.message import Message
from StringIO import StringIO
from gpsdatetime import GPSDateTime
from ftplib import FTP

try:
    import requests
except:
    requests = False
try:
    import paramiko
except ImportError:
    paramiko = False
try:
    import MySQLdb as sql
except ImportError:
    sql = False
try:
    from lxml import etree
except ImportError:
    etree = False

PI = math.pi
rad2deg = 180./PI
deg2rad = PI/180.

_path_ = os.path.dirname(os.path.abspath(__file__))

scinda_vars = {
    'west'    : (2, '-', (0, 1.2), 'S4', 'UHF West Antenna'),
    'west_r'  : (3, '-', (0, 1.2), 'S4', 'UHF West Antenna Redundant'),
    'L'       : (5, '-', (0, 1.2), 'S4', 'L-Band Antenna'),
    'east'    : (6, '-', (0, 1.2), 'S4', 'UHF East Antenna'),
    'east_r'  : (7, '-', (0, 1.2), 'S4', 'UHF East Antenna Redundant'),
    'drift'   : (14, '.r', (-100, 300), 'Velocity m/s',
                 'Drift for West Channels')}

class Notification(object):
    '''
    Class to handle and send email notifications
    '''
   
    def __init__(self, alert=''):
        self.alert = alert
        self.params = Config().notification()
           
    def send(self, receiver, subject, content=''):
        '''
        '''

        message = Message()
        message['To']      = receiver
        message['From']    = self.params['from']
        message['Subject'] = 'LISN %s: %s' % (self.alert, subject)
        message.set_payload(content)

        try:
            print 'Sending Mail to: %s' % receiver
            if self.params['server']=='localhost':
                smtp = smtplib.SMTP('localhost')
                smtp.sendmail(self.params['from'], receiver.split(','), message.as_string())
                smtp.quit()
            else:
                smtp = smtplib.SMTP(self.params['server'], self.params['port'])
                smtp.ehlo()
                smtp.starttls()
                smtp.ehlo()
                smtp.login(self.params['user'], self.params['password'])
                smtp.sendmail(self.params['from'], receiver.split(','), message.as_string())
                smtp.quit()                
            print "Successfully sent email"
        except Exception as e:
            print "Error: %s" % e
            return e
        return True

class LoginError(Exception):
    pass

class LISNRequest(object):
    '''
    Class to access django api information, uses requests module
    '''
    
    def __init__(self, debug=False, config=None):
        if requests is False:
            raise ImportError('Install requests module')
        self.debug = debug
        self.params = Config(filename=config)
        self.client = requests.session()
        self.login()
    
    def login(self):
        self.client.get(self.params.url_login)
        csrftoken = self.client.cookies['csrftoken']
        login_data = {'username': self.params.username, 'password': self.params.password, 'next': '/',
                      'this_is_the_login_form': '1', 'csrfmiddlewaretoken': csrftoken}
        req = self.client.post(self.params.url_login, data=login_data)
        #TODO check for successful login
        if '<a href="/accounts/logout">[logout]' in req.text:
            self.ok = True
        else: 
            self.ok = False
    
    def close(self):
        self.client.get(self.params.url_logout)
            
    def stations(self, **kwargs):
        url = os.path.join(self.params.url, 'utils/stations')        
        req = self.client.get(url, params=kwargs)
        if self.debug:
            print req.url
        return req.json()
    
    def station(self, code, instrument):
        url = os.path.join(self.params.url, 'utils/station')
        params = {'code':code, 'instrument':instrument}
        req = self.client.get(url, params=params)
        if self.debug:
            print req.url
        return req.json()
    
    def update_station(self, code, instrument, **kwargs):
        url = os.path.join(self.params.url, 'utils/update', instrument, code)
        req = self.client.get(url, params=kwargs)
        if self.debug:
            print req.url
        return req.json()
    
    def update_date(self, code, instrument, datatype, **kwargs):
        url = os.path.join(self.params.url, 'data/update', instrument, code, datatype)
        req = self.client.get(url, params=kwargs)
        if self.debug:
            print req.url
        return req.json()

class Data(object):
    
    TYPE = ''
    BREAK_TIME = 60    
    
    def __len__(self):
        return len(self.data)

    def __str__(self):
        '''
        '''
        return '%sData[site=%s, records=%s]' % (self.TYPE, self.station['code'], len(self))
    
    def __getitem__(self, key):
        '''
        '''
        if isinstance(key, slice):
            return [self.data[self.dates[x]] for x in xrange(*key.indices(len(self.dates)))]
        elif isinstance(key, int):
            return self.data[self.dates[key]]
        elif isinstance(key, GPSDateTime):
            return self.data[key]
    
    def get_records(self, prn=None, dt1=None, dt2=None):
        '''
        '''
        if prn is None:
            if dt1 is None:
                return [self.data[dt] for dt in self.dates]
            else:
                return [rec for rec in self if dt1<=rec.epoch<=dt2]
        else:
            if dt1 is None:
                return [self.data[dt] for dt in self.dates if prn in self.data[dt]]
            else:
                return [rec for rec in self if prn in rec and dt1<=rec.epoch<=dt2]
    
    def get_max(self, key, prns=None, values=None):
        '''
        '''
        if prns is None:
            prns = self.prns
        if values is None:                          
            return max([rec[prn][key] for rec in self for prn in prns \
                        if prn in rec and key in rec[prn]])
        else:
            return max([rec[prn][key] for rec in self for prn in prns \
                        if prn in rec and key in rec[prn] and rec[prn][values[0]]>values[1]])
    
    def get_min(self, key, values=None, prns=None):
        '''
        '''
        if prns is None:
            prns = self.prns
        if values is None:                          
            return min([rec[prn][key] for rec in self for prn in prns \
                        if prn in rec and key in rec[prn]])
        else:
            return min([rec[prn][key] for rec in self for prn in prns \
                        if prn in rec and key in rec[prn] and rec[prn][values[0]]>values[1]])
    
    def get_array(self, prn, keys=None, dt1=None, dt2=None, values=None, attr='doys'):
        '''
        '''
        
        def check(rec, prn, key, values):
            dum = prn in rec and key in rec[prn]
            if values:
                dum = dum and values[0] in rec[prn] and rec[prn][values[0]]>values[1]
            if dum:
                if key=='epoch':
                    return getattr(rec.epoch, attr)
                else:
                    return rec[prn][key]
            return np.nan
        
        if isinstance(keys, (list, tuple, set, dict)):
            if not keys:
                keys = None
            elif isinstance(keys, (tuple, set, dict)):
                keys = [o for o in keys]
        elif isinstance(keys, str):
            keys = [keys]
        if keys is None:
            keys = ['epoch']+self.KEYS
        
        if dt1 is not None and dt2 is not None:
            array = [[check(rec, prn, key, values) for key in keys] for rec in self if dt1<=rec.epoch<=dt2]
        else:
            array = [[check(rec, prn, key, values) for key in keys] for rec in self]
            
        if array:
            return np.array(array, dtype=float).T
        else:
            return np.array([[] for key in keys])

    def breakphase(self, prn):
        '''
        Begin new phase-connected-arc for satellite prn.
        '''
        if isinstance(prn, (list, tuple, set, dict)):
            [self.breakphase(p) for p in prn]
        elif prn not in self.arcs:
            self.arcs[prn] = [[len(self.data)-1, None]]
        else:
            self.endphase(prn)
            self.arcs[prn] += [[len(self.data)-1, None]]
    
    def endphase(self, prn, last=False):
        '''
        End current phase-connected-arc, if any, for satellite prn.
        Ends arc just before the current record.
        '''
        if isinstance(prn, (list, tuple, set, dict)):
            [self.endphase(p, last) for p in prn]
        elif not last and prn in self.arcs and self.arcs[prn][-1][1] is None:
            self.arcs[prn][-1][1] = len(self.data)-2
        elif last and prn in self.arcs and self.arcs[prn][-1][1] is None:
            self.arcs[prn][-1][1] = len(self.data)-1

class MySQL(object):
    '''
    '''
    def __new__(cls, host, user, passwd, db, verbose=False):
        
        if sql is False:
            raise ImportError('Install MySQL-python module')
        mysql = object.__new__(cls)
        mysql.conn = sql.connect(host   = host,
                                    user   = user,
                                    passwd = passwd,
                                    db     = db)
        mysql.cursor = mysql.conn.cursor()
        mysql.verbose = verbose
        return mysql

    def query(self, query):
        '''
        '''
        result = self.cursor.execute(query)
        if self.verbose:
            print '>>>%s -> %s' % (query,result)
        self.commit()
        return result

    def select(self, table, fields, filter='1', distinct=False):
        '''
        '''
        ret    = []
        fields = ','.join([s for s in fields])
        if distinct:
            query  = 'SELECT DISTINCT %s FROM %s WHERE %s' % (fields, table, filter)
        else:
            query  = 'SELECT %s FROM %s WHERE %s' % (fields, table, filter)
        

        if self.verbose:
            print '>>>%s -> ' % query,

        self.cursor.execute(query)
        result = self.cursor.fetchall()
        if self.verbose:
            print result
        for tup in result:
            if '*' in fields:
                ret.append(list(tup))
            else:
                ret.append(dict(zip(fields.split(','), list(tup))))
        return ret

    def update(self, table, fields, values, filter='1'):
        '''
        '''
        vars = ""
        for field, value in zip(fields, values):
            if type(value) in (float,):
                vars += "%s=%f," % (field, value)
            elif type(value) in (int, long):
                vars += "%s=%d," % (field, value)
            elif value in ('NULL', 'null', 'Null', 'False'):
                vars += "%s=NULL," % field
            else:
                vars += "%s='%s'," % (field, value)
        vars = vars[:-1]
        query = '''UPDATE %s SET %s WHERE %s''' % (table, vars, filter)
        
        if self.verbose:
            print '>>>%s -> ' % query,
        result = self.cursor.execute(query)
        if self.verbose:
            print result
        self.commit()
        return result

    def insert(self, table, fields, values):
        '''
        '''
        
        fields = ','.join(fields)
        vars   = ""
        for var in values:
            if type(var) in (float,):
                vars += "%f," % var
            elif type(var) in (int, long):
                vars += "%d," % var
            elif var in ('NULL', 'null', 'Null', False, 'False'):
                vars += "NULL,"
            else:
                vars += "'%s'," % var
        vars = vars[:-1]
        query  = '''INSERT INTO %s (%s) VALUES (%s)''' % (table, fields, vars)
        result = self.cursor.execute(query)
        if self.verbose:
            print '>>>%s -> %s' % (query,result)
        self.commit()
        return result

    def select_one(self, table, field, fields, values, like=False):
        '''
        '''
        
        filter = ""
        for fild, value in zip(fields, values):
            if type(value) in (float,):
                filter += "%s=%f AND " % (fild, value)
            elif type(value) in (int, long):
                filter += "%s=%d AND " % (fild, value)
            elif value in ('NULL', 'null', 'Null'):
                filter += "%s=NULL AND " % fild
            elif value.startswith('@str:'):
                filter += "%s = '%s' AND " % (fild, value[5:])
            elif value.startswith('@var:'):
                filter += "%s=%s AND " % (fild, value[5:])
            elif like:
                filter += "%s LIKE '%%%s%%' AND " % (fild, value)
            else:
                filter += "%s = '%s' AND " % (fild, value)
        filter = filter[:-4]

        query = '''SELECT %s FROM %s WHERE %s''' % (field, table, filter)

        self.cursor.execute(query)
        result = self.cursor.fetchall()

        if self.verbose:
            print '>>>%s -> %s' % (query,result)

        if len(result)==1:
            return result[0][0]
        elif len(result)>1:
            return [tup[0] for tup in result]
        else:
            return False

    def commit(self):
        self.conn.commit()

    def close(self):
        self.conn.commit()
        self.cursor.close()
        self.conn.close()

class MyXML(object):
    '''
    '''

    def __init__(self, s=''):
        '''
        '''
        if etree is False:
            raise ImportError('Install lxml module')
        if s:
            self.data  = etree.fromstring(s, )            
        else:
            self.data  = etree.Element('data')

    def __str__(self):
        '''
        '''
        
        return self.tostring()

    def update_attr(self, attr, value):
        '''
        '''
        self.data.set(attr, value)

    def delete_attr(self, attr):
        '''
        '''
        #print dir(self.doc)
        if attr in self.data.keys():
            del self.data.attrib[attr]

    def find_text(self, element, tag, text):
        '''
        '''
        for el in element:
            if el.tag.strip()==tag.strip() and el.text.strip()==text.strip():
                return el
        return False

    def add_date(self, year, month, day, dt=None):
        '''
        Append or add new childs for given year, month and day to XML data
        '''    
        if dt:
            year,month,day = dt.year,dt.month,dt.day
        # add year
        el_year = self.find_text(self.data, 'y', '%s' % year)
        if el_year is False :
            el_year = etree.SubElement(self.data, 'y')
            el_year.text = '%s' % year
        # add month
        el_month = self.find_text(el_year, 'm', '%s' % month)
        if el_month is False:
            el_month = etree.SubElement(el_year, 'm')
            el_month.text = '%s' % month
        # add day
        el_day = self.find_text(el_month, 'd', '%s' % day)
        if el_day is False:
            el_day = etree.SubElement(el_month, 'd')
            el_day.text = '%s' % day
            return True
        else:
            return False

    def delete_data(self):
        '''
        '''
        etree.strip_elements(self.data, 'y', 'm', 'd')

    def tostring(self, pretty_print=True, encoding='utf-8'):
        '''
        '''
        s = etree.tostring(self.data, pretty_print=pretty_print,
                           encoding=encoding, xml_declaration=True)
        return s.replace("'",'"')

class MyFTP(object):
    '''
    Class to manage sftp and ftp connections with servers, functions supported:
    listdir, download, getfileutime
    '''

    def __init__(self, host, username, password, secure=False, havekey=False):
        if secure:
            self.secure = secure
            if paramiko is False:
                raise ImportError('Install paramiko module')
        else:
            self.secure = False
        self.transport = 0
        self.ftp = 0
        self.error = False
        try:
            print 'Connecting to ' + host
            if secure:
                self.transport = paramiko.Transport((host, 22))
                if havekey:
                    sha = False
                    try:
                        print 'Trying SHA key...'
                        privatekeyfile = os.path.expanduser('~/.ssh/id_rsa')
                        mykey = paramiko.RSAKey.from_private_key_file(privatekeyfile)
                        self.transport.connect(username=username, pkey=mykey)
                        sha = True
                    except:
                        pass
                    if not sha:
                        try:
                            print 'Trying DSA key...'
                            privatekeyfile = os.path.expanduser('~/.ssh/id_dsa')
                            mykey = paramiko.DSSKey.from_private_key_file(privatekeyfile)
                            self.transport.connect(username=username, pkey=mykey)
                        except:
                            print 'Error connecting with key'
                            self.error = True
                else:
                    self.transport.connect(username=username, password=password)
                self.ftp = paramiko.SFTPClient.from_transport(self.transport)
            else:
                self.ftp = FTP(host)
                self.ftp.login(username, password)
        except:
            traceback.print_exc(2)
            self.error = True

    def listdir(self, path='.'):
        if self.secure:
            return self.ftp.listdir(path)
        else:
            return [x.split('/')[-1] for x in self.ftp.nlst(path)]

    def download(self, remotepath, localpath):
        print 'Downloading: %s' % remotepath
        if self.secure:
            self.ftp.get(remotepath, localpath)
        else:
            self.ftp.retrbinary('RETR '+remotepath, open(localpath, 'wb').write)

    def upload(self, localpath, remotepath):
        print 'Uploading: %s' % localpath
        if self.secure:
            self.ftp.put(localpath, remotepath)
        else:
            self.ftp.storbinary('STOR '+remotepath, open(localpath, 'rb'))

    def remove(self, path):
        print 'Deleting: %s' % path
        if self.secure:
            self.ftp.remove(path)
        else:
            self.ftp.delete(path)
    
    def exists(self, filename):
        path = os.path.dirname(filename)
        name = filename.split('/')[-1]
        return name in self.listdir(path)
    
    def utime(self, path):
        if self.secure:
            return self.ftp.stat(path).st_mtime
        else:
            return int(time.mktime(time.strptime(str(self.ftp.sendcmd('MDTM '+ path)[-14:]),
                                                 '%Y%m%d%H%M%S')))
    
    def open(self, filename, mode='rb'):        
        if self.exists(filename):
            if self.secure:
                return self.ftp.open(filename, mode)
            else:
                so = StringIO()
                so.name = filename
                self.ftp.retrbinary('RETR '+filename, so.write)
                so.seek(0)
                return so
        else:
            print 'File not found in server: %s' % filename
            return None
    
    def close(self):
        if self.ftp <> 0: self.ftp.close()
        if self.transport <> 0: self.transport.close()

class scindaData(list):
    '''
    class to ...
    '''

    def __init__(self, filename, vars='all'):
        '''
        '''
        if vars=='all':
            vars = ['west', 'west_r', 'east', 'east_r', 'L', 'drift']
        elif vars=='west':
            vars = ['west', 'west_r']
        elif vars=='east':
            vars = ['east', 'east_r']
        else:
            vars = [var for var in vars.split(',') if var in scinda_vars]
        self.vars = [scinda_vars[var] for var in vars]

        fo = open_file(filename)
        while True:
            try:
                line = fo.next()
            except StopIteration:
                break
            if line[0] == '#' or len(line)==0: continue
            cols = [x.strip() for x in line.split(' ') if len(x)<>0]
            dt = GPSDateTime(cols[0]+cols[1])
            cols = [float(x.strip()) for x in line.split(' ') if len(x)<>0]
            self.append((dt,)+tuple(cols[2:]))
        fo.close()    

def copy_file(filei, fileo, move=False):
    '''
    '''
    output_path = os.path.dirname(fileo)
    if not os.path.exists(output_path):
        os.makedirs(output_path)
    shutil.copy(filei, fileo)
    if move:
        os.remove(filei)

def rglob(src, pattern, recursive=True):
    '''
    Recursive glob
    '''
    if not recursive:
        return glob.glob(os.path.join(src, pattern))
    ret = []
    for root, dirs, files in os.walk(src):
        ret.extend([os.path.join(root, file) for file in files \
                   if fnmatch.fnmatch(file, pattern)])
    return ret

class AttrDict(dict):
    '''
    '''
    def __init__(self, *args, **kwargs):
        super(AttrDict, self).__init__(*args, **kwargs)
        self.__dict__ = self

    def from_keys(self, iterable):
        new = {}
        for key in iterable:
            new[key] = self[key]
        return AttrDict(new)

class EpochDict(dict):
    '''
    '''
    def __contains__(self, index):
        if index=='epoch':
            return True
        return dict.__contains__(self, index)
        
def path_size(path, name='*', human=False):
    '''
    Return the size of the files that math with path
    '''
    size = 0
    cnt = 0
    for filename in glob.glob(os.path.join(path, name)):
        size += os.path.getsize(filename)
        cnt += 1
    if human:
        return human_size(size)
    else:
        return (size, cnt)

def readable_filesize(size, fmt='%3.2f'):
    '''
    Get size in human readable format

    return (string, factor, unit)
    '''
    f = 1
    fmt += ' %s'
    for unit in ['b','Kb','Mb','Gb', 'Tb']:
        if size<1024.0:
            return (fmt % (size, unit), f, unit)
        size /= 1024.0
        f *= 1024.0

class open_file(object):
    '''
    Wrap "sufficiently file-like objects" (ie those with readline())
    in an iterable which counts line numbers, strips newlines, and raises
    StopIteration at EOF.
    '''
    def __new__(cls, file):
        '''Create a open_file object.

        Input can be filename string, file descriptor number, or any object
        with `readline'.
        '''

        if isinstance(file, open_file):
            file.reset()
            return file
        of = object.__new__(cls)
        if isinstance(file, (str, unicode)):
            if file.endswith('.gz'):
                of.fid = gzip.open(file)
            else:
                of.fid = open(file, 'rb')
        elif isinstance(file, int):
            of.fid = os.fdopen(file)
        elif hasattr(file, 'readline'):
            of.fid = file
        else:
            raise IOError('Type ' + `type(file)` + ' is not supported.')

        if hasattr(of.fid, 'name'):
            of.name = of.fid.name
        elif hasattr(of.fid, 'filename'):
            of.name = of.fid.filename
        else:
            of.name = str(file)
        try:
            of.fid.read()
        except:
            raise IOError('Error reading file: %s' % of.name)
        of.reset()
        return of

    def next(self):
        '''
        Return the next line, also incrementing `lineno'.
        '''
        self.last_pos = self.fid.tell()
        self.line = self.fid.readline()
        if not self.line:
            raise StopIteration()
        self.lineno += 1
        return self.line.rstrip('\r\n')

    def readline(self):
        '''
        A synonym for next() which doesn't strip newlines or raise
        StopIteration.
        '''
        self.last_pos = self.fid.tell()
        self.line = self.fid.readline()
        if self.line:
            self.lineno += 1
        return self.line

    def backline(self):
        '''
        Set the file's current position to the previous line.
        '''
        self.fid.seek(self.last_pos)

    def values(self, char=' ', fmt=None):
        '''
        Return the line as python list.
        the line is splited by char and formatted as fmt
        fmt could be int, float or a list of types 
        '''
        
        L = [s.strip() for s in self.next().split(char) if s]

        if isinstance(fmt, type):
            return [fmt(s) for s in L]
        elif isinstance(fmt, (list, tuple)):            
            return [f(s) for f,s in zip(fmt, L)]
        else:
            return L                    

    def read(self, n=None):
        self.lineno = None
        if n is None:
            return self.fid.read()
        else:
            return self.fid.read(n)

    def write(self, s):
        return self.fid.write(s)

    def find(self, s):
        return self.fid.find(s)

    def seek(self, pos, start=0):
        self.lineno = None
        self.fid.seek(pos)#, start)

    def tell(self):
        return self.fid.tell()

    def __iter__(self):
        return self

    def reset(self):
        '''Go back to the beginning if possible. Set lineno to 0 regardless.'''
        if hasattr(self.fid, 'seek'):
            try:
                self.fid.seek(0)
            except IOError:
                pass
        self.lineno = 0

    def close(self):
        '''Close the file.  A closed file cannot be used for further I/O.'''
        if hasattr(self.fid, 'fileno') and self.fid.fileno() < 3:
            # Closing stdin, stdout, stderr can be bad
            return
        if hasattr(self.fid, 'close'):
            try:
                self.fid.close()
            except IOError:
                pass

class Config(object):
    '''
    Read configuration file
    '''
    def __init__(self, filename=None):
        if filename is None:
            filename = '/etc/lisn.conf'
        if not os.path.exists(filename):
            raise IOError('Configuration file does not exists')
        self.all = json.load(open(filename))
       
    def username(self):
        return self.all['user']['username']
    username = property(username)
    
    def password(self):
        return self.all['user']['password']
    password = property(password)
    
    def url(self):
        return self.all['server']['url']
    url = property(url)
    
    def url_login(self):
        return os.path.join(self.url, self.all['server']['url_login'])
    url_login = property(url_login)
    
    def url_logout(self):
        return os.path.join(self.url, self.all['server']['url_logout'])
    url_logout = property(url_logout)
    
    def mag_server(self, server=None):
        if server is not None:
            return self.all['mag_server'][server]
        else:
            return self.all['mag_server']
    
    def gps_server(self, server=None):
        if server is not None:
            return self.all['gps_server'][server]
        else:
            return self.all['gps_server']
        
    def rinex_server(self, server=None):
        if server is not None:
            return self.all['rinex_server'][server]
        else:
            return self.all['rinex_server']
    
    def notification(self):
        return self.all['notification']

if __name__=='__main__':
    pass