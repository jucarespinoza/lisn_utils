# -*- coding: utf-8 -*-
'''
Constants
'''

GPS_GRAVITY_CONSTANT                       =  3.986005e14
GPS_CLOCK_CORRECTION_RELATIVISTIC_CONSTANT = -4.442807633e-10
GPS_WGS84_EARTH_ROTATION_RATE              =  7.2921151467e-05
LIGHT_SPEED                                =  299792458.0

SECONDS_IN_WEEK = 7*86400
#TECUns          = 2.854       # TECU/ns according to GPS-Scinda, Charles Carrano, 4-7-08
TECUns        = 2.852       # TECU/ns according to TECalc_rinex, Pat Doherty, 2-21-94
F1              = 1.57542     # L1 Frequency (GHz)
F2              = 1.22760     # L2 Frequency (GHz)

gps_names_new = {
             '.nvd'         :('%s.%s.nvd',    'nvd', 'binary', ''),
             '.nvd.gz'      :('%s.%s.nvd.gz', 'nvd', 'binary', ''),
             '.lb2'         :('%s.%s.lb2',    'lb2', 'binary', ''),
             'gps.dat.gz'   :('%s.%s.s4.gz',  's4',  'scint',  ''),
             's4.txt'       :('%s.%s.txt',    's4',  'scint',  ''),
             'pos.dat.gz'   :('%s.%s.pos.gz', 'pos', 'posit',  ''),
             'scn.dat.gz'   :('%s.%s.scn.gz', 'scn', 'stats',  ''),
             '.obs'         :('%s.%s.obs',    'obs', 'binary', ''),
             '.obs.gz'      :('%s.%s.obs.gz', 'obs', 'binary', ''),
             '.rnx'         :('%s.%s.rnx',    'rnx', 'rinex',  ''),
             '.tec'         :('%s.tec',       'tec', 'tec',    ''),
             }

gps_names_old = {
             '.nvd'         :('%s.%s.nvd',    'nvd', 'binary', ''),
             '.lb2'         :('%s.%s.lb2',    'lb2', 'binary', ''),
             'gps.dat.gz'   :('%ss4.%s.gz',   's4',  'scint',  's4'),
             's4.txt'       :('%ss4.%s.txt',  's4',  'scint',  ''),
             'pos.dat.gz'   :('%s.%s.pos.gz', 'pos', 'posit' , ''),
             '.obs'         :('%s.obs.%s.gz', 'obs', 'binary', 'obs'),
                }

rnx_names = {
             'd.tar.gz'   : 'lisn',
             '1.zip'      : 'rbmc',
             'd.Z'        : 'unavco'
             }


lisn_header = {
               'agency'   : 'Low-Latitude Ionospheric Sensor Network',
               'run_by'   : 'LISN',
               'observer' : 'Juan C. Espinoza',
               }

lisn_comment = \
'''
For permission to use this data, contact:
Cesar Valladares (cesar.valladares@bc.edu)
Boston College, Institute for Scientific Research
COMM: 617-552-8789, 617-552-4328
                                                            
For technical details about the data or software contact:
Juan C. Espinoza (juan.espinoza@jro.igp.gob.pe)
Jicamarca Radio Observatory, LISN-Engineering
COMM: 511-317-2313

S1, if present, is the SNR for the C/A data stream on L1.
SNR is mapped to RINEX snr flag value [2-9]
L1&L2:  = 26dBHz -> 1; 26-28dBHz -> 2; 28-32dBHz -> 3
       32-36dBHz -> 4; 36-39dBHz -> 5; 39-42dBHz -> 6
       42-45dBHz -> 7; 45-49dBHz -> 8;  >=49dBHz -> 9

'''

STATION_VARS = ['station_id', 'name', 'country', 'latitude', 'longitude', 'altitude', 'instrument_id',
                'lisn_code', 'code', 'itype_id', 'network', 'server', 'status',
                'last_file', 'last_file_date', 'receiver', 'antenna', 'receiver_number', 'antenna_number',
                'receiver_firmware']

PLOT_LABELS = {'lat': 'Latitude',
               'ele': 'Elevation',
               'azi': 'Azimuth',
               'lon': 'Longitude',
               'sTEC': 'Slant TEC',
               'TEC': 'Uncalibrated TEC',
               'eqTEC': 'Equivalent TEC',
               'aTEC': 'TEC from Code',
               'rTEC': 'TEC from Phase',
               's4': 'S4 index',
               'epoch': 'Time'}

PLOT_UNITS = {'lat': 'Degrees [deg]',
              'ele': 'Degrees [deg]',
              'azi': 'Degrees [deg]',
              'lon': 'Degrees [deg]',
              'sTEC': 'TEC [UNITS]',
              'TEC': 'TEC [UNITS]',
              'eqTEC': 'TEC [UNITS]',
              'aTEC': 'TEC [UNITS]',
              'rTEC': 'TEC [UNITS]',
              's4': 'S4 index',
              'epoch': 'Time'}
