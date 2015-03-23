'''
Module containing a UTC-based datetime class with additional gps time constants.
'''

import os
import re
import time
from datetime import datetime, timedelta, date
from math import fmod

GPS_EPOCH_JDAY = 2444245
TIMESTAMP0 = datetime(1970, 1, 1, 0, 0)
GPSTIMESTAMP0 = datetime(1980, 1, 6, 0, 0)

class GPSDateTime(object):
    '''
    A UTC-based datetime object.

    This datetime class is based on the POSIX time, a system for describing
    instants in time, defined as the number of seconds elapsed since midnight
    Coordinated Universal Time (UTC) of Thursday, January 1, 1970. Using a
    single float timestamp allows higher precision as the default Python
    :class:`datetime.datetime` class. It features additional string patterns 
    during object initialization and GPS time constants.
    
    '''
    
    timestamp = 0.0
    DEFAULT_PRECISION = 6
    
    def __init__(self, *args, **kwargs):
        '''
        Creates a new GPSDateTime object.
        '''
        self.precision = kwargs.pop('precision', self.DEFAULT_PRECISION)
        filename = kwargs.pop('filename', False)
        if len(args) == 0 and len(kwargs) == 0:
            self.timestamp = time.time()
            return
        elif len(args) == 1:
            value = args[0]
            # check types
            try:
                # got a timestamp
                self.timestamp = value.__float__()
                return
            except:
                pass
            if isinstance(value, datetime):
                return self.from_datetime(value)
            elif isinstance(value, date):
                dt = datetime(value.year, value.month, value.day)
                return self.from_datetime(dt)
            elif isinstance(value, basestring):
                value = value.strip()
                if filename is True:
                    value = ''.join([n for n in re.findall(r'\d+', value) if len(n)>1])
                elif isinstance(filename, basestring) and filename.lower() == 'rinex':
                    value = value[4:7]+value[9:11]
                # try to apply some standard patterns
                value = value.replace('T', ' ')
                value = value.replace('_', ' ')
                value = value.replace('-', ' ')
                value = value.replace(':', ' ')
                value = value.replace('/', ' ')
                value = value.replace(',', ' ')
                value = value.replace('Z', ' ')
                value = value.replace('W', ' ')
                # check for ordinal date (julian date)
                parts = [val for val in value.split(' ') if val]
                # check for patterns                
                if isinstance(filename, basestring) and filename.lower() == 'rinex':
                    pattern = "%j%y"
                elif len(parts) == 1 and len(value) == 7 and value.isdigit():
                    pattern = "%Y%j"
                elif len(parts) == 2 and len(parts[1]) == 3 and parts[1].isdigit():                    
                    value = ''.join(parts)
                    pattern = "%Y%j"
                elif len(parts) == 2 and len(parts[0]) == 7 and len(parts[1]) == 6:                    
                    value = ''.join(parts)
                    pattern = "%Y%j%H%M%S"
                elif len(parts) == 3 and len(parts[1]) == 3 and len(parts[2]) == 6:                    
                    value = ''.join(parts)
                    pattern = "%Y%j%H%M%S"
                else:
                    # some parts should have 2 digits
                    for i in range(1, min(len(parts), 6)):
                        if len(parts[i]) == 1:
                            parts[i] = '0' + parts[i]
                    # standard date string
                    value = ''.join(parts)
                    if len(value) == 6:
                        pattern = "%y%m%d"
                    elif len(value) == 10:
                        pattern = "%y%m%d%H%M"
                    elif len(value) == 12:
                        pattern = "%y%m%d%H%M%S"
                    elif len(value) == 14:
                        pattern = "%Y%m%d%H%M%S"
                    else:
                        pattern = "%Y%m%d"
                ms = 0
                if '.' in value:
                    parts = value.split('.')
                    value = parts[0].strip()
                    try:
                        ms = float('.' + parts[1].strip())
                    except:
                        pass
                # all parts should be digits now - here we filter unknown
                # patterns and pass it directly to Python's  datetime
                if not ''.join(parts).isdigit():
                    dt = datetime(*args, **kwargs)
                    self.from_datetime(dt)
                    return                
                
                dt = datetime.strptime(value, pattern)
                self.from_datetime(dt, ms)
                return
        # check new kwargs
        if 'julday' in kwargs:
            jd     = kwargs.pop('julday')
            a      = int(jd+0.5)
            b      = a+1537
            c      = int((b-122.1)/365.25)
            d      = int(365.25*c)
            e      = int((b-d)/30.6001)
            td     = b-d-int(30.6001*e)+fmod(jd+0.5, 1.0)
            day    = int(td)
            td    -= day
            td    *= 24.
            hour   = int(td)
            td    -= hour
            td    *= 60.
            minute = int(td)
            td    -= minute
            td    *= 60.
            second = int(td)
            month  = (e-1-12*int(e/14))
            year   = (c-4715-int((7+month)/10.))
            
            kwargs['year']   = year
            kwargs['month']  = month
            kwargs['day']    = day
            kwargs['hour']   = hour
            kwargs['minute'] = minute
            kwargs['second'] = second

        if 'doy' in kwargs and 'year' in kwargs:
            year  = kwargs['year']
            day   = kwargs.pop('doy')
            month = 0
            while day>0:
                nd     = GPSDateTime().days_in_month(year, month+1)
                day   -= nd
                month += 1
            day += nd
            kwargs['month'] = month
            kwargs['day']   = day

        if 'GPSweek' in kwargs:
            week = kwargs.pop('GPSweek')
            sow  = kwargs.pop('GPSsow', 0)
            dt   = GPSTIMESTAMP0 + timedelta(days=7*week, seconds=sow)
            self.from_datetime(dt)
            return
        # check if seconds are given as float value
        if len(args) == 6 and isinstance(args[5], float):
            kwargs['microsecond'] = int(args[5] % 1 * 1000000)
            kwargs['second']      = int(args[5])
            args                  = args[0:5]
        dt = datetime(*args, **kwargs)
        self.from_datetime(dt)        

    def from_datetime(self, dt, ms=0.0):
        '''
        Set timestamp from a datetime object
        '''
        try:
            td = (dt - TIMESTAMP0)
        except TypeError:
            td = (dt.replace(tzinfo=None) - dt.utcoffset()) - TIMESTAMP0
        self.timestamp = (td.microseconds + (td.seconds + td.days * 86400) *
                          1000000) / 1000000.0 +ms    
    
    def _set(self, **kwargs):
        '''
        Sets current timestamp using kwargs.
        '''
        year = kwargs.get('year', self.year)
        month = kwargs.get('month', self.month)
        day = kwargs.get('day', self.day)
        hour = kwargs.get('hour', self.hour)
        minute = kwargs.get('minute', self.minute)
        second = kwargs.get('second', self.second)
        microsecond = kwargs.get('microsecond', self.microsecond)
        doy = kwargs.get('doy', None)
        if doy:
            self.timestamp = GPSDateTime(year=year, doy=doy, hour=hour,
                                         minute=minute, second=second,
                                         microsecond=microsecond).timestamp
        else:
            self.timestamp = GPSDateTime(year, month, day, hour, minute,
                                         second, microsecond).timestamp
    
    def get_timestamp(self):
        return self.timestamp


    def get_datetime(self):
        return TIMESTAMP0 + timedelta(seconds=self.timestamp)

    datetime = property(get_datetime)

    def get_year(self):
        return self.get_datetime().year
    
    def set_year(self, value):
        self._set(year=value)

    year = property(get_year, set_year)
    
    def get_month(self):
        return self.get_datetime().month
    
    def set_month(self, value):
        self._set(month=value)

    month = property(get_month, set_month)
    
    def get_day(self):
        return self.get_datetime().day
    
    def set_day(self, value):
        self._set(day=value)

    day = property(get_day, set_day)
    
    def get_hour(self):
        return self.get_datetime().hour
    
    def set_hour(self, value):
        self._set(hour=value)

    hour = property(get_hour, set_hour)
    
    def get_minute(self):
        return self.get_datetime().minute
    
    def set_minute(self, value):
        self._set(minute=value)

    minute = property(get_minute, set_minute)
    
    def get_second(self):
        return self.get_datetime().second
    
    def set_second(self, value):
        self._set(second=value)

    second = property(get_second, set_second)
    
    def get_microsecond(self):
        return self.get_datetime().microsecond
    
    def set_microsecond(self, value):
        self._set(microsecond=value)

    microsecond = property(get_microsecond, set_microsecond)
    
    def get_doy(self):
        '''
        Returns the day of the year of the current GPSDateTime object.
        '''        
        return self.utctimetuple().tm_yday
    
    def set_doy(self, value):
        '''
        Set the day of the year of the current GPSDateTime object.
        '''
        self._set(doy=value)
    
    doy = property(get_doy, set_doy)

    def get_julday(self):
        '''
        Returns the Julian day of the current GPSDateTime object.
        '''
        yy, mm, dd     = self.day, self.month, self.year
        hour, minute, second = self.hour, self.minute, self.second
        
        
        a = (14 - self.month)//12
        y = self.year + 4800 - a
        m = self.month + 12*a - 3
        return self.day + ((153*m + 2)//5) + 365*y + y//4 - y//100 + y//400 - 32045

        
        
        return dd-32075+1461*(yy+4800+(mm-14)/12)/4+ \
                367*(mm-2-(mm-14)/12*12)/12-3*((yy+4900+(mm-14)/12)/100)/4
        

        a = (14 - self.month)//12
        y = self.year + 4800 - a
        m = self.month + 12*a - 3
        jd = self.day + ((153*m + 2)//5) + 365*y + y//4 - y//100 + y//400 - 32045

        if self.month<= 2:
            y = self.year - 1
            m = self.month + 12
        else:
            y = self.year
            m = self.month
        jd1 = int(365.25*y)+int(30.6001*(m+1.))+self.day+self.hour/24.+self.minute/1440.+self.second/86400.+1720981.5

        return jd1

    julday = property(get_julday)

    def get_GPSweek(self):
        """
        Returns the GPS Week of the current GPSDateTime object.
        """
        if self >= GPSTIMESTAMP0:
            return (self.datetime-GPSTIMESTAMP0).days/7
        else:
            return 0
        

    GPSweek = property(get_GPSweek)

    def get_GPSdow(self):
        """
        Returns the day of the Week number of the current GPSDateTime object.
        """
        
        return int((self.julday%7+1)%7)

    GPSdow  = property(get_GPSdow)

    def get_GPStow(self):
        """
        Returns the seconds of the Week of the current GPSDateTime object.
        """
        
        return self-GPSDateTime(GPSweek=self.GPSweek)

    GPStow  = property(get_GPStow)

    def get_doys(self):
        """
        Returns a decimal day of the year of the current GPSDateTime object.
        """
        return self.doy+self.hours/24.

    doys = property(get_doys)

    def get_hours(self):
        """
        Returns a decimal hours of the current GPSDateTime object.
        """
        return self.hour+self.minute/60.+self.second/3600.

    hours = property(get_hours)

    def get_seconds(self):
        """
        Returns total of seconds of the current GPSDateTime object.
        """
        return self.hour*3600+self.minute*60+self.second

    seconds = property(get_seconds)

    def __add__(self, value):
        '''
        Adds seconds and microseconds to current GPSDateTime object.
        '''
        if isinstance(value, timedelta):
            # see timedelta.total_seconds
            value = (value.microseconds + (value.seconds + value.days *
                     86400) * 1000000) / 1000000.0
        return GPSDateTime(self.timestamp + value)

    def __sub__(self, value):
        '''
        Subtracts seconds and microseconds from current GPSDateTime object.
        '''
        if isinstance(value, (GPSDateTime, datetime)):
            return round(self.timestamp - GPSDateTime(value).timestamp, self.precision)
        elif isinstance(value, timedelta):
            # see timedelta.total_seconds
            value = (value.microseconds + (value.seconds + value.days *
                     86400) * 1000000) / 1000000.0
        return GPSDateTime(self.timestamp - value)

    def __float__(self):
        return self.timestamp
    
    def __int__(self):
        return int(self.timestamp)
    
    def __eq__(self, other):
        '''
        Rich comparison operator '=='
        '''
        try:
            return round(self.timestamp - GPSDateTime(other).timestamp, self.precision) == 0
        except (TypeError, ValueError):
            return False

    def __ne__(self, other):
        '''
        Rich comparison operator '!='
        '''
        return not self.__eq__(other)

    def __lt__(self, other):
        '''
        Rich comparison operator '<'.
        '''
        try:
            return round(self.timestamp - GPSDateTime(other).timestamp, self.precision) < 0
        except (TypeError, ValueError):
            return False

    def __le__(self, other):
        '''
        Rich comparison operator '<='.
        '''
        try:
            return round(self.timestamp - GPSDateTime(other).timestamp, self.precision) <= 0
        except (TypeError, ValueError):
            return False

    def __gt__(self, other):
        '''
        Rich comparison operator '>'
        '''
        try:
            return round(self.timestamp - GPSDateTime(other).timestamp, self.precision) > 0
        except (TypeError, ValueError):
            return False

    def __ge__(self, other):
        '''
        Rich comparison operator '>='
        '''
        try:
            return round(self.timestamp - GPSDateTime(other).timestamp, self.precision) >= 0
        except (TypeError, ValueError):
            return False
    
    def __unicode__(self):
        return str(self.__str__())
    
    def __str__(self):
        return self.strftime('%Y-%m-%d %H:%M:%S.%fUTC')    
    
    def __repr__(self):
        return 'GPSDateTime' + self.get_datetime().__repr__()[17:]
    
    def __hash__(self):
        return hash(self.timestamp)
    
    def __abs__(self):
        return abs(self.timestamp)
    
    def replace(self, **kwargs):
        return GPSDateTime(self.datetime.replace(**kwargs))
    
    def date(self):
        return GPSDateTime(datetime(self.year, self.month, self.day))
    
    def strftime(self, fmt):
        return self.get_datetime().strftime(fmt)
    
    @staticmethod
    def strptime(s, fmt):
        return GPSDateTime(datetime.strptime(s, fmt))
    
    def timetz(self):        
        return self.get_datetime().timetz()
    
    def utcoffset(self):
        return self.get_datetime().utcoffset()
    
    def dst(self):
        return self.get_datetime().dst()
    
    def tzname(self):
        return self.get_datetime().tzname()
    
    def ctime(self):
        return self.get_datetime().ctime()
    
    def isoweekday(self):
        return self.get_datetime().isoweekday()
    
    def isocalendar(self):
        return self.get_datetime().isocalendar()
    
    def isoformat(self, sep="T"):
        return self.get_datetime().isoformat(sep=sep)
    
    def timetuple(self):
        return self.get_datetime().timetuple()

    def utctimetuple(self):        
        return self.get_datetime().utctimetuple()
    
    @staticmethod
    def today():
        dt = datetime.now()
        return GPSDateTime(dt.year, dt.month, dt.day)
    
    @staticmethod    
    def now():       
        return GPSDateTime()
        
    @staticmethod
    def utcnow():      
        return GPSDateTime()
    
    def toordinal(self):       
        return self.get_datetime().toordinal()

    def get_str(self, all=False):
        """
        Returns common string representation.
        """
        if all:
            return self.strftime('%y%m%d_%H%M%S')
        else:
            return self.strftime('%y%m%d')

    def is_leap_year(self, year=None):
        '''
        Returns True if year is a leap year, False otherwise.
        '''
        if not year: year = self.year
        if year%4 == 0:
            if year%100 == 0:
                if year%400 == 0:
                    return True
                else:
                    return False
            else:
                return True
        else:
            return False

    def days_in_month(self, year=None, month=None):
        '''
        Return the number of days in month for the given month and year
        '''
        if not year : year  = self.year
        if not month: month = self.month
        if month == 2:
            if GPSDateTime().is_leap_year(year):
                return 29
            else:
                return 28
        elif month in (1, 3, 5, 7, 8, 10, 12):
            return 31
        else:
            return 30

if __name__ == '__main__':
    dt = GPSDateTime(2014,12,31)
    print 'Date           :', dt
    print 'Julday         :', dt.julday
    print 'Day of the year:', dt.doy
    print 'GPS week       :', dt.GPSweek
    print 'Day of week    :', dt.GPSdow
    print 'GPS seconds    :', dt.GPStow
    print 'Timestamp      :', dt.timestamp
    