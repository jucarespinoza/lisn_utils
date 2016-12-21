# -*- coding: utf-8 -*-
'''
Concatenation and decimation of magnetometer files

@author: Juan C. Espinoza
@contact: jucar.espinoza@gmail.com

'''

import gzip
from utility import open_file
from gpsdatetime import GPSDateTime
from plotter import plot_mag

class MAGData(list):
    '''
    '''
    vars =  ('DD', 'MM', 'YYYY', 'hh', 'mm', 'D(deg)', 'H(nT)', 
             'Z(nT)', 'I(deg)', 'F(nT)', 'T1(deg)', 'T2(ºC)')

    def __init__(self, url, station={}, date=None, format=None):
        '''
        '''
        self.arcs = []
        self.station = station        
        self.header = ''
                
        try:              
            print 'Opening: %s' % url
            fo = open_file(url)
            if format is None:
                self.header = fo.next()
        except BaseException as e:
            print e
            return
        if 'code' not in station:
            self.station['code'] = fo.name.split('/')[-1][:4]
        while True:
            try:
                if format is None:            
                    values = fo.values(fmt=(int,int,int,int,int,float,float,float,float,float))
                else:
                    values = fo.values(char=',', fmt=(int,int,int,int,float,float,float,float,float, float, float))
            except StopIteration:
                break
            except (ValueError, TypeError):
                continue

            if not values or len(values)not in (10, 11): continue
            if format is None:
                try:
                    rec = MAGRecord(values)
                except:
                    continue
            else:
                rec = MAGRecordCVS(values)
            if date is not None and date<>rec.epoch.date():
                continue
            if len(self)>0 and rec.epoch<=self[-1].epoch:
                continue
            self.append(rec)
            self.checkbreak()
        fo.close()
        if len(self)>0:
            self.date = self[-1].epoch.date()
        else:
            self.date = None
    
    def __str__(self):
        '''
        '''
        return 'MAGData[site=%s, records=%s]' % (self.station['code'], len(self))
    
    def checkbreak(self):
        '''
        '''

        if not self.arcs:
            self.arcs.append([len(self)-1, len(self)-1])
        elif self[-1].epoch-self[self.arcs[-1][-1]].epoch>60:
            self.arcs.append([len(self)-1, len(self)-1])
        else:
            self.arcs[-1][-1] = len(self)-1
        
    def iterlist(self, keys=['H']):
        '''
        '''
        if len(keys)==1:
            return [rec[keys[0]] for rec in self]
        else:
            return [[rec[key] for key in keys]for rec in self]
    
    def save(self, filename, date=[], gz=True):
        '''
        '''
        if isinstance(date, (list, tuple)):
            date = [dt.date() for dt in date]
        elif isinstance(date, (datetime, GPSDateTime)):
            date = [GPSDateTime(date).date()]
        else:
            date = []            
        
        if gz:
            if not filename.endswith('.gz'):
                filename += '.gz'
            fo = gzip.open(filename, 'wb')
        else:
            fo = open(filename, 'wb')
        
        vars = self[0].keys()
        if len(vars)==5:
            fmt = '%s%9.4f%9.1f%9.1f%9.4f%9.1f\n'
            values = ('D', 'H', 'Z', 'I', 'F')
        elif len(vars)==6:
            fmt = '%s%9.4f%9.1f%9.1f%9.4f%9.1f%9.2f\n'
            values = ('D', 'H', 'Z', 'I', 'F', 'T1')
        elif len(vars)==7:
            fmt = '%s%9.4f%9.1f%9.1f%9.4f%9.1f%9.2f%9.2f\n'
            values = ('D', 'H', 'Z', 'I', 'F', 'T1', 'T2')
        
        pos = self.header.find('<')
        if pos>=0:
            if date:
                dt = date[0]
            else:
                dt=self[0].epoch
            self.header = self.header.replace(self.header[pos:pos+5], '<%03d>' % dt.doy)
        
        fo.write('%s\n\n' % self.header)
        fo.write('%3s%3s%5s%3s%3s%9s%9s%9s%9s%9s%9s%9s\n\n' % self.vars)
        cnt = 0
        for rec in self:
            if date and rec.epoch.date() not in date:
                continue
            allvalues = (rec.epoch.strftime(' %d %m %Y %H %M'),)+tuple([rec[x] for x in values])
            fo.write(fmt % allvalues)
            cnt += 1
        fo.close()
        
        print '%s of %s records saved at %s' % (cnt, len(self), filename)
        
    def get_records(self, dt1, dt2=None):
        '''
        '''
        if dt2 is None:
            return [rec for rec in self if rec.epoch.date()==dt1.date()]
        else:
            return [rec for rec in self if dt1<=rec.epoch<=dt2]
    
    def merge(self, data_i):
        '''
        Merge current magdata object with another.
        '''

        if isinstance(data_i, str):
            data_i = MAGData(data_i)
        data_j = self[:]
        self.arcs = []
        
        del self[:]
        
        i = 0
        j = 0
        imax = len(data_i)
        jmax = len(data_j)

        while True:
            if i==imax or j==jmax:
                break          
            if data_i[i]==data_j[j]:
                self.append(data_j[j])
                i+=1
                j+=1
            elif data_i[i]>data_j[j]:
                self.append(data_j[j])
                j+=1
            else:
                self.append(data_i[i])
                i+=1 
            self.checkbreak()            
            
        if i==imax:
            temp = data_j[j:]
        else:
            temp = data_i[i:]
        
        for rec in temp:
            self.append(rec)
            self.checkbreak()
        
        self.header = data_i.header            

        print 'Merging %d + %d -> %d records' % (jmax, imax, len(self))
        
    def plot(self, **kwargs):
        '''
        '''
        return plot_mag(self, station=self.station, **kwargs)
    
class MAGRecord(dict):
    '''
    '''
    def __init__(self, values):
        '''
        '''  
        self.epoch = GPSDateTime(values[2],values[1],values[0],values[3],values[4])
        self['D'] = values[5]
        self['H'] = values[6]
        self['Z'] = values[7]
        self['I'] = values[8]
        self['F'] = values[9]
        try:
            self['T1'] = values[10]
            self['T2'] = values[11]
        except IndexError:
            pass

    def __getitem__(self, index):
        '''
        '''
        if index == 'epoch':
            return self.epoch
        else:
            return dict.__getitem__(self, index)
        
    def __eq__(self, other):
        '''
        Equal comparison reference to the epoch of the record
        '''        
        return self.epoch==other.epoch

    def __ne__(self, other):
        '''
        Not equal comparison reference to the epoch of the record
        '''
        return self.epoch!=other.epoch

    def __lt__(self, other):
        '''
        Comparison reference to the epoch of the record
        '''
        return self.epoch<other.epoch

    def __le__(self, other):
        '''
        Comparison reference to the epoch of the record
        '''
        return self.epoch<=other.epoch

    def __gt__(self, other):
        '''
        Comparison reference to the epoch of the record
        '''
        return self.epoch>other.epoch

    def __ge__(self, other):
        '''
        Comparison reference to the epoch of the record
        '''
        return self.epoch>=other.epoch

class MAGRecordCVS(MAGRecord):
    '''
    '''
    def __init__(self, values):
        '''
        '''  
        self.epoch = GPSDateTime(values[0],values[1],values[2])+values[3]*60
        self['D'] = values[6]/100
        self['H'] = values[5]
        self['Z'] = values[9]
        self['I'] = values[8]
        self['F'] = values[9]
        try:
            self['T1'] = values[10]
            self['T2'] = values[11]
        except IndexError:
            pass


if __name__=='__main__':
    pass
