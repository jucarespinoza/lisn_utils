# -*- coding: utf-8 -*-
'''
Utilities to plot gps data

@author: Juan C. Espinoza
@contact: jucar.espinoza@gmail.com

'''

__version__ = 1.0

import os, time
import numpy as np
from StringIO import StringIO
from pkg_resources import resource_string
import utility as utl
import gps
from gpsdatetime import GPSDateTime, timedelta
from headers import PLOT_LABELS, PLOT_UNITS

deg2rad = utl.deg2rad
PI      = utl.PI
PLT     = False
MPL     = False
Basemap = False

def setup(environ=None):
    '''
    '''
    global PLT, MPL, Basemap
    if environ=='web':
        os.environ['HOME'] = '/var/www/mpl'
    import matplotlib as MPL
    if isinstance(environ, basestring):
        MPL.use('Agg')
    import matplotlib.pyplot as PLT
    try:
        from mpl_toolkits.basemap import Basemap
    except ImportError:
        Basemap = False

class MyFigure(object):
    '''
    '''
    def __init__(self, nrows=1, ncols=1, figname=None, **kwargs):
        '''
        '''
        global PLT, MPL

        if not MPL:
            setup(environ=figname)
        
        self.LC = MPL.collections.LineCollection
        self.CM = MPL.cm
        self.py = PLT
        self.mpl = MPL
        self.figure  = PLT.figure(**kwargs)
        self.figname = figname
        self.axes    = []
        self.nrows   = int(nrows)
        self.ncols   = int(ncols)
        self.rowcnt  = 1
        self.colcnt  = 1
        jetcmap = self.CM.get_cmap("spectral", 20) 
        jet_vals = jetcmap(np.arange(20))[3:19]
        self.cmap = self.mpl.colors.LinearSegmentedColormap.from_list("new", jet_vals, 100)

    def adjust(self, left=None, bottom=None, right=None, top=None, wspace=None,
               hspace=None):
        '''
        Subplots adjust in figure
        '''
        width  = self.figure.get_figwidth()
        height = self.figure.get_figheight()
        if not bottom:  bottom  = 0.5*height**-1
        if not top:     top     = 0.74*height**0.11
        if not left:    left    = 0.125
        if not right:   right   = 0.9
        if hspace is None:  hspace  = 0.2
        if wspace is None:  wspace  = 0.2
        self.figure.subplots_adjust(left, bottom, right, top, wspace, hspace)

    def suptitle(self, s, **kwargs):
        '''
        '''
        x = kwargs.pop('x', 0.5)
        y = kwargs.pop('y', 1-0.15/self.figure.get_figheight())
        if ('horizontalalignment' not in kwargs) and ('ha' not in kwargs):
            kwargs['horizontalalignment'] = 'center'

        if ('verticalalignment' not in kwargs) and ('va' not in kwargs):
            kwargs['verticalalignment'] = 'top'

        t = self.figure.text(x, y, s, **kwargs)
        return t

    def text(self, x, y, s, **kwargs):
        '''
        '''
        t = self.figure.text(x, y, s, **kwargs)
        return t

    def show(self, watermark=True):
        '''
        '''
        if watermark:
            self.figure.text(1, 0, 'LISN\nLow-Latitude Ionospheric Sensor Network',
                             fontsize=8, color='gray', ha='right', va='bottom',
                             alpha=0.5)
        if self.figname == 'web':
            imgdata = StringIO()
            self.figure.savefig(imgdata)
            self.py.close(self.figure)
            imgdata.seek(0)
            return imgdata.read()
        elif self.figname:
            if not os.path.exists(os.path.dirname(self.figname)):
                os.makedirs(os.path.dirname(self.figname))
            self.figure.savefig(self.figname)
            self.py.close(self.figure)
            print 'figure saved: %s' % self.figname
        else:
            PLT.show()

    def add_axes(self, shared=False, **kwargs):
        '''
        '''
        polar = kwargs.pop('polar', False)
        axisbg = kwargs.pop('axisbg', None)
        if self.colcnt%(self.ncols+1)==0:
            self.colcnt = 1
            self.rowcnt += 1
        
        ax = MyAxes(**kwargs)
        pos = self.colcnt+((self.ncols)*(self.rowcnt-1))
        ax.n = pos
        ax.ax = self.figure.add_subplot(self.nrows, self.ncols, pos, polar=polar, axisbg=axisbg)
        if shared:
            ax.shax = self.figure.add_axes(ax.ax.get_position(), sharex=ax.ax,
                                           frameon=False)
            ax.shax.yaxis.tick_right()
            ax.shax.yaxis.set_label_position('right')
            ax.ax.yaxis.tick_left()
        self.axes.append(ax)
        self.colcnt += 1
        return ax

    def set_labels(self, **kwargs):
        '''
        '''
        for ax in self.axes:
            ax.set_labels(**kwargs)

    def set_xticks(self, hidden=[], **kwargs):
        '''
        '''
        for ax in self.axes:
            if ax.n in hidden:
                ax.set_xticks(True, **kwargs)
            else:
                ax.set_xticks(**kwargs)

    def set_yticks(self, hide_labels=False, **kwargs):
        '''
        '''
        for ax in self.axes:
            if hide_labels:
                ax.set_yticks(self.ncols, **kwargs)
            else:
                ax.set_yticks(**kwargs)

class MyAxes(object):
    '''
    '''
    def __init__(self, **kwargs):
        self.n            = 1
        self.title        = kwargs.get('title', '')
        self.text         = kwargs.get('text', '')
        self.xlabel       = kwargs.get('xlabel', '')
        self.ylabel       = kwargs.get('ylabel', '')
        self.y2label      = kwargs.get('y2label', '')
        self.xticks       = kwargs.get('xticks', [])
        self.xticklabels  = kwargs.get('xticklabels', [])
        self.yticks       = kwargs.get('yticks', [])
        self.yticklabels  = kwargs.get('yticklabels', [])
        self.y2ticks      = kwargs.get('y2ticks', [])
        self.y2ticklabels = kwargs.get('y2ticklabels', [])
        self.ylim         = kwargs.get('ylim', None)
        self.y2lim        = kwargs.get('y2lim', None)
        self.xlim         = kwargs.get('xlim', None)
        self.xminmax      = set()
        self.yminmax      = set()
        self.y2minmax     = set()
        self.ax           = None
        self.shax         = None

    def set_labels(self,  **kwargs):
        '''
        '''
        if 'size' not in kwargs or 'fontsize' not in kwargs:
            kwargs['fontsize'] = 10
                
        if self.ylabel: self.ax.set_ylabel(self.ylabel, **kwargs)
        if self.xlabel: self.ax.set_xlabel(self.xlabel, **kwargs)
        if self.title : self.ax.set_title(self.title, **kwargs)
        
        if self.shax:                        
            if self.y2label: self.shax.set_ylabel(self.y2label, **kwargs)
            
        if self.text:
            xlim = self.ax.get_xlim()
            ylim = self.ax.get_ylim()
            self.ax.text(xlim[0]+(xlim[1]-xlim[0])*0.95, ylim[1]*0.8,
                         self.text, color='0.4', size=8, ha='right', va='top')

    def set_xticks(self, hidden=False, **kwargs):
        '''
        '''
        
        if self.xticks:
            self.ax.set_xticks(self.xticks)
        if hidden:
            ticklabels = self.ax.get_xticklabels()
            PLT.setp(ticklabels, visible=False)
            self.ax.set_xlabel('')
        else:
            if self.xticklabels:
                self.ax.set_xticklabels(self.xticklabels, **kwargs)
            elif self.xticks:
                self.ax.set_xticklabels(self.xticks, **kwargs)

        if self.xlim  : self.ax.set_xlim(*self.xlim)

        if self.shax:
            ticklabels = self.shax.get_xticklabels()
            PLT.setp(ticklabels, visible=False)

    def set_yticks(self, ncols=None, **kwargs):
        '''
        '''
        if self.yticks:
            self.ax.set_yticks(self.yticks)

        if self.yticklabels:
            self.ax.set_yticklabels(self.yticklabels, **kwargs)
        elif self.yticks:
            self.ax.set_yticklabels(self.yticks, **kwargs)

        if self.ylim  : self.ax.set_ylim(*self.ylim)
        
        if self.shax:
            if self.y2ticks:
                self.shax.set_yticks(self.y2ticks)
            if self.y2ticklabels:
                self.shax.set_yticklabels(self.y2ticklabels, **kwargs)
            elif self.y2ticks:
                self.shax.set_yticklabels(self.y2ticks, **kwargs)

            if self.y2lim  : self.shax.set_ylim(*self.y2lim)

        if ncols:
            cols = range(ncols)
            if self.col in cols[:-1] and self.shax:
                self.shax.set_ylabel('')
                ticklabels = self.shax.get_yticklabels()
                PLT.setp(ticklabels, visible=False)
            if self.col in cols[1:]:
                self.ax.set_ylabel('')
                ticklabels = self.ax.get_yticklabels()
                PLT.setp(ticklabels, visible=False)

    def plot(self, *args, **kwargs):
        '''
        '''

        shared = kwargs.pop('shared', False)

        if shared:
            return self.shax.plot(*args, **kwargs)
        else:
            return self.ax.plot(*args, **kwargs)

def plot_magfield(ax, xpos, dH=None):
    '''
    '''
    X = []
    Y = []
    magfile = resource_string(__name__, 'data/field.dat')
    for line in StringIO(magfile):
        vals = [s.strip() for s in line.split(' ') if s]
        X.append(float(vals[1]))
        Y.append(float(vals[0]))            
    dum = [int(x)==xpos for x in X]
    if True in dum:
        pos = dum.index(True)
    else:
        pos = 0
    xpos = X[pos]
    ypos = Y[pos]
    ax.plot(X,Y, color='k', lw=2)
    if dH is not None:
        ax.plot(X, np.array(Y)-dH, color='k', lw=1)
        ax.plot(X, np.array(Y)+dH, color='k', lw=1)
        ax.text(xpos+0.5, ypos+dH, '+%s deg' % dH, ha='left', 
                va='bottom', size=7, rotation=-10, weight='bold')
    ax.text(xpos+0.6, ypos-0.5, 'Magnetic\n\nEquator', ha='left', 
            va='center', size=7, rotation=-10, weight='bold')

def plot_map(title, magfield=True, colorbar=True, figsize=(6, 6),
             figname=None, lat_limit=(-60, 20), lon_limit=(-90, -30),
             deltaH=None):
    '''
    Create Figure and SA map with coastlines, countries paralels and meridians    
    '''
    
    fig = MyFigure(figsize=figsize, figname=figname)
    
    if Basemap is False:
        raise ImportError('Install matplotlib toolkits (basemap)')
    
    if colorbar:
        axm = fig.figure.add_axes([0.04,0.1,0.9,0.8])
        axc = fig.figure.add_axes([0.87, 0.2, 0.03, 0.6])
    else:    
        axm = fig.figure.add_subplot(1,1,1)
        axc = None
    
    m = Basemap(llcrnrlon=lon_limit[0], llcrnrlat=lat_limit[0],
                  urcrnrlon=lon_limit[1], urcrnrlat=lat_limit[1], 
                  projection='cyl', resolution='l', area_thresh=2500, ax=axm)
    m.drawcoastlines(linewidth=0.5)
    m.drawcountries(linewidth=0.5)
    m.drawmeridians(np.arange(lon_limit[0]+10,lon_limit[1],10), linewidth=0.25, 
                      labels=[1,0,0,1], labelstyle='+/-')
    m.drawparallels(np.arange(lat_limit[0]+10,lat_limit[1],10), linewidth=0.25, 
                      labels=[1,0,0,1], labelstyle='+/-')
    if magfield:
        plot_magfield(axm, lon_limit[0], deltaH)
    tx = axm.get_xlim()[1]
    ty = axm.get_ylim()[0]
    axm.text(tx, ty, 'LISN\nLow-Latitude Ionospheric Sensor Network', 
             size=8, color='black', ha='right', va='bottom', alpha=0.8)
    axm.set_xlabel('Geographic Longitude [deg]', labelpad=20)
    axm.set_ylabel('Geographic Latitude [deg]', labelpad=30)     
    fig.suptitle(title, size=12, weight='bold')
    
    return fig, axm, axc

def plot_s4_map(S4, date, delta=30, figname=None, scale=(0, 1), 
                prns=None, min_s4=0.15, min_ele=30, threshold=1.5):
    '''
    '''
    dt0 = date
    dt1 = date+delta*60
    fig,axm,axc = plot_map('S4 index Map', figname=figname)
    norm = fig.mpl.colors.Normalize(*scale)
    fig.mpl.colorbar.ColorbarBase(axc, cmap=fig.cmap, norm=norm, orientation='vertical')
    axm.set_title('%s     %s - %s (UT)    S4>%3.2f' % (dt0.strftime('%Y/%m/%d'), 
                                                       dt0.strftime('%H:%M'), 
                                                       dt1.strftime('%H:%M'), 
                                                       min_s4),
                  size=11)
    
    for data in S4:
        for prn in data.arcs:
            if (prns is not None and prn not in prns) or int(prn[1:])>32: continue
            x, y, c = data.get_array(prn, ('lon', 'lat', 's4'), dt0, dt1, ('ele', min_ele))
            if len(c)>0 and not np.all(np.isnan(c)):
                s = data.get_array(prn, 's4', values=('ele', min_ele))[0]
                noise = np.nanmean(s)+threshold*np.nanstd(s)
                noise = noise if noise > min_s4 else min_s4
                cm = np.ma.masked_where(((c <= noise)|(~np.isfinite(c))), c)
                cm = cm*(1./(scale[1]-scale[0]))-(scale[0]/(scale[1]-scale[0]))
                axm.scatter(x, y, c=cm, cmap=fig.cmap, norm=norm, s=20, edgecolor='none')
    
    fig.show(watermark=False)

def plot_tec_map(TEC, date, ptype, delta=60, figname=None, scale=(0, 80), 
                prns=None, bad_prns=None, min_ele=30):
    '''
    #TODO
    '''
    dt1 = date
    dt2 = date+delta*60
    ###create figure
    fig,axm,axc = plot_map('Total Electron Content (TEC)', figsize=(7,6), 
                           lat_limit=(-60,40), lon_limit=(-120,-30),
                           figname=figname, deltaH=20)
    norm = fig.mpl.colors.Normalize(*scale)
    ###colorbar    
    cb = fig.mpl.colorbar.ColorbarBase(axc, cmap=fig.cmap, norm=norm, orientation='vertical')
    cb.set_label('TECU', size=10)
    cb.set_ticks(range(0, scale[1]+1, 10))
    ###subtitle
    axm.set_title('%s%30s%s - %s (UT)' % (date.strftime('%Y/%m/%d'), ' ', 
                                          dt1.strftime('%H:%M'), 
                                          dt2.strftime('%H:%M')), 
                   size=11)
    
    if ptype=='measured':    
        for data in TEC:
            for prn in data.arcs:            
                x,y,c = data.get_array(prn, ('lon', 'lat', 'eqTEC'), dt1, dt2, ('ele', min_ele))
                axm.scatter(x, y, c=c, cmap=fig.cmap, norm=norm, s=5, edgecolor='none')
    else:
        lat = 60        
        lon = 120        
        mTEC = np.zeros((90,100))
        iTEC = np.zeros((90,100))
        iMAX = np.ones((90,100))*-99
        #mTEC = np.ma.masked_where((mTEC==0), mTEC)
        for data in TEC:
            for prn in data.arcs:
                for rec in data.get_records(prn, dt1, dt2):
                    if rec[prn]['ele']<min_ele:
                        continue
                    x=int(rec[prn]['lon']+lon+0.5)
                    y=int(rec[prn]['lat']+lat+0.5)
                    mTEC[x,y] = mTEC[x,y]+rec[prn]['eqTEC']
                    iTEC[x,y] = iTEC[x,y]+1
                    if rec[prn]['eqTEC']>iMAX[x,y]: iMAX[x,y] = rec[prn]['eqTEC']
        
        mTEC = np.ma.masked_where((iTEC==0), mTEC)
        Z = mTEC/iTEC
        
        x = np.arange(-5, 7, 1)
        y = np.arange(5, -7, -1)
        
        x, y = np.meshgrid(x, y)
        
        z = Z[70:82,40:52]
        ave = np.nanmean(z)
        std = np.nanstd(z)
        
        p = [z, z*x, z*y, z*x**2, z*x*y, z*y**2, z*x**3, z*y*x**2, z*x*y**2, z*y**3]
        p = [np.nansum(a) for a in p]
        print p
        print ave
        print std
        
                
        #P = (1+x+y+x**2+y**2+x*y+x**2*y+x*y**2+x**3+y**3)*t        
        
        X = np.arange(-lon, -29)
        Y = np.arange(-lat, 41)        
        X, Y = np.meshgrid(X, Y)        
        axm.pcolor(X, Y, Z.T, cmap=fig.cmap, norm=norm)
        
    fig.show(watermark=False)

def plot_bubbles_map(bub, date=None, delta=60, figname=None, prns=None, codes=None,
                     scale=(0, 20)):
    '''
    '''
    if date is None:
        date = bub[0].epoch.date()
    dt1 = date
    dt2 = date+delta*60-1
    if codes is not None:
        codes = [code for code in codes if code in bub.codes]
    else:
        codes = bub.codes
    if prns is not None:
        prns = [prn for prn in prns if prn in bub.prns]
    else:
        prns = bub.prns
    ###create figure
    fig,axm,axc = plot_map('TEC Depletions', figsize=(7,6), 
                           lat_limit=(-60,40), lon_limit=(-120,-30),
                           figname=figname, deltaH=20)
    norm = fig.mpl.colors.Normalize(*scale)
    ###colorbar    
    cb = fig.mpl.colorbar.ColorbarBase(axc, cmap=fig.cmap, norm=norm, orientation='vertical')
    cb.set_label('DEP', size=10)
    cb.set_ticks(range(0, scale[1]+1, 10))    
    ###Plot Data
    segments = []
    dTEC = []
    cnt = 0
    for rec in bub:
        if rec.prn not in prns:
            continue
        if rec.code not in codes:
            continue
        if dt1>rec.epoch or rec.epoch>dt2:
            continue        
        lat1 = rec['lat'][0]
        lat2 = rec['lat'][-1]
        lon1 = rec['lon'][0]
        lon2 = rec['lon'][-1]
        if gps.GTODDBL(lat1, lon1, lat2, lon2, 350)>650:            
            continue
        dtec = abs(min(rec['dTEC']))
        segments.append(((lon1, lat1), (lon2, lat2)))
        dTEC.append(dtec)
        cnt += 1 
    
    LC = fig.LC(segments, cmap=fig.cmap, norm=norm)
    LC.set_linewidth(3)
    
    LC.set_array(np.array(dTEC))
    axm.add_collection(LC)    
    axm.set_title('%s     %s - %s (UT)    #bubbles=%s' % (dt1.strftime('%Y/%m/%d'), 
                                                          dt1.strftime('%H:%M'), 
                                                          dt2.strftime('%H:%M'), 
                                                          cnt),
                  size=11)
    
    if len(codes)==len(bub.codes):
        chcodes = 'all'
    else:
        chcodes = ','.join(codes)
    if len(prns)==len(bub.prns):
        chprns = 'all'
    else:
        chprns = ','.join(prns)
    fig.text(0.82, 0.92,
             "Stations:%s\nPRN's:%s" % (chcodes, chprns),
             size=8, ha='left')
    fig.show(watermark=False)
    
def plot_data_vars(gdo, X, Y, figname=None, localtime=False, prns=None, min_ele=30, 
                  colormap=None, marks=False, legend=False, s4legend=False, 
                  xscale='auto', yscale='auto', figsize='auto', pkwargs={}):
    '''
    '''
    if isinstance(Y, basestring):
        Ys = [Y]
    elif hasattr(Y, '__iter__'):
        Ys = list(Y)
    else:
        raise TypeError('Invalid type:%s for Y' % type(Y))
    
    for var in Ys+[X]:
        if var not in gdo.vars:
            raise RuntimeError('Variable: %s does not exists, choose from %s' % (var, list(gdo.vars)))
    
    if colormap and colormap not in gdo.vars:
        raise RuntimeError('Variable: %s does not exists, choose from %s' % (colormap, gdo.vars))    
    
    if localtime:
        timezone = timedelta(gdo.tzoffset/24.)
    
    if prns:
        prns = [x for x in gdo.prns if x in gdo.arcs and x not in gdo.bad_prn and x in prns]
    else:
        prns = [x for x in gdo.prns if x in gdo.arcs and x not in gdo.bad_prn]
    prns.sort()
    
    if figsize=='auto':
        figsize = (9, 6) if legend else (8, 6)
    
    fig = MyFigure(figname=figname, figsize=figsize)
    ax = fig.figure.add_subplot(1,1,1)
    cmap = fig.CM.get_cmap("jet", len(Ys)*len(prns))
    
    if colormap:
        maxcolor = gdo.get_max(colormap, values=('ele', min_ele))
        mincolor = gdo.get_min(colormap, values=('ele', min_ele))
        norm = fig.mpl.colors.Normalize(mincolor, maxcolor)
    
    maxisY = []
    minisY = []
    L = {}
    for Y in Ys:        
        for prn in prns:
            if colormap is None:
                x, y = gdo.get_array(prn, (X, Y), values=('ele', min_ele))
                if len(y)>0 and not np.all(np.isnan(y)):
                    maxisY.append(np.nanmax(y))
                    minisY.append(np.nanmin(y))
                    if marks: 
                        ax.text(np.nanmin(x), y[np.nonzero(x==np.nanmin(x))[0][0]], prn, size=8)
                    P = ax.plot(x, y, **pkwargs)
                    if s4legend and np.nanmax(y)>0.35:
                        L['%s' % prn] = P[0]
                    elif legend:
                        L['%s(%s)' % (Y, prn)] = P[0]
            else:
                x, y, c = gdo.get_array(prn, (X, Y, colormap), values=('ele', min_ele))
                if len(y)>0 and not np.all(np.isnan(y)):
                    maxisY.append(np.nanmax(y))
                    minisY.append(np.nanmin(y))
                    if marks:
                        ax.text(np.nanmin(x), y[np.nonzero(x==np.nanmin(x))[0][0]], prn, size=8)
                    P = ax.scatter(x, y, c=c, s=2, edgecolor='none', cmap=fig.cmap, norm=norm)
            
    if colormap:
        legend = False
        cbar = fig.figure.colorbar(P, orientation='vertical', fraction=0.05,
                                   shrink=0.5, cmap=fig.cmap, norm=norm)
        cbar.set_label(PLOT_LABELS[colormap], size=10)
        fig.adjust(right=0.92)
        for label in cbar.ax.get_yticklabels():
            label.set_fontsize(8)

    dt1 = gdo[-1].epoch
    dt0 = gdo[0].epoch
   
    if xscale=='auto' and X=='epoch':
        if dt1-dt0<20*60*60 or dt1-dt0>24*60*60:
            dt0 = dt1 - timedelta(1)
        ax.set_xlim(dt0.doys, dt1.doys)
    elif xscale!='auto' and X=='epoch':
        dt0 = dt0.replace(hour=xscale[0], minute=0, second=0)
        dt1 = dt1.replace(hour=xscale[1], minute=0, second=0)
        ax.set_xlim(dt0.doys, dt1.doys)
        print dt0.doys, dt1.doys
    elif xscale!='auto':
        ax.set_xlim(*xscale)
        
    if yscale=='auto':
        if Y=='s4':
            ax.set_ylim(0, 1)
            ax.yaxis.grid(True)
        elif Y=='eqTEC':
            ax.set_ylim(0, (int(max(maxisY))+10)/10*10)
        else:
            ax.set_ylim((int(min(minisY))-5)/5*5, (int(max(maxisY))+10)/10*10)
    else:
        ax.set_ylim(*yscale)
    
    if X=='epoch':
        dt0 = GPSDateTime(dt0.year, dt0.month, dt0.day, dt0.hour)
        dt1 = GPSDateTime(dt1.year, dt1.month, dt1.day, dt1.hour)
        if xscale=='auto':
            tk_list = [dt0+timedelta(x/24.) for x in range(24)]
        else:
            tk_list = [dt0+timedelta(x/24.) for x in range(xscale[1]-xscale[0])]
        ax.set_xticks([dt.doys for dt in tk_list])
    
    xlabel = PLOT_UNITS[X]
    if X=='epoch' and localtime:
        ax.set_xticklabels([`(dt+timezone).hour` for dt in tk_list])
        xlabel += ' [LT]'
    elif X=='epoch':
        ax.set_xticklabels([`dt.hour` for dt in tk_list])
        xlabel += ' [UT]'
    
    ax.set_xlabel(xlabel)
    ylabel = set()
    for var in Ys:
        ylabel.add(PLOT_UNITS[var])
    ax.set_ylabel('/'.join(ylabel))

    if legend or s4legend:
        if len(L)>0:
            fig.figure.legend(L.values(), L.keys(), prop={'size':10}, loc=7, 
                              numpoints=3)
        fig.adjust(right=0.84)
    
    if 'fullname' in gdo.station:
        ax.set_title(gdo.station['fullname'], size=12)
    else:
        ax.set_title('Station: %s' % gdo.station['code'].upper(), size=12)
    
    if gdo[0].epoch.date()==gdo[-1].epoch.date():
        fig.figure.text(0.05, 0.95, '%s' % dt1.strftime('%Y/%m/%d'),
                        size=10, ha='left')
    else:
        fig.figure.text(0.05, 0.95, '%s - %s' % (dt0.strftime('%Y/%m/%d'),
                        dt1.strftime('%Y/%m/%d')), size=10, ha='left')
    
    info = ''
    ypos = 0.95
    if 'eqTEC' in Ys and gdo.rec_bias!=0:
        info += 'rx_bias = %.1f\n' % gdo.rec_bias
        ypos = 0.93
    
    info += 'elevation > %s' % min_ele
    fig.figure.text(0.95, ypos, info, size=10, ha='right')
    return fig.show()

def plot_rnx(gdo, figname=None, localtime=False):
    '''
    '''
    if localtime:
        timezone = timedelta(gdo.tzoffset/24.)
    prns = [x for x in gdo.prns if x in gdo.arcs]
    prns.sort()
    fig = MyFigure(figname=figname, figsize=(7,7))
    ax = fig.figure.add_subplot(1,1,1)

    dt1 = gdo[-1].epoch
    dt0 = dt1 - timedelta(1)
    data = {}
    all = []
    for rec in gdo:
        time = rec.epoch.doys        
        for prn in prns:
            if prn not in data:
                data[prn] = []
            if prn in rec:
                data[prn].append(time)
        all.append(time)

    prns = data.keys()
    prns.sort()
    i = 2
    for prn in prns:
        ax.plot(data[prn], [i for x in data[prn]], '|', ms=4, color='blue')
        i+=1

    ax.plot(all, [1 for x in all], '|', ms=3, color='green')
    ax.grid(True)

    ax.set_xlim(dt0.doys, dt1.doys)
    ax.set_ylim(0, len(prns))
    ax.set_yticks(range(i+1))
    ax.set_yticklabels(['', 'All']+prns)
    dt0 = GPSDateTime(dt0.year, dt0.month, dt0.day, dt0.hour)
    dt1 = GPSDateTime(dt1.year, dt1.month, dt1.day, dt1.hour)
    tk_list = [dt0+timedelta(x/24.) for x in range(1,25)]
    ax.set_xticks([dt.doys for dt in tk_list])
    if localtime:
        ax.set_xticklabels([`(dt+timezone).hour` for dt in tk_list])
        ax.set_xlabel('Time [LT]')
    else:
        ax.set_xticklabels([`dt.hour` for dt in tk_list])
        ax.set_xlabel('Time [UT] - %s' % gdo.date.strftime('%Y/%m/%d'))
    ax.set_ylabel('PRN')
    ax.set_title('GPS Data', size=12)    
    fig.adjust()
    return fig.show()

def plot_mag(mdo, figname=None, station=None, figure=None, scale=None):
    '''
    '''
    if not figure:
        fig = MyFigure(figname=figname, figsize=(6,3))
    else:
        fig = figure
    
    H = [x-mdo[0]['H'] for x in mdo.iterlist(['H'])]
    Z = [x-mdo[0]['Z'] for x in mdo.iterlist(['Z'])]
    D = [(x-mdo[0]['D'])*180/3.1415 for x in mdo.iterlist('D')]
    T = [x.hours for x in mdo.iterlist(['epoch'])]
    
    if len(H)>0:
        max_value = max(abs(max(H)),abs(min(H)))
    else:
        max_value = 0
    
    if not scale: 
        scale = int(max_value+20)/10*10
        if scale<200:
            scale = 200
    else: 
        scale = int(scale)

    [y1max, y1min, y2max, y2min] = [scale, -scale, scale/10., -scale/10.]
    
    ax  = fig.add_axes(shared=True,
                       xlim=(0,24), ylim=(y1min,y1max), y2lim=(y2min,y2max),                       
                       xlabel='Hour [UTC]',
                       ylabel='H & Z [nT]',
                       y2label='D [min]',
                       xticks=range(0,24,2),
                       yticks=[y1min+(y1max-y1min)*x/4 for x in range(5)],
                       y2ticks=[y2min+(y2max-y2min)*x/4 for x in range(5)],
                       )
    msg = '\nNo Data Available'
    if len(mdo)>0:
        for arc in mdo.arcs:
            i,j = arc
            ax.plot(T[i:j], H[i:j], 'b')
            ax.plot(T[i:j], Z[i:j], 'g')
            ax.plot(T[i:j], D[i:j], 'r', shared=True)    
        ax.ax.text(23.5, y1max-0.2*y1max, mdo.date.strftime('%d/%m/%Y'), size=8, ha='right', va='top')
        msg = ''
    
    fig.set_xticks(size=8)
    fig.set_yticks(size=8)
    fig.set_labels()
    fig.suptitle('H---', x=0.4, size=10, color = 'b', ha='center', va='bottom')
    fig.suptitle('Z---', x=0.5, size=10, color = 'g', ha='center', va='bottom')
    fig.suptitle('D---', x=0.6, size=10, color = 'r', ha='center', va='bottom')
    
    if station:
        ax.ax.text(0.5, y1max-0.2*y1max, station['fullname']+msg, size=8, ha='left', va='top')
    ax.ax.grid(color='0.8')
    fig.adjust(right=0.85)

    if not figure:
        return fig.show()
    else:
        return ax

def plot_s4(sdo, plot_type='0', figname=None, prns=None, 
            min_ele=30, localtime=False):
    '''
    '''
    if 'fullname' in sdo.station:
        title = sdo.station['fullname']
    else:
        title = sdo.station['code']
    if not prns:
        prns = [prn for prn in sdo.prns]
        prns.sort()

    if localtime and 'longitude' in sdo.station:
        timezone = timedelta((sdo.station['longitude']/15)/24.)
    else:
        localtime = False

    dt1 = sdo[-1].epoch
    dt0 = sdo[0].epoch
   
    if dt1-dt0<20*60*60 or dt1-dt0>24*60*60:
        dt0 = dt1 - timedelta(1) 
        
    tk_list = [dt0.replace(minute=0, second=0)+timedelta(x/24.) for x in range(0,25)]
    
    if localtime:
        if (dt1+timezone).doy==(dt0+timezone).doy:
            date = (dt0+timezone).strftime('%Y/%m/%d')
        else:
            date = '%s-%s' % ((dt0+timezone).strftime('%Y/%m/%d'), (dt1+timezone).strftime('%Y/%m/%d'))
        xticklabels = [`(dt+timezone).hour` for dt in tk_list]
        xlabel= 'Time [LT]'
    else:
        if dt1.doy==dt0.doy:
            date = dt0.strftime('%Y/%m/%d')
        else:
            date = '%s-%s' % (dt0.strftime('%Y/%m/%d'), dt1.strftime('%Y/%m/%d'))
        xticklabels = [`dt.hour` for dt in tk_list]
        xlabel = 'Time [UT]'

    if plot_type == '0':
        plots = {}      
        fig  = MyFigure(figname=figname, figsize=(6, 3))
        cmap = fig.CM.get_cmap("jet", 32)
        ax  = fig.add_axes(xlim=((sdo[-1].epoch-timedelta(1)).doys, sdo[-1].epoch.doys), 
                           ylim=(0,1),
                           xlabel=xlabel,
                           ylabel='Scintillation index S4',
                           xticks=[dt.doys for dt in tk_list],
                           xticklabels = xticklabels,
                           yticks=[0.2, 0.4, 0.6, 0.8, 1.0],
                           )
        for prn in prns:
            if int(prn[1:])>32: continue              
            x, y = sdo.get_array(prn, ('epoch', 's4'), values=('ele', min_ele))
            mask = np.isfinite(y)
            p = ax.plot(x[mask], y[mask], '.', ms=2, color=cmap(int(prn[1:])-1), label='PRN %s' % prn)
            if np.nanmax(y)>0.5 and prn not in plots:
                plots[prn] = p[0]
        
        ax.ax.grid(True)    
        ax.ax.xaxis.grid(False)
        ax.ax.text(tk_list[-1].doys,0.85, r'$S4\ values\ for\ elevation\ angles\ >\ %s\degree$' % min_ele, 
                   ha='right', size=8, color='0.5')
        fig.set_xticks(size=8)
        fig.set_yticks(size=8)
        fig.set_labels()
        fig.adjust(bottom=0.2,left=0.1, right=0.84)
        ax.ax.set_title('%s' % title, size=10)
        fig.figure.text(0.05, 0.9, date, size=8, ha='left')
        if plots:
            labels = ['PRN %s' % s for s in plots.keys()]
            leg=fig.figure.legend(plots.values(), labels, loc=7, prop={'size':8}, borderpad=0.5)
            leg.set_title("S4>0.5")

    elif plot_type == '1':
        ###create fig and axes
        fig = MyFigure(figname=figname, figsize=(5, 5.5), facecolor='white')
        ax  = fig.add_axes(xlim=(-100,100), ylim=(-100,100))
        ax.ax.axison = False
        ###plot custom grid
        ax.plot([-100, 100], [0, 0], color='0.6', lw=0.5)
        ax.plot([0, 0], [-100, 100], color='0.6', lw=0.5)
        for loc in [30,60,0]:
            th = np.linspace(0, 2*utl.PI, 100)
            xr = 100*np.cos(loc*utl.deg2rad)*np.cos(th)
            yr = 100*np.cos(loc*utl.deg2rad)*np.sin(th)
            ax.plot(xr,yr,color='0.6', lw=0.5)
        ###plot coloured line for time arc for prn in data
        for prn in prns:
            if prn>32: continue
            for arc in sdo.arcs[prn]:
                s4_array = np.array([rec[prn]['s4'] \
                                    for rec in sdo[arc[0]:arc[1]]]+[0,1])
                el_array = np.array([100*np.cos(rec[prn]['ele']*utl.deg2rad) \
                                    for rec in sdo[arc[0]:arc[1]]])
                az_array = np.array([utl.PI/2-rec[prn]['azi']*utl.deg2rad \
                                    for rec in sdo[arc[0]:arc[1]]])
                X = el_array*np.cos(az_array)
                Y = el_array*np.sin(az_array)
                points = zip(X,Y)
                segments = zip(points[:-1], points[1:])
                LC = fig.LC(segments)
                LC.set_linewidth(2)
                LC.set_array(s4_array)
                ax.ax.add_collection(LC)
                ax.ax.text(X[-1],Y[-1], prn, ha='center', va='center',
                           fontsize=9, color='0.5')
        ###add colorbar and labels to figure
        cbar = fig.figure.colorbar(LC, orientation='horizontal', ax=ax.ax,
                                   shrink=0.8, pad=0.08, fraction=0.05)
        cbar.set_label('GPS S4', size=10)
        ticklabels = cbar.ax.get_xticklabels()
        fig.py.setp(ticklabels, size=10)
        ax.set_labels()
        ax.ax.text(102, 0,'E', va='center')
        ax.ax.text(-102,0,'W', va='center', ha='right')
        ax.ax.text(0, 102,'N', ha='center')
        ax.ax.text(0,-102,'S', ha='center', va='top')
        fig.suptitle('GPS S4 in Skymap at %s\n%s ' % (sdo.station['fullname'],
                                                      sdo.date))

    return fig.show()