import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import cubehelix
from temperature import *

aurora = cubehelix.cmap(start=0.5,rot=-0.75,minLight=0.05,gamma=1.2)
aurora_r = cubehelix.cmap(start=0.5,rot=-0.75,reverse=True,minLight=0.05,gamma=1.2)

def read_ncdf(file,vars=[]):
    """ Reads a netCDF 3 file and returns a dict with its variables
    
    Parameters
    ----------
    file : string
        The netCDF file to be read
    vars : string list, optional
        List of variables to be read from file

    Returns
    -------
    d : A python dict containing the files variables

    Examples
    --------
    >>> d=read_ncdf('test.cdf')
    >>> d1=read_ncdf('test.cdf',vars=['var1'])
    >>> d2=read_ncdf('test.cdf',vars=['var1','var2'])

    """
    from scipy.io import netcdf

    try:
        f=netcdf.netcdf_file(file,'r',mmap=False)
    except IOError:
        print('Error: Cannot open file'+file)
        return 0

    if vars == []: vars=f.variables.keys()
    d=dict((v,f.variables[v].data) for v in vars)
    f.close()
    return d


from numpy import ma
from matplotlib import cbook
from matplotlib.colors import Normalize

class MidPointNorm(Normalize):    
    def __init__(self, midpoint=0, vmin=None, vmax=None, clip=False):
        Normalize.__init__(self,vmin, vmax, clip)
        self.midpoint = midpoint

    def __call__(self, value, clip=None):
        if clip is None:
            clip = self.clip

        result, is_scalar = self.process_value(value)

        self.autoscale_None(result)
        vmin, vmax, midpoint = self.vmin, self.vmax, self.midpoint

        if not (vmin < midpoint < vmax):
            raise ValueError("midpoint must be between maxvalue and minvalue.")       
        elif vmin == vmax:
            result.fill(0) # Or should it be all masked? Or 0.5?
        elif vmin > vmax:
            raise ValueError("maxvalue must be bigger than minvalue")
        else:
            vmin = float(vmin)
            vmax = float(vmax)
            if clip:
                mask = ma.getmask(result)
                result = ma.array(np.clip(result.filled(vmax), vmin, vmax),
                                  mask=mask)

            # ma division is very slow; we can take a shortcut
            resdat = result.data

            #First scale to -1 to 1 range, than to from 0 to 1.
            resdat -= midpoint            
            resdat[resdat>0] /= abs(vmax - midpoint)            
            resdat[resdat<0] /= abs(vmin - midpoint)

            resdat /= 2.
            resdat += 0.5
            result = ma.array(resdat, mask=result.mask, copy=False)                

        if is_scalar:
            result = result[0]            
        return result

    def inverse(self, value):
        if not self.scaled():
            raise ValueError("Not invertible until scaled")
        vmin, vmax, midpoint = self.vmin, self.vmax, self.midpoint

        if mpl.cbook.iterable(value):
            val = ma.asarray(value)
            val = 2 * (val-0.5)  
            val[val>0]  *= abs(vmax - midpoint)
            val[val<0] *= abs(vmin - midpoint)
            val += midpoint
            return val
        else:
            val = 2 * (val - 0.5)
            if val < 0: 
                return  val*abs(vmin-midpoint) + midpoint
            else:
                return  val*abs(vmax-midpoint) + midpoint

__author__="Paul H, Horea Christian"
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid

def remappedColorMap(cmap, start=0, midpoint=0.5, stop=1.0, 
name='shiftedcmap'):
    '''
    Function to offset the median value of a colormap, and scale the
    remaining color range. Useful for data with a negative minimum and
    positive maximum where you want the middle of the colormap's dynamic
    range to be at zero.

    Input
    -----
      cmap : The matplotlib colormap to be altered
      start : Offset from lowest point in the colormap's range.
          Defaults to 0.0 (no lower ofset). Should be between
          0.0 and 0.5; if your dataset mean is negative you should leave 
          this at 0.0, otherwise to (vmax-abs(vmin))/(2*vmax) 
      midpoint : The new center of the colormap. Defaults to 
          0.5 (no shift). Should be between 0.0 and 1.0; usually the
          optimal value is abs(vmin)/(vmax+abs(vmin)) 
      stop : Offset from highets point in the colormap's range.
          Defaults to 1.0 (no upper ofset). Should be between
          0.5 and 1.0; if your dataset mean is positive you should leave 
          this at 1.0, otherwise to (abs(vmin)-vmax)/(2*abs(vmin)) 
    '''
    cdict = {
        'red': [],
        'green': [],
        'blue': [],
        'alpha': []
    }

    # regular index to compute the colors
    reg_index = np.hstack([
        np.linspace(start, 0.5, 128, endpoint=False), 
        np.linspace(0.5, stop, 129)
    ])

    # shifted index to match the data
    shift_index = np.hstack([
        np.linspace(0.0, midpoint, 128, endpoint=False), 
        np.linspace(midpoint, 1.0, 129)
    ])

    for ri, si in zip(reg_index, shift_index):
        r, g, b, a = cmap(ri)

        cdict['red'].append((si, r, r))
        cdict['green'].append((si, g, g))
        cdict['blue'].append((si, b, b))
        cdict['alpha'].append((si, a, a))

    newcmap = matplotlib.colors.LinearSegmentedColormap(name, cdict)
    plt.register_cmap(cmap=newcmap)

    return newcmap

def make_ep_plot(energy,pitch,x,cmap=temperature,cb_range= None,cb_ticks=None, filename= None,cb_midpoint=None):
    

    if isinstance(cmap,str):
        cm = plt.cm.get_cmap(cmap)
        cm.set_under("w")
    else:
        cm=cmap
        cm.set_under("w")

    if cb_range == None:
        v_min, v_max = np.amin(x), np.amax(x)
    else:
        if any([x == None for x in cb_range]):
            if cb_range[0] == None:
                v_min, v_max = np.amin(x),cb_range[1]
            if cb_range[1] == None:
                v_min, v_max = cb_range[0], np.amax(x)
        else:
            v_min, v_max = cb_range

    fig,ax = plt.subplots()
    fig.set_size_inches((5,3.5))
    
    if cb_midpoint != None:
        #norm = MidPointNorm(midpoint = cb_midpoint,vmin=v_min,vmax=v_max,clip=False)
        #positive_ticks = np.linspace(0,v_max,5)
        #negative_ticks = np.linspace(v_min,0,4,endpoint=False)
        #cb_ticks = np.concatenate((negative_ticks, positive_ticks))
        #cb_ticks = [np.around(x,-int(np.floor(np.log10(np.abs(x))-1))) if x != 0.0 else 0.0 for x in cb_ticks]
        if np.abs(v_max) > np.abs(v_min):
            start = 1.0 - (np.abs(v_max) - v_min)/(2*np.abs(v_max))
            fin = 1.0
        else:
            start = 0.0
            fin = 1.0-(abs(v_min) - v_max)/(2*np.abs(v_min))
        
        cm = remappedColorMap(cm,start=start,stop=fin,midpoint=(np.abs(v_min)/(np.abs(v_max)+np.abs(v_min))))
        
        norm = None
    else:
        norm = None

    p = ax.pcolor(energy,pitch,x,cmap = cm,vmin=v_min,vmax=v_max,lw=0,rasterized=True,norm=norm)
    ax.set_xlim((20.0,100.0))
    ax.set_yticks([-1.0,-0.5,0.0,0.5,1.0])
    ax.set_xticks([20.0,40.0,60.0,80.0,100.0])
    ax.tick_params(labelsize="medium")
    ax.set_xlabel("Energy [keV]",labelpad=-0.05,fontsize="16")
    ax.set_ylabel("Pitch",labelpad=-5,fontsize="16")
    cb = fig.colorbar(p,pad=0)
    cb.set_ticks(cb_ticks)
    cb.ax.tick_params(labelsize="medium")
    cb.ax.yaxis.get_offset_text().set_size("small")
    cb.ax.yaxis.get_offset_text().set_position((-1.2,1))
    cb.solids.set_edgecolor("face")
    fig.tight_layout()

    if filename != None:
        fig.savefig(filename,bbox_inches="tight")
        fig.clear()
        plt.close(fig)

    return

def make_stacked_plot(energy,pitch,x,err,dist,cb_range=None, cmap=aurora,filename=None):

    if isinstance(cmap,str):
        cm = plt.cm.get_cmap(cmap)
        cm.set_under("w")
    else:
        cm=cmap
        cm.set_under("w")

    nr,nc = np.shape(x)
    inds = np.arange(2,nc,2)
    nrow= len(inds)

    if cb_range == None:
        v_min, v_max = np.amin(x), np.amax(x)
    else:
        if any([x == None for x in cb_range]):
            if cb_range[0] == None:
                v_min, v_max = np.amin(x),cb_range[1]
            if cb_range[1] == None:
                v_min, v_max = cb_range[0], np.amax(x)
        else:
            v_min, v_max = cb_range

    fig, ax = plt.subplots(nrows=nrow,sharex=True)
    fig.set_size_inches((5,5))    
    for i,ind in enumerate(inds):
        j = -(i+1)
        ax[j].set_frame_on(False)
        ax[j].get_xaxis().set_ticks_position("none")
        ax[j].fill_between(pitch,x[:,ind] - err[:,ind],x[:,ind] + err[:,ind],facecolor="lightgray",color="lightgray")
        p = ax[j].scatter(pitch,dist[:,ind],c=dist[:,ind],vmin=v_min,vmax=v_max,cmap=cm,marker="o",lw=0.5)
        ax[j].set_xticks([-1.0,-0.5,0.0,0.5,1.0])
        ax[j].set_yticks([])
        ax[j].text(-1.0,dist[1,ind],"{0:d} keV -".format(int(energy[ind])),fontsize="small",va="center",ha="right")

    ax[-1].set_xlabel("Pitch",fontsize="small")
    ax[-1].tick_params(labelsize="small")

    fig.tight_layout(h_pad=0.0)
    cax,aw = mpl.colorbar.make_axes(fig.axes,aspect=40,shrink=0.95,pad=0.0)
    cb = fig.colorbar(p,cax=cax,ticks=[])
    cb.solids.set_edgecolor("face")

    if filename != None:
        fig.savefig(filename,bbox_inches="tight")
        fig.clear()
        plt.close(fig)

    return

