"""
Time series highlights

Module to help process time series arrays and, in particular, to help
visualise areas of change in time series data. 

Routines available presently:
    - simple differences between timesteps
    - calculate natural breaks (jenks) between time steps to optimize
      visualisation of regions of change
    - colouring array(s) using natural breaks to highlight areas of change  
"""

import numpy as np
import rasterio
#from jenks import jenks
from jenkspy import jenks_breaks # faster
import warnings
from . import colours as cl
from functools import reduce 
from operator import mul 

def nearest_rc(array, value):
    """ return the row, col in an array closest to a given value """
    A = np.ma.abs(array - value)
    return(np.unravel_index(A.argmin(), A.shape))

def nels(arr):
    """ returns number of elements in an ndim array """
    return(reduce(mul, arr.shape))

class TimeSeries(object):
    """
    TimeSeries object
    """
    def __init__(self, arrays, affine = None, nodata = None, 
                 deltas = None, breaks = None):
        """
        Params
        ------
        arrays: list of filepaths, np.array, np.ma.array
        affine: supply if arrays is an array
        nodata: supply if arrays is an array
        deltas: difference arrays between timesteps
        breaks: list of break points (value classification)
        """
        if isinstance(arrays, list):
            ## assume filepaths
            self.src = arrays
            self.arr = None
            self.affine = None
            self.nodata = None
            self.deltas = None
            self.breaks = None
        
        elif isinstance(arrays, np.ndarray):
            self.array = arrays
            self.src = None
            if affine is None:
                warnings.warn('No affine transform supplied - using default')
                self.affine = 'TODO'
            else: 
                self.affine = affine
            if nodata is None:
                warnings.warn('No nodata value supplied - using -9999')
                self.nodata = -9999
            self.array = np.ma.array(self.array, mask = (self.array == self.nodata))
            self.deltas = None
            self.breaks = None
            
        elif isinstance(arrays, np.ma.core.MaskedArray):
            self.array = arrays
            self.src = None
            if affine is None:
                warnings.warn('No affine transform supplied - using default')
                self.affine = 'TODO'
            else: 
                self.affine = affine
            if nodata is None:
                ## use fill value
                warnings.warn('Using masked fill value as nodata')
                self.nodata = self.array.fill_value
            self.deltas = None
            self.breaks = None
            
    def from_files(self):
        """ Read rasters from files into a 3D masked array """
        if self.src is not None:
            this = np.dstack([Raster(i) for i in self.src])
            self.arr = np.ma.array(this.array, mask = (this.array == this.nodata))
            self.affine = this.affine
            self.nodata = this.nodata
            return(self)
    
    def calc_deltas(self):
        """ 
        Calculate timestep delta arrays. 
        
        Notes
        -----
        Assumes 3rd dimension is time.
        """
        temp = np.ma.zeros((self.array.shape))
        temp = temp[:, :, :-1] # rm last slot of time dim
        for a in range(temp.shape[-1]):
            temp[:, :, a] = self.array[:, :, a + 1] - self.array[:, :, a]
        self.deltas = temp
        return(self)
    
    def calc_breaks(self, K, sample = 100000):
        """
        Calculate jenks breaks
        
        Params
        ------
        self: TimeSeries object
        K: number of breaks to calculate
        sample: breaks take ~ 30s to calculate on an array of shape (100000, ).
                sample arg will randomly subsample the arrays by specified number
                of elements. Default 100000. If array contains < sample elements, 
                all elements will be considered. When a subsample is performed it
                includes the max and min values.
        """
        ## deal with sampling
        n_els = self.array[:, :, 0].count()
        if n_els < sample:
            samp = np.arange(0, n_els, 1)
        else:
            pos = np.random.choice(n_els, size = sample, replace = False)
            rs, cs = np.take((~x.array[:, :, 0].mask).nonzero(), pos, axis=1)
            ## make sure min, max are included: first calc for time 0
            rmn, cmx = np.unravel_index(
                self.array[:, :, 0].argmin(), (self.array[:, :, 0].shape))
            rmx, cmx = np.unravel_index(
                self.array[:, :, 0].argmax(), (self.array[:, :, 0].shape))
            rs = np.hstack([rmn, rs, rmx])
            cs = np.hstack([rmx, cs, cmx])
            ## then for each delta
            for a in range(self.deltas.shape[-1]):
                rmn, cmx = np.unravel_index(
                    self.deltas[:, :, a].argmin(), (self.array[:, :, 0].shape))
                rmx, cmx = np.unravel_index(
                    self.deltas[:, :, a].argmax(), (self.array[:, :, 0].shape))
                rs = np.hstack([rmn, rs, rmx])
                cs = np.hstack([rmx, cs, cmx])
        jbreaks = []
        ## do array 0
        jbreaks.append(jenks_breaks(self.array[rs, cs, 0].flatten(), K))
        ## iterate over the rest
        for a in range(self.deltas.shape[-1]):
            jbreak_a = jenks_breaks(self.deltas[rs, cs, a].flatten(), K)
            ## find elements where these breaks are in orig array and get vals
            n_breaks = []
            for n in jbreak_a:
                nr, nc = nearest_rc(self.deltas[:, :, a], n)
                n_breaks.append(self.array[nr, nc, a + 1])
            jbreaks.append(n_breaks)
        self.breaks = jbreaks
        return(self)
    
    def colour(self, cols, N):
        """
        Classify raster using jenks breaks and apply colour map
        
        Params
        ------
        cols: list of (hexcodes only?) colours or matplotlib cmap (?)
        N: is length of requested colour steps
        
        Returns
        -------
        list: 4D masked arrays of length time steps (RGBA, uint8)
        """
        rgb = []
        for i in range(len(self.breaks)):
            cm, norm = cl.pal(col_list).seq(N, self.breaks[i])
            rgb.append(cm(norm(self.array[:, :, i]), bytes = True))
        return(rgb)
    
def highlights(obj, K, N, cols, sample = 100000):
    """ 
    Params:
    obj: one of a numpy array, masked array, or list of filepaths to rasters
    """
    x = TimeSeries(obj)
    x.calc_deltas()
    x.calc_breaks(K, sample)
    return(x.colour(cols, N))