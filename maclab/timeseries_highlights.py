"""
Time series highlights

Module to help process time series arrays and, in particular, to help
visualise areas of change in time series data. 

Routines available:
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
from . import spatial as sp
from functools import reduce 
from operator import mul 

def nearest_rc(array, value):
    """Return the row, col in an array closest to a given value."""
    A = np.ma.abs(array - value)
    return(np.unravel_index(A.argmin(), A.shape))

def nels(arr):
    """Returns number of elements in an ndim array."""
    if isinstance(arr, np.ndarray):
        return(reduce(mul, arr.shape))
    elif isinstance(arrays, np.ma.core.MaskedArray):
        return(arrays.count())
    else:
        print('Object passed was not recognised')
        return
        

class TimeSeries(object):
    """
    Methods for processing and visualising a time series of rasters
    
    Parameters
    ------
    arrays: list of filepaths, np.array, or np.ma.array, required
    affine: Affine object, optional
        supply if arrays is an array in memory. Defaults to globe.
    nodata: int, float, optional
        supply if arrays is an array in memory. Default: -9999
    deltas: arrays, optional
        difference arrays between timesteps.
    breaks: list, optional
        break points (array classification)
    """
    def __init__(self, arrays, affine = None, nodata = None, 
                 deltas = None, breaks = None):
    
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
                self.affine = sp.default_transform(arrays)
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
                self.affine = sp.default_transform(arrays)
            else: 
                self.affine = affine
            if nodata is None:
                ## use fill value
                warnings.warn('Using masked fill value as nodata')
                self.nodata = self.array.fill_value
            self.deltas = None
            self.breaks = None
            
        elif arrays.__class__.__name__:
            if not isinstance(arrays.array, np.ma.core.MaskedArray):
                self.array = np.ma.array(self.array, mask = \
                                         (self.array == self.nodata))
            else:
                self.array = arrays.array
            self.src = None
            self.affine = arrays.affine
            self.nodata = arrays.nodata
            self.deltas = None
            self.breaks = None
            
    def from_files(self):
        """ Read rasters from files into a 3D masked array """
        if self.src is not None:
            this = np.dstack([Raster(i) for i in self.src])
            self.array = np.ma.array(this.array, mask = (this.array == this.nodata))
            self.affine = this.affine
            self.nodata = this.nodata
            return(self)
    
    @property
    def ts_range(self):
        if self.array is not None:
            return(np.min(self.array), np.max(self.array))
    
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
       
    def calc_breaks(self, K, sample = 100000, include_bounds = False, 
                    force_bounds = None):
        """
        Calculate jenks breaks
        
        Params
        ------
        self: TimeSeries object
        K: int, required
            number of breaks to calculate
        sample: int, default 100000
            breaks take ~ 30s to calculate on an array of shape (100000, ).
            sample arg will randomly subsample the arrays by specified number
            of elements. Default 100000. If array contains < sample elements, 
            all elements will be considered. When a subsample is performed it
            includes the max and min values.
        include_bounds: boolean, default False
            if True, includes the min and max of array to be inlcuded as break 
            points
        force_bounds: tuple, optional
            if provided, use these as the upper and lower break points
        
        TODO
        ----
        Should allow the specification of min and max bounds as well as both
        """
        # deal with sampling
        n_els = self.array[:, :, 0].count()
        if n_els < sample:
            samp = np.arange(0, n_els, 1)
        else:
            pos = np.random.choice(n_els, size = sample, replace = False)
            rs, cs = np.take((~self.array[:, :, 0].mask).nonzero(), pos, axis=1)
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
        
        # calc breaks
        jbreaks = []
        # iterate over deltas
        for a in range(self.deltas.shape[-1]):
            jbreak_a = jenks_breaks(self.deltas[rs, cs, a].flatten(), K)
            # find elements where these breaks are in orig array and get vals
            n_breaks = []
            for n in jbreak_a:
                nr, nc = nearest_rc(self.deltas[:, :, a], n)
                n_breaks.append(self.array[nr, nc, a + 1])
            if include_bounds:
                mn = np.min(self.array[:, :, a + 1])
                mx = np.max(self.array[:, :, a + 1])
                n_breaks = [mn] + n_breaks + [mx]
            if force_bounds is not None:
                n_breaks = [force_bounds[0]] + n_breaks + \
                    [force_bounds[1]]
            jbreaks.append(n_breaks)
        # double breaks 0 for array 0
        jbreaks = [jbreaks[0]] + jbreaks
        self.breaks = jbreaks
        return(self)
    
    def colour(self, cols, N):
        """
        Classify raster using jenks breaks and apply colour map
        
        Parameters
        ----------
        cols: list, required
            (hexcodes only?) colours or matplotlib cmap (?)
        N: int, required 
            length of requested colour steps
        
        Returns
        -------
        list: 4D masked arrays of length time steps (RGBA, uint8)
        """
        rgb = []
        for i in range(len(self.breaks)):
            cm, norm = cl.pal(cols).seq(N, self.breaks[i])
            rgb.append(cm(norm(self.array[:, :, i]), bytes = True))
        return(rgb)
    
def highlights(obj, K, N, cols, sample = 100000):
    """ 
    Parameters:
    -----------
    obj: numpy array, masked array, or filepath list, required
    K: int, required
        number of breaks to calculate
    N: int, required
        length of requested colour steps. Should probably match K
    cols: list, required
        (hexcodes only?) colours or matplotlib cmap (?)
    sample: int, default: 100000
        number of subsampled elements to randomly select in calculating
        break points. 
    """
    x = TimeSeries(obj)
    x.calc_deltas()
    x.calc_breaks(K, sample)
    return(x.colour(cols, N))