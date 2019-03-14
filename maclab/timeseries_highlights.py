"""
Time series highlights

Module to help process time series arrays an dvisualise areas of  
change in time series data. 

Routines available:
    - simple differences between timesteps
    - calculate jenks breaks between time steps to optimize
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
from . import utils as ut
from functools import reduce 
from operator import mul 
import pandas as pd

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
            # assume filepaths
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

    def calc_breaks(self, K, k = 20, background_breaks = 'linear', 
                    sample = 100000, include_bounds = False, 
                    force_bounds = None):
        """
        Calculate jenks breaks
        
        Params
        ------
        self: TimeSeries object
        K: int, required
            number of breaks to calculate on deltas
        k: int, optional
            number of colour breaks to apply as background breaks on input 
            arrays. Default: 20. These are modified by K. 
        background_breaks: str, optional.
            Options are 'linear' (default) or 'jenks.' This arg describes the 
            types of breaks calculated for the input arrays. These are then 
            updated by the breaks calculated for the delta arrays to highlight 
            differences.
        quantile: float, optional.
            Quantile from which to calculate breaks on the delta arrays. 
            By default, these will be calculated on the upper 5% of delta
            values (e.g. quantile = 0.95, representing the greatest difference 
            in values between two time points).
        sample: int, default 100000
            breaks take ~ 30s to calculate on an array of shape (100000, ).
            sample arg will randomly subsample the arrays by specified number
            of elements. Default 100000. If array contains < sample elements, 
            all elements will be considered. When a subsample is performed it
            includes the max and min values.
         force_bounds: tuple, optional
            if provided, use these as the upper and lower break points.
        
        """
        # deal with sampling
        n_els = self.array[:, :, 0].count()
        if n_els < sample:
            samp = np.arange(0, n_els, 1)
            rs, cs = np.ma.where(self.array[:, :, 0] != self.nodata)
        else:
            pos = np.random.choice(n_els, size = sample, replace = False)
            rs, cs = np.take((~self.array[:, :, 0].mask).nonzero(), \
                             pos, axis=1)
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
                    self.deltas[:, :, a].argmin(), \
                    (self.array[:, :, 0].shape))
                rmx, cmx = np.unravel_index(
                    self.deltas[:, :, a].argmax(), \
                            (self.array[:, :, 0].shape))
                rs = np.hstack([rmn, rs, rmx])
                cs = np.hstack([rmx, cs, cmx])
          
        # calc breaks between time points
        delta_breaks = []
        # iterate over deltas
        for a in range(self.deltas.shape[-1]):
            ordered = np.sort(self.deltas[rs, cs, a].flatten())
            # pandas because np.quantile is not in installed numpy version
            df = pd.DataFrame(ordered) 
            ordered = ordered[ordered > df.quantile(quantile)[0]]
            if len(ordered) > K:
                # drop min so len(dbreaks) == K
                dbreaks = jenks_breaks(ordered, K)[1:] 
            else: 
                dbreaks = ordered
            # find elements where these breaks are in data array and get vals
            t0_breaks = []
            t1_breaks = []
            for n in dbreaks:
                nr, nc = nearest_rc(self.deltas[:, :, a], n)
                t0_breaks.append(self.array[nr, nc, a])
                t1_breaks.append(self.array[nr, nc, a + 1])
            delta_breaks.append((t0_breaks, t1_breaks))

        # calc breaks for time points (i.e. the data)
        back_breaks = []
        # loop over the time points-1 as each time point uses breaks 
        # from the previous
        for a in range(self.array.shape[-1]-1):
            if background_breaks == 'linear':
                if force_bounds:
                    back_breaks.append(\
                        np.linspace(force_bounds[0], force_bounds[1], \
                                    k+1).tolist())
                else:
                    mn = np.ma.min(self.array[:, :, a])
                    mx = np.ma.max(self.array[:, :, a])
                    back_breaks.append(np.linspace(mn, mx, k+1).tolist())
                    
            if background_breaks == 'jenks':
                jbreaks = jenks_breaks(self.array[rs, cs, a].flatten(), k)
                if force_bounds is not None:
                    jbreaks = [force_bounds[0]] + jbreaks  + [force_bounds[1]]
                    # remove possible duplicate bounds
                    jbreaks = list(set(jbreaks))
                back_breaks.append(jbreaks)
        # duplicate breaks for final time point
        back_breaks = back_breaks + [back_breaks[-1]]
               
        # insert delta breaks into background breaks
        combined_breaks = [back_breaks[0]] 
        for b in back_breaks[1:]: # index 0 isn't changed
            for d in delta_breaks:
                idx0 = [(ut.which(b, '==', ut.closest_val(b, val))[0]) \
                        for val in d[0]]
                idx1 = [(ut.which(b, '==', ut.closest_val(b, val))[0]) \
                        for val in d[1]]
                time_idx = (idx0, idx1)
            new_breaks = b.copy()
            for idx in range(len(time_idx[0])):
                if time_idx[0][idx] == time_idx[1][idx]:
                    new_breaks[time_idx[1][idx]] = d[1][idx]
            c_breaks = b + new_breaks
            c_breaks.sort()
            combined_breaks.append(c_breaks)
                  
        self.breaks = combined_breaks
        return(self)
    
    def colour(self, cols, N = None):
        """
        Classify raster using jenks breaks and apply colour map
        
        Parameters
        ----------
        cols: list, required
            (hexcodes only?) colours or matplotlib cmap (?)
        N: int, optional 
            length of requested colour steps. Default is len(cols) +1
        
        Returns
        -------
        list: 4D masked arrays of length time steps (RGBA, uint8)
        """
        rgb = []
        for i in range(len(self.breaks)):
            if N is None:
                _N = len(self.breaks[i])-1
            else:
                _N = N
            cm, norm = cl.pal(cols).seq(_N, self.breaks[i])
            rgb.append(cm(norm(self.array[:, :, i]), bytes = True))
        return(rgb)
    
def highlights(obj, K, cols, k = 20, N = None, background_breaks = 'linear', 
               quantile = 0.95, sample = 100000, force_bounds = True):
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
    if force_bounds:
        bounds = x.ts_range
        x.calc_breaks(K, sample, force_bounds = bounds)
    else:
        x.calc_breaks(K, sample)
    return(x.colour(cols, N))