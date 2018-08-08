"""Methods for reading, writing, and processing spatial datasets

Largely convennience wrappers around the following excellent packages:
    * `rasterio <https://github.com/mapbox/rasterio>`_; and
    * `rasterstats <https://github.com/perrygeo/python-rasterstats routines>`_. 

Example:
    ::
    from maclab import spatial as sp
    src = sp.Raster('/home/usr/dir/eg.tif')
    src.meta
    src.affine

Todo:
    * transform, crop, and extract methods
    * link to vector geometries / rasterize
"""
import numpy as np
import rasterio, warnings
import rasterio.transform
from rasterio.crs import CRS

class Raster(object):
    """ 
    Read raster files into 2/3D arrays
   
    Params
    ------
    raster: (str, array) 2/3D array-like data source, required 
        Filepath to rasterio (GDAL?) supported raster formats or
        numpy arrays with spatial context (just Affine transform?).
    affine: Affine object
        Maps row/col to coordinate reference system
        required if raster is ndarray
    nodata: (int, float) nodata value, optional (default: -9999)
        Overrides the datasource's internal nodata if specified
    band: raster band number, optional (default: 1)
    
    Methods
    -------
    - read
    - write
    - xy_rc
    - crop 
    - extract (TODO)
    - transform (TODO)
    
    """

    def __init__(self, raster, affine=None, nodata=None, band=None, driver = None, crs = None):
        self.array = None
        self.src = None
        self.driver = None
        self.crs = None

        if isinstance(raster, np.ndarray):
            if affine is None:
                # supply default affine transform:
                # (x res, row rotation, xmin, column rotation, y res, ymax)
                aff = default_transform(raster)
                warnings.warn('No affine transform supplied - using default')
            self.array = raster
            self.affine = affine
            self.shape = raster.shape
            if nodata is None:
                nodata = -9999
                warnings.warn('No no data value supplied - using -9999')
            self.nodata = nodata
            if len(raster.shape) > 2:
                self.nbands = raster.shape[-1] # assumes 3rd dim is bands
            else:
                self.nbands = 1
            self.band = band
            self.dtype = raster.dtype
                        
        else:
            self.src = rasterio.open(raster, 'r')
            self.affine = guard_transform(self.src.transform)
            self.shape = (self.src.height, self.src.width)
            if nodata is not None:
                # override with specified nodata
                self.nodata = float(nodata)
            else:
                self.nodata = self.src.nodata
            self.nbands = self.src.count
            self.band = band
            
    @property
    def meta(self, compression = None):
        if self.src is not None:
            if hasattr(self.src, 'profile'):
                m = self.src.profile
        else:
            # update to gen from self.__dict__ 
            m = {}
            m['transform'] = self.affine
            m['dtype'] = self.dtype
            if self.driver is None:
                m['driver'] = 'GTiff'
            else:
                m['driver'] = driver
            if self.crs is None:
                m['crs'] = CRS(init='epsg:4326')
            else:
                m['crs'] = self.crs
            m['height'] = self.shape[0]
            m['width'] = self.shape[1]
            m['count'] = self.nbands
            if comression is not None:
                m['compression'] = None
        return(m)


    def xy_rc(self, x, y):
        """ Given (x, y) in crs, return the (row, column) on the raster"""
        
        # TODO: Generalise to work on arrays so the operation is vectorised
        col, row = [math.floor(a) for a in (~self.affine * (x, y))]
        return row, col

    def read(self, extent=None, window=None, masked=False):
        """ 
        Reads raster into array
        
        Params
        ----------
        bounds: (tuple) bounding box
            in w, s, e, n order, iterable, optional
        window: (tuple) row, col window, optional
        masked: boolean
            return a masked numpy array, default: False
            
        Returns
        -------
        (array, Affine object) Raster object with updated affine and array
        """
        # bounds and window opts
        if extent and window:
            raise ValueError("Specify either bounds or window")

        if extent:
            win = bounds_window(extent, self.affine)
        if window:
            win = window
        
        if extent or window:
            c, _, _, f = window_bounds(win, self.affine)  # c ~ west, f ~ north
            a, b, _, d, e, _, _, _, _ = tuple(self.affine)
            aff = Affine(a, b, c, d, e, f)

        elif self.src:
            # It's an open rasterio dataset
            if window:
                if self.band:
                    new_array = self.src.read(
                        self.band, window=win, masked = masked)
                    aff = self.src.transform
                elif self.src.count == 1:
                    new_array = self.src.read(
                        1, window=win, masked = masked)
                    aff = self.src.transform
                else:
                    bands = []
                    for b in range(self.src.count):
                        bands.append(self.src.read(
                            b + 1, window = win, masked = masked))
                    new_array = np.dstack(bands)
                    aff = self.src.transform
            else:
                if self.band:
                    new_array = self.src.read(
                        self.band, masked = masked)
                    aff = self.src.transform
                elif self.src.count == 1:
                    new_array = self.src.read(masked = masked)
                    aff = self.src.transform
                else:
                    bands = []
                    for b in range(self.src.count):
                        bands.append(self.src.read(
                            b + 1, masked = masked))
                    new_array = np.dstack(bands)
                    aff = self.src.transform

        return Raster(new_array, aff, self.nodata)
    
    def write(self, filepath, tags = None):
        """Generalise for ESRI flt grids etc"""
        with rasterio.open(filepath, 'w', **self.meta) as output:           
            if tags is not None:
                d = {}
                for key, value in tags.items():
                    d[key] = value
                output.tags = d
            if len(self.array.shape) > 2 and self.nbands > 1:
                for b in range(self.nbands):
                    output.write(self.array[:, :, b + 1].\
                                astype(self.meta['dtype']), b + 1)
            else:
                output.write(self.array.astype(self.meta['dtype']), 1)
                                
    
    def crop(self, mask = None, extent = None):
        """Crop a raster to the extent of another raster or vector.
        
        Parameters
        ----------
        mask: str, Raster, Fiona obj, optional
            One of a filepath to a spatial dataset (raster and shp for now...), 
            existing Raster class object, or Fiona object
        extent: tuple, optional
            The extent to crop to (xmin, xmax, ymin, ymax). 
            One of mask or extent must be supplied
        
        Returns
        -------
        Raster object
        """
        
        if not mask or extent:
            print('A mask or extent must be supplied')
            return
        
        if mask is None:
            if isinstance(mask, str):
                try:
                    # replace with Vector.read() at some point
                    mask_src = fiona.open(mask)
                except 'IOError':
                    try:
                        mask_src = rasterio.open(mask)
                    except 'IOError':
                        print('Do not recognise the mask data source', \
                               'as a valid vector or raster type')
                        return

        """
        TODO: block to deal with mask files in memory
        """
        if extent is None: 
            ext = mask_src.bounds
        else:
            ext = extent
        llx, lly = ext[0], ext[1]
        col, row = ~aff * (llx, lly)
        rows = round(row) - mask.shape[0] 
        cols = round(col) + mask.shape[1]
        row, col = round(row), round(col)

        fitted = self.array[rows:row, col:cols]
        
        # handle bounds that go > 180 | < -180
        if cols > self.shape[1]:
            # probably > 180E: take chunk from western side
            westside = self.array[rows:row, 0:cols-self.shape[1]]
            fitted = np.hstack([fitted, westside])

        if col < 0:
            # probably < 180W: take chunk from eastern side
            eastside = self.array[rows:row, self.shape[1]-abs(col):self.shape[1]]
            fitted = np.hstack([eastside, fitted])

        return(fitted)
    
    def __enter__(self):
        return self

    def __exit__(self, *args):
        if self.src is not None:
            # close the rasterio reader
            self.src.close()
            
def bounds_window(bounds, affine):
    """Create window to read from file or array."""
    w, s, e, n = bounds
    row_start, col_start = rowcol(w, n, affine)
    row_stop, col_stop = rowcol(e, s, affine, op=math.ceil)
    return (row_start, row_stop), (col_start, col_stop)

def window_bounds(window, affine):
    """Create bounds from window."""
    (row_start, row_stop), (col_start, col_stop) = window
    w, s = (col_start, row_stop) * affine
    e, n = (col_stop, row_start) * affine
    return w, s, e, n

def stack(filepaths):
    """Simple wrapper to stack consistent data sources."""
    objs = [sp.Raster(i) for i in filepaths]
    # cursory checks
    check_aff = all([i.affine for i in objs])
    check_shp = all([i.shape for i in objs])
    check_na = all([i.nodata for i in objs])
    if check_aff and check_na and check_na:
        dat = np.dstack([i.read().array for i in objs])
        return(sp.Raster(dat, objs[0].affine, objs[0].nodata))
    else:
        print('Not all files are consistent')
        
def default_transform(arr):
    """Generate a default global affine from an array."""
    cols, rows = arr.shape[0], arr.shape[1]
    return(rasterio.transform.from_bounds(-180, -90, 180, 90, cols, rows))

def write_meta(filepath):
    """Writes a metadata file for convenient inspection and provenance tracking
    
    Notes
    -----
    Uses tags to embed extra metadata - not sure how this works beyond GTiff.
    """
    f = rasterio.open(filepath)
    meta = f.meta
    tags = f.tags()
    output = open('%s.meta' %filepath, 'w')
    for key, value in meta.items():
        output.write('%s: %s\n' %(key, value))
    for key, value in tags.items():
        output.write('%s: %s\n' %(key, value))
    output.close
    f.close()