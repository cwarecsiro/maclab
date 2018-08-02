maclab
--------

To install:: 
        
    >>> pinstall git+https://github.com/cwarecsiro/maclab.git

To use, something like::  

    >>> from maclab import timeseries_hightlights as ts
    >>> rgb_arrays = ts.highlights(nd_array, 10, 5, 'YlGn')