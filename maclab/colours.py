import matplotlib
from matplotlib import pyplot as plt

class pal(object):
    """ 
    Colour stuff
    """
    def __init__(self, cols, name = None):
        """
        Params
        ------
        cols: either a list of (hex codes only?) colours or
              a matplotlib cmap.
        """
        self.cols = cols
        if name is None:
            self.name = 'pal'
        else:
            self.name = name
    
    def seq(self, N, breaks = None):
        cm = matplotlib.colors.LinearSegmentedColormap.from_list(\
               self.name, self.cols, N = N)
        if breaks is not None:
            norm = matplotlib.colors.BoundaryNorm(\
               breaks, len(breaks))
            return(cm, norm)
        else:
            return(cm)    
        