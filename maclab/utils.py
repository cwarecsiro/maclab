import re, operator, sys, os

digits = re.compile(r'(\d+)')
def tokenize(filename):
    return tuple(int(token) if match else token
                 for token, match in
                 ((fragment, digits.search(fragment))
                  for fragment in digits.split(filename)))

def readfile(filepath):
    if os.path.exists(filepath):
        output = subprocess.run(["cat", filepath],stdout=subprocess.PIPE).stdout.decode('utf-8')
        print(output)
    else:
        print('The filepath does not exist')

def which(x, op, val):
    mappings = {'<': operator.lt, '<=': operator.le,
                '==': operator.eq, '>': operator.gt,
               '>=': operator.ge} 
    return([i for i, x in enumerate(x) if mappings[op](x, val) ])

def asnumber(number, dtype = 'float32'):
    formtype = '{0:.%sf}' %dtype[-2:]
    print(formtype.format(number))
    
def get(x):
    """Equiavlent of R get? Get name of variable as string"""
    return([ k for k,v in locals().items() if v is x][0])

def closest_val(l, val):
    """Find closets value to a value in a list. 
    
    Parameters
    ----------
    l: list, required
    val: int, float, required
    
    Reference
    ---------
    https://stackoverflow.com/questions/12141150/from-list-of-integers-get-number-closest-to-a-given-value
    """
    return(min(l, key=lambda l:abs(l-val)))


import os, re, rasterio, numpy as np, subprocess

def listdirs(dir):
    """ http://stackoverflow.com/questions/141291/how-to-list-only-top-level-directories-in-python """
    return [os.path.join(os.path.join(dir, x)) for x in os.listdir(dir) 
        if os.path.isdir(os.path.join(dir, x))]

def listfiles(wd, pattern = None, invert = False, recursive = False):
    files = []
    for file in os.listdir(wd):
        if pattern:
            if invert:
                if re.search(pattern, file):
                    pass
                else:
                    files.append(os.path.join(wd, file))
            else:
                if re.search(pattern, file):
                    files.append(os.path.join(wd, file))
        else:
            files.append(os.path.join(wd, file))
    return(files)


def grep(obj, pattern = None, invert = False):
    """ < doc string to come - add more instances... """
    if pattern is None:
            print('Error - need a pattern to search for')
            return
    if isinstance(obj, list):
        if invert:
            new_obj = []
            for i in obj:
                if re.search(pattern, i):
                    pass
                else:
                    new_obj.append(i)
        else:
            new_obj = []
            for i in obj:
                if re.search(pattern, i):
                    new_obj.append(i)
                    
    if isinstance(obj, str):
        n_match = len(re.findall(pattern, obj))
        first_match = re.search(pattern, obj).span()
        if n_match != 0:
            new_obj = ['found ' + str(n_match) + ' matches ' + 'the first between ', first_match]
        else:
            new_obj = 'found no matches'
        
    #if isinstance(obj, DataFrame)
    #if isinstance(obj, array)
                    
    return(new_obj)

    