import os

def downloadfile(url, dst = None, auth = False, overwrite = False, 
                 verbose = False):
    """Download file using wget.
    
    Parameters
    ----------
    url: str, required
        url of the file to download
    dst: str, optional
        Dir to save file to. Default: tempdir
    auth: tuple, optional
        Provide username and password for authentification
    overwrite: boolean, optional
        If datasets exists, should it be written over?
    verbose: boolean, optional
    
    Returns
    -------
    Filepath to downloaded file
    """
    
    ## config dst
    if dst is not None:
        if not os.path.isdir(dst):
            raise IOError('dst is not a dir')
    else:
        dst = tempfile.gettempdir()
    
    filename = os.path.join(dst, os.path.basename(url))
    if os.path.exists(filename):
        if not overwrite:
            raise IOError('File already exists and overwrite is set to False')
        else:
            os.remove(filename)
    
    ## config call
    dst = '--directory-prefix=' + dst
    if auth:
        us, pw = auth[0], auth[1] 
        call = ' '.join(['wget', '--user', us, '--password', pw, dst, url]) 
    else:
        call = ' '.join(['wget', dst, url])
    
    if verbose:
        print(call)
    
    ## make call
    response = os.system(call)
    
    ## check
    if response == 0:
        if os.path.exists(filename):
            return(filename)
        else:
            raise FileNotFoundError('Response code was 0, but the file does not exist')
    else:
        return(response)