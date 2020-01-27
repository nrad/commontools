"""
Module with 'useful' and general (non-PyROOT) python functions and tricks


"""

import uuid
import hashlib


def dict_func ( func , d):
    """
    creates a new dictionary with the same structure and depth as the input dictionary
    but the final values are determined by func(val)
    """
    new_dict = {}
    for k in d.keys():
        v = d.get(k)
        if type(v)==dict:
            ret = dict_function( v , func)
        else:
            ret = func(v)
        new_dict[k] = ret
    return new_dict



def dict_function ( d,  func ):
    raise Exception("Sorry but use dict_func instead :(")
    

##
## Hashing Tools
##

def uniqueHash():
    #return hashlib.md5("%s"%time.time()).hexdigest()    
    return "tmp" + str( uuid.uuid4() ).replace("-","")

def uniqueName(prefix="tmp"):
    return "%s_%s"%(prefix, uuid.uuid4().hex.replace("-","") )


def hashString( string ):
    import base64
    s = base64.b64encode( hashlib.md5( string ).digest() )
    return s


def hashObj( obj ):
    hsh = hashlib.md5()
    hsh.update( str(obj).encode('utf-8') )
    return hsh.hexdigest()

def getFileChecksum(fname):
    """
    https://stackoverflow.com/questions/3431825/generating-an-md5-checksum-of-a-file
    """
    import hashlib
    hash_md5 = hashlib.md5()
    with open(fname, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest()


###

#def runFuncInParal( func, args , nProc = 15 ):
#    import multiprocessing
#    #nProc=1
#    if nProc >1:
#        pool         =   multiprocessing.Pool( processes = nProc )
#        results      =   pool.map( func , args)
#        pool.close()
#        pool.join()
#    else:
#        results = map(func,args)
#    return results

def runFuncInParal( func, args , nProc = 15 , starmap=False):
    import multiprocessing
    #nProc=1
    if nProc >1:
        pool         =   multiprocessing.Pool( processes = nProc )
        if starmap:
            results      =   pool.starmap( func , args)
        else:
            results      =   pool.map( func , args)
        pool.close()
        pool.join()
    else:
        if starmap:
            from itertools import starmap
            results = starmap(func,args)
        else:
            results = map(func,args)
    return results

