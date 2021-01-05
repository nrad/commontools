"""
Module with 'useful' and general (non-PyROOT) python functions and tricks


"""

import uuid
import hashlib
import os 
import pickle


class PickleDB():
    """
        little class to pickle and update dictionaries
        results = PickleDB('results.pkl')
        results.update(**{'param':123, 'histo':ROOT.TH1D()})
        

    """
    def __init__(self, pkl_path): 
        self.path     = os.path.expanduser(pkl_path)
        self.dir      = os.path.dirname(self.path)
        self.filename = os.path.basename(self.path)

    @property
    def db(self):
        return self.read(create=True)

    def write(self, path=None, db=None):
        path = path if path else self.path
        db   = db if not db == None  else getattr(self,'db',None)
        if db == None:
            print ("Database is missing... was asked to write %s"%db)
            return
        with open(path,'wb') as f:
            pickle.dump(db,f)
        

    def exists(self):
        isfile = os.path.isfile(self.path)
        return isfile 

    def read(self, create=False):
        exists = self.exists()
        if not exists:
            if create:
                if not os.path.isdir( self.dir ):
                    print("making dir {self.dir}")
                    os.makedirs(self.dir)
                print(f"Creating File: {self.path}")
                self.write( path=self.path, db={})
            else:
                print(f"Expected file at {self.path} but didnt find it")
                return False
        with open(self.path,'rb') as f:
            db  = pickle.load(f)
        return db

    def update(self, **kwargs):
        db = self.db
        db.update(**kwargs)
        self.write(db=db)
        return self.db 




#def dict_func ( func , d):
#    """
#    creates a new dictionary with the same structure and depth as the input dictionary
#    but the final values are determined by func(val)
#    """
#    new_dict = {}
#    for k in d.keys():
#        v = d.get(k)
#        if type(v)==dict:
#            ret = dict_func(func, v)
#        else:
#            ret = func(v)
#        new_dict[k] = ret
#    return new_dict

def dict_func ( func , d):
    """
    creates a new dictionary with the same structure and depth as the input dictionary
    but the final values are determined by func(val)
    depth is the maximum depth we will dive into
    """
    new_dict = {}
    for k in d.keys():
        v = d.get(k)
        #if type(v)==dict:
        #    ret = dict_func(func, v)
        if hasattr(v, 'items') and (depth!=0):
            ret = dict_func(func, v, depth=depth-1)
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


def runFuncInParalNative( func, args , nProc = 15 , starmap=False):
    import  multiprocessing
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
            results = list(map(func,args))
    return results

def runFuncInParal( func, args , nProc = 15 , starmap=False, native=False):
    """
        uses pathos.multiprocessing instead of native multiprocessing one
        to avoid problems with pickling functions
        
        if it fails, the pool might need to get restarted again (pool.terminate(), pool.restart() )

    """

    try:
        import pathos
        useNativeMP = False
    except ImportError:
        useNativeMP = True

    if useNativeMP or nProc==1 or native:
        print("unable to import pathos... multiprocessing might crash because of pickling...")
        return runFuncInParalNative( func, args, nProc, starmap )

    if starmap:
        pool   = pathos.helpers.mp.Pool( processes=nProc )
        mp_map = pool.starmap
    else:
        pool = pathos.pools.ProcessPool( nodes=nProc)
        mp_map = pool.map
    results = mp_map( func, args )
    pool.close()
    pool.join()
    pool.terminate()
    if hasattr(pool, 'restart'):
        pool.restart()
    return results 



