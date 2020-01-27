import os, shutil,re
#RootToolsPath = re.search( ".*/RootTools/" , RootTools.__file__).group()
def copyIndexPHP( directory ):
    ''' Copy index.php to directory
    '''
    index_php = os.path.join( directory, 'index.php' )
    if not os.path.exists( directory ): os.makedirs( directory )
    #shutil.copyfile( os.path.expandvars( '$CMSSW_BASE/src/RootTools/plot/php/index.php' ), index_php )
    indexPhpPath = os.path.expandvars( '$COMMONTOOLS_BASE/RootTools/plot/php/index.php' )
    try:
        shutil.copyfile( indexPhpPath, index_php )
    except:
        print("Failed to copy index.php... looked for it in:\n%s"%indexPhpPath)
    
    #if not directory[-1] == '/': directory = directory+'/'
    #subdirs = directory.split('/')
    #for i in range(1,len(subdirs)):
    #  p = '/'.join(subdirs[:-i])
    #  index_php = os.path.join( p, 'index.php' )
    #  if os.path.exists( index_php ): break
    #  else:
    #    shutil.copyfile( os.path.expandvars( '$CMSSW_BASE/src/RootTools/plot/php/index.php' ), index_php )
