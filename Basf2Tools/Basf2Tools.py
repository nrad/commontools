import basf2



import sys
#sys.path.append("/home/belle2/prazcyri/repositories")
#sys.path.append("/home/belle2/nrad/_old/from_cryll/")
#from my_utils.basf2_helpers import PathDefiner, remove_module

tf2.TestBit(ROOT.TFile.kRecovered)




class PathDefiner:

    def __init__(self, global_tags=[], release=3):

        compatible_releases = [2,3,4,5]
        if not release in compatible_releases:
            raise ValueError("Compatible releases: {}".format(compatible_releases))

        self.global_tags = global_tags
        self.release = release

    def setup_database(self):

        if self.global_tags:

            if self.release in [2,3]:
                basf2.reset_database()
                basf2.use_database_chain()

                for global_tag in self.global_tags.split():
                    if global_tag.startswith("/"):
                        basf2.use_local_database(global_tag)
                    else:
                        basf2.use_central_database(global_tag)

            elif self.release in [4]:
                basf2.conditions.override_globaltags()
                basf2.conditions.expert_settings(usable_globaltag_states={"TESTING", "VALIDATED",
                                                                          "PUBLISHED", "RUNNING", "OPEN"})

                for global_tag in self.global_tags.split():
                    if global_tag.startswith("/"):
                        basf2.conditions.append_testing_payloads(global_tag)
                    else:
                        basf2.conditions.append_globaltag(global_tag)

            elif self.release in [5]:
                #basf2.conditions.override_globaltags()
                #basf2.conditions.expert_settings(usable_globaltag_states={"TESTING", "VALIDATED",
                #                                                          "PUBLISHED", "RUNNING", "OPEN"})
                #basf2.conditions.reset()
                basf2.conditions.expert_settings(usable_globaltag_states={"TESTING", "VALIDATED",
                                                                          "PUBLISHED", "RUNNING", "OPEN"})
                for global_tag in self.global_tags.split():
                    if global_tag.startswith("/"):
                        basf2.conditions.append_testing_payloads(global_tag)
                    else:
                        basf2.conditions.prepend_globaltag(global_tag)

        


    def define_path(self): # pure virtual method
        raise NotImplementedError



class Basf2Path(PathDefiner):
    
    def define_path(self):
        #global_tags = global_tags if isinstance(global_tags, (list,tuple)) else global_tags.split()
        self.setup_database()
        self.path = basf2.create_path()


    def track_reconstruction(self, components=["PXD","SVD","CDC"] ):
        assert False
        import rawdata 
        from tracking import add_hit_preparation_modules, add_track_finding, add_track_fit_and_track_creator

        self.path.add_module("Gearbox")
        self.path.add_module("Geometry", useDB=True)
        rawdata.add_unpackers(self.path)
        add_hit_preparation_modules(self.path)                                                   
        self.path.add_module('SetupGenfitExtrapolation',                                         
                        energyLossBrems=False, noiseBrems=False)                            
        add_track_finding(self.path, prune_temporary_tracks=False, components=components)        
        add_track_fit_and_track_creator(self.path, components=components)

    def add_reconstruction(self, **kwargs):
        import reconstruction
        reconstruction.add_reconstruction(self.path, **kwargs)

    def add_raw_data_reconstruction(self, **kwargs):
        import rawdata 
        import reconstruction
        self.path.add_module("Gearbox")
        self.path.add_module("Geometry", useDB=True)
        rawdata.add_unpackers(self.path)
        reconstruction.add_reconstruction(self.path, **kwargs)

        path = self.path
        basf2.set_module_parameters(path, "PXDPostErrorChecker", CriticalErrorMask=0)
        basf2.set_module_parameters(path, "PXDUnpacker", logLevel=basf2.LogLevel.ERROR, SuppressErrorMask=0xFFFFFFFFFFFFFFFF)
        self.path = remove_module(path, module="TRGCDCT3DUnpacker")


    def etc(self):
            # Patch by Maiko Takahashi - BIIDP-1366
            basf2.set_module_parameters(self.path, "PXDPostErrorChecker", CriticalErrorMask=0)
            # Silence PXDUnpacker warnings - BIIDP-1531
            #basf2.set_module_parameters(self.path, "PXDUnpacker", logLevel=basf2.LogLevel.ERROR, SuppressErrorMask=0xFFFFFFFFFFFFFFFF)
            # Patch by Koga-san + Ritter-san:
            # exclude the TRGCDCT3DUnpackerModule from the path.
            #self.path = remove_module(self.path, module="TRGCDCT3DUnpacker")
            #return path


    def get_track_collection(self, track_selection='p>1.0 and abs(dz)<2.0 and dr<0.5 and nSVDHits>=6', event_selection='nParticlesInList(Upsilon(4S):vertexKFit)==1'):
        modularAnalysis.fillParticleList('mu+:vertexKFit', selection, path=self.path)
        modularAnalysis.reconstructDecay('Upsilon(4S):vertexKFit -> mu+:vertexKFit mu-:vertexKFit', '9.5<M<11.5', path=self.path)    
        #vertex.vertexKFit('Upsilon(4S):vertexKFit', conf_level=0.0, path=path)
        modularAnalysis.applyCuts('Upsilon(4S):vertexKFit', 'nParticlesInList(Upsilon(4S):vertexKFit)==1', path=self.path)


  
     



    def _set_module_params__old(self, module_name, **new_params ):
        path = self.path
        mod_list = path.modules()
        for mod in mod_list:
            if not mod.name() == module_name:
                continue
            print (" found module: ", mod)
            for param_name, param_val in new_params.items():
                mod.param(param_name, param_val)
    
    def set_module_params(self, name=None, type=None, recursive=False, **new_params): 
        return set_module_params( path=self.path, name=name, type=type, recursive=recursive, strict=True, **new_params)

    def find_module(self, name=None, type=None):
        return find_module(self.path, name=name, type=type)

    def add_simulation_to_path(self,  nEvents = 10000, experiment = 0 , run = 1, vertex=None, covVertex=None, mode = None, output_cdst = None, bkg_files=None, bkg_overlay=True, path=None, **kwargs):
        from generators import add_kkmc_generator, add_babayaganlo_generator, add_evtgen_generator
        from simulation import add_simulation
        from reconstruction import add_reconstruction, add_cdst_output
        from beamparameters import add_beamparameters
        # specify number of events to be generated

        path = path if path else self.path

        path.add_module("EventInfoSetter",
                        expList=experiment,
                        runList=run,
                        evtNumList=nEvents)
    
        basf2.set_random_seed(run)
    
        # beam parameters
        #beamparameters = add_beamparameters(path, "Y4S")
        if vertex:
            beamparameters.param("vertex", vertex)
        if covVertex:
            beamparameters.param("covVertex", covVertex)
        if mode=='ee':
            add_babayaganlo_generator(path, finalstate='ee')
        elif mode in ['mumu', 'mu+mu-', 'mu-mu+', 'dimu', 'dimuon']:
            # generate events with kkmc
            add_kkmc_generator(path, finalstate='mu-mu+')
        elif mode in ['charged', 'mixed', 'signal']:
            add_evtgen_generator(path, finalstate=mode)
        elif mode in ['mugun']:
            path.add_module('ParticleGun',
                             pdgCodes=[13],
                             nTracks=1,
                             momentumGeneration='uniform',
                             momentumParams=[0.50, 3.00],
                             thetaGeneration='uniform',
                             thetaParams=[17, 150],
                             #thetaParams=[40, 90],
                             phiGeneration='uniform',
                             #phiParams=[0, 360],
                             phiParams=[91.5, 137.5],
                             vertexGeneration='uniform',
                             xVertexParams=[0.0, 0.0],
                             yVertexParams=[0.0, 0.0],
                             zVertexParams=[0.0, 0.0]
                             )



        else:
            raise Exception("only ee, mumu, charged, mixed, signal modes. Instead got: %s"%mode)
    
        # activate simulation of dead/masked pixel and reproduce detector gain, which will be
        #   applied at reconstruction level when the data GT is present in the DB chain
        path.add_module("ActivatePXDPixelMasker")
        path.add_module("ActivatePXDGainCalibrator")
    
        # detector simulation
        print('bkg overlay:', bkg_overlay, 'bkg files:',  bkg_files)
        #add_simulation(path, usePXDDataReduction=False, bkgfiles=bkg_files, bkgOverlay=bkg_overlay, **kwargs)
        add_simulation(path, bkgfiles=bkg_files, bkgOverlay=bkg_overlay, **kwargs)
    
        if output_cdst:
            add_cdst_output(path, filename=output_cdst)
    
        #return path



#def setup_database(path, global_tags=None, release=4):
#
#    global_tags = global_tags if not global_tags or isinstance(global_tags, (list,tuple)) else global_tags.split()
#
#    if global_tags:
#        if release in [2,3]:
#            basf2.reset_database()
#            basf2.use_database_chain()
#
#            for global_tag in global_tags:
#                if global_tag.startswith("/"):
#                    basf2.use_local_database(global_tag)
#                else:
#                    basf2.use_central_database(global_tag)
#
#        elif release in [4]:
#            basf2.conditions.override_globaltags()
#            basf2.conditions.expert_settings(usable_globaltag_states={"TESTING", "VALIDATED", 
#                                                                      "PUBLISHED", "RUNNING", "OPEN"})
#
#            for global_tag in global_tags:
#                if global_tag.startswith("/"):
#                    basf2.conditions.append_testing_payloads(global_tag)
#                else:
#                    basf2.conditions.append_globaltag(global_tag)
#        else:
#            raise ValueError(f"Selected release {release} is not implemented")
#


def remove_module(path, module):
    """
    Hack to remove a module from the input path.
    See https://questions.belle2.org/question/127/how-to-remove-a-module-from-a-path-in-basf2/

    Args:
        path (basf2.Path): the input path.
        module (str): the name of the module to be removed.
    Returns:
        basf2.Path: the path w/o the module.
    """

    new_path = basf2.create_path()
    found = False
    for m in path.modules():
        if m.name() != module:
            new_path.add_module(m)
        else:
            found = True
    if not found:
        basf2.B2ERROR(f"Could not find module {module} in path, cannot remove")
    return new_path



def set_module_params(path, name=None, type=None, recursive=False, strict=True, **kwargs):
    """Similar to basf2.set_module_parameters, with slightly different (maybe more correct?)
    recusrive behavior. The basf2 one raises exception if module is found ONLY in the subpaths, and
    not in the main one.

    Set the given set of parameters for all `modules <Module>` in a path which
    have the given ``name`` (see `Module.set_name`)

    Usage is similar to `register_module()` but this function will not create
    new modules but just adjust parameters for modules already in a `Path`

    >>> set_module_parameters(path, "Geometry", components=["PXD"], logLevel=LogLevel.WARNING)

    Parameters:
      path (basf2.Path): The path to search for the modules
      name (str): Then name of the module to set parameters for
      type (str): The type of the module to set parameters for.
      recursive (bool): if True also look in paths connected by conditions or `Path.for_each()`
      kwargs: Named parameters to be set for the module, see  `register_module()`
    """

    if name is None and type is None:
        raise ValueError("At least one of name or type has to be given")

    if not kwargs:
        raise ValueError("no module parameters given")

    #found = False
    founds = []
    for module in path.modules():
        if (name is None or module.name() == name) and (type is None or module.type() == type):
            # use register_module as this automatically takes care of logLevel
            # and debugLevel parameters
            basf2.register_module(module, **kwargs)
            #found = True
            founds.append(1)

        if recursive:
            if module.has_condition():
                for condition_path in module.get_all_condition_paths():
                    found = set_module_params(condition_path, name, type, recursive, strict=False, **kwargs)
                    founds.append(found)
            if module.type() == "SubEvent":
                for subpath in [p.values for p in module.available_params() if p.name == "path"]:
                    found = set_module_params(subpath, name, type, recursive, strict=False, **kwargs)
                    founds.append(found)

    found = sum(founds)
    if strict:
        if not found:
            raise KeyError("No module with given name found anywhere in the path")
        else:
            print(f"Parameter values changes in { found } place(s)")
    return found






def find_module(path, name=None, type=None):
    found = []
    for m in path.modules():
        #print(m.name(),m.type())
        if name and m.name()==name:
            found.append(m)
        elif type and m.type()==type:
            found.append(m)
        if hasattr(m, 'get_all_condition_paths'):
            for p in m.get_all_condition_paths():
                new = find_module(p, name, type)
                for m_ in new:
                    found.append(m_)
    return found



def get_module_params(module):
    return {x.name:x.values for x in module.available_params()}


