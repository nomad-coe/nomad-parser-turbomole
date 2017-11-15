

from builtins import object
import setup_paths
import numpy as np
from datetime import datetime
import nomadcore.ActivateLogging
from nomadcore.caching_backend import CachingLevel
from nomadcore.simple_parser import AncillaryParser, mainFunction
from nomadcore.simple_parser import SimpleMatcher as SM
from TurbomoleCommon import get_metaInfo, RE_FLOAT, RE_DATE, RE_TIME
import logging, os
import TurbomoleCommon as Common
from GradientParser import GradientParser
from MethodParser import MethodParser
from OrbitalParser import OrbitalParser
from SystemParser import SystemParser
from ESCFparser import ESCFparser
from DSCFparser import DSCFparser
from GRADparser import GRADparser
from RIDFTparser import RIDFTparser
from RICC2parser import RICC2parser

############################################################
# This is the parser for the main file of turbomole.
############################################################

logger = logging.getLogger("nomad.turbomoleParser")


class TurbomoleParserContext(object):

    def __init__(self):
        self.__data = dict()
        self.functionals = []
        self.generic = False
        self.__general_info = {
            "module": list(),
            "version": list(),
            "node": list(),
            "clean_end": list(),
            "start_time": list(),
            "end_time": list()
        }

    def __getitem__(self, item):
        return self.__data[item]

    def __setitem__(self, key, value):
        if key in self.__data:
            raise Exception("sub data storage '%s' already declared!" % key)
        self.__data[key] = value

    def __iter__(self):
        for key, sub_parser in self.__data.items():
            yield key, sub_parser

    def purge_subparsers(self):
        for sub_parser in self.__data.values():
            sub_parser.purge_data()

    def get_module_invocation(self, module_name):
        return r"\s*("+module_name+")\s*\(([^\)]+)\)\s*\:\s*TURBOMOLE\s+([a-zA-Z0-9.]+)"

    def process_module_invocation(self, backend, groups):
        self.purge_subparsers()
        self.generic = False
        self.__general_info["version"].append(groups[2])
        self.__general_info["node"].append(groups[1])
        self.__general_info["module"].append(groups[0])

    def build_start_time_matcher(self):

        def set_start_time(backend, groups):
            utc_time = datetime.strptime("%sT%sZ" % (groups[0], groups[1]), "%Y-%m-%dT%H:%M:%S.%fZ")
            epoch_time = (utc_time - datetime(1970, 1, 1)).total_seconds()
            self.__general_info["start_time"].append(epoch_time)

        return SM(r"\s*("+RE_DATE+r")\s+("+RE_TIME+r")\s*$",
                  name="start timestamp",
                  startReAction=set_start_time
                  )

    def build_end_time_matcher(self, module_name):

        def get_run_time(backend, groups):
            time = float(groups[3])
            if groups[2]:
                time += 60.0 * float(groups[2])
            if groups[1]:
                time += 3600.0 * float(groups[1])
            if groups[0]:
                time += 86400.0 * float(groups[0])
            backend.addRealValue("time_calculation", time)

        def set_clean_end(backend, groups):
            self.__general_info["clean_end"].append(True)

        def set_end_time(backend, groups):
            utc_time = datetime.strptime("%sT%sZ" % (groups[0], groups[1]), "%Y-%m-%dT%H:%M:%S.%fZ")
            epoch_time = (utc_time - datetime(1970, 1, 1)).total_seconds()
            self.__general_info["end_time"].append(epoch_time)

        walltime = SM(r"\s*total\s+wall-time\s*:(?:\s*("+RE_FLOAT+")\s+days)?"
                      r"(?:\s*("+RE_FLOAT+")\s+hours)?"
                      r"(?:\s*("+RE_FLOAT+")\s+minutes\s+and)?"
                      r"(?:\s*("+RE_FLOAT+")\s+seconds)\s*$",
                      name="wall time",
                      startReAction=get_run_time
                      )
        clean_end = SM("\s*\*{4}\s*"+module_name+"\s*:\s*all\s+done\s*\*{4}\s*$",
                       name="clean end",
                       startReAction=set_clean_end
                       )
        end_date = SM(r"\s*("+RE_DATE+r")\s+("+RE_TIME+r")\s*$",
                      name="end timestamp",
                      startReAction=set_end_time
                      )

        return SM(r"\s*total\s+cpu-time\s*:(?:\s*("+RE_FLOAT+")\s+days)?"
                  r"(?:\s*("+RE_FLOAT+")\s+hours)?"
                  r"(?:\s*("+RE_FLOAT+")\s+minutes\s+and)?"
                  r"(?:\s*("+RE_FLOAT+")\s+seconds)\s*$",
                  name="cpu time",
                  subMatchers=[
                      walltime,
                      SM(r"\s*-{20,}\s*$", name="<format>", coverageIgnore=True),
                      clean_end,
                      end_date
                  ]
                  )

    def initialize_values(self):
        """Initializes the values of certain variables.

        This allows a consistent setting and resetting of the variables,
        when the parsing starts and when a section_run closes.
        """
        self.secMethodIndex = None
        self.secSystemDescriptionIndex = None
        self.lastCalculationGIndex = None
        self.singleConfCalcs = []
        self.geoConvergence = None

    def startedParsing(self, fInName, parser):
        """Function is called when the parsing starts.

        Get compiled parser, filename and metadata.

        Args:
            fInName: The file name on which the current parser is running.
            parser: The compiled parser. Is an object of the class SimpleParser in nomadcore.simple_parser.py.
        """
        self.parser = parser
        self.fName = fInName
        # save metadata
        self.metaInfoEnv = self.parser.parserBuilder.metaInfoEnv
        # allows to reset values if the same superContext is used to parse different files
        self.initialize_values()

    def onClose_section_run(self, backend, gIndex, section):
        if len(self.__general_info["version"]) > 0:
            if len(set(self.__general_info["version"])) > 1:
                logger.warning("found inconsistent code version information: %s" %
                               self.__general_info["version"])
            backend.addValue("program_version", self.__general_info["version"][0])
            if len(self.__general_info["clean_end"]) != len(self.__general_info["module"]):
                backend.addValue("run_clean_end", False)
            else:
                backend.addValue("run_clean_end", False in self.__general_info["clean_end"])
        if len(self.__general_info["node"]) > 0:
            if len(set(self.__general_info["node"])) > 1:
                logger.warning("found inconsistent host node information: %s" %
                               self.__general_info["node"])
            backend.addValue("x_turbomole_nodename", self.__general_info["node"][0])
        if len(self.__general_info["start_time"]) > 0:
            backend.addRealValue("time_run_date_start", min(self.__general_info["start_time"]))
        if len(self.__general_info["end_time"]) > 0:
            backend.addRealValue("time_run_date_end", min(self.__general_info["end_time"]))

        if self.geoConvergence is not None:
            backend.addValue('x_turbomole_geometry_optimization_converged', self.geoConvergence)
            if self.geoConvergence is True:
                sampling_method = "geometry_optimization"
            elif len(self.singleConfCalcs) > 1:
                pass
            else:
                return
            samplingGIndex = backend.openSection("section_sampling_method")
            backend.addValue("sampling_method", sampling_method)
            backend.closeSection("section_sampling_method", samplingGIndex)
            frameSequenceGIndex = backend.openSection("section_frame_sequence")
            backend.addValue("frame_sequence_to_sampling_ref", samplingGIndex)
            backend.addArrayValues("frame_sequence_local_frames_ref", np.asarray(self.singleConfCalcs))
            backend.closeSection("section_frame_sequence", frameSequenceGIndex)

    def onOpen_section_method(self, backend, gIndex, section):
        # keep track of the latest method section
        self.secMethodIndex = gIndex

    ###################################################################
    # (3.4) onClose for geometry and force (section_system)
    # todo: maybe we can move the force to onClose_section_single_configuration_calculation in the future. 
    ###################################################################

    def onOpen_section_system(self, backend, gIndex, section):
        # keep track of the latest system description section
        self.secSystemDescriptionIndex = gIndex

    def onOpen_section_single_configuration_calculation(self, backend, gIndex, section):
        self.singleConfCalcs.append(gIndex)

    def onClose_section_single_configuration_calculation(self, backend, gIndex, section):
        if section['x_turbomole_geometry_optimization_converged'] is not None:
            if section['x_turbomole_geometry_optimization_converged'][-1] == 'FULFILLED':
                self.geoConvergence = True
            else:
                self.geoConvergence = False

        if self.generic:
            if self.secMethodIndex:
                backend.addValue('single_configuration_to_calculation_method_ref',
                                 self.secMethodIndex)
            if self.secSystemDescriptionIndex:
                backend.addValue('single_configuration_calculation_to_system_ref',
                                 self.secSystemDescriptionIndex)

    def setStartingPointCalculation(self, parser):
        backend = parser.backend
        backend.openSection('section_calculation_to_calculation_refs')
        if self.lastCalculationGIndex:
            backend.addValue('calculation_to_calculation_ref', self.lastCalculationGIndex)
        backend.addValue('calculation_to_calculation_kind', 'pertubative GW')
        #        backend.closeSection('section_calculation_to_calculation_refs')
        return None

#############################################################
#################[2] MAIN PARSER STARTS HERE  ###############
#############################################################

def build_root_parser(context):
    """Builds the SimpleMatcher to parse the main file of turbomole.
    Matches for subsections of the output file are generated in dedicated
    functions below.

    Returns:
       SimpleMatcher that parses main file of Turbomole.
    """

    def set_backends(backend, gIndex, section):
        for key, sub_parser in context:
            sub_parser.set_backend(backend)

    # shared subparsers created here are automatically stored in the context
    GradientParser(context)
    SystemParser(context)
    OrbitalParser(context)
    MethodParser(context)

    def set_generic(backend, groups):
        context.process_module_invocation(backend, groups)
        context.generic = True

    def finalize_system_data(backend, gIndex, section):
        context["geo"].finalize_sections()
        context["method"].close_method_section()

    generic_modules = r"aoforce|cosmoprep|egrad|evib|frog|gradsammel|" \
                      r"hesssammel|moloch|odft|relax|rirpa|sdg|thirdsammel|" \
                      r"vibration|atbandbta|define|eigerf|fdetools|grad|haga|" \
                      r"intense|mpgrad|proper|ricc2|rimp2|ruecker|statpt|tm2molden|" \
                      r"woelfling|bsseenergy|freeh|gradruecker|riper|hessruecker|" \
                      r"mdprep|mpshift|rdgrad|ricctools|rimp2prep|sammler|thirdruecker|uff"

    # matches only those subprograms without dedicated parser
    generic = SM(context.get_module_invocation(generic_modules),
                 name="NewRun",
                 repeats=True,
                 sections=["section_single_configuration_calculation"],
                 onClose={"section_single_configuration_calculation": finalize_system_data},
                 startReAction=set_generic,
                 subMatchers=[
                     SM(name="general info",
                         startReStr=r"\s*Copyright \(C\) ",
                         subMatchers=[
                             context.build_start_time_matcher(),
                             context["geo"].build_qm_geometry_matcher(simple_mode=True),
                             context["geo"].build_orbital_basis_matcher(),
                             context["method"].build_dft_functional_matcher(simple_mode=True)
                         ]),
                     SM(r"\s*1e\-*integrals will be neglected if expon",
                        name="Single Config",
                        subMatchers=[
                            SM(r"\s*\|\s*EMBEDDING IN PERIODIC POINT CHARGES\s*\|",
                               name = "Embedding",
                               subMatchers=[
                                   #SmearingOccupation,
                                   context["geo"].build_embedding_matcher(),
                               ]
                               ),
                            SM(name='TotalEnergyForEachScfCycle',
                                startReStr = r"\s*scf convergence criterion",
                                subMatchers=[
                                    #SmearingOccupation,
                                    Common.build_total_energy_matcher()
                                ]),
                            context["orbitals"].build_eigenstate_matcher(),
                        ]),
                     SM(r"\s*Energy of reference wave function is",
                        name="PostHFTotalEnergies",
                        subMatchers=[
                            build_total_energy_coupled_cluster_matcher()
                        ]
                        ),
                     context["gradient"].build_gradient_matcher(),
                     context.build_end_time_matcher("(?:"+generic_modules+")")
                 ]
                 )

    return SM(name="Root",
              startReStr="",
              forwardMatch=True,
              sections=["section_run"],
              onOpen={"section_run": set_backends},
              subFlags=SM.SubFlags.Unordered,
              fixedStartValues={
                  'program_name': 'turbomole',
                  'program_basis_set_type': 'GTOs'
              },
              subMatchers=[
                  ESCFparser(context).build_parser(),
                  DSCFparser(context).build_parser(),
                  GRADparser(context).build_parser(),
                  RIDFTparser(context).build_parser(),
                  RICC2parser(context).build_parser(),
                  generic
                  # build_relaxation_matcher()
              ]
              )

def build_occupation_smearing_matcher():
    return SM (name = "smearing",
               startReStr = r"\s*and increment of one",
               sections = ["section_method"],
               subMatchers = [
                   SM (r"\s*(?P<smearing_kind>[a-zA-Z]+)\s*smearing switched on"),
                   SM (r"\s*Final electron temperature\:\s*(?P<smearing_width>[0-9.eEdD]+)")
               ])

def build_total_energy_coupled_cluster_matcher():
    return SM (name = 'TotalEnergyCC',
               #startReStr = r"\s*Calculate\s*integrals\s*\(*ia\|*jb\)*\s*for MP2 start guess",
               startReStr = r"\s*=========================================================================",
               forwardMatch = True,
               subMatchers = [
                   #SM (r"\s*\*\s*RHF  energy\s*\:\s*(?P<x_turbomole_HF_total_energy_final__eV>[-+0-9.eEdD]+)"),
                   SM (r"\s*\*\s*RHF  energy\s*\:\s*(?P<energy_total__eV>[-+0-9.eEdD]+)"),
                   #SM (r"\s*\*\s*UHF  energy\s*\:\s*(?P<x_turbomole_HF_total_energy_final__eV>[-+0-9.eEdD]+)"),
                   SM (r"\s*\*\s*UHF  energy\s*\:\s*(?P<energy_total__eV>[-+0-9.eEdD]+)"),
                   SM (r"\s*\*\s*Final MP2 energy\s*\:\s*(?P<x_turbomole_MP2_total_energy_final__eV>[-+0-9.eEdD]+)"),
                   SM (r"\s*\*\s*Final CCSD energy\s*\:\s*(?P<x_turbomole_CCSD_total_energy_final__eV>[-+0-9.eEdD]+)"),
                   SM (r"\s*\*\s*Final CC2 energy\s*\:\s*(?P<x_turbomole_CC2_total_energy_final__eV>[-+0-9.eEdD]+)"),
                   SM (r"\s*\*\s*Final CCSD\(T\) energy\s*\:\s*(?P<x_turbomole_CCSDparT_total_energy_final__eV>[-+0-9.eEdD]+)"),
                   SM (r"\s*\*\s*D1 diagnostic \(CCSD\)\s*\:\s*(?P<x_turbomole_D1_diagnostic>[-+0-9.eEdD]+)")

               ])

def build_relaxation_matcher():
    return SM (name = "relaxation",
               #sections = ["section_single_configuration_calculation"],
               sections = ["section_single_configuration_calculation"],
               startReStr = r"\s*CONVERGENCY CRITERIA (?P<x_turbomole_geometry_optimization_converged>FULFILLED) IN CYCLE",
               subMatchers = [])

def get_cachingLevelForMetaName(metaInfoEnv):
    """Sets the caching level for the metadata.

    Args:
        metaInfoEnv: metadata which is an object of the class InfoKindEnv in nomadcore.local_meta_info.py.

    Returns:
        Dictionary with metaname as key and caching level as value. 
    """
    # manually adjust caching of metadata
    cachingLevelForMetaName = {
        'eigenvalues_eigenvalues': CachingLevel.Cache,
        'eigenvalues_kpoints':CachingLevel.Cache,
        'x_turbomole_geometry_optimization_converged': CachingLevel.Cache
    }

    # Set caching for temparary storage variables
    for name in metaInfoEnv.infoKinds:
        if (   name.startswith('x_turbomole_controlInOut')
               or name.startswith('x_turbomole_geometry')
               or name.startswith('x_turbomole_embed')):
            cachingLevelForMetaName[name] = CachingLevel.Cache
    return cachingLevelForMetaName


def main():
    """Main function.

    Set up everything for the parsing of the turbomole main file and run the parsing.
    """
    # get main file description
    context = TurbomoleParserContext()
    # loading metadata from nomad-meta-info/meta_info/nomad_meta_info/turbomole.nomadmetainfo.json
    path = "../../../../nomad-meta-info/meta_info/nomad_meta_info/turbomole.nomadmetainfo.json"
    metaInfoPath = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)), path))
    metaInfoEnv = get_metaInfo(metaInfoPath)
    # set parser info
    parserInfo = {'name': 'turbomole-parser', 'version': '1.0'}
    # get caching level for metadata
    cachingLevelForMetaName = get_cachingLevelForMetaName(metaInfoEnv)
    # start parsing
    mainFunction(mainFileDescription=build_root_parser(context),
                 metaInfoEnv=metaInfoEnv,
                 parserInfo=parserInfo,
                 cachingLevelForMetaName=cachingLevelForMetaName,
                 superContext=context)


if __name__ == "__main__":
    main()
