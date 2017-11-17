from future.utils import raise_

from builtins import object
import setup_paths
import numpy as np
from datetime import datetime
import nomadcore.ActivateLogging
from nomadcore.caching_backend import CachingLevel
from nomadcore.simple_parser import AncillaryParser, mainFunction, SkipFileException
from nomadcore.simple_parser import SimpleMatcher as SM
from TurbomoleCommon import get_metaInfo, RE_FLOAT, RE_DATE, RE_TIME
import logging, os
import TurbomoleCommon as Common
from EmbeddingParser import EmbeddingParser
from GradientParser import GradientParser
from MethodParser import MethodParser
from OrbitalParser import OrbitalParser
from SystemParser import SystemParser
from DSCFparser import DSCFparser
from ESCFparser import ESCFparser
from GRADparser import GRADparser
from RIDFTparser import RIDFTparser
from RICC2parser import RICC2parser
from STATPTparser import STATPTparser

############################################################
# This is the parser for the main file of turbomole.
############################################################

logger = logging.getLogger("nomad.turbomoleParser")

# the modules in this list don't generate any output useful for Nomad and are thus not parsed
IGNORED_MODULES = [
    "dscfserver",
    "eiger",
    "gradserver",
]


class _GeneralInfo(object):

    def __init__(self):
        self.module = None
        self.version = None
        self.node = None
        self.clean_end = None
        self.start_time = None
        self.end_time = None
        self.index_config = None
        self.index_geo = None
        self.index_method = None
        self.kinetic_energy = None
        self.potential_energy = None


class TurbomoleParserContext(object):

    def __init__(self):
        self.__data = dict()
        self.__invocations = list()
        self.__sampling_mode_section = None

    def __getitem__(self, item):
        return self.__data[item]

    def __setitem__(self, key, value):
        if key in self.__data:
            raise Exception("sub data storage '%s' already declared!" % key)
        self.__data[key] = value

    def __iter__(self):
        for key, sub_parser in self.__data.items():
            yield key, sub_parser

    def index_configuration(self):
        return self.__invocations[-1].index_config

    def index_method(self):
        return self.__invocations[-1].index_method

    def index_system(self):
        return self.__invocations[-1].index_geo

    def set_sampling_mode_section(self, index_settings):
        self.__sampling_mode_section = index_settings

    def purge_subparsers(self):
        for sub_parser in self.__data.values():
            sub_parser.purge_data()

    def build_module_matcher(self, module_name, subMatchers, generic=False):
        def open_section_config(backend, gIndex, section):
            self.__invocations.append(_GeneralInfo())
            self.__invocations[-1].index_config = gIndex

        def open_section_method(backend, gIndex, section):
            self.__invocations[-1].index_method = gIndex

        def open_section_geo(backend, gIndex, section):
            self.__invocations[-1].index_geo = gIndex

        def close_section_method(backend, gIndex, section):
            self["geo"].write_basis_set_mapping(self.index_configuration(), gIndex)
            self.__invocations[-1].kinetic_energy = self["method"].get_energy_kinetic()
            self.__invocations[-1].potential_energy = self["method"].get_energy_potential()
            backend.addValue("single_configuration_to_calculation_method_ref", gIndex)

        def close_section_system(backend, gIndex, section):
            backend.addValue("single_configuration_calculation_to_system_ref", gIndex)

        def process_module_invocation(backend, groups):
            if generic:
                logger.error("Turbomole module without dedicated parser found: %s" % groups[0])
            self.purge_subparsers()
            self.__invocations[-1].version = groups[2]
            self.__invocations[-1].node = groups[1]
            self.__invocations[-1].module = groups[0]

        return SM(r"\s*("+module_name+")\s*\(([^\)]+)\)\s*\:\s*TURBOMOLE\s+([a-zA-Z0-9.]+)",
                  name=(module_name if not generic else "unknown") + " module",
                  sections=[
                      "section_single_configuration_calculation",
                      "section_system",
                      "section_method"
                  ],
                  onOpen={
                      "section_single_configuration_calculation": open_section_config,
                      "section_system": open_section_geo,
                      "section_method": open_section_method
                  },
                  onClose={
                      "section_system": close_section_system,
                      "section_method": close_section_method
                  },
                  startReAction=process_module_invocation,
                  subMatchers=subMatchers
                  )

    def build_start_time_matcher(self):

        def set_start_time(backend, groups):
            utc_time = datetime.strptime("%sT%sZ" % (groups[0], groups[1]), "%Y-%m-%dT%H:%M:%S.%fZ")
            epoch_time = (utc_time - datetime(1970, 1, 1)).total_seconds()
            self.__invocations[-1].start_time = epoch_time

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
            self.__invocations[-1].clean_end = True

        def set_end_time(backend, groups):
            utc_time = datetime.strptime("%sT%sZ" % (groups[0], groups[1]), "%Y-%m-%dT%H:%M:%S.%fZ")
            epoch_time = (utc_time - datetime(1970, 1, 1)).total_seconds()
            self.__invocations[-1].end_time = epoch_time

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

    def onClose_section_run(self, backend, gIndex, section):
        if len(self.__invocations) > 0:
            versions = list(set(x.version for x in self.__invocations))
            if len(versions) > 1:
                logger.warning("found inconsistent code version information: %s" % str(versions))
            backend.addValue("program_version", versions[0])
            backend.addValue("run_clean_end", all(x.clean_end for x in self.__invocations))
            nodes = list(set(x.node for x in self.__invocations))
            if len(nodes) > 1:
                logger.warning("found inconsistent host node information: %s" % str(nodes))
            backend.addValue("x_turbomole_nodename", nodes[0])
            backend.addRealValue("time_run_date_start",
                                 min(x.start_time for x in self.__invocations))
            backend.addRealValue("time_run_date_end", max(x.end_time for x in self.__invocations))

        if self.__sampling_mode_section is not None:
            index = backend.openSection("section_frame_sequence")
            backend.addValue("frame_sequence_to_sampling_ref", self.__sampling_mode_section, index)
            frames_all = np.asarray([x.index_config for x in self.__invocations], dtype=int)
            frames_kinetic = np.asarray([x.index_config for x in self.__invocations
                                         if x.kinetic_energy], dtype=int)
            frames_potential = np.asarray([x.index_config for x in self.__invocations
                                         if x.potential_energy], dtype=int)
            energies_kinetic = np.asarray([x.kinetic_energy for x in self.__invocations
                                if x.kinetic_energy], dtype=float)
            energies_potential = np.asarray([x.potential_energy for x in self.__invocations
                                             if x.potential_energy], dtype=float)
            backend.addValue("number_of_frames_in_sequence", len(self.__invocations), index)
            backend.addValue("number_of_kinetic_energies_in_sequence", len(energies_kinetic), index)
            backend.addValue("number_of_potential_energies_in_sequence",
                             len(energies_potential), index)
            backend.addArrayValues("frame_sequence_local_frames_ref", frames_all, index)
            backend.addArrayValues("frame_sequence_kinetic_energy_frames", frames_kinetic, index)
            backend.addArrayValues("frame_sequence_potential_energy_frames", frames_potential,
                                   index)
            backend.addArrayValues("frame_sequence_kinetic_energy", energies_kinetic, index,
                                   unit="hartree")
            backend.addArrayValues("frame_sequence_potential_energy", energies_potential, index,
                                   unit="hartree")
            backend.closeSection("section_frame_sequence", index)


def build_root_parser(context):
    """Builds the SimpleMatcher to parse the main files of turbomole.
    Matches for subsections of the output file are generated in dedicated
    sub parsers provided as custom classes.

    Returns:
       SimpleMatcher that parses the main files of Turbomole.
    """

    def set_backends(backend, gIndex, section):
        for key, sub_parser in context:
            sub_parser.set_backend(backend)

    # shared subparsers created here are automatically stored in the context
    EmbeddingParser(context)
    GradientParser(context)
    MethodParser(context)
    OrbitalParser(context)
    SystemParser(context)

    # matches only those subprograms without dedicated parser

    sub_matcher = [
        context.build_start_time_matcher(),
        context["geo"].build_qm_geometry_matcher(),
        context["geo"].build_orbital_basis_matcher(),
        context["method"].build_dft_functional_matcher(),
        context["embedding"].build_embedding_matcher(),
        context["method"].build_dftd3_vdw_matcher(),
        context["method"].build_total_energy_matcher(),
        context["orbitals"].build_eigenstate_matcher(),
        context["gradient"].build_gradient_matcher(),
        context.build_end_time_matcher("[A-z0-9]+")
    ]
    generic = context.build_module_matcher("[A-z0-9]+", sub_matcher, generic=True)

    def ignore_mpi_slaves(backend, groups):
        if int(groups[0]) > 1:  # RICC2 module starts indexing processes from 0, others at 1
            backend.closeSection("section_run", 0)
            raise_(SkipFileException, "MPI slave output only")
    mpi_slaves = SM(r"\s*this\s+is\s+node-proc.\s+number\s+([0-9]+)\s+running\s+on\s+node\s+"
                    r"([^ ].*[^ ])\s*$",
                    name="MPI rank",
                    startReAction=ignore_mpi_slaves)

    def skip_ignored_modules(backend, groups):
        backend.closeSection("section_run", 0)
        raise_(SkipFileException, "MPI slave output only")
    ignored_modules = "|".join(IGNORED_MODULES)
    skip_modules = SM(r"\s*("+ignored_modules+")\s*\(([^\)]+)\)\s*\:\s*TURBOMOLE\s+([a-zA-Z0-9.]+)",
                      name="ignored module",
                      startReAction=skip_ignored_modules)

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
                  mpi_slaves,
                  skip_modules,
                  ESCFparser(context).build_parser(),
                  DSCFparser(context).build_parser(),
                  GRADparser(context).build_parser(),
                  RIDFTparser(context).build_parser(),
                  RICC2parser(context).build_parser(),
                  STATPTparser(context).build_parser(),
                  generic
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
