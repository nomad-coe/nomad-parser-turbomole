# Copyright 2016-2018 Arvid Conrad Ihrig, Aliaksei Mazheika
#                     Fritz-Haber-Institut der Max-Planck-Gesellschaft
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from builtins import object
import setup_paths
import numpy as np
import nomadcore.ActivateLogging
from nomadcore.caching_backend import CachingLevel
from nomadcore.simple_parser import mainFunction
from nomadcore.simple_parser import SimpleMatcher as SM
from turbomoleparser.TurbomoleCommon import get_metaInfo
#from TurbomoleCommon import get_metaInfo, write_controlIn, write_k_grid, write_xc_functional
import logging, os, re, sys

############################################################
# This is the parser for the control.in file of FHI-aims.
############################################################

logger = logging.getLogger("nomad.TurbomoleControlInParser")

class TurbomoleControlInParserContext(object):
    """Context for parsing Turbomole control file.

    Attributes:
        sectionRun: Stores the parsed value/sections found in section_run.

    The onClose_ functions allow processing and writing of cached values after a section is closed.
    They take the following arguments:
        backend: Class that takes care of writing and caching of metadata.
        gIndex: Index of the section that is closed.
        section: The cached values and sections that were found in the section that is closed.
    """
    def __init__(self, writeSectionRun = True):
        """Args:
            writeSectionRun: Deteremines if metadata is written on close of section_run
                or stored in sectionRun.
        """
        self.writeSectionRun = writeSectionRun

    def startedParsing(self, fInName, parser):
        """Function is called when the parsing starts and the compiled parser is obtained.

        Args:
            fInName: The file name on which the current parser is running.
            parser: The compiled parser. Is an object of the class SimpleParser in nomadcore.simple_parser.py.
        """
        self.parser = parser
        # save metadata
        self.metaInfoEnv = parser.parserBuilder.metaInfoEnv
        # allows to reset values if the same superContext is used to parse different files
        self.sectionRun = None


def build_TurbomoleControlInKeywordsSimpleMatchers():
    """Builds the list of SimpleMatchers to parse the control keywords of Turbomole.

    SimpleMatchers are called with 'SM (' as this string has length 4,
    which allows nice formating of nested SimpleMatchers in python.

    Returns:
       List of SimpleMatchers that parses control.in keywords of Turbomole.
    """
    # Now follows the list to match the keywords from the control.in.
    # Explicitly add ^ to ensure that the keyword is not within a comment.
    # Repating occurrences of the same keywords are captured.
    # List the matchers in alphabetical order according to keyword name.
    #
    return [
        SM (r"\s*\$operating system\s+(?P<x_turbomole_controlIn_operating_system>[a-zA-Z\s]+)", repeats = True),
        SM (r"\s*\$symmetry\s+(?P<x_turbomole_controlIn_symmetry>[a-zA-Z0-9])", repeats = True),
        SM (startReStr = r"\s*\$atoms\s+",
            forwardMatch = True,
            #repeats = True,
            subMatchers = [
            SM (r"\s*(?P<x_turbomole_controlIn_atom_label>[a-zA-Z]+)\s+(?P<x_turbomole_controlIn_atom_number>[0-9]+)\s*(?: \\)", repeats = True)#,
            #SM (r"\s*(?P<x_turbomole_controlIn_atom_label>[a-zA-Z]+)\s+(?P<x_turbomole_controlIn_atom_number>[-0-9]+)\s*(?: \\)", repeats = True)
            ]),
        SM (r"\s*\$pople\s+(?P<x_turbomole_controlIn_pople_kind>[a-zA-Z]+)", repeats = True),
        SM (r"\s*\$uhfmo_alpha\s*file\=alpha"),
        SM (r"\s*\$uhfmo_beta\s*file\=beta"),
        SM (r"\s*\$scfiterlimit\s+(?P<x_turbomole_controlIn_scf_iter_limit>[0-9]+)"),
        SM (r"\s*\$scfconv\s+(?P<x_turbomole_controlIn_scf_conv>[0-9]+)"),
        #SM (r"\s*\$thize\s+(?P<x_turbomole_controlIn_number_of_integral_stored>[.0-9.eEdD]+)"),
        SM (r"\s*\$thime\s+(?P<x_turbomole_controlIn_time_for_integral_calc>[0-9]+)"),
        SM (r"\s*\$scfdamp\s+start\=(?P<x_turbomole_controlIn_damping_parameter_start>[.0-9]+)\s*step\=\s*(?P<x_turbomole_controlIn_damping_parameter_step>[.0-9]+)\s*min\=\s*(?P<x_turbomole_controlIn_damping_parameter_min>[.0-9]+)"),
        SM (startReStr = r"\s*\$scfintunit",
            forwardMatch = True,
            subMatchers = [
            SM (r"\s*unit\=(?P<x_turbomole_controlIn_scfintunit_unit>[0-9]+)\s*size\=\s*(?P<x_turbomole_controlIn_scfintunit_size>[0-9]+)\s*file\=\s*(?P<x_turbomole_controlIn_scfintunit_file>[a-zA-Z]+)")#,
            ]),
        SM (startReStr = r"\s*\$drvopt",
            forwardMatch = True,
            subMatchers = [
            SM (r"\s*(?:cartesian)\s+(?P<x_turbomole_controlIn_cartesian_status>[a-zA-Z]+)"),
            SM (r"\s*(?:basis)\s+(?P<x_turbomole_controlIn_basis_status>[a-z-A-Z]+)"),
            SM (r"\s*(?:global)\s+(?P<x_turbomole_controlIn_global_status>[a-zA-Z]+)"),
            SM (r"\s*(?:hessian)\s+(?P<x_turbomole_controlIn_hessian_status>[a-z-A-Z]+)"),
            SM (r"\s*(?:dipole)\s+(?P<x_turbomole_controlIn_dipole_status>[a-zA-Z]+)"),
            SM (r"\s*(?:\$interconversion)\s+(?P<x_turbomole_controlIn_interconversion_status>[a-zA-Z]+)")
            ])
        ]

def build_TurbomoleControlInFileSimpleMatcher():
    """Builds the SimpleMatcher to parse the control.in file of Turbomole.

    SimpleMatchers are called with 'SM (' as this string has length 4,
    which allows nice formating of nested SimpleMatchers in python.

    Returns:
       SimpleMatcher that parses control.in file of Turbomole.
    """
    return SM (name = 'Root1',
        startReStr = "",
        sections = ['section_run'],
        forwardMatch = True,
        weak = True,
        subMatchers = [
        SM (name = 'Root2',
            startReStr = "",
            sections = ['section_method'],
            forwardMatch = True,
            weak = True,
            # The search is done unordered since the keywords do not appear in a specific order.
            subFlags = SM.SubFlags.Unordered,
            subMatchers = build_TurbomoleControlInKeywordsSimpleMatchers()
            )
        ])

def get_cachingLevelForMetaName(metaInfoEnv, CachingLvl):
    """Sets the caching level for the metadata.

    Args:
        metaInfoEnv: metadata which is an object of the class InfoKindEnv in nomadcore.local_meta_info.py.
        CachingLvl: Sets the CachingLevel for the sections method and run. This allows to run the parser
            without opening new sections.

    Returns:
        Dictionary with metaname as key and caching level as value.
    """
    # manually adjust caching of metadata
    cachingLevelForMetaName = {
                               'section_method': CachingLvl,
                               'section_run': CachingLvl,
                              }
    # Set all controlIn metadata to Cache to capture multiple occurrences of keywords and
    # their last value is then written by the onClose routine in the TurbomoleControlInParserContext.
    for name in metaInfoEnv.infoKinds:
        if (name.startswith('x_turbomole_controlIn_') and not name.startswith('x_turbomole_controlIn_basis_')
            and not name.startswith('x_turbomole_controlIn_species_')):
            cachingLevelForMetaName[name] = CachingLevel.Cache
    return cachingLevelForMetaName

def main(CachingLvl):
    """Main function.

    Set up everything for the parsing of the Turbomole control.in file and run the parsing.

    Args:
        CachingLvl: Sets the CachingLevel for the sections method and run. This allows to run the parser
            without opening new sections
    """
    # get control.in file description
    TurbomoleControlInSimpleMatcher = build_TurbomoleControlInFileSimpleMatcher()
    # loading metadata from nomad-meta-info/meta_info/nomad_meta_info/turbomole.nomadmetainfo.json
    metaInfoPath = os.path.normpath(os.path.join(os.path.dirname(os.path.abspath(__file__)),"../../../../nomad-meta-info/meta_info/nomad_meta_info/turbomole.nomadmetainfo.json"))
    metaInfoEnv = get_metaInfo(metaInfoPath)
    # set parser info
    parserInfo = {'name':'x_turbomole-control-in-parser', 'version': '1.0'}
    # get caching level for metadata
    cachingLevelForMetaName = get_cachingLevelForMetaName(metaInfoEnv, CachingLvl)
    # start parsing
    mainFunction(mainFileDescription = TurbomoleControlInSimpleMatcher,
                 metaInfoEnv = metaInfoEnv,
                 parserInfo = parserInfo,
                 cachingLevelForMetaName = cachingLevelForMetaName,
                 superContext = TurbomoleControlInParserContext())

if __name__ == "__main__":
    main(CachingLevel.Forward)

