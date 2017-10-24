import logging
import numpy as np
from nomadcore.simple_parser import SimpleMatcher as SM

logger = logging.getLogger("nomad.turbomoleParser")


class MethodParser(object):

    __functional_map = {
        "S-VWN":  ["LDA_X", "LDA_C_VWN_3"],
        "PWLDA":  ["LDA_X", "LDA_C_PW"],
        "B-VWN":  ["GGA_X_B88", "LDA_C_VWN"],
        "B-LYP":  ["GGA_X_B88", "GGA_C_LYP"],
        "B-P":    ["GGA_X_B88", "GGA_C_P86"],
        "B-P86":  ["GGA_X_B88", "GGA_C_P86"],
        "PBE":    ["GGA_X_PBE", "GGA_C_PBE"],
        "TPSS":   ["MGGA_X_TPSS", "MGGA_C_TPSS"],
        "M06":    ["MGGA_X_M06", "MGGA_C_M06"],
        "BH-LYP": ["HYB_GGA_XC_BHANDHLYP"],
        "B3-LYP": ["HYB_GGA_XC_B3LYP"],
        "PBE0":   ["HYB_GGA_XC_PBEH"],
        "TPSSh":  ["HYB_MGGA_XC_TPSSH"],
        "M06-2X": ["MGGA_X_M06_2X", "MGGA_C_M06_2X"],
        "B2-PLYP":["HYB_GGA_XC_B2PLYP"]
    }

    def __init__(self, context, key="method"):
        context[key] = self
        self.__context = context
        self.__backend = None
        self.__spin_channels = 1
        self.__index_method = -1
        self.__method = None
        self.__functional = None

    def set_backend(self, backend):
        self.__backend = backend

    # getter methods

    def spin_channels(self):
        return self.__spin_channels

    def index_method(self):
        return self.__index_method

    # matcher generation methods

    def close_method_section(self):
        if self.__index_method != -1:
            self.__backend.closeSection("section_method", self.__index_method)

    def add_default_functional(self):
        self.__index_method = self.__backend.openSection("section_method")
        self.__backend.addValue("electronic_structure_method", "DFT")
        self.__backend.addValue("calculation_method_kind", "absolute")
        index = self.__backend.openSection("section_XC_functionals")
        self.__backend.addValue('XC_functional_name', "HF_X")
        self.__backend.closeSection("section_XC_functionals", index)

    def build_uhf_matcher(self):

        def set_spin_polarized(backend, groups):
            self.__spin_channels = 2

        return SM(r"\s*UHF mode switched on !",
                  name="UHF switch",
                  startReAction=set_spin_polarized
                  )

    # TODO: add support for remaining XC-functionals in Turbomole + custom mixing combinations
    def build_dft_functional_matcher(self):
        exchange = SM(r"\s*exchange:\s*(?P<x_turbomole_functional_type_exchange>.+)")
        correlation = SM(r"\s*correlation:\s*(?P<x_turbomole_functional_type_correlation>.+)")

        def set_functional(backend, groups):
            if groups[0] in self.__functional_map:
                self.__functional = self.__functional_map[groups[0]]
            else:
                self.__functional = ["UNKNOWN"]
                logger.warning("XC-functional '%s' not known!" % groups[0])
            backend.addValue("electronic_structure_method", "DFT")
            backend.addValue("calculation_method_kind", "absolute")
            for component in self.__functional:
                index = backend.openSection("section_XC_functionals")
                backend.addValue('XC_functional_name', component)
                backend.closeSection("section_XC_functionals", index)

        def start_section(backend, groups):
            self.__index_method = backend.openSection("section_method")

        return SM(r"\s*density functional\s*$",
                  name="DFT functional",
                  startReAction=start_section,
                  subMatchers=[
                      SM(r"\s*-{5,}\s*$",
                         name="<format>",
                         coverageIgnore=True
                         ),
                      SM(r"\s*([A-z0-9-]+)\s+functional",
                         name="XC Functional",
                         startReAction=set_functional
                         ),
                      SM(r"\s*([A-z0-9-]+)\s+meta-GGA functional",
                         name="XC Functional",
                         startReAction=set_functional
                         ),
                      SM(r".*functional\s*:\s*([A-z0-9-]+)",
                         name="XC Functional",
                         startReAction=set_functional
                         ),
                      exchange,
                      correlation
                  ]
                  )

