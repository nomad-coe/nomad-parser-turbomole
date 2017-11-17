import logging
import numpy as np
from nomadcore.simple_parser import SimpleMatcher as SM
from TurbomoleCommon import RE_FLOAT

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

    # TODO: verify mapping of approximations to general Electronic Structure Method class
    __wavefunction_models_map = {
        "MP2": "MP2",
        "CCS": "CCS",
        "CIS": "CIS",
        "CIS(D)": "CISD",
        "CIS(Dinf)": "CISD",
        "ADC(2)": "MP2",  # TODO: check paper to verify this mapping
        "CC2": "CCSD"
    }


    def __init__(self, context, key="method"):
        context[key] = self
        self.__context = context
        self.__backend = None
        self.spin_channels = 1
        self.k_points = 1
        self.__method = None
        self.__functional = None

    def purge_data(self):
        self.spin_channels = 1
        self.k_points = 1
        self.__method = None
        self.__functional = None

    def set_backend(self, backend):
        self.__backend = backend

    # matcher generation methods

    def add_default_functional(self, backend, groups):
        if not self.__method:
            self.__method = "DFT"
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
            self.__method = "DFT"
            backend.addValue("electronic_structure_method", "DFT")
            backend.addValue("calculation_method_kind", "absolute")
            for component in self.__functional:
                index = backend.openSection("section_XC_functionals")
                backend.addValue('XC_functional_name', component)
                backend.closeSection("section_XC_functionals", index)

        return SM(r"\s*density functional\s*$",
                  name="DFT functional",
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

    def build_wave_function_model_matcher(self):

        def determine_spin(backend, groups):
            if groups[0] == "restricted closed":
                pass
            elif groups[0] == "restricted open":
                self.__context["method"].spin_channels = 2
            elif groups[0] == "unrestricted open":
                self.__context["method"].spin_channels = 2
            else:
                logger.error("found unknown spin configuration in : %s" %  groups[0])

        def extract_wf_method(backend, groups):
            method = self.__wavefunction_models_map.get(groups[0], None)
            self.__method = method
            if method:
                backend.addValue("electronic_structure_method", method)
            else:
                logger.error("unknown wave-function model encountered: %s - %s" % groups)

        method_matcher = SM(r"\s*([^\s].*[^\s])\s*-\s*([^\s].*[^\s])\s*$",
                            name="WF model",
                            required=True,
                            startReAction=extract_wf_method
                            )

        # TODO: extract further RICC2 parameters?
        return SM(r"\s*([^ ].+)\s+shell\s+calculation\s+for\s+the\s+wavefunction\s+models:\s*$",
                  startReAction=determine_spin,
                  name="spin treatment",
                  required=True,
                  subMatchers=[
                      method_matcher
                  ]
                  )

    def build_dftd3_vdw_matcher(self):

        energy_matcher = SM(r"\s*Edisp\s+/kcal,\s*au\s*:\s*"+RE_FLOAT+"\s+"
                            "(?P<energy_van_der_Waals__hartree>"+RE_FLOAT+")\s*$",
                            name="vdW energy"
                            )

        return SM(r"\s*\|\s*DFTD3\s+(?P<x_turbomole_dft_d3_version>"
                  r"V[0-9]+.[0-9]+\s+Rev\s+[0-9]+)\s*\|\s*$",
                  name="DFT-D3 version",
                  subMatchers=[
                      energy_matcher
                  ],
                  fixedStartValues={"van_der_Waals_method": "DFT-D3"}
                  )
