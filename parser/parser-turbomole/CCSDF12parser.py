"""This module constructs the parser for the CCSDF12 module from TurboMole"""

import logging
import re
from nomadcore.simple_parser import SimpleMatcher as SM
from TurbomoleCommon import RE_FLOAT
import TurbomoleCommon as Common

logger = logging.getLogger("nomad.turbomoleParser")


class CCSDF12parser(object):

    def __init__(self, context, key="ccsdf12"):
        context[key] = self
        self.__context = context
        self.__backend = None
        self.__previous_energy = None
        self.__correlation_indices = dict()

    def purge_data(self):
        self.__previous_energy = None
        self.__correlation_indices = dict()

    def set_backend(self, backend):
        self.__backend = backend

    def build_parser(self):
        references = SM(r"\s{5,}[^+ ]+",
                        name="references",
                        coverageIgnore=True,
                        repeats=True,
                        )
        header = SM(r"\s*\*\s*C\s*C\s*S\s*D\s*F\s*1\s*2\s+P\s*R\s*O\s*G\s*R\s*A\s*M\s*\*\s*$",
                    name="header",
                    coverageIgnore=True,
                    subMatchers=[references],
                    endReStr=r"\s*\+-+\+"
                    )

        sub_matchers = [
            self.__context.build_start_time_matcher(),
            header,
            self.__context["geo"].build_qm_geometry_matcher(),
            self.__context["geo"].build_orbital_basis_matcher(),
            self.__context["method"].build_wave_function_model_matcher(),
            self.__context["geo"].build_auxiliary_basis_matcher(),
            self.__build_mp2_starting_point_matcher(),
            self.__build_cc_scf_matcher(),
            self.__context["method"].build_correlation_energy_matcher(),
            self.__context.build_end_time_matcher("ccsdf12")
        ]

        return self.__context.build_module_matcher("ccsdf12", sub_matchers, "CCSDF12")

    def __create_auxiliary_sections(self, backend, method):
        if method == self.__context["method"].get_main_method():
            # main method invoked, write into the primary single_configuration_calculation
            index_config = self.__context.index_configuration()
            index_method = self.__context.index_method()
        else:
            # intermediate calculation, create extra method/config sections for results
            index_config = backend.openSection("section_single_configuration_calculation")
            index_method = backend.openSection("section_method")
            backend.addValue("single_configuration_to_calculation_method_ref", index_method,
                             index_config)
            backend.addValue("single_configuration_calculation_to_system_ref",
                             self.__context.index_system(), index_config)
            backend.addValue("electronic_structure_method", method, index_method)
        return index_config, index_method

    def __close_auxiliary_sections(self, backend, index_config, index_method):
        if index_config != self.__context.index_configuration():
            backend.closeSection("section_single_configuration_calculation", index_config)
        if index_method != self.__context.index_method():
            backend.closeSection("section_method", index_method)

    def __write_correlation_data(self, backend, data, index_config):
        self.__correlation_indices[data["method"]] = index_config
        backend.addRealValue("energy_current", data["e_corr"], index_config,
                             unit="hartree")
        backend.addRealValue("energy_total", data["e_total"], index_config,
                             unit="hartree")

    def __build_mp2_starting_point_matcher(self):

        def write_data(backend, gIndex, section):
            data = self.__context["method"].get_correlation_method_data()
            index_config, index_method = self.__create_auxiliary_sections(backend, data["method"])
            self.__write_correlation_data(backend, data, index_config)
            self.__close_auxiliary_sections(backend, index_config, index_method)

        return SM(r"\s*Calculate\s+integrals\s+\(ia\|jb\)\s+for\s+MP2\s+start\s+guess\s*$",
                  name="MP2 starting point",
                  subMatchers=[
                      self.__context["method"].build_correlation_energy_matcher()
                  ],
                  onClose={None: write_data}
                  )


    def __build_cc_scf_matcher(self):
        iterations = list()

        def extract_iteration_data(backend, groups):
            e_tot = float(groups[0])
            if self.__previous_energy:
                iterations.append({"e_tot": e_tot, "e_change": e_tot - self.__previous_energy})
            else:
                iterations.append({"e_tot": e_tot})
            self.__previous_energy = e_tot

        def convergence(backend, groups):
            if len(iterations) != int(groups[0]):
                logger.error("number of CCSD-iteration doesn't match the collected data! "
                             "(found %i, but there should be %i)" %
                             (len(iterations), int(groups[0])))

        def write_data(backend, gIndex, section):
            data = self.__context["method"].get_correlation_method_data()
            index_config, index_method = self.__create_auxiliary_sections(backend, data["method"])
            # TODO: write SCF data
            # backend.addRealValue("energy_change_scf_iteration",
            #                      float(groups[0])-self.__previous_energy, unit="hartree")
            # energy_total_scf_iteration
            self.__write_correlation_data(backend, data, index_config)
            self.__close_auxiliary_sections(backend, index_config, index_method)

        scf_iteration = SM(r"\s*[0-9]+\s+("+RE_FLOAT+")" + 5 * (r"\s+"+RE_FLOAT) + "\s*$",
                           name="iteration",
                           repeats=True,
                           startReAction=extract_iteration_data
                           )
        scf_converged = SM(r"\s*CC\s+equations\s+converged\s+in\s+([0-9]+)\s+iterations\.\s*$",
                           name="CC SCF converged",
                           startReAction=convergence)

        scf_cycle = SM(r"\s*Iter\.\s+CCSD\s+energy\s+Norm\(Omega\)\s+"
                       r"Norm\(t1\)\s+Norm\(t2\)\s+cpu\s+wall",
                       name="CCSD SCF iterations",
                       subMatchers=[
                           scf_iteration,
                           scf_converged
                       ],
                       )

        return SM(r"\s*\*\s*OPTIMIZATION\s+OF\s+THE\s+GROUND\s+"
                  r"STATE\s+CLUSTER\s+AMPLITUDES\s*\*\s*$",
                  name="CC ground state",
                  subMatchers=[
                      scf_cycle,
                      self.__context["method"].build_correlation_energy_matcher()
                  ],
                  onClose={None: write_data}
                  )
