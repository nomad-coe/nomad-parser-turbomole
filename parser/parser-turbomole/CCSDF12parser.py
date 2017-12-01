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

    def purge_data(self):
        self.__previous_energy = None

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
            self.__build_cc_scf_matcher(),
            self.__context.build_end_time_matcher("ccsdf12")
        ]

        return self.__context.build_module_matcher("ccsdf12", sub_matchers, "CCSDF12")

    def __build_cc_scf_matcher(self):
        correlation = SM(r"\s*\*\s*correlation\s+energy\s*:\s*"
                         r"(?P<energy_current__hartree>"+RE_FLOAT+r")\s*\*\s*",
                         name="correlation energy"
                         )
        total_energy = SM(r"\s*\*\s*Final\s+[^\s]+\s+energy\s*:\s*"
                          r"(?P<energy_total__hartree>"+RE_FLOAT+r")\s*\*\s*",
                          name="total energy"
                          )
        reference_energy = SM(r"\s*\*\s*(?:RHF|UHF|ROHF)\s+energy\s*:\s*"+RE_FLOAT+r"\s*\*\s*$",
                              name="HF energy",
                              )

        def compute_energy_change(backend, groups):
            if self.__previous_energy:
                backend.addRealValue("energy_change_scf_iteration",
                                     float(groups[0])-self.__previous_energy, unit="hartree")
            self.__previous_energy = float(groups[0])
        scf_iteration = SM(r"\s*[0-9]+\s+(?P<energy_total_scf_iteration__hartree>"+RE_FLOAT+")"
                           + 5 * (r"\s+"+RE_FLOAT) + "\s*$",
                           name="iteration",
                           repeats=True,
                           sections=["section_scf_iteration"],
                           startReAction=compute_energy_change
                           )
        scf_converged = SM(r"\s*CC\s+equations\s+converged\s+in\s+"
                           r"(?P<number_of_scf_iterations>[0-9]+)\s+iterations\.\s*$",
                           name="CC SCF converged")

        scf_cycle = SM(r"\s*Iter\.\s+CCSD\s+energy\s+Norm\(Omega\)\s+"
                       r"Norm\(t1\)\s+Norm\(t2\)\s+cpu\s+wall",
                       name="CCSD SCF iterations",
                       subMatchers=[
                           scf_iteration,
                           scf_converged
                       ])

        return SM(r"\s*\*\s*OPTIMIZATION\s+OF\s+THE\s+GROUND\s+"
                  r"STATE\s+CLUSTER\s+AMPLITUDES\s*\*\s*$",
                  name="CC ground state",
                  subMatchers=[
                      scf_cycle,
                      reference_energy,
                      correlation,
                      total_energy
                  ]
                  )