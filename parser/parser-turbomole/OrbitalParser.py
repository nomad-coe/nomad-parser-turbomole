import logging
import numpy as np
from nomadcore.simple_parser import SimpleMatcher as SM
from TurbomoleCommon import RE_FLOAT

logger = logging.getLogger("nomad.turbomoleParser")


class OrbitalParser(object):

    def __init__(self, context, key="orbitals"):
        context[key] = self
        self.__context = context
        self.__backend = None
        self.__num_orbitals = 0
        self.__current_spin_channel = 0

    def set_backend(self, backend):
        self.__backend = backend

    def build_ir_rep_matcher(self):

        def sum_orbitals(backend, groups):
            self.__num_orbitals += int(groups[1])

        # irreducible representation name, num orbitals, num occupied
        ir_rep = SM(r"\s*([a-z]+[0-9]*)\s+([0-9]+)\s+([0-9]+)",
                    name="irRep info",
                    repeats=True,
                    startReAction=sum_orbitals
                    )
        return SM(r"\s*mo occupation :",
                  name="MO distribution",
                  subMatchers=[
                      SM(r"\s*irrep\s+mo's\s+occupied",
                         name="MO distribution",
                         required=True),
                      ir_rep
                  ]
                  )

    def build_eigenstate_matcher(self):

        def process_mo_file(backend, groups):
            pass

        def extract_states(backend, groups):
            pass

        def extract_eigenvalues(backend, groups):
            pass


        def extract_occupation(backend, groups):
            pass

        def next_spin_channel(backend, groups):
            self.__current_spin_channel += 1

        def states():
            eigenvals_hartree = SM(r"\s*eigenvalues H\s+("+RE_FLOAT+")"
                                   + 4 * (r"(?:\s+("+RE_FLOAT+"))?"),
                                   name="Hartree eigenvalues",
                                   startReAction=extract_eigenvalues,
                                   required=True
                                   )
            eigenvals_ev = SM(r"\s*eV\s+("+RE_FLOAT+")" + 4 * (r"(?:\s+("+RE_FLOAT+"))?"),
                              name="eV eigenvalues",
                              coverageIgnore=True
                              )
            occupations = SM(r"\s*occupation\s+("+RE_FLOAT+")" + 4 * (r"(?:\s+("+RE_FLOAT+"))?"),
                             name="occupation",
                             startReAction=extract_occupation
                             )
            return SM(r"\s*irrep\s+([0-9]+[a-z][1-9'\"]?)" + 4 * r"(?:\s+([0-9]+[a-z][1-9'\"]?))?",
                      name="irRep list",
                      required=True,
                      repeats=True,
                      subMatchers=[
                          eigenvals_hartree,
                          eigenvals_ev,
                          occupations
                      ]
                      )

        return SM(r"\s*orbitals\s+\$([A-z]+)\s+will be written to file ([A-z]+)",
                  name="MOs to file",
                  repeats=True,
                  startReAction=process_mo_file,
                  subMatchers=[
                      SM(r"\s*alpha\s*:\s*$",
                         name="alpha spin"
                         ),
                      states(),
                      SM(r"\s*beta\s*:\s*$",
                         name="beta spin",
                         startReAction=next_spin_channel
                         ),
                      states()
                  ]
                  )
