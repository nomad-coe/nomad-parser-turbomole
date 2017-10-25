import logging
import numpy as np
from nomadcore.simple_parser import SimpleMatcher as SM
from TurbomoleCommon import RE_FLOAT

logger = logging.getLogger("nomad.turbomoleParser")


class _EigenValue(object):

    def __init__(self, index):
        self.index = index
        self.eigenvalue = float("nan")
        self.occupation = 0.0


class OrbitalParser(object):

    def __init__(self, context, key="orbitals"):
        context[key] = self
        self.__context = context
        self.__backend = None
        self.__index_eigenvalues = 0
        self.__num_orbitals = 0
        self.__current_spin_channel = 0
        self.__last_row_offset = -5
        self.__eigenstates = {(0, 0): list()}

    def set_backend(self, backend):
        self.__backend = backend

    # getter functions

    def index_eigenvalues(self):
        return self.__index_eigenvalues

    # parsing functions

    def __write_eigenvalues(self, backend, gIndex, section):
        num_eigenvalues = max(len(x) for x in self.__eigenstates.values())
        num_kpoints = self.__context["method"].num_kpoints()
        num_spin = self.__context["method"].spin_channels()
        eigenvalue_type = "normal" if num_eigenvalues == self.__num_orbitals else "partial"
        eigenvalues = np.ndarray(shape=(num_spin, num_kpoints, num_eigenvalues), dtype=float)
        occupation = np.ndarray(shape=(num_spin, num_kpoints, num_eigenvalues), dtype=float)
        if num_kpoints > 1:
            logger.error("no support for multi k-point eigenvalues implemented, skipping!")
            return
        for (spin, k_point), state_list in self.__eigenstates.items():
            for i, state in enumerate(state_list):
                eigenvalues[spin, k_point, i] = state.eigenvalue
                occupation[spin, k_point, i] = state.occupation
        # TODO: add irreducible representation information once publicly available
        self.__index_eigenvalues = self.__backend.openSection("section_eigenvalues")
        self.__backend.addValue("number_of_eigenvalues", num_eigenvalues)
        self.__backend.addValue("number_of_eigenvalues_kpoints", num_kpoints)
        self.__backend.addValue("eigenvalues_kind", eigenvalue_type)
        self.__backend.addArrayValues("eigenvalues_values", eigenvalues, unit="hartree")
        self.__backend.addArrayValues("eigenvalues_occupation", occupation)
        self.__backend.closeSection("section_eigenvalues", self.__index_eigenvalues)

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
            eigenstates = self.__eigenstates[self.__current_spin_channel, 0]
            for state in groups:
                if state:
                    eigenstates.append(_EigenValue(state))
            self.__last_row_offset += 5

        def extract_eigenvalues(backend, groups):
            eigenstates = self.__eigenstates[self.__current_spin_channel, 0]
            for i, eigenvalue in enumerate(groups):
                if eigenvalue:
                    eigenstates[self.__last_row_offset+i].eigenvalue = float(eigenvalue)

        def extract_occupation(backend, groups):
            eigenstates = self.__eigenstates[self.__current_spin_channel, 0]
            for i, occupation in enumerate(groups):
                if occupation:
                    eigenstates[self.__last_row_offset+i].occupation = float(occupation)

        def next_spin_channel(backend, groups):
            self.__current_spin_channel += 1
            self.__eigenstates[self.__current_spin_channel, 0] = list()
            self.__last_row_offset = -5

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
                      startReAction=extract_states,
                      subMatchers=[
                          eigenvals_hartree,
                          eigenvals_ev,
                          occupations
                      ]
                      )

        subfiles = SM(r"\s*orbitals\s+\$([A-z]+)\s+will be written to file ([A-z]+)",
                      name="MOs to file",
                      repeats=True,
                      startReAction=process_mo_file
                      )

        return SM(r"\s*orbitals\s+\$([A-z]+)\s+will be written to file ([A-z]+)",
                  name="MOs to file",
                  startReAction=process_mo_file,
                  sections=["section_eigenvalues"],
                  onClose={"section_eigenvalues": self.__write_eigenvalues},
                  subMatchers=[
                      subfiles,
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
