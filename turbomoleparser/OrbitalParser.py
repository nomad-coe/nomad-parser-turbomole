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

import logging
import numpy as np
import os.path
import re
from functools import total_ordering
from nomadcore.simple_parser import SimpleMatcher as SM
from nomadcore.simple_parser import AncillaryParser
from turbomoleparser.TurbomoleCommon import RE_FLOAT

logger = logging.getLogger("nomad.turbomoleParser")


@total_ordering
class _EigenState(object):

    def __init__(self, index, ir_rep, spin, k_point):
        self.index = index
        self.eigenvalue = float("nan")
        self.occupation = -1.0
        self.spin_channel = spin
        self.k_point = k_point
        self.irreducible_representation = ir_rep

    def __eq__(self, other):
        if isinstance(other, _EigenState):
            return self.irreducible_representation == other.irreducible_representation \
                   and self.spin_channel == other.spin_channel \
                   and self.k_point == other.k_point and self.index == other.index
        else:
            return False

    def __lt__(self, other):
        if isinstance(other, _EigenState):
            if self.irreducible_representation < other.irreducible_representation:
                return True
            elif self.spin_channel < other.spin_channel:
                return True
            elif self.k_point < other.k_point:
                return True
            elif self.index < other.index:
                return True
        return False


class OrbitalParser(object):

    def __init__(self, context, key="orbitals"):
        context[key] = self
        self.__context = context
        self.__backend = None
        self.__index_eigenvalues = 0
        self.__num_orbitals = 0
        self.__current_spin_channel = -1
        self.__current_k_point = 0
        self.__eigenstates = list()
        self.__eigenstates_staging = list()

    def purge_data(self):
        self.__index_eigenvalues = 0
        self.__num_orbitals = 0
        self.__current_spin_channel = -1
        self.__current_k_point = 0
        self.__eigenstates = list()
        self.__eigenstates_staging = list()

    def set_backend(self, backend):
        self.__backend = backend

    # getter functions

    def index_eigenvalues(self):
        return self.__index_eigenvalues

    # parsing functions

    def __write_eigenvalues(self, backend, gIndex, section):
        num_k_points = self.__context["method"].k_points
        num_spin = self.__context["method"].spin_channels
        num_eigenvalues = 0
        for spin in range(0, num_spin):
            for k_point in range(0, num_k_points):
                subset = [state for state in self.__eigenstates if state.spin_channel == spin and
                          state.k_point == k_point]
                num_eigenvalues = max(num_eigenvalues, len(subset))
        eigenvalue_type = "normal" if num_eigenvalues == self.__num_orbitals else "partial"
        eigenvalues = np.ndarray(shape=(num_spin, num_k_points, num_eigenvalues), dtype=float)
        eigenvalues[:, :, :] = 0.0
        occupation = np.ndarray(shape=(num_spin, num_k_points, num_eigenvalues), dtype=float)
        occupation[:, :, :] = 0.0
        representation = np.ndarray(shape=(num_spin, num_k_points, num_eigenvalues), dtype="U2")
        representation[:, :, :] = "  "
        if num_k_points > 1:
            logger.error("no support for multi k-point eigenvalues implemented, skipping!")
            return
        for spin in range(0, num_spin):
            for k_point in range(0, num_k_points):
                subset = list(state for state in self.__eigenstates if state.spin_channel == spin
                              and state.k_point == k_point)
                ir_reps = set(state.irreducible_representation for state in subset)
                offset = 0
                for ir_rep in sorted(ir_reps):
                    ir_rep_states = sorted([state for state in subset
                                            if state.irreducible_representation == ir_rep])
                    highest_occupied = max(x.index for x in ir_rep_states if x.occupation > 0.0)
                    max_occupation = max(x.occupation for x in ir_rep_states if x.occupation > 0.0)
                    for i, state in enumerate(ir_rep_states):
                        eigenvalues[spin, k_point, offset + i] = state.eigenvalue
                        if state.occupation <= 0.0:
                            state.occupation = max_occupation if state.index <= highest_occupied \
                                else 0.0
                        occupation[spin, k_point, offset + i] = state.occupation
                        representation[spin, k_point, offset + i] = state.irreducible_representation
                    offset += len(ir_rep_states)
        # TODO: move irreducible representation information to public meta data, once available
        self.__index_eigenvalues = self.__backend.openSection("section_eigenvalues")
        self.__backend.addValue("number_of_eigenvalues", num_eigenvalues)
        self.__backend.addValue("number_of_eigenvalues_kpoints", num_k_points)
        self.__backend.addValue("eigenvalues_kind", eigenvalue_type)
        self.__backend.addArrayValues("eigenvalues_values", eigenvalues, unit="hartree")
        self.__backend.addArrayValues("eigenvalues_occupation", occupation)
        self.__backend.addArrayValues("x_turbomole_eigenvalues_irreducible_representation",
                                      representation)
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

        re_state_split = re.compile(r"([0-9]+)([a-z][0-9'\"]?)")

        def process_mo_file(backend, groups):
            self.__current_spin_channel += 1
            states_parser = self.__build_eigenstate_file_parser()
            dir_name = os.path.dirname(os.path.abspath(self.__context.fName))
            file_name = os.path.normpath(os.path.join(dir_name, groups[1]))
            try:
                with open(file_name) as file_in:
                    states_parser.parseFile(file_in)
            except IOError:
                logger.warning("Could not find auxiliary file '%s' in directory '%s'." % (
                    groups[1], dir_name))

        def reset_spin_count(backend, groups):
            self.__current_spin_channel = 0

        def extract_states(backend, groups):
            self.__eigenstates_staging = list()
            for state in groups:
                if state:
                    match = re_state_split.match(state)
                    index = int(match.group(1))
                    ir_rep = match.group(2)
                    to_add = _EigenState(index, ir_rep, self.__current_spin_channel,
                                         self.__current_k_point)
                    if to_add in self.__eigenstates:
                        to_update = self.__eigenstates[self.__eigenstates.index(to_add)]
                        self.__eigenstates_staging.append(to_update)
                    else:
                        self.__eigenstates_staging.append(to_add)
                        self.__eigenstates.append(to_add)

        def extract_eigenvalues(backend, groups):
            for i, eigenvalue in enumerate(groups):
                if eigenvalue:
                    self.__eigenstates_staging[i].eigenvalue = float(eigenvalue)

        def extract_occupation(backend, groups):
            for i, occupation in enumerate(groups):
                if occupation:
                    self.__eigenstates_staging[i].occupation = float(occupation)

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
            return SM(r"\s*irrep\s+([0-9]+[a-z][0-9'\"]?)" + 4 * r"(?:\s+([0-9]+[a-z][0-9'\"]?))?",
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
                         name="alpha spin",
                         startReAction=reset_spin_count
                         ),
                      states(),
                      SM(r"\s*beta\s*:\s*$",
                         name="beta spin",
                         startReAction=next_spin_channel
                         ),
                      states()
                  ]
                  )

    def __build_eigenstate_file_parser(self):
        def add_eigenstate(backend, groups):
            to_add = _EigenState(int(groups[0]), groups[1], self.__current_spin_channel,
                                 self.__current_k_point)
            to_add.eigenvalue = float(groups[2].replace("D", "E").replace("d", "e"))
            if to_add in self.__eigenstates:
                logger.warning("eigenstate from file already known: " + str(to_add.__dict__))
            else:
                self.__eigenstates.append(to_add)

        coefficients = SM(r"\s*"+4*("(?:"+RE_FLOAT+")?")+r"\s*$",
                          name="coefficients",
                          repeats=True,
                          coverageIgnore=True)
        state = SM(r"\s*([0-9]+)\s+([a-z'\"])\s+eigenvalue=("+RE_FLOAT+")\s+nsaos=([0-9]+)\s*$",
                   name="eigenstate",
                   repeats=True,
                   startReAction=add_eigenstate,
                   subMatchers=[
                       coefficients
                   ]
                   )

        return AncillaryParser(fileDescription=state,
                               parser=self.__context.parser,
                               cachingLevelForMetaName=dict(),
                               superContext=self.__context
                               )
