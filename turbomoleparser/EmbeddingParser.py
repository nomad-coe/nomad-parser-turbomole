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
from nomadcore.simple_parser import SimpleMatcher as SM
import nomadcore.elements as elements
from turbomoleparser.TurbomoleCommon import RE_FLOAT

logger = logging.getLogger("nomad.turbomoleParser")


class Atom(object):

    def __init__(self, x, y, z, elem, charge):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        self.elem = elem.capitalize()
        self.charge = float(charge)
        self.label = self.elem


class EmbeddingParser(object):

    def __init__(self, context, key="embedding"):
        context[key] = self
        self.__context = context
        self.__backend = None

    def purge_data(self):
        pass

    def set_backend(self, backend):
        self.__backend = backend

    # match builders

    def build_embedding_matcher(self):

        def store_max_multipole(backend, groups):
            index = self.__context.index_system()
            backend.addValue("x_turbomole_pceem_max_multipole", int(groups[0]), index)

        def store_multipole_precision(backend, groups):
            index = self.__context.index_system()
            backend.addValue("x_turbomole_pceem_multipole_precision", float(groups[0]), index)

        def store_cell_separation(backend, groups):
            index = self.__context.index_system()
            backend.addValue("x_turbomole_pceem_min_separation_cells", float(groups[0]), index)

        max_multipole = SM(r"\s*Maximum multipole moment used\s*:\s*([0-9]+)\s*$",
                           name="max multipole",
                           startReAction=store_max_multipole
                           )
        multipole_precision = SM(r"\s*Multipole precision parameter\s*:\s*("+RE_FLOAT+")\s*$",
                                 name="multipole precision",
                                 startReAction=store_multipole_precision
                                 )
        cell_separation = SM(r"\s*Minimum separation between cells\s*:\s*("+RE_FLOAT+")\s*$",
                             name="cell seperation",
                             startReAction=store_cell_separation
                             )
        header = SM(r"\s*\+-+\s*Parameters\s*-*\+\s*$",
                    name="PCEEM parameters",
                    subMatchers=[
                        max_multipole,
                        multipole_precision,
                        cell_separation
                    ]
                    )

        return SM(r"\s*\|\s*EMBEDDING IN PERIODIC POINT CHARGES\s*\|\s*$",
                  name="embedding (PEECM)",
                  subMatchers=[
                      SM(r"\s*\|\s*M. Sierka and A. Burow\s*\|\s*$", name="credits"),
                      header,
                      self.__build_pc_cell_matcher(),
                      self.__build_removed_pc_cluster_matcher(),
                      self.__build_shifted_qm_cluster_matcher()
                  ]
                  )

    def __write_embedding_atom_data(self, gIndex, atoms, linkd_description, add_charges):
        references = {"section_system": self.__context.index_system()}
        index = self.__backend.openSection("section_system_to_system_refs")
        self.__backend.addValue("system_to_system_kind", linkd_description)
        self.__backend.addValue("system_to_system_ref", gIndex)
        self.__backend.setSectionInfo("section_system_to_system_refs", index, references)
        self.__backend.closeSection("section_system_to_system_refs", index)

        pos = np.ndarray(shape=(len(atoms), 3), dtype=float)
        labels = list()
        self.__backend.addValue("number_of_atoms", len(atoms))
        atom_numbers = np.ndarray(shape=(len(atoms),), dtype=float)
        charges = np.ndarray(shape=(len(atoms),), dtype=float)
        for i, atom in enumerate(atoms):
            pos[i, 0:3] = (atom.x, atom.y, atom.z)
            labels.append(atom.elem)
            atom_numbers[i] = elements.get_atom_number(atom.elem)
            charges[i] = atom.charge
        self.__backend.addArrayValues("atom_positions", pos, unit="angstrom")
        self.__backend.addArrayValues("atom_labels", np.asarray(labels, dtype=str))
        self.__backend.addArrayValues("atom_atom_number", atom_numbers)
        if add_charges:
            self.__backend.addArrayValues("x_turbomole_pceem_charges", charges)

    def __build_pc_cell_matcher(self):
        lattice_vectors = list()
        embedding_atoms = list()

        def add_lattice_vector(backend, groups):
            lattice_vectors.append(tuple(float(x) for x in groups[0:3]))

        def add_point_charge(backend, groups):
            new_atom = Atom(x=groups[1], y=groups[2], z=groups[3], elem=groups[0], charge=groups[4])
            embedding_atoms.append(new_atom)

        def finalize(backend, gIndex, section):
            if len(lattice_vectors) != 3:
                logger.error("didn't get expected 3 lattice vectors for PCEEM embedding cell:"
                             + str(len(lattice_vectors)))
            else:
                lattice = np.ndarray(shape=(3, 3))
                for i, vector in enumerate(lattice_vectors[0:3]):
                    lattice[i, :] = vector[:]
                backend.addArrayValues("lattice_vectors", lattice, gIndex, unit="bohr")
            self.__write_embedding_atom_data(gIndex, embedding_atoms,
                                             "periodic point charges for embedding", True)

        lattice_vector = SM(r"\s*("+RE_FLOAT+")"+ 2 * ("\s+("+RE_FLOAT+")")+"\s*$",
                            name="PCEEM lattice",
                            repeats=True,
                            startReAction=add_lattice_vector,
                            required=True
                            )
        point_charge_in_unit_cell = SM(r"\s*([A-z]+)"+4*("\s+("+RE_FLOAT+")")+"\s*$",
                                       name="point charge",
                                       repeats=True,
                                       startReAction=add_point_charge,
                                       required=True)
        point_charge_cell = SM(r"\s*Redefined unit cell content \(au\):\s*$",
                               name="header",
                               subMatchers=[
                                   SM(r"\s*Label\s+Cartesian\s+Coordinates\s+Charge\s*$",
                                      name="header"),
                                   point_charge_in_unit_cell
                               ]
                               )
        return SM(r"\s*Cell vectors \(au\):\s*$",
                  name="PCEEM lattice",
                  sections=["section_system"],
                  required=True,
                  subMatchers=[
                      lattice_vector,
                      point_charge_cell
                  ],
                  onClose={"section_system": finalize}
                  )

    def __build_removed_pc_cluster_matcher(self):
        embedding_atoms = list()

        def add_removed_point_charge(backend, groups):
            new_atom = Atom(x=groups[1], y=groups[2], z=groups[3], elem=groups[0], charge=0)
            embedding_atoms.append(new_atom)

        def finalize(backend, gIndex, section):
            self.__write_embedding_atom_data(gIndex, embedding_atoms,
                                             "removed point charge cluster", False)

        removed_point_charge = SM(r"\s*([A-z]+)"+6*("\s+("+RE_FLOAT+")")+"\s*$",
                                  name="removed point charge",
                                  repeats=True,
                                  startReAction=add_removed_point_charge,
                                  required=True)
        return SM(r"\s*PC cluster transformed to the center of cell 0 \(au\):\s*$",
                  name="header",
                  sections=["section_system"],
                  subMatchers=[
                      SM(r"\s*Label\s+Cartesian\s+Coordinates\s+Cell\s+Indices\s*$", name="header"),
                      removed_point_charge
                  ],
                  onClose={"section_system": finalize}
                  )

    def __build_shifted_qm_cluster_matcher(self):
        embedding_atoms = list()

        def add_shifted_qm_atom(backend, groups):
            new_atom = Atom(x=groups[1], y=groups[2], z=groups[3], elem=groups[0], charge=0)
            embedding_atoms.append(new_atom)

        def finalize(backend, gIndex, section):
            self.__write_embedding_atom_data(gIndex, embedding_atoms,
                                             "shifted embedded QM cluster", False)

        shifted_qm_atom = SM(r"\s*([A-z]+)"+3*("\s+("+RE_FLOAT+")")+"\s*$",
                             name="shifted QM atom",
                             repeats=True,
                             startReAction=add_shifted_qm_atom,
                             required=True)
        return SM(r"\s*QM cluster transformed to the center of cell 0 \(au\):\s*$",
                  name="header",
                  sections=["section_system"],
                  subMatchers=[
                      SM(r"\s*Atom\s+Cartesian\s+Coordinates\s*$", name="header"),
                      shifted_qm_atom
                  ],
                  onClose={"section_system": finalize}
                  )
