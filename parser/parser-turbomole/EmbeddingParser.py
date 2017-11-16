import logging
import numpy as np
from nomadcore.simple_parser import SimpleMatcher as SM
import nomadcore.elements as elements
from TurbomoleCommon import RE_FLOAT

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
        self.__index_peecm_unit_cell = -1
        self.__index_peecm_pc_cluster = -1
        self.__index_peecm_qm_cluster = -1
        self.__pceem_parameters = dict()

    def purge_data(self):
        self.__index_peecm_unit_cell = -1
        self.__index_peecm_pc_cluster = -1
        self.__index_peecm_qm_cluster = -1
        self.__pceem_parameters = dict()

    def set_backend(self, backend):
        self.__backend = backend

    # match builders

    def link_embedding_systems_to_qm(self, qm_geo_index):
        references = {"section_system": qm_geo_index}
        if self.__index_peecm_unit_cell != -1:
            index = self.__backend.openSection("section_system_to_system_refs")
            self.__backend.addValue("system_to_system_kind", "periodic point charges for embedding")
            self.__backend.addValue("system_to_system_ref", self.__index_peecm_unit_cell)
            self.__backend.setSectionInfo("section_system_to_system_refs", index, references)
            self.__backend.closeSection("section_system_to_system_refs", index)
        if self.__index_peecm_pc_cluster != -1:
            index = self.__backend.openSection("section_system_to_system_refs")
            self.__backend.addValue("system_to_system_kind", "removed point charge cluster")
            self.__backend.addValue("system_to_system_ref", self.__index_peecm_pc_cluster)
            self.__backend.setSectionInfo("section_system_to_system_refs", index, references)
            self.__backend.closeSection("section_system_to_system_refs", index)
        if self.__index_peecm_qm_cluster != -1:
            index = self.__backend.openSection("section_system_to_system_refs")
            self.__backend.addValue("system_to_system_kind", "shifted embedded QM cluster")
            self.__backend.addValue("system_to_system_ref", self.__index_peecm_qm_cluster)
            self.__backend.setSectionInfo("section_system_to_system_refs", index, references)
            self.__backend.closeSection("section_system_to_system_refs", index)
        if self.__index_peecm_unit_cell != -1:
            self.__backend.addValue("embedded_system", True, qm_geo_index)
        for key, value in self.__pceem_parameters.items():
            self.__backend.addValue(key, value, qm_geo_index)

    def build_embedding_matcher(self):
        embedding_atoms = list()
        lattice_vectors = list()

        def store_pc_cell_index(backend, gIndex, section):
            del embedding_atoms[:]
            self.__index_peecm_unit_cell = gIndex
            if len(lattice_vectors) != 3:
                logger.error("didn't get expected 3 lattice vectors for PCEEM embedding cell:"
                             + str(len(lattice_vectors)))
            else:
                lattice = np.ndarray(shape=(3, 3))
                for i, vector in enumerate(lattice_vectors[0:3]):
                    lattice[i, :] = vector[:]
                backend.addArrayValues("lattice_vectors", lattice, gIndex, unit="bohr")

        def store_pc_cluster_index(backend, gIndex, section):
            del embedding_atoms[:]
            self.__index_peecm_pc_cluster = gIndex

        def store_qm_cluster_index(backend, gIndex, section):
            del embedding_atoms[:]
            self.__index_peecm_qm_cluster = gIndex

        def write_data(backend, gIndex, section):
            pos = np.ndarray(shape=(len(embedding_atoms), 3), dtype=float)
            labels = list()
            backend.addValue("number_of_atoms", len(embedding_atoms))
            atom_numbers = np.ndarray(shape=(len(embedding_atoms),), dtype=float)
            charges = np.ndarray(shape=(len(embedding_atoms),), dtype=float)
            for i, atom in enumerate(embedding_atoms):
                pos[i, 0:3] = (atom.x, atom.y, atom.z)
                labels.append(atom.elem)
                atom_numbers[i] = elements.get_atom_number(atom.elem)
                charges[i] = atom.charge
            backend.addArrayValues("atom_positions", pos, unit="angstrom")
            backend.addArrayValues("atom_labels", np.asarray(labels, dtype=str))
            backend.addArrayValues("atom_atom_number", atom_numbers)
            if self.__index_peecm_unit_cell == gIndex:
                backend.addArrayValues("x_turbomole_pceem_charges", charges)

        def add_point_charge(backend, groups):
            new_atom = Atom(x=groups[1], y=groups[2], z=groups[3], elem=groups[0], charge=groups[4])
            embedding_atoms.append(new_atom)

        def add_removed_point_charge(backend, groups):
            new_atom = Atom(x=groups[1], y=groups[2], z=groups[3], elem=groups[0], charge=0)
            embedding_atoms.append(new_atom)

        def add_shifted_qm_atom(backend, groups):
            new_atom = Atom(x=groups[1], y=groups[2], z=groups[3], elem=groups[0], charge=0)
            embedding_atoms.append(new_atom)

        def add_lattice_vector(backend, groups):
            lattice_vectors.append(tuple(float(x) for x in groups[0:3]))

        lattice_vector = SM(r"\s*("+RE_FLOAT+")"+ 2 * ("\s+("+RE_FLOAT+")")+"\s*$",
                            name="PCEEM lattice",
                            repeats=True,
                            startReAction=add_lattice_vector,
                            required=True
                            )

        lattice_setup = SM(r"\s*Cell vectors \(au\):\s*$",
                           name="PCEEM lattice",
                           required=True,
                           subMatchers=[
                               lattice_vector
                           ]
                           )

        point_charge_in_unit_cell = SM(r"\s*([A-z]+)"+4*("\s+("+RE_FLOAT+")")+"\s*$",
                                       name="point charge",
                                       repeats=True,
                                       startReAction=add_point_charge,
                                       required=True)
        point_charge_cell = SM(r"\s*Redefined unit cell content \(au\):\s*$",
                               name="header",
                               sections=["section_system"],
                               subMatchers=[
                                   SM(r"\s*Label\s+Cartesian\s+Coordinates\s+Charge\s*$",
                                      name="header"),
                                   point_charge_in_unit_cell
                               ],
                               onOpen={"section_system": store_pc_cell_index},
                               onClose={"section_system": write_data}
                               )

        removed_point_charge = SM(r"\s*([A-z]+)"+6*("\s+("+RE_FLOAT+")")+"\s*$",
                                  name="removed point charge",
                                  repeats=True,
                                  startReAction=add_removed_point_charge,
                                  required=True)
        point_charge_cluster = SM(r"\s*PC cluster transformed to the center of cell 0 \(au\):\s*$",
                                  name="header",
                                  sections=["section_system"],
                                  subMatchers=[
                                      SM(r"\s*Label\s+Cartesian\s+Coordinates\s+Cell\s+Indices\s*$",
                                         name="header"),
                                      removed_point_charge
                                  ],
                                  onOpen={"section_system": store_pc_cluster_index},
                                  onClose={"section_system": write_data}
                                  )

        shifted_qm_atom = SM(r"\s*([A-z]+)"+3*("\s+("+RE_FLOAT+")")+"\s*$",
                             name="shifted QM atom",
                             repeats=True,
                             startReAction=add_shifted_qm_atom,
                             required=True)
        shifted_qm_cluster = SM(r"\s*QM cluster transformed to the center of cell 0 \(au\):\s*$",
                                name="header",
                                sections=["section_system"],
                                subMatchers=[
                                    SM(r"\s*Atom\s+Cartesian\s+Coordinates\s*$",
                                       name="header"),
                                    shifted_qm_atom
                                ],
                                onOpen={"section_system": store_qm_cluster_index},
                                onClose={"section_system": write_data}
                                )

        def store_max_multipole(backend, groups):
            self.__pceem_parameters["x_turbomole_pceem_max_multipole"] = int(groups[0])

        def store_multipole_precision(backend, groups):
            self.__pceem_parameters["x_turbomole_pceem_multipole_precision"] = float(groups[0])

        def store_cell_separation(backend, groups):
            self.__pceem_parameters["x_turbomole_pceem_min_separation_cells"] = float(groups[0])

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
                      lattice_setup,
                      point_charge_cell,
                      point_charge_cluster,
                      shifted_qm_cluster
                  ]
                  )