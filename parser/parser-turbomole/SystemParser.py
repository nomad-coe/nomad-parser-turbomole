import logging
import numpy as np
from nomadcore.simple_parser import SimpleMatcher as SM
import nomadcore.elements as elements
from TurbomoleCommon import RE_FLOAT

logger = logging.getLogger("nomad.turbomoleParser")


class Atom(object):

    def __init__(self, x, y, z, elem, charge, shells, isotope, pseudo):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        self.elem = elem.capitalize()
        self.charge = float(charge)
        self.shells = int(shells) if shells else -1
        self.isotope = int(isotope)
        if pseudo:
            self.is_pseudo = True if pseudo == "1" else False
        else:
            self.is_pseudo = None  # effective core potential not specified
        self.label = self.elem if self.shells > 0 else self.elem+"_1"


class BasisSet(object):

    def __init__(self, name, index, cartesian, num_atoms):
        self.name = name
        self.index = index
        self.cartesian = cartesian
        self.num_atoms = int(num_atoms)


class SystemParser(object):

    def __init__(self, context, key="geo"):
        context[key] = self
        self.__context = context
        self.__backend = None
        self.__index_peecm_unit_cell = -1
        self.__index_peecm_pc_cluster = -1
        self.__index_peecm_qm_cluster = -1
        self.__index_basis_set = -1
        self.__atoms = list()
        self.__basis_sets = dict()
        self.__auxbasis_sets = dict()
        self.__pceem_parameters = dict()

    def purge_data(self):
        self.__index_peecm_unit_cell = -1
        self.__index_peecm_pc_cluster = -1
        self.__index_peecm_qm_cluster = -1
        self.__index_basis_set = -1
        self.__atoms = list()
        self.__basis_sets = dict()
        self.__auxbasis_sets = dict()
        self.__pceem_parameters = dict()

    def set_backend(self, backend):
        self.__backend = backend

    # getter methods

    def index_basis_set(self):
        return self.__index_basis_set

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

    def write_basis_set_mapping(self):
        """the caller is responsible for opening the enclosing
        section_single_configuration_calculation"""
        showed_warning = False
        self.__index_basis_set = self.__backend.openSection("section_basis_set")
        mapping = np.ndarray(shape=(len(self.__atoms),), dtype=int)
        total_basis_atoms = sum(x.num_atoms for x in self.__basis_sets.values())
        index_pseudo = None
        if total_basis_atoms < len(self.__atoms):
            # some atoms have no basis set assigned, need to create an empty set
            index_pseudo = self.__backend.openSection("section_basis_set_atom_centered")
            self.__backend.addValue("basis_set_atom_centered_short_name", "empty", index_pseudo)
            self.__backend.addValue("basis_set_atom_number", 0, index_pseudo)
            self.__backend.addValue("number_of_basis_functions_in_basis_set_atom_centered", 0,
                                    index_pseudo)
            self.__backend.closeSection("section_basis_set_atom_centered", index_pseudo)
        kind_counts = dict()
        # TODO: if the geometry data doesn't list shells, a simple ordering assumption is used...
        for i, atom in enumerate(self.__atoms):
            if atom.shells == 0:
                mapping[i] = index_pseudo
            elif atom.elem in kind_counts:
                if kind_counts[atom.elem] >= self.__basis_sets[atom.elem].num_atoms:
                    if not showed_warning:
                        logger.warning("guessing basis set to atom map from ordering!")
                        showed_warning = True
                    mapping[i] = index_pseudo
                else:
                    kind_counts[atom.elem] = kind_counts[atom.elem] + 1
                    mapping[i] = self.__basis_sets[atom.elem].index
            else:
                kind_counts[atom.elem] = 1
                mapping[i] = self.__basis_sets[atom.elem].index
        self.__backend.addArrayValues("mapping_section_basis_set_atom_centered", mapping)
        self.__backend.closeSection("section_basis_set", self.__index_basis_set)

    def write_method_basis_set_mapping(self):
        pass

    def build_qm_geometry_matcher(self):

        def add_atom(backend, groups):
            self.__atoms.append(Atom(x=groups[0], y=groups[1], z=groups[2], elem=groups[3],
                                     charge=groups[5], isotope=groups[7], shells=groups[4],
                                     pseudo=groups[6]
                                     )
                                )

        def finalize_data(backend, groups):
            pos = np.ndarray(shape=(len(self.__atoms), 3), dtype=float)
            labels = list()
            atom_numbers = np.ndarray(shape=(len(self.__atoms),), dtype=float)
            backend.addValue("number_of_atoms", len(self.__atoms))
            for i, atom in enumerate(self.__atoms):
                pos[i, 0:3] = (atom.x, atom.y, atom.z)
                labels.append(atom.elem)
                atom_numbers[i] = elements.get_atom_number(atom.elem)
            backend.addArrayValues("atom_positions", pos, unit="angstrom")
            backend.addArrayValues("atom_labels", np.asarray(labels, dtype=str))
            backend.addArrayValues("atom_atom_number", atom_numbers)

        # x, y, z, element, (shells), charge, (pseudo), isotope
        atom_re = r"\s*([-+]?[0-9]+\.?[0-9]*)" \
                  r"\s+([-+]?[0-9]+\.?[0-9]*)" \
                  r"\s+([-+]?[0-9]+\.?[0-9]*)" \
                  r"\s+([a-zA-Z]+)" \
                  r"(\s+[0-9]+)?" \
                  r"\s+([-+]?[0-9]+\.?[0-9]*)" \
                  r"(\s+[-+]?[0-9]+)?" \
                  r"\s+([-+]?[0-9]+)"
        atom = SM(atom_re, repeats=True, name="single atom", startReAction=add_atom)
        header_re = r"\s*atomic\s+coordinates\s+atom(?:\s+shells)?\s+charge(?:\s+pseudo)?\s+isotop"

        return SM(name="geometry",
                  startReStr=r"\s*\|\s+Atomic coordinate, charge and isotope? information\s+\|",
                  subMatchers=[
                      SM(r"\s*-{20}-*", name="<format>", coverageIgnore=True),
                      SM(header_re, name="atom list", subMatchers=[atom]),
                      SM("\s*center of nuclear mass", startReAction=finalize_data)
                  ]
                  )

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
            embedding_atoms.append(Atom(x=groups[1], y=groups[2], z=groups[3], elem=groups[0],
                                        charge=groups[4], isotope=0, shells=None, pseudo=None
                                        )
                                   )

        def add_removed_point_charge(backend, groups):
            embedding_atoms.append(Atom(x=groups[1], y=groups[2], z=groups[3], elem=groups[0],
                                        charge=0, isotope=0, shells=None, pseudo=None
                                        )
                                   )

        def add_shifted_qm_atom(backend, groups):
            embedding_atoms.append(Atom(x=groups[1], y=groups[2], z=groups[3], elem=groups[0],
                                        charge=0, isotope=0, shells=None, pseudo=None
                                        )
                                   )

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

    def build_orbital_basis_matcher(self, title_regex=None):
        if not title_regex:
            title_regex = r"\s*\|\s*basis\s+set\s+information\s*\|"
        return self.__build_basis_matcher(False, title_regex)

    def build_auxiliary_basis_matcher(self, title_regex=None):
        if not title_regex:
            title_regex = r"\s*\|\s*Auxiliary\s+basis\s+set\s+information\s*\|"
        return self.__build_basis_matcher(True, title_regex)

    def __build_basis_matcher(self, is_auxbasis, title_regex):
        basis_type = "Auxiliary Basis" if is_auxbasis else "Orbital Basis"

        class LocalBasisData(object):
            spherical = False

        def set_spherical_basis(backend, groups):
            LocalBasisData.spherical = True

        def add_basis_set(backend, groups):
            # TODO: support assignment of different basis sets to atoms of the same element
            index = backend.openSection("section_basis_set_atom_centered")
            key = groups[0].capitalize()
            basis_set = BasisSet(name=groups[4], index=index,
                                 cartesian=not LocalBasisData.spherical, num_atoms=groups[1])
            if is_auxbasis:
                self.__auxbasis_sets[key] = basis_set
            else:
                self.__basis_sets[key] = basis_set
            atom_number = elements.get_atom_number(groups[0].capitalize())
            backend.addValue("basis_set_atom_centered_short_name", groups[4], index)
            backend.addValue("basis_set_atom_number", atom_number, index)
            if LocalBasisData.spherical:
                backend.addValue("number_of_basis_functions_in_basis_set_atom_centered",
                                 int(groups[3]), index)
            else:
                logger.warning("no basis function count information for cartesian Gaussians!")
            backend.closeSection("section_basis_set_atom_centered", index)

        # elem, num atoms, prim gauss, contracted gauss, name, contracted details, prim details
        basis_re = r"\s*([A-z]+)\s+([0-9]+)\s+" \
                   r"([0-9]+)\s+([0-9]+)\s+" \
                   r"([A-z0-9-]+)\s+" \
                   r"\[((?:[0-9]+[spdfghij])+)\|((?:[0-9]+[spdfghij])+)\]"
        basis = SM(basis_re, repeats=True, name="basis assignment", startReAction=add_basis_set)
        header_re = r"\s*we\s+will\s+work\s+with\s+the\s+1s\s+3p\s+5d\s+7f\s+9g\s+...\s+basis\s+set"
        gauss_type_spherical = SM(header_re,
                                  name="spherical Gaussians",
                                  subMatchers=[
                                      SM(r"\s*...i.e. with spherical basis functions...",
                                         name="spherical Gaussians",
                                         required=True)
                                  ],
                                  startReAction=set_spherical_basis
                                  )

        return SM(title_regex,
                  name=basis_type,
                  subMatchers=[
                      SM(r"\s*\+----*\+", name="<format>", coverageIgnore=True),
                      gauss_type_spherical,
                      SM(r"\s*type\s+atoms\s+prim\s+cont\s+basis",
                         name=basis_type),
                      SM(r"\s*----*", name="<format>", coverageIgnore=True),
                      basis,
                      SM(r"\s*----*", name="<format>", coverageIgnore=True),
                      SM(r"\s*total:\s*([0-9]+)\s+([0-9]+)\s+([0-9]+)",
                         name=basis_type),
                      SM(r"\s*----*", name="<format>", coverageIgnore=True),
                      SM(r"\s*total number of primitive shells\s*:\s*([0-9]+)",
                         name=basis_type),
                      SM(r"\s*total number of contracted shells\s*:\s*([0-9]+)",
                         name=basis_type),
                      SM(r"\s*total number of cartesian basis functions\s*:\s*([0-9]+)",
                         name=basis_type),
                      SM(r"\s*total number of SCF-basis functions\s*:\s*([0-9]+)",
                         name=basis_type),
                  ]
                  )
