import logging
import numpy as np
from nomadcore.simple_parser import SimpleMatcher as SM
import nomadcore.elements as elements

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
        self.is_pseudo = True if pseudo == "1" else False
        self.label = self.elem if self.shells > 0 else self.elem+"_1"

class BasisSet(object):

    def __init__(self, name, index, cartesian):
        self.name = name
        self.index = index
        self.cartesian = cartesian

class SystemParser(object):

    def __init__(self, context, key="geo"):
        context[key] = self
        self.__context = context
        self.__backend = None
        self.__index_qm_geo = -1
        self.__atoms = list()
        self.__basis_sets = dict()

    def set_backend(self, backend):
        self.__backend = backend

    def finalize_sections(self):
        if self.__index_qm_geo >= -1:
            self.__backend.closeSection("section_system", self.__index_qm_geo)
            self.__index_qm_geo = -2

    def build_qm_geometry_matcher(self):

        def open_section(backend, groups):
            self.__index_qm_geo = backend.openSection("section_system")

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
                  ],
                  startReAction=open_section
                  )


    def build_orbital_basis_matcher(self):

        class LocalBasisData(object):
            spherical = False

        def set_spherical_basis(backend, groups):
            LocalBasisData.spherical = True

        def add_basis_set(backend, groups):
            index = backend.openSection("section_basis_set_atom_centered")
            self.__basis_sets[groups[0]] = BasisSet(name=groups[4], index=index,
                                                    cartesian=not LocalBasisData.spherical)
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

        return SM(name="Orbital Basis",
                  startReStr=r"\s*\|\s*basis set information\s*\|",
                  subMatchers=[
                      SM(r"\s*\+----*\+", name="<format>", coverageIgnore=True),
                      gauss_type_spherical,
                      SM(r"\s*type\s+atoms\s+prim\s+cont\s+basis", name="Orbital Basis"),
                      SM(r"\s*----*", name="<format>", coverageIgnore=True),
                      basis,
                      SM(r"\s*----*", name="<format>", coverageIgnore=True),
                      SM(r"\s*total:\s*([0-9]+)\s+([0-9]+)\s+([0-9]+)",
                         name="Orbital Basis"),
                      SM(r"\s*----*", name="<format>", coverageIgnore=True),
                      SM(r"\s*total number of primitive shells\s*:\s*([0-9]+)",
                         name="Orbital Basis"),
                      SM(r"\s*total number of contracted shells\s*:\s*([0-9]+)",
                         name="Orbital Basis"),
                      SM(r"\s*total number of cartesian basis functions\s*:\s*([0-9]+)",
                         name="Orbital Basis"),
                      SM(r"\s*total number of SCF-basis functions\s*:\s*([0-9]+)",
                         name="Orbital Basis"),
                  ]
                  )
