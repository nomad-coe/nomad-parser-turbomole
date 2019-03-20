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


class _Atom(object):

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


class _AtomKind(object):

    def __init__(self, index, elem, charge, isotope, pseudo):
        self.index = int(index)
        self.elem = elem.capitalize()
        self.charge = float(charge)
        self.isotope = int(isotope)
        if pseudo:
            self.is_pseudo = True if pseudo == "1" else False
        else:
            self.is_pseudo = None  # effective core potential not specified


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
        self.__index_basis_set = -1
        self.__atoms = list()
        self.__atom_kinds = dict()
        self.__basis_sets = dict()
        self.__basis_one_to_one = True
        self.__auxbasis_sets = dict()
        self.__auxbasis_one_to_one = True

    def purge_data(self):
        self.__index_basis_set = -1
        self.__atoms = list()
        self.__atom_kinds = dict()
        self.__basis_sets = dict()
        self.__basis_one_to_one = True
        self.__auxbasis_sets = dict()
        self.__auxbasis_one_to_one = True

    def set_backend(self, backend):
        self.__backend = backend

    # match builders

    def num_atoms(self):
        return len(self.__atoms)

    def write_atomic_data(self, system_index):
        pos = np.ndarray(shape=(len(self.__atoms), 3), dtype=float)
        labels = list()
        atom_numbers = np.ndarray(shape=(len(self.__atoms),), dtype=float)
        self.__backend.addValue("number_of_atoms", len(self.__atoms), system_index)
        for i, atom in enumerate(self.__atoms):
            pos[i, 0:3] = (atom.x, atom.y, atom.z)
            labels.append(atom.elem)
            atom_numbers[i] = elements.get_atom_number(atom.elem)
        self.__backend.addArrayValues("atom_positions", pos, gIndex=system_index, unit="bohr")
        self.__backend.addArrayValues("atom_labels", np.asarray(labels, dtype=str),
                                      gIndex=system_index)
        self.__backend.addArrayValues("atom_atom_number", atom_numbers, gIndex=system_index)

    def write_basis_set_mapping(self, config_index, method_index):
        # TODO: if these conditions are not true, we need to read the files basis and coords instead
        if self.__basis_one_to_one and len(self.__basis_sets) > 0:
            self.__write_simple_method_basis_set_mapping(method_index, False, self.__basis_sets)
        if self.__auxbasis_one_to_one and len(self.__auxbasis_sets) > 0:
            self.__write_simple_method_basis_set_mapping(method_index, True, self.__auxbasis_sets)
        if sum(x.num_atoms for x in self.__basis_sets.values()) == len(self.__atoms):
            self.__write_simple_atomic_basis_set(config_index, False, self.__basis_sets)
        if sum(x.num_atoms for x in self.__auxbasis_sets.values()) == len(self.__atoms):
            self.__write_simple_atomic_basis_set(config_index, True, self.__auxbasis_sets)

    def __write_simple_atomic_basis_set(self, config_index, is_auxiliary, basis_sets):
        references = {"section_single_configuration_calculation": config_index}
        index = self.__backend.openSection("section_basis_set")
        self.__backend.setSectionInfo("section_basis_set", index, references)
        mapping = np.ndarray(shape=(len(self.__atoms),), dtype=int)
        for i, atom in enumerate(self.__atoms):
            mapping[i] = basis_sets[atom.elem].index
        self.__backend.addArrayValues("mapping_section_basis_set_atom_centered", mapping)
        self.__backend.addValue("basis_set_kind", "density" if is_auxiliary else "wavefunction",
                                index)
        self.__backend.closeSection("section_basis_set", index)

    def __write_simple_method_basis_set_mapping(self, method_index, is_auxiliary, basis_sets):
        references = {"section_method": method_index}
        index = self.__backend.openSection("section_method_basis_set")
        self.__backend.setSectionInfo("section_method_basis_set", index, references)
        self.__backend.addValue("number_of_basis_sets_atom_centered", len(self.__atom_kinds), index)
        refs = np.ndarray(shape=(len(self.__atom_kinds),2), dtype=int)
        refs[:, :] = -1
        for i, elem in enumerate(sorted(self.__atom_kinds)):
            kind = self.__atom_kinds[elem]
            refs[i, 0] = kind.index
            refs[i, 1] = basis_sets[kind.elem].index
        self.__backend.addValue("method_basis_set_kind",
                                "density" if is_auxiliary else "wavefunction", index)
        self.__backend.addArrayValues("mapping_section_method_basis_set_atom_centered", refs, index)
        self.__backend.closeSection("section_method_basis_set", index)

    def build_qm_geometry_matcher(self):

        def add_atom(backend, groups):
            new_atom = _Atom(x=groups[0], y=groups[1], z=groups[2], elem=groups[3],charge=groups[5],
                             isotope=groups[7], shells=groups[4], pseudo=groups[6])
            self.__atoms.append(new_atom)
            if groups[3].capitalize() not in self.__atom_kinds:
                index = backend.openSection("section_method_atom_kind")
                kind = _AtomKind(index=index, elem=groups[3],charge=groups[5],
                                 isotope=groups[7], pseudo=groups[6])
                atom_number = elements.get_atom_number(kind.elem)
                backend.addValue("method_atom_kind_atom_number", atom_number, index)
                # backend.addValue("method_atom_kind_explicit_electrons", ..., index)
                backend.addValue("method_atom_kind_label", kind.elem, index)
                # backend.addValue("method_atom_kind_mass", ..., index)
                # backend.addValue("method_atom_kind_pseudopotential_name", ..., index)
                backend.closeSection("section_method_atom_kind", index)
                self.__atom_kinds[kind.elem] = kind

        # x, y, z, element, (shells), charge, (pseudo), isotope
        atom_re = r"\s*" + 3 * ("("+RE_FLOAT+")\s+") + r"([a-zA-Z]+)(\s+[0-9]+)?" \
                  r"\s+("+RE_FLOAT+")(\s+[-+]?[0-9]+)?\s+([-+]?[0-9]+)"
        atom = SM(atom_re, repeats=True, name="single atom", startReAction=add_atom)
        header_re = r"\s*atomic\s+coordinates\s+atom(?:\s+shells)?\s+charge(?:\s+pseudo)?\s+isotop"

        return SM(name="geometry",
                  startReStr=r"\s*\|\s+Atomic coordinate, charge and isotope? information\s+\|",
                  subMatchers=[
                      SM(r"\s*-{20}-*", name="<format>", coverageIgnore=True),
                      SM(header_re, name="atom list", subMatchers=[atom])
                  ]
                  )

    def build_qm_geometry_matcher_statpt(self):

        def add_atom(backend, groups):
            new_atom = _Atom(x=groups[1], y=groups[2], z=groups[3], elem=groups[0], charge=0,
                             isotope=0, shells=None, pseudo=None)
            self.__atoms.append(new_atom)

        atom = SM(r"\s*[0-9]+\s+([A-z]+)"+ 3 * (r"\s+("+RE_FLOAT+")") + "\s*$",
                  name="single atom",
                  repeats=True,
                  startReAction=add_atom)

        return SM(name="geometry",
                  startReStr=r"\s*ATOM\s+CARTESIAN\s+COORDINATES\s*$",
                  subMatchers=[
                      atom
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
                if key in self.__auxbasis_sets:
                    self.__auxbasis_one_to_one = False
                else:
                    self.__auxbasis_sets[key] = basis_set
            else:
                if key in self.__basis_sets:
                    self.__basis_one_to_one = False
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
                   r"([A-z0-9-\(\)]+)\s+" \
                   r"\[((?:[0-9]+[spdfghij])+)\|((?:[0-9]+[spdfghij])+)\]"
        basis = SM(basis_re, repeats=True, name="basis assignment", startReAction=add_basis_set)
        basis_regex = r"1s\s+" + \
            "".join((r"(?:"+str(2*(i+1)+1)+x+r"\s+)?" for i,x in enumerate("pdfghijkl")))
        basis_type = "(\s+auxiliary)?\s+basis\s+set"
        header_re = r"\s*we\s+will\s+work\s+with\s+the\s+"+basis_regex+"\s*\.\.\." + basis_type
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
