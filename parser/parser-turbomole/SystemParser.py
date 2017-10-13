import numpy as np
from nomadcore.simple_parser import SimpleMatcher as SM


class Atom(object):

    def __init__(self, x, y, z, elem, charge, shells, isotope, pseudo):
        self.x = float(x)
        self.y = float(y)
        self.z = float(z)
        self.elem = elem
        self.charge = int(float(charge))
        self.shells = int(shells) if shells else -1
        self.isotope = int(isotope)
        self.is_pseudo = True if pseudo == "1" else False


class SystemParser(object):

    def __init__(self, context, key="geo"):
        context[key] = self
        self.__context = context
        self.__backend = None
        self.__index_qm_geo = -1
        self.__atoms = list()

    def finalize_sections(self):
        if self.__index_qm_geo >= -1:
            self.__backend.closeSection("section_system", self.__index_qm_geo)
            self.__index_qm_geo = -2

    def build_qm_geometry_matcher(self):

        def open_section(backend, groups):
            self.__backend = backend
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
            charges = np.ndarray(shape=(len(self.__atoms),), dtype=float)
            for i, atom in enumerate(self.__atoms):
                pos[i, 0:3] = (atom.x, atom.y, atom.z)
                labels.append(atom.elem.capitalize())
                charges[i] = atom.charge
            backend.addArrayValues("atom_positions", pos, unit="angstrom")
            backend.addArrayValues("atom_labels", np.asarray(labels, dtype=str))
            backend.addArrayValues("atom_atom_number", charges)

        # x, y, z, element, (shells), charge, (pseudo), isotope
        atom_re = r"\s*([-+]?[0-9]+\.?[0-9]*)" \
                  r"(\s+[-+]?[0-9]+\.?[0-9]*)" \
                  r"(\s+[-+]?[0-9]+\.?[0-9]*)" \
                  r"(\s+[a-zA-Z]+)" \
                  r"(\s+[0-9]+)?" \
                  r"(\s+[-+]?[0-9]+\.?[0-9]*)" \
                  r"(\s+[-+]?[0-9]+)?" \
                  r"(\s+[-+]?[0-9]+)"
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
