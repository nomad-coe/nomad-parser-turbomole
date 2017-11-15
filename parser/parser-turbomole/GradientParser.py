import logging
import math
import re
from functools import total_ordering
import numpy as np
from nomadcore.simple_parser import SimpleMatcher as SM
from TurbomoleCommon import RE_FLOAT

logger = logging.getLogger("nomad.turbomoleParser")


@total_ordering
class _Forces(object):

    def __init__(self, elem, index):
        self.x = math.nan
        self.y = math.nan
        self.z = math.nan
        self.elem = elem.capitalize()
        self.index = int(index)

    def __eq__(self, other):
        if isinstance(other, _Forces):
            return self.elem == other.elem and self.index == other.index
        else:
            return False

    def __lt__(self, other):
        if isinstance(other, _Forces):
            if self.index < other.index:
                return True
            elif self.elem < other.elem:
                return True
        return False


class GradientParser(object):

    __re_atoms = re.compile(r"\s*([0-9]+)\s+([A-z]+)")
    __re_force = re.compile(r"\s*("+RE_FLOAT+")")

    def __init__(self, context, key="gradient"):
        context[key] = self
        self.__context = context
        self.__backend = None
        self.__forces = list()
        self.__forces_staging = list()

    def purge_data(self):
        self.__forces = list()
        self.__forces_staging = list()

    def set_backend(self, backend):
        self.__backend = backend

    # matcher generation methods

    # improve: check if Turbomole allows constraints in force calculations
    def __write_forces(self, backend, groups):
        num_atoms = len(self.__forces)
        forces = np.ndarray(shape=(num_atoms, 3), dtype=float)
        forces[:, :] = 0.0
        for i, atom in enumerate(sorted(self.__forces)):
            forces[i, 0:3] = atom.x, atom.y, atom.z
        self.__backend.addArrayValues("atom_forces_raw", forces, unit="forceAu")

    def build_gradient_matcher(self):

        def extract_atoms(backend, groups):
            self.__forces_staging = list()
            match = self.__re_atoms.match(groups[0], 0)
            while match:
                index = int(match.group(1))
                element = match.group(0)
                to_add = _Forces(element, index)
                if to_add in self.__forces:
                    to_update = self.__forces[self.__forces.index(to_add)]
                    self.__forces_staging.append(to_update)
                else:
                    self.__forces_staging.append(to_add)
                    self.__forces.append(to_add)
                match = self.__re_atoms.match(groups[0], match.endpos)

        def extract_derivative(backend, groups):
            match = self.__re_force.match(groups[1], 0)
            axis = groups[0]
            index = 0
            while match:
                gradient_component = float(match.group(1).replace("D", "E").replace("d", "e"))
                setattr(self.__forces_staging[index], axis, gradient_component)
                match = self.__re_force.match(groups[1], match.endpos)
                index += 1

        def gradients():
            def derivative():
                return SM(r"dE/d([xyz])((?:\s+"+RE_FLOAT+")+)\s*$",
                          name="gradients",
                          startReAction=extract_derivative,
                          required=True
                          )
            return SM(r"\s*ATOM((?:\s+[0-9]+\s+[A-z]+)+)\s*$",
                      name="atoms list",
                      required=True,
                      repeats=True,
                      startReAction=extract_atoms,
                      subMatchers=[
                          derivative(),
                          derivative(),
                          derivative()
                      ]
                      )
        max_component = SM(r"\s*\|maximum\s+component\s+of\s+gradient\|\s+:\s+("
                           +RE_FLOAT+")\s+\(atom\s+([0-9]+)\s+([A-z]+)\s*\)\s*$",
                           name="max gradient component",
                           startReAction=self.__write_forces,
                           required=True)

        return SM(r"\s*cartesian\s+gradient\s+of\s+the\s+energy\s+\(hartree/bohr\)\s*$",
                  name="gradients header",
                  subMatchers=[
                      gradients(),
                      max_component
                  ]
                  )
