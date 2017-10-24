import logging
import numpy as np
from nomadcore.simple_parser import SimpleMatcher as SM

logger = logging.getLogger("nomad.turbomoleParser")


class OrbitalParser(object):

    def __init__(self, context, key="orbitals"):
        context[key] = self
        self.__context = context
        self.__backend = None
        self.__num_orbitals = 0

    def set_backend(self, backend):
        self.__backend = backend

    def build_state_matcher(self):

        def sum_orbitals(backend, groups):
            self.__num_orbitals += int(groups[1])

        ir_rep = SM(r"\s*([a-z]+[0-9]+)\s+([0-9]+)\s+([0-9]+)",
                    # irreducible representation name, num orbitals, num occupied
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
