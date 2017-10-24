import logging
import numpy as np
from nomadcore.simple_parser import SimpleMatcher as SM

logger = logging.getLogger("nomad.turbomoleParser")


class MethodParser(object):

    def __init__(self, context, key="method"):
        context[key] = self
        self.__context = context
        self.__backend = None
        self.spin_channels = 1

    def set_backend(self, backend):
        self.__backend = backend

    def build_uhf_matcher(self):

        def set_spin_polarized(backend, groups):
            self.spin_channels = 2

        return SM(r"\s*UHF mode switched on !",
                  name="UHF switch",
                  startReAction=set_spin_polarized
                  )
