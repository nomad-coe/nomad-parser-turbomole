"""This module constructs the parser for the STATPT module from TurboMole"""

import logging
import re
from nomadcore.simple_parser import SimpleMatcher as SM
from TurbomoleCommon import RE_FLOAT
import TurbomoleCommon as Common

logger = logging.getLogger("nomad.turbomoleParser")


class STATPTparser(object):

    def __init__(self, context, key="statpt"):
        context[key] = self
        self.__context = context
        self.__backend = None

    def purge_data(self):
        pass

    def set_backend(self, backend):
        self.__backend = backend

    def build_parser(self):
        references = SM(r"\s{15,}[^+ ]+",
                        name="references",
                        coverageIgnore=True,
                        repeats=True,
                        )
        header = SM(name="Credits",
                    startReStr=r"\s*this is S T A T P T\s*$",
                    coverageIgnore=True,
                    subMatchers=[references],
                    endReStr=r"\s*\*+\s*Stationary\s*point\s*options\s*\*+\s*$"
                    )

        sub_matchers = [
            self.__context.build_start_time_matcher(),
            header,
            self.__context.build_end_time_matcher("statpt")
        ]

        return self.__context.build_module_matcher("statpt", sub_matchers)
