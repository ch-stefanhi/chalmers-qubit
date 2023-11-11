from abc import ABC
import numpy as np

from qutip_qip import Instruction


class Compiler(ABC):
    """
    The abstract base class that defines the basic methods and their call
    signature for any compiler that derives from this class. The compiler
    defines the rules to translate a circuit into pulses.

    .. note::

        This is an abstract base class to outline the general API to define
        new compilers but cannot be used by itself. It provides a series of
        methods and attributes expected while developing a new compiler.

    """

    def __init__(
        self,
        model,
    ):
        self.model = model
        self.global_phase = 0.0

        # Compilers are functions mapped by the name of the Gate to the
        # instructions on how to translate them to an Instruction
        self.gate_compilers = {"GLOBALPHASE": self.global_phase}

    def global_phase(self):
        """Compiles global phase operations"""
        pass
