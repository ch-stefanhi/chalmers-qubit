from abc import ABC
import numpy as np


class Compiler(ABC):
    """
    The abstract base class that defines the basic methods and their call
    signature for any compiler that derives from this class. The compiler
    defines the rules to translate a circuit into pulses with qutip.qobjevo

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
