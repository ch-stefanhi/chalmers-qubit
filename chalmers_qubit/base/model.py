from abc import ABC
from abc import abstractmethod


class Model(ABC):
    """
    Template class for a physical model representing quantum hardware.
    The concrete model class does not have to inherit from this,
    as long as the following methods are defined.

    .. note::

        This is an abstract base class to outline the general API to define
        new compilers but cannot be used by itself. It provides a series of
        methods and attributes expected while developing a new compiler.

    Parameters
    ----------
    num_The number of qubits
        The number of qubits.
    dims : list, optional
        The dimension of each component system.
        Default value is a qubit system of ``dim=[2,2,2,...,2]``.
    **params :
        Hardware parameters for the model.

    Attributes
    ----------
    num_The number of qubits
        The number of qubits.
    dims : list, optional
        The dimension of each component system.
    params : dict
        Hardware parameters for the model.
    """
    def __init__(self, num_qubits, dims=None, **params):
        if dims is None:
            dims = num_qubits * [3]
        self.num_qubits = num_qubits
        self.dims = dims

        self.controls = {}
        self.drift = {}

    @abstractmethod
    def setup_controls(self):
        """Creates the control Hamiltonian terms and adds it to the model
        """
        # num_qubits = self.num_qubits
        # dims = self.dims
        # controls = {}

        # for m in range(num_qubits):
        #     destroy_op = destroy(dims[m])
        #     op = destroy_op + destroy_op.dag()
        #     self.controls["sx" + str(m)] = (2 * np.pi / 2 * op, [m])
        pass

    @abstractmethod
    def setup_drift(self):
        """Creates the drift Hamiltonian terms and adds it to the model
        """
        # for m in range(self.num_qubits):
        #     destroy_op = destroy(self.dims[m])
        #     coeff = 2 * np.pi * self.params["alpha"][m] / 2.0
        #     self.drift["adag_a"+str(m)] = (coeff * destroy_op.dag() ** 2 * destroy_op**2, [m]))
        pass
