from abc import ABC, abstractmethod


class Model(ABC):
    """
    Abstract base class representing a physical model for quantum hardware simulation.

    This class serves as a template for defining the Hamiltonians (both control and drift)
    associated with quantum hardware models. Subclasses should implement the methods
    `setup_controls` and `setup_drift` to define the specific Hamiltonians for their
    quantum hardware.

    Parameters
    ----------
    num_qubits : int
        The number of qubits in the quantum hardware model.
    dims : list, optional
        The dimensions of each qubit system. For instance, [2, 2, 2] for a three-qubit system.
        Default is [3, 3, 3, ..., 3] for `num_qubits` qubits.
    **params : dict
        Additional parameters specific to the quantum hardware model, such as frequencies,
        anharmonicities, coupling strengths, etc.

    Attributes
    ----------
    num_qubits : int
        The number of qubits.
    dims : list
        The dimensions of each component in the quantum system.
    controls : dict of {str: (Qobj, int)}
        A dictionary to store control Hamiltonian where the key defines a
        string representation of the Hamiltonian and the values are
        tuples of (Qobj, list) representing the operator and the qubit
        e.g., {"sx": (sigmax(), 0)}.
    drift : dict of {str: (Qobj, int)}
        A dictionary to store drift Hamiltonian terms where the key defines a
        string representation of the Hamiltonian and the values are
        tuples of (Qobj, list), e.g., {"a": (destroy(self.dims[0]), 0)}.
    params : dict
        Hardware-specific parameters.

    Examples
    --------
    A subclass can implement the Model class like this:

    class CustomQubitModel(Model):
        def __init__(self, num_qubits, **params):
            super().__init__(num_qubits, **params)
            # Additional initialization if needed

        def setup_controls(self):
            # Implementation of control Hamiltonians
            pass

        def setup_drift(self):
            # Implementation of drift Hamiltonians
            pass
    """

    def __init__(self, num_qubits, dims=None, **params):
        if dims is None:
            dims = [3] * num_qubits
        self.num_qubits = num_qubits
        self.dims = dims
        self.params = params
        self._controls = {}
        self._drift = {}

    @abstractmethod
    def setup_controls(self):
        """
        Abstract method to create and add control Hamiltonian terms to the model.

        Subclasses should override this method to define the control Hamiltonians
        specific to their quantum hardware model.

        Example:
        --------
        def setup_controls(self):
            for qubit in range(self.num_qubits):
                # Define control Hamiltonian for each qubit
                pass
        """
        # self._controls = ...
        pass

    @abstractmethod
    def setup_drift(self):
        """
        Abstract method to create and add drift Hamiltonian terms to the model.

        Subclasses should override this method to define the drift Hamiltonians
        specific to their quantum hardware model.

        Example:
        --------
        def setup_drift(self):
            for qubit in range(self.num_qubits):
                # Define drift Hamiltonian for each qubit
                pass
        """
        # self._drift = ...
        pass

    @property
    def controls(self):
        """
        Property to get the control Hamiltonian terms.

        Returns
        -------
        dict
            A dictionary containing the control Hamiltonian terms.
        """
        return self._controls

    @property
    def drift(self):
        """
        Property to get the drift Hamiltonian terms.

        Returns
        -------
        dict
            A dictionary containing the drift Hamiltonian terms.
        """
        return self._drift
