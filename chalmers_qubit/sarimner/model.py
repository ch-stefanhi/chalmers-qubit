import numpy as np
from typing import Optional
from qutip import destroy, tensor
from qutip_qip.device import Model

__all__ = ["SarimnerModel"]
class SarimnerModel(Model):
    """
    Initializes a new quantum system simulation configuration.

    This method sets up the essential parameters and defaults needed for simulating a quantum system
    with specified qubit characteristics and interactions. It also initializes the internal state
    required for managing the system's dynamics, such as drift and controls, and prepares an empty noise
    model.

    Parameters
    ----------
    qubit_frequencies : list of float
        Frequencies of each qubit in GHz, defining the energy level spacings.
    anharmonicities : list of float
        Anharmonicities of each qubit in GHz, indicating the deviation from harmonic oscillator behavior.
    rotating_frame_frequencies : list of float, optional
        Frequencies defining the rotating frame for each qubit. Defaults to the frequencies of the qubits
        themselves if not provided.
    cpl_matrix : np.ndarray or int, optional
        Coupling matrix between qubits. If an integer is provided, it initializes a matrix filled with this
        integer in the upper triangular part. If not provided, the coupling effect is considered absent.
    dims : list of int, optional
        Dimensions for the state space of each qubit, defaulting to three levels (qutrits) per qubit if not specified.

    Raises
    ------
    ValueError
        If the lengths of `anharmonicities` does not match the number of qubits.
        If `cpl_matrix` is provided but is neither an integer nor a numpy.ndarray.

    Attributes
    ----------
    num_qubits : int
        Number of qubits.
    qubit_frequencies : list of float
        Qubit frequencies stored.
    anharmonicities : list of float
        Stored anharmonicities of each qubit.
    rotating_frame_frequencies : list of float
        Rotating frame frequencies used.
    cpl_matrix : np.ndarray
        Coupling matrix used for the simulation.
    dims : list of int
        Dimensions of each qubit's state space.
    params : dict
        Dictionary holding system parameters for easy access.
    _drift : list
        Internal representation of the system's drift.
    _controls : dict
        Internal setup for system controls.
    _noise : list
        List initialized for adding noise models.
    """
    def __init__(self, 
                 qubit_frequencies: list, 
                 anharmonicities: list, 
                 rotating_frame_frequencies: Optional[list] = None,
                 cpl_matrix: Optional[np.ndarray] = None,
                 dims: Optional[list] = None):

        # number of qubits
        num_qubits = len(qubit_frequencies)

        if len(anharmonicities) != num_qubits:
            raise ValueError("The length of anharmonicities must be the same as num_qubits.")

        if isinstance(cpl_matrix, int):
            # Create an n x n matrix filled with zeros
            matrix = np.zeros((num_qubits, num_qubits))

            # Fill the upper triangular part of the matrix with x
            # NOT SURE IF THIS IS CORRECT
            for i in range(num_qubits):
                for j in range(i+1, num_qubits):
                    matrix[i, j] = cpl_matrix
            cpl_matrix = matrix

        elif isinstance(cpl_matrix, np.ndarray) is False and cpl_matrix is not None: 
            raise ValueError("cpl_matrix should be type int or numpy.ndarray.")

        # Initialize class variables if all checks pass
        self.num_qubits = num_qubits
        self.qubit_frequencies = qubit_frequencies
        self.anharmonicity = anharmonicities
        self.cpl_matrix = cpl_matrix
        self.dims = dims if dims is not None else [3] * num_qubits

        if rotating_frame_frequencies is None:
            self.rotating_frame_frequencies = self.qubit_frequencies
        else: 
            self.rotating_frame_frequencies = rotating_frame_frequencies

        self.params = {
            "wq": self.qubit_frequencies,
            "alpha": self.anharmonicity,
            "wr": self.rotating_frame_frequencies,
            "cpl_matrix": self.cpl_matrix
        }

        # setup drift, controls an noise
        self._drift = self._set_up_drift()
        self._controls = self._set_up_controls()
        self._noise = []

    def _set_up_drift(self):
        drift = []
        for m in range(self.num_qubits):
            destroy_op = destroy(self.dims[m])
            alpha = self.anharmonicity[m] / 2.0
            omega = self.qubit_frequencies[m]
            omega_rot = self.rotating_frame_frequencies[m]
            drift.append(
                ((omega - omega_rot) * destroy_op.dag() * destroy_op
                 + alpha * destroy_op.dag()**2 * destroy_op**2, [m])
            )
        return drift

    def _set_up_controls(self):
        """
        Generate the Hamiltonians and save them in the attribute `controls`.
        """
        num_qubits = self.num_qubits
        dims = self.dims
        controls = {}

        for m in range(num_qubits):
            destroy_op = destroy(dims[m])
            controls["x" + str(m)] = (destroy_op.dag() + destroy_op, [m])
            controls["y" + str(m)] = (1j*(destroy_op.dag() - destroy_op), [m])

        if self.cpl_matrix is not None:
            # Looping through non-zero elements of the coupling matrix
            for (m, n), value in np.ndenumerate(self.cpl_matrix):
                if value != 0:
                    destroy_op1 = destroy(dims[m])
                    destroy_op2 = destroy(dims[n])
                    op1 = tensor(destroy_op1.dag(), destroy_op2)
                    op2 = tensor(destroy_op1, destroy_op2.dag())
                    controls["(xx+yy)" + str(m) + str(n)] = (op1+op2, [m, n])
                    controls["(yx-xy)" + str(m) + str(n)] = (1j*(op1-op2), [m, n])

        return controls

    def get_control_latex(self):
        """
        Get the labels for each Hamiltonian.
        It is used in the method method :meth:`.Processor.plot_pulses`.
        It is a 2-d nested list, in the plot,
        a different color will be used for each sublist.
        """
        num_qubits = self.num_qubits
        labels = [
            {f"x{n}": r"$a_{" + f"{n}" + r"}^\dagger + a_{"
                + f"{n}" + r"}$" for n in range(num_qubits)},
            {f"y{n}": r"$i(a_{" + f"{n}" + r"}^\dagger - a_{"
             + f"{n}" + r"}$)" for n in range(num_qubits)},
        ]
        label_zz = {}

        for m in range(num_qubits - 1):
            for n in range(m + 1, num_qubits):
                label_zz[f"(xx+yy){m}{n}"] = r"$a^\dagger_{"+ f"{m}" + r"}a_{" + f"{n}" + r"} + a^\dagger_{" + f"{n}" + r"}a_{" + f"{m}" + r"}$"
                label_zz[f"(yx-xy){m}{n}"] = r"$i(a^\dagger_{"+ f"{m}" + r"}a_{" + f"{n}" + r"} - a^\dagger_{" + f"{n}" + r"}a_{" + f"{m}" + r"})$"

        labels.append(label_zz)
        return labels
