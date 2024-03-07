import numpy as np
from qutip import destroy, tensor
from chalmers_qubit.base.model import Model
from chalmers_qubit.sarimner.noise import DecoherenceNoise

from typing import List, Dict, Any, Optional, Tuple, Hashable

__all__ = ["SarimnerModel"]
class SarimnerModel(Model):
    """
    The processor based on the physical implementation of
    a Transmon qubit.

    Parameters
    ----------
    N: int
        The number of qubits
    t1: float-array
        T1 times
    t2: float-array
        T2 time

    Attributes
    ----------
    params: dict
        A Python dictionary contains the name and the value of the parameters
        in the physical realization, such as laser frequency, detuning etc.
    """

    def __init__(self, 
                 num_qubits: int,
                 qubit_frequencies: list, 
                 anharmonicities: list, 
                 rotating_frame_frequencies: list = None,
                 cpl_matrix: np.ndarray = None,
                 dims: list = None):

        # Check if num_qubits is a positive integer greater than 0
        if not isinstance(num_qubits, int) or num_qubits <= 0:
            raise ValueError("num_qubits must be a positive integer greater than 0.")

        # Check if the length of qubit_frequencies is the same as num_qubits
        if len(qubit_frequencies) != num_qubits:
            raise ValueError("The length of qubit_frequencies must be the same as num_qubits.")

         # Check if the length of anharmonicities is the same as num_qubits
        if len(anharmonicities) != num_qubits:
            raise ValueError("The length of anharmonicities must be the same as num_qubits.")
        
        if isinstance(cpl_matrix, int):
            # Create an n x n matrix filled with zeros
            matrix = np.zeros((num_qubits, num_qubits))
            
            # Fill the upper triangular part of the matrix with x
            # NOT SURE IF THIS IS CORRECT
            for i in range(num_qubits):
                for j in range(i, num_qubits):
                    matrix[i, j] = cpl_matrix
            cpl_matrix = matrix
            
        elif isinstance(cpl_matrix, np.ndarray) is False and cpl_matrix is not None: 
            raise ValueError("cpl_matrix should be type int or numpy.ndarray.")
        
        # Initialize class variables if all checks pass
        self.num_qubits = num_qubits
        self.qubit_frequencies = qubit_frequencies # Qubit frequency in (GHz)
        self.anharmonicity = anharmonicities # Anharmonicity in (GHz)
        self.cpl_matrix = cpl_matrix # coupling matrix
        self.dims = dims if dims is not None else [3] * num_qubits

        # Choose rotating frame frequency as the qubit freq
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
        # init empty noise list
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
        Generate the Hamiltonians and save them in the attribute `ctrls`.
        """
        num_qubits = self.num_qubits
        dims = self.dims
        controls = {}

        for m in range(num_qubits):
            destroy_op = destroy(dims[m])
            controls["X" + str(m)] = (destroy_op.dag() + destroy_op, [m])
            controls["Y" + str(m)] = (1j*(destroy_op.dag() - destroy_op), [m])

        if self.cpl_matrix is not None:
            # Looping through non-zero elements of the coupling matrix
            for (m, n), value in np.ndenumerate(self.cpl_matrix):
                if value != 0:
                    destroy_op1 = destroy(dims[m])
                    destroy_op2 = destroy(dims[n])
                    op1 = tensor(destroy_op1.dag(), destroy_op2)
                    op2 = tensor(destroy_op1, destroy_op2.dag())
                    controls["ab" + str(m) + str(n)] = (op1+op2, [m, n])
                    controls["ba" + str(m) + str(n)] = (1j*(op1-op2), [m, n])

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
            {f"X{n}": r"$a_{" + f"{n}" + r"}^\dagger + a_{"
                + f"{n}" + r"}$" for n in range(num_qubits)},
            {f"Y{n}": r"$i(a_{" + f"{n}" + r"}^\dagger - a_{"
             + f"{n}" + r"}$)" for n in range(num_qubits)},
        ]
        label_zz = {}

        for m in range(num_qubits - 1):
            for n in range(m + 1, num_qubits):
                label_zz[f"ab{m}{n}"] = r"$a^\dagger_{"+ f"{m}" + r"}a_{" + f"{n}" + r"} + a^\dagger_{" + f"{n}" + r"}a_{" + f"{m}" + r"}$"
                label_zz[f"ba{m}{n}"] = r"$i(a^\dagger_{"+ f"{m}" + r"}a_{" + f"{n}" + r"} - a^\dagger_{" + f"{n}" + r"}a_{" + f"{m}" + r"})$"

        labels.append(label_zz)
        return labels
