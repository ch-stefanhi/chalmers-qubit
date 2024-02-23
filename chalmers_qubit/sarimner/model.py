import numpy as np
from qutip import destroy, tensor
from chalmers_qubit.base.model import Model
from chalmers_qubit.sarimner.noise import DecoherenceNoise

from typing import List, Dict, Any, Optional, Tuple, Hashable

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
                 dims: list = None,
                 t1: list = None, 
                 t2: list = None):
        
        self.num_qubits = num_qubits
        self.dims = dims if dims is not None else [3] * num_qubits
        # Qubit frequency in (GHz)
        self.qubit_frequencies = qubit_frequencies
        # Choose rotating frame frequency as the qubit freq
        if rotating_frame_frequencies is None:
            self.rotating_frame_frequencies = self.qubit_frequencies
        # Anharmonicity in (GHz)
        self.anharmonicity = anharmonicities
        # Decoherence in (ns)
        self.t1 = t1
        self.t2 = t2

        self.params = {
            "wq": self.qubit_frequencies,
            "alpha": self.anharmonicity,
            "wr": self.rotating_frame_frequencies,
        }

        self._drift = self._set_up_drift()
        self._controls = self._set_up_controls()

        if t1 is not None and t2 is not None:
            self._noise = [DecoherenceNoise(self.num_qubits, self.dims, t1=t1, t2=t2)]

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

        for m in range(self.num_qubits - 1):
            for n in range(m + 1, self.num_qubits):
                d1 = dims[m]
                d2 = dims[n]
                destroy_op1 = destroy(d1)
                destroy_op2 = destroy(d2)
                op1 = tensor(destroy_op1.dag(), destroy_op2)
                op2 = tensor(destroy_op1, destroy_op2.dag())
                controls["ab" + str(m) + str(n)] = (op1, [m, n])
                controls["ba" + str(m) + str(n)] = (op2, [m, n])

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
                label_zz[f"ab{m}{n}"] = r"$a^\dagger_{" + \
                    f"{m}" + r"}a_{" + f"{n}" + r"}$"
                label_zz[f"ba{m}{n}"] = r"$a^\dagger_{" + \
                    f"{n}" + r"}a_{" + f"{m}" + r"}$"

        labels.append(label_zz)
        return labels
