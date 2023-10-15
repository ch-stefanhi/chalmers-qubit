import numpy as np
from qutip import destroy, tensor, basis, qeye
from copy import deepcopy

from qutip_qip.device import Model
from qutip_qip.noise import RelaxationNoise


class SarimnerModel(Model):
    """Model definitions for the physical system defining the processor.

    The model requires the number of qubits and stores the drift
    and control Hamiltonian terms used in the processor. It also stores
    various parameters that the processor will use for the simulations, e.g.,
    t1, t2 values etc.
    """

    def __init__(
        self,
        num_qubits: int,
        wq: list,
        alpha: list,
        wr: list,
        t1: list,
        t2: list,
        zz_crosstalk_static: dict = None,
        dims: int = None,
    ):
        if dims is None:
            dims = num_qubits * [3]

        self.num_qubits = num_qubits
        self.wq = wq
        self.alpha = alpha
        self.wr = wr
        self.t1 = t1
        self.t2 = t2
        self.zz_crosstalk_static = zz_crosstalk_static
        self.dims = dims
        self._drift = self.get_drift_hamiltonian_terms()
        self._controls = self.get_control_hamiltonian_terms()
        self.params = {"t1": t1, "t2": t2}

    def get_drift_hamiltonian_terms(self):
        drift = []
        for m in range(self.num_qubits):
            destroy_op = destroy(self.dims[m])
            alpha = self.alpha[m] / 2.0
            omega = self.wq[m]
            omega_rot = self.wr[m]
            drift.append(
                (
                    (omega - omega_rot) * destroy_op.dag() * destroy_op
                    + alpha * destroy_op.dag() ** 2 * destroy_op**2,
                    [m],
                )
            )
        return drift

    def get_control_hamiltonian_terms(self):
        """
        Generate the Hamiltonians and save them in the attribute `ctrls`.
        """
        num_qubits = self.num_qubits
        dims = self.dims
        controls = {}

        for m in range(num_qubits):
            destroy_op = destroy(dims[m])
            controls["sx" + str(m)] = (destroy_op.dag() + destroy_op, [m])
            controls["sy" + str(m)] = (1j * (destroy_op.dag() - destroy_op), [m])

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
