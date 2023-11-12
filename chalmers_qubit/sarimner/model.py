import numpy as np
from qutip import destroy, tensor, basis, qeye, Qobj
from copy import deepcopy

from ..base.model import Model


class SarimnerModel(Model):
    """
    A quantum hardware model representing a processor based on specific qubit characteristics.
    
    This model is designed to simulate the physical system of a quantum processor, incorporating
    both drift and control Hamiltonian terms based on qubit frequencies, anharmonicities,
    and resonance frequencies.

    Parameters
    ----------
    num_qubits : int
        The number of qubits in the model.
    dims : list, optional
        The dimensions of each qubit system. Default is [3, 3, 3, ..., 3] for `num_qubits` qubits.
    qubit_frequencies : list
        List of qubit frequencies.
    anharmonicities : list
        List of anharmonicities for each qubit.
    resonance_frequencies : list
        List of resonance frequencies for the qubits.
    zz_crosstalk_static : dict of {(m, n): (Qobj, int)}, optional
        Dictionary containing static ZZ crosstalk terms between qubit m and n
        (default is None).

    Attributes
    ----------
    wq : list
        Qubit frequencies.
    alpha : list
        Anharmonicities for each qubit.
    wr : list
        Resonance frequencies for the qubits.
    zz_crosstalk_static : dict of {(m, n): (Qobj, int)}
        Static ZZ crosstalk terms.
    controls: dict[Qobj]
        Control Hamiltoninan terms.
    drift: dict[Qobj]
        Drift Hamiltoninan terms.
    """

    def __init__(
        self,
        num_qubits: int,
        qubit_frequencies: list,
        anharmonicities: list,
        resonance_frequencies: list = None,
        zz_crosstalk_static: dict = None,
        dims: list = None,
    ):
        super().__init__(num_qubits, dims=dims)
        self.qubit_frequencies = qubit_frequencies
        self.anharmonicities = anharmonicities
        if resonance_frequencies is None:
            resonance_frequencies = qubit_frequencies
        self.resonance_frequencies = resonance_frequencies
        self.zz_crosstalk = (
            zz_crosstalk_static if zz_crosstalk_static is not None else {}
        )

        self._drift = self.setup_drift()
        self._controls = self.setup_controls()

    def setup_drift(self):
        """
        Sets up the drift Hamiltonian terms based on qubit frequencies, anharmonicities, and
        resonance frequencies. Each term is stored in the drift dictionary with keys indicating
        the nature of the term and associated qubit number.

        Returns
        -------
        drift: dict of {str: (Qobj, [int])}
            The drift Hamiltonian terms.
        """
        drift = {}
        for m in range(self.num_qubits):
            destroy_op = destroy(self.dims[m])
            alpha = self.anharmonicities[m] / 2.0
            omega = self.qubit_frequencies[m]
            omega_rot = self.resonance_frequencies[m]

            # Frequency difference term
            drift[(r"\Delta \omega a.dag a", str(m))] = ((omega - omega_rot)*destroy_op.dag() * destroy_op, [m])

            # Anharmonicity term
            drift[(r"$\alpha a.dag^2 a^2$", str(m))] = (
                alpha * destroy_op.dag() ** 2 * destroy_op ** 2,
                [m],
            )

        # Adding static ZZ coupling terms
        for (m, n), coupling_strength in self.zz_crosstalk.items():
            sigma_z1 = tensor(
                [
                    qeye(self.dims[m])
                    if i != m
                    else destroy(self.dims[m]).dag() * destroy(self.dims[m])
                    for i in range(self.num_qubits)
                ]
            )
            sigma_z2 = tensor(
                [
                    qeye(self.dims[n])
                    if i != n
                    else destroy(self.dims[n]).dag() * destroy(self.dims[n])
                    for i in range(self.num_qubits)
                ]
            )

            zz_term = coupling_strength * sigma_z1 * sigma_z2
            drift[("ZZ", "".join([m, n]))] = (zz_term, [m, n])

        return drift

    def setup_controls(self):
        """
        Sets up the control Hamiltonian terms based on the system's qubit dimensions. These terms
        are used to manipulate the qubit states. The control terms are stored in a dictionary with
        keys given by a tuple of () representing the control action and
        the qubit number.

        Returns
        -------
        controls: dict of {str: (Qobj, int)}
            The drift Hamiltonian terms.
        """
        controls = {}
        for m in range(self.num_qubits):
            destroy_op = destroy(self.dims[m])
            controls[("X", str(m))] = (destroy_op.dag() + destroy_op, [m])
            controls["Y", str(m)] = (1j * (destroy_op.dag() - destroy_op), [m])

        # Inter-qubit coupling terms
        for m in range(self.num_qubits - 1):
            for n in range(m + 1, self.num_qubits):
                d1 = self.dims[m]
                d2 = self.dims[n]
                destroy_op1 = destroy(d1)
                destroy_op2 = destroy(d2)
                op1 = tensor(destroy_op1.dag(), destroy_op2)
                op2 = tensor(destroy_op1, destroy_op2.dag())
                controls[("ab", "".join([m, n]))] = (op1, [m, n])
                controls[("ba", "".join([m, n]))] = (op2, [m, n])

        return controls
