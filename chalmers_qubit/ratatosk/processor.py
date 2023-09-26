import numpy as np
from qutip import basis, destroy, qeye, tensor, propagator
from qutip_qip.device import ModelProcessor, Model
from qutip_qip.transpiler import to_chain_structure
from .compiler import ChalmersCompiler
from .model import ChalmersQubitsModel

__all__ = ["ChalmersQubits"]


class ChalmersQubits(ModelProcessor):
    """
    Description goes here

    Parameters
    ----------
    num_qubits: int
        The number of qubits in the system.
    dims: list, optional
        The dimension of each component system.
        Default value is a qubit system of ``dim=[2,2,2,...,2]``.
    """

    def __init__(self, num_qubits, dims=None, **params):
        if dims is None:
            dims = [3] * num_qubits
        model = ChalmersQubitsModel(
            num_qubits=num_qubits,
            dims=dims,
            **params,
        )
        super(ChalmersQubits, self).__init__(model=model)
        self.native_gates = None
        self._default_compiler = ChalmersCompiler
        self.pulse_mode = "discrete"

    def run_propagator(self, qc=None, noisy=False, **kwargs):
        """

        Parameters
        ----------
        qc: :class:`qutip.qip.QubitCircuit`, optional
            A quantum circuit. If given, it first calls the ``load_circuit``
            and then calculate the evolution.

        states: :class:`qutip.Qobj`, optional
         Old API, same as init_state.

        **kwargs
           Keyword arguments for the qutip solver.

        Returns
        -------
        evo_result: :class:`qutip.Result`
            If ``analytical`` is False,  an instance of the class
            :class:`qutip.Result` will be returned.

            If ``analytical`` is True, a list of matrices representation
            is returned.
        """
        if qc is not None:
            self.load_circuit(qc)
        # construct qobjevo for unitary evolution
        noisy_qobjevo, c_ops = self.get_qobjevo(noisy=noisy)

        # time steps
        tlist = noisy_qobjevo.tlist
        H = noisy_qobjevo.to_list()

        # Compute drift Hamiltonians
        H_drift = 0
        drift = self._get_drift_obj()
        for drift_ham in drift.drift_hamiltonians:
            H_drift += drift_ham.get_qobj(self.dims)
        H[0] = H_drift

        # compute the propagator
        evo_result = propagator(H=H, t=tlist, c_op_list=c_ops, **kwargs)
        return evo_result
