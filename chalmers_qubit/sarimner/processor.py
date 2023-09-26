from qutip_qip.device.modelprocessor import ModelProcessor
from qutip_qip.transpiler import to_chain_structure

from chalmers_qubit.sarimner.model import SarimnerModel
from chalmers_qubit.sarimner.compiler import SarimnerCompiler


class SarimnerProcessor(ModelProcessor):
    """
    Sarimner processor

    Parameters
    ----------
    num_qubits: int
        The number of qubits in the system.
    dims: list, optional
        The dimension of each component system.
        Default value is a qutrit system of ``dim=[3, 3, ..., 3]``.
    """

    def __init__(self, num_qubits, dims=None, zz_crosstalk=False, **params):
        if dims is None:
            dims = [3] * num_qubits
        model = SarimnerModel(
            num_qubits=num_qubits,
            dims=dims,
            zz_crosstalk=zz_crosstalk,
            **params,
        )
        super(SarimnerProcessor, self).__init__(model=model)
        self.native_gates = ["RX", "RY", "CNOT"]
        self._default_compiler = SarimnerCompiler
        self.pulse_mode = "continuous"
