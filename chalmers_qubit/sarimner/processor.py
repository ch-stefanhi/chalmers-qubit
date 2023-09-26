import numpy as np

from qutip_qip.device.processor import Processor
from qutip_qip.transpiler import to_chain_structure

from chalmers_qubit.sarimner.model import SarimnerModel
from chalmers_qubit.sarimner.compiler import SarimnerCompiler


class SarimnerProcessor(Processor):
    """Sarimner processor

    Args:
        Processor (Processor): A Processor subclass that inherits its methods
    """
    def __init__(self,
                 num_qubits: int,
                 dims: int = None,
                 wq: list = None,
                 alpha: list = None,
                 wr: list = None,
                 t1: list = None,
                 t2: list = None,
                 zz_crosstalk_static: dict = None,
                 gate_params: dict = None,
                 pulse_mode: str = "continuous"):
        """_summary_

        Args:
            num_qubits (int): _description_
            dims (int, optional): _description_. Defaults to None.
            wq (list, optional): _description_. Defaults to None.
            alpha (list, optional): _description_. Defaults to None.
            wr (list, optional): _description_. Defaults to None.
            t1 (list, optional): _description_. Defaults to None.
            t2 (list, optional): _description_. Defaults to None.
            zz_crosstalk_static (dict, optional): _description_. Defaults to None.
            gate_params (dict, optional): _description_. Defaults to None.
            pulse_mode (str, optional): _description_. Defaults to "continous".
        """
        # Construct a model for the processor. The model is expected to have
        # certain methods defined. We choose some default parameters if not
        # specified while initializing the Processor
        if dims is None:
            dims = [3] * num_qubits

        if wq is None:
            wq = [2 * np.pi * 5.0] + [2 * np.pi * 5.4] * (num_qubits - 1)

        # Choose rotating frame frequency as the qubit freq
        wr = wq

        # Anharmonicity in (GHz)
        if alpha is None:
            alpha = [-2 * np.pi * 0.3] * num_qubits
        
        model = SarimnerModel(
            num_qubits=num_qubits,
            dims=dims,
            wq=wq,
            alpha=alpha,
            wr=wr,
            t1=t1,
            t2=t2,
            zz_crosstalk_static=zz_crosstalk_static)
        super(SarimnerProcessor, self).__init__(model=model)
        self.pulse_mode = pulse_mode
