import numpy as np

from qutip_qip.circuit import QubitCircuit


from ..base.processor import Processor


from chalmers_qubit.sarimner.model import SarimnerModel
from chalmers_qubit.sarimner.compiler import SarimnerCompiler
from chalmers_qubit.sarimner.scheduler import SarimnerScheduler


class SarimnerProcessor(Processor):
    """Sarimner processor

    Args:
        Processor (Processor): A Processor subclass that inherits its methods
    """

    def __init__(
        self,
        model: SarimnerModel,
        compiler: SarimnerCompiler,
        scheduler=None,
        noise=None,
        transpiler=None,
    ):
        if scheduler is None:
            scheduler = SarimnerScheduler()

        super().__init__(
            model=model,
            compiler=compiler,
            scheduler=scheduler,
            noise=noise,
            transpiler=transpiler,
        )

    def load_circuit(self, qc: QubitCircuit):
        return None

    def run_state(self, init_state=None, solver="mesolve", **kwargs):
        return None
