"""Tests the basic functionalities for the processor"""

import numpy as np
from qutip import basis, tensor
from qutip_qip.circuit import QubitCircuit

from chalmers_qubit.sarimner.model import SarimnerModel
from chalmers_qubit.sarimner.processor import SarimnerProcessor
from chalmers_qubit.sarimner.compiler import SarimnerCompiler


def test_default_processor_initialization():
    """Tests if a processor is correctly initialized with default values"""

    # Define Model
    model = SarimnerModel(
        2,
        wq=[
            1.0,
            2.0,
        ],
        wr=[2.0, 3.0],
        alpha=[100, 200],
        t1=[1, 2],
        t2=[20, 30],
        zz_crosstalk_static=None,
    )

    compiler = SarimnerCompiler(model, g=2.0)

    # Initialize and check processor attributes
    processor = SarimnerProcessor(model=model, compiler=compiler)

    assert processor.num_qubits == 2
    assert processor.dims == 2 * [3]  # Default
    assert processor.model.__class__.__name__ == "SarimnerModel"
    assert processor.compiler.__class__.__name__ == "SarimnerCompiler"
    assert processor.compiler.g == 2
