"""Tests the basic functionalities for the processor"""

import numpy as np
from qutip import basis, tensor
from qutip_qip.circuit import QubitCircuit
from chalmers_qubit.sarimner.processor import SarimnerProcessor
from chalmers_qubit.ratatosk.processor import ChalmersQubits


def test_basic_circuit():
    """Tests if a simple circuit runs with the processor"""
    # Define a circuit
    qc = QubitCircuit(2)
    qc.add_gate("X", targets=[0])

    # Run gate-level simulation
    init_state = tensor(basis(3, 0), basis(3, 0))

    # Run pulse-level simulation
    processor = SarimnerProcessor(num_qubits=2)
    processor.load_circuit(qc)
    tlist = np.linspace(0, 20, 300)
    result = processor.run_state(init_state, tlist=tlist)
    assert result.states[-1].dims == [[3, 3], [1, 1]]
