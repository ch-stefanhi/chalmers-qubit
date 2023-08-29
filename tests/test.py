import unittest
import numpy
from qutip import Qobj, basis, tensor, fidelity
from qutip_qip.circuit import QubitCircuit
from chalmers_qubit.processor import ChalmersProcessor

import numpy as np
from qutip import Qobj


def project_on_qubit(rho: Qobj) -> Qobj:
    """
    Project a given quantum object (density matrix or state vector) onto the computational basis of qubits.

    Parameters
    ----------
    rho : Qobj
        The input quantum object, either a density matrix ('oper' type) or a state vector.

    Returns
    -------
    qubit_state : Qobj
        The projected quantum object onto the computational basis (0 and 1 states) of qubits.
    """
    num_qubits = len(rho.dims[0])
    # Generate all states for N levels in base 3
    base_N = [np.base_repr(i, base=3) for i in range(3**num_qubits)]
    # Generate all computational basis states in base 2
    base_2 = [np.base_repr(i, base=2) for i in range(2**num_qubits)]
    # Find indices of computational basis states in base N list
    l = [base_N.index(i) for i in base_2]
    if rho.type == 'oper':
        # If rho is a density matrix, project it onto the computational basis
        qubit_state = Qobj([[rho.data[j, i] for i in l] for j in l], dims=[
                           [2] * num_qubits, [2] * num_qubits])
    else:
        # If rho is a state vector, project it onto the computational basis
        qubit_state = Qobj([rho.full()[j] for j in l], dims=[
                           [2] * num_qubits, [1] * num_qubits])
    return qubit_state


class TestSingleQubitGates(unittest.TestCase):
    def setUp(self) -> None:
        self.num_qubits = 1
        self.myprocessor = ChalmersProcessor(self.num_qubits)
        return super().setUp()

    def tearDown(self) -> None:
        return super().tearDown()

    def test_single_qubit_gate(self):
        initial_state = basis(3, 0)
        # Create a circuit
        qc = QubitCircuit(self.num_qubits)
        qc.add_gate("RX", targets=0, arg_value=numpy.pi/2)
        qc.add_gate("RZ", targets=0, arg_value=5*numpy.pi/2)
        qc.add_gate("RX", targets=0, arg_value=numpy.pi/2)
        qc.add_gate("RZ", targets=0, arg_value=numpy.pi/2)
        qc.add_gate("RX", targets=0, arg_value=-numpy.pi/2)
        # Simulate
        result = self.myprocessor.run_state(initial_state, qc=qc)
        final_state = result.states[-1]
        # Project final state to the qubit subspace
        qubit_state = Qobj(final_state.data[:2])
        qubit_initial_state = Qobj(initial_state.full()[:2])
        # Get the ideal target state
        U = qc.compute_unitary()
        target_state = U * qubit_initial_state
        # Compute the fidelity
        f = fidelity(target_state, qubit_state)
        self.assertAlmostEqual(1, f, places=3)


class TestTwoQubitGates(unittest.TestCase):
    def setUp(self) -> None:
        self.num_qubits = 2
        self.myprocessor = ChalmersProcessor(self.num_qubits)
        return super().setUp()

    def tearDown(self) -> None:
        return super().tearDown()

    def test_single_qubit_gate(self):
        initial_state = tensor(basis(3, 1), basis(3, 1))
        # Create a circuit
        qc = QubitCircuit(self.num_qubits)
        qc.add_gate("CZ", controls=0, targets=1)
        # Simulate
        result = self.myprocessor.run_state(initial_state, qc=qc)
        final_state = result.states[-1]
        # Project final state to the qubit subspace
        qubit_state = project_on_qubit(final_state)
        # Compute the fidelity
        phase = float(numpy.angle(qubit_state[-1]))
        self.assertAlmostEqual(phase, numpy.pi, places=1)


if __name__ == '__main__':
    unittest.main()
