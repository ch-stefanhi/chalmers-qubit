import unittest
import numpy
from qutip import Qobj, basis, tensor, fidelity
from qutip_qip.circuit import QubitCircuit
from device import ChalmersQubits
from device.operations import project_on_qubit


class TestSingleQubitGates(unittest.TestCase):
    def setUp(self) -> None:
        self.num_qubits = 1
        self.myprocessor = ChalmersQubits(self.num_qubits)
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
        self.myprocessor = ChalmersQubits(self.num_qubits)
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
