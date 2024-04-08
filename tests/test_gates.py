import unittest
import numpy as np
from qutip import Qobj, basis, tensor, fidelity, average_gate_fidelity
from qutip_qip.circuit import QubitCircuit
from chalmers_qubit.sarimner import SarimnerModel, SarimnerProcessor
from chalmers_qubit.base.operations import project_on_qubit
from chalmers_qubit.base.gates import cczs


class TestSingleQubitGates(unittest.TestCase):
    def setUp(self) -> None:
        self.num_qubits = 1
        self.model = SarimnerModel(
            qubit_frequencies=[2 * np.pi * 5.0],
            anharmonicities=[-2 * np.pi * 0.3],
        )
        self.processor = SarimnerProcessor(model=self.model)
        return super().setUp()

    def tearDown(self) -> None:
        return super().tearDown()

    def test_global_phase_gate(self):
        # Create a circuit
        qc = QubitCircuit(self.num_qubits)
        arg_value = np.pi / 2
        qc.add_gate("RX", targets=0, arg_value=arg_value)
        qc.add_gate("GLOBALPHASE", targets=0, arg_value=arg_value)
        qc.add_gate("RZ", targets=0, arg_value=arg_value)
        # Load circuit onto processor
        self.processor.load_circuit(qc)
        # Extract global phase
        global_phase = self.processor.global_phase
        self.assertEqual(global_phase, arg_value)

    def test_rx_gate(self):
        qc = QubitCircuit(self.num_qubits)
        qc.add_gate("RX", targets=0, arg_value=np.pi / 2)
        # Load circuit onto processor
        coeffs, tlist = self.processor.load_circuit(qc)
        # Compute the propagator
        result = self.processor.run_propagator()
        # Get final propagator and project onto qubit subspace
        prop = project_on_qubit(result[-1])
        # Get the ideal unitary
        U = qc.compute_unitary()
        # Compute the average gate fidelity
        f = average_gate_fidelity(U, prop)
        # We want 99.99% Average gate fidelity
        self.assertAlmostEqual(1, f, places=3, msg="Precision of RX gate failed.")

    def test_ry_gate(self):
        qc = QubitCircuit(self.num_qubits)
        qc.add_gate("RY", targets=0, arg_value=np.pi / 2)
        # Load circuit onto processor
        coeffs, tlist = self.processor.load_circuit(qc)
        # Compute the propagator
        result = self.processor.run_propagator()
        # Get final propagator and project onto qubit subspace
        prop = project_on_qubit(result[-1])
        # Get the ideal unitary
        U = qc.compute_unitary()
        # Compute the average gate fidelity
        f = average_gate_fidelity(U, prop)
        # We want 99.99% Average gate fidelity
        self.assertAlmostEqual(1, f, places=3, msg="Precision of RY gate failed.")

    def test_rz_gate(self):
        qc = QubitCircuit(self.num_qubits)
        qc.add_gate("RZ", targets=0, arg_value=np.pi / 2)
        qc.add_gate("RX", targets=0, arg_value=np.pi / 2)
        # Load circuit onto processor
        coeffs, tlist = self.processor.load_circuit(qc)
        # Compute the propagator
        init_state = (basis(3, 0) + basis(3, 1)) / np.sqrt(2)
        result = self.processor.run_state(init_state)
        final_state = project_on_qubit(result.states[-1])
        ideal_state = basis(2, 0)
        f = fidelity(final_state, ideal_state)
        # We want 99.99% average state fidelity
        self.assertAlmostEqual(1, f, places=3, msg="Precision of RZ gate failed")

    def test_x_gate(self):
        # Create a quantum circuit with an X-gate
        qc = QubitCircuit(self.num_qubits)
        qc.add_gate("X", targets=0)
        # Load circuit onto processor
        coeffs, tlist = self.processor.load_circuit(qc)
        # Compute the propagator
        result = self.processor.run_propagator()
        # Extract the global phase
        global_phase = self.processor.global_phase
        # Project onto qubit subspace and correct the global phase
        final_prop = np.exp(1j*global_phase)*project_on_qubit(result[-1])
        # Compute ideal unitary
        U_ideal = qc.compute_unitary()
        # Calculate the average gate fidelity
        f = average_gate_fidelity(final_prop, U_ideal)
        # We want 99.99% average state fidelity
        self.assertAlmostEqual(1, f, places=3)

    def test_hadamard_gate(self):
        qc = QubitCircuit(self.num_qubits)
        qc.add_gate("H", targets=0)
        # Load circuit onto processor
        coeffs, tlist = self.processor.load_circuit(qc)
        # Prepare "plus" state as initial state
        init_state = (basis(3, 0) + basis(3, 1)) / np.sqrt(2)
        result = self.processor.run_state(init_state)
        final_state = project_on_qubit(result.states[-1])
        # This should be the ideal state
        ideal_state = basis(2, 0)
        f = fidelity(final_state, ideal_state)
        # We want 99.99% average state fidelity
        self.assertAlmostEqual(1, f, places=3, msg="Precision of H-gate failed")
        # Prepare "1" as initial state
        init_state = basis(3, 1) 
        result = self.processor.run_state(init_state)
        final_state = project_on_qubit(result.states[-1])
        # This should be the ideal state up to a rz(pi) rotation
        ideal_state = (basis(2, 0) + basis(2, 1)) / np.sqrt(2)
        f = fidelity(final_state, ideal_state)
        # We want 99.99% average state fidelity
        self.assertAlmostEqual(1, f, places=3, msg="Precision of H-gate failed")

    def test_idling_gate(self):
        t_total = 1e3 # time in (ns)
        qc = QubitCircuit(self.num_qubits)
        qc.add_gate("IDLE", targets=0, arg_value=t_total)
        # Load circuit onto processor
        coeffs, tlist = self.processor.load_circuit(qc)
        # Simulate
        init_state = (basis(3, 0) + basis(3, 1)).unit()
        result = self.processor.run_state(init_state)
        final_state = result.states[-1]
        # Compute the state fidelity
        f = fidelity(final_state, init_state)
        self.assertAlmostEqual(1, f, places=7)

    def test_single_qubit_circuit(self):
        # Create a somewhat random qubit circuit
        qc = QubitCircuit(self.num_qubits)
        qc.add_gate("RX", targets=0, arg_value=np.pi/2)
        qc.add_gate("RZ", targets=0, arg_value=np.pi/2)
        qc.add_gate("RY", targets=0, arg_value=np.pi/2)
        qc.add_gate("RZ", targets=0, arg_value=-np.pi/2)
        qc.add_gate("RX", targets=0, arg_value=np.pi)
        # Load circuit onto processor
        coeffs, tlist = self.processor.load_circuit(qc)
        # Simulate
        initial_state = basis(3, 0)
        result = self.processor.run_state(initial_state)
        final_state = result.states[-1]
        # Project final state to the qubit subspace
        qubit_state = Qobj(final_state.data[:2])
        qubit_initial_state = Qobj(initial_state.full()[:2])
        # Get the ideal target state
        U = qc.compute_unitary()
        target_state = U * qubit_initial_state
        print(np.abs(qubit_state.full()))
        # Compute the fidelity
        print(np.abs(target_state.full()))
        f = fidelity(target_state, qubit_state)
        self.assertAlmostEqual(1, f, places=3)


class TestTwoQubitGates(unittest.TestCase):
    def setUp(self) -> None:
        self.num_qubits = 2
        self.model = SarimnerModel(
            qubit_frequencies=[2 * np.pi * 5.0, 2 * np.pi * 5.4],
            anharmonicities=[-2 * np.pi * 0.3, -2 * np.pi * 0.3],
            cpl_matrix=np.array([[0, 2* np.pi * 1e-3], 
                                 [0, 0]]),
        )
        self.processor = SarimnerProcessor(model=self.model)
        return super().setUp()

    def tearDown(self) -> None:
        return super().tearDown()

    def test_cz_gate(self):
        # Create a circuit
        qc = QubitCircuit(self.num_qubits)
        qc.add_gate("CZ", controls=0, targets=1)
        # Load circuit onto processor
        coeffs, tlist = self.processor.load_circuit(qc)
        # Simulate
        init_state = tensor(basis(3,1), basis(3,1))
        result = self.processor.run_state(init_state)
        # Project final state onto the qsubit subspace
        qubit_state = project_on_qubit(result.states[-1])
        # Compute the phase and make sure that it is -1
        phase = np.angle(qubit_state[-1])[0][0]
        self.assertAlmostEqual(phase, np.pi, places=1)

    def test_iswap_gate(self):
        # Create a circuit
        qc = QubitCircuit(self.num_qubits)
        qc.add_gate("ISWAP", controls=0, targets=1)
        # Load circuit onto processor
        coeffs, tlist = self.processor.load_circuit(qc)
        # Simulate
        init_state = tensor(basis(3, 0), basis(3, 1))
        result = self.processor.run_state(init_state)
        # Project final state onto state 10
        final_state = project_on_qubit(result.states[-1])[2]
        # Compute the phase
        phase_1 = np.angle(final_state)[0][0]

        # Simulate
        init_state = tensor(basis(3, 1), basis(3, 0))
        result = self.processor.run_state(init_state)
        # Project final state onto state 01
        final_state = project_on_qubit(result.states[-1])[1]
        # Compute the phase
        phase_2 = np.angle(final_state)[0][0]
        # Assert that the phase is -i
        self.assertAlmostEqual(phase_1, -np.pi/2, places=1)
        self.assertAlmostEqual(phase_2, -np.pi / 2, places=1)


class TestThreeQubitGates(unittest.TestCase):
    def setUp(self) -> None:
        self.num_qubits = 3
        self.model = SarimnerModel(
            qubit_frequencies=[2 * np.pi * 5.0, 2 * np.pi * 5.4, 2 * np.pi * 5.2],
            anharmonicities=[-2 * np.pi * 0.3, -2 * np.pi * 0.3, -2 * np.pi * 0.3],
            cpl_matrix=np.array([[0, 2 * np.pi * 1e-3, 2 * np.pi * 1e-3], 
                                 [0, 0, 0], 
                                 [0, 0, 0]]),
        )
        self.processor = SarimnerProcessor(model=self.model)
        return super().setUp()

    def tearDown(self) -> None:
        return super().tearDown()

    def test_cczs_gate(self):
        # Create a circuit
        qc = QubitCircuit(self.num_qubits)
        qc.user_gates = {"CCZS": cczs}
        qc.add_gate("CCZS", targets=[0,1,2], arg_value=[np.pi/2,np.pi,0])
        # Load circuit onto processor
        coeffs, tlist = self.processor.load_circuit(qc)
        # Simulate
        initial_state = tensor(basis(3, 1), basis(3, 1), basis(3, 1))
        result = self.processor.run_state(initial_state)
        # Project final state onto the qsubit subspace
        qubit_state = project_on_qubit(result.states[-1]) # type: ignore
        # Compute the phase
        phase = np.angle(qubit_state[-1])[0][0]
        self.assertAlmostEqual(phase, np.pi, places=1)


if __name__ == "__main__":
    unittest.main()
