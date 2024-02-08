import numpy as np

import numpy as np
from qutip_qip.operations.gateclass import Gate
from qutip_qip.circuit import QubitCircuit
from qutip_qip.pulse import Pulse

from ..base.compiler import Compiler
from ..base.model import Model


class SarimnerCompiler(Compiler):
    """
    Compiler for the Sarimner processor.

    This compiler translates quantum gates into control pulses specific
    to the Sarimner processor's architecture.

    Parameters
    ----------
    model : Model
        The quantum hardware model for the Sarimner processor.
    amp : float
        Amplitude for the control pulses.
    scale : float
        Scaling factor for the control pulses.
    """

    def __init__(
        self,
        model: Model,
        amplitude: float,
        scale: float,
        drive_freq: float,
        detuning: float,
    ):
        super().__init__(model)
        self.amplitude = amplitude
        self.scale = scale
        self.drive_freq = drive_freq
        self.detuning = (detuning,)

        self.phase = [0] * model.num_qubits
        self.gate_compilers = {
            "GLOBALPHASE": self.global_phase_compiler,
            "RX": self.rx_compiler,
        }

    # Implement the gate compiler methods (e.g., rx_compiler, ry_compiler, etc.)
    # These methods will use self.amp, self.scale, and other attributes from self.model
    # Similarly implement other gate compilers

    def compile_gate(self, gate: Gate):
        """
        Compiles a single quantum gate into hardware instructions.
        """
        # Use gate compiler based on gate name
        if gate.name in self.gate_compilers:
            return self.gate_compilers[gate.name](gate)
        else:
            raise NotImplementedError(
                f"Gate {gate.name} is not supported by this compiler."
            )

    def compile_circuit(self, circuit: QubitCircuit):
        """
        Compiles a quantum circuit into a sequence of hardware instructions.
        """
        instructions = []
        for gate in circuit.gates:
            instructions.extend(self.compile_gate(gate))
        return instructions

    def global_phase_compiler(self, gate: Gate):
        """
        Compiles the global phase operation.
        """
        self.global_phase += gate.arg_value

    def compile_rotation(self, gate, gate_duration=1.0, num_points: int = 100):
        """
        Compiles a single rotation Gate object into a Pulse.

        This is a simplified example assuming single-qubit rotation gates
        and Gaussian-shaped control pulses.

        Parameters
        ----------
        gate : Gate
            A qutip Gate object to be compiled.

        gate_duration : float, optional
            Duration of the gate (and thus the pulse).

        Returns
        -------
        Pulse
            A Pulse object representing the control pulse for the gate.
        """
        if gate.name not in ["RX", "RY"]:
            raise NotImplementedError("Gate type not supported for compilation.")

        # Define the basic pulse shape (Gaussian for simplicity)
        tlist = np.linspace(0, gate_duration, num_points)
        coeff = np.exp(-(tlist - gate_duration/2)**2 / (2 * (gate_duration/4)**2))
        
        # Scale the pulse amplitude based on the gate rotation angle
        coeff *= gate.arg_value / np.pi

        # Create the qobj representing the Hamiltonian of the pulse
        # For simplicity, assuming sigmax() for RX and sigmay() for RY
        if gate.name == "RX":
            # instead of sigmax, use the operator defined for the particular
            # dim in the model
            qobj = self.model.controls[("X", gate.targets[0])]

        elif gate.name == "RY":
            qobj = self.model.controls[("Y", gate.targets[0])]

        # Create the Pulse object
        pulse = Pulse(qobj, gate.targets, tlist, coeff)

        return pulse

    def rx_compiler(self, gate: Gate):
        """Pulses for an RX gate

        Args:
            gate (Gate): _description_
        """
        return self.compile_rotation(self, gate)
