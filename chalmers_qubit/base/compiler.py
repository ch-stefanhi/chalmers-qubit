from abc import ABC, abstractmethod

from qutip_qip.operations.gateclass import Gate
from qutip_qip.circuit import QubitCircuit
from qutip_qip.pulse import Pulse


class Compiler(ABC):
    """
    Abstract base class for compilers that translate quantum circuits into 
    control pulses for quantum hardware.

    Attributes
    ----------
    model : Model
        The quantum hardware model for which the compiler is designed.
    gate_compilers : dict
        Dictionary mapping gate names to their compilation methods.
    _pulses : dict
        Dictionary storing pulses for each target qubit.
    _tlist : list
        List of time points for the control sequence.

    Methods
    -------
    compile_gate(gate)
        Compiles a single quantum gate into a Pulse.
    compile(circuit, schedule_mode)
        Compiles a quantum circuit into a control pulse sequence.
    """

    def __init__(self, model, **params):
        self.model = model
        self.params = params
        self.gate_compilers = {"RX": self.rx_compiler}
        self._pulses = {}
        self._tlist = []
    
        @property
    def tlist(self):
        """Returns the full tlist for the computation

        Args:
            gate_time (_type_): _description_

        Raises:
            ValueError: _description_

        Returns:
            _type_: _description_
        """
        return self._tlist
    
    @property
    def pulses(self)->dict:
        """Returns the dictionary of Pulses for all qubits

        Returns:
            _type_: _description_
        """
        return self._pulses

    def rx_compiler(self, gate):
        """
        Placeholder method for compiling RX gates.

        Parameters
        ----------
        gate : Gate
            Quantum gate to be compiled.

        Returns
        -------
        Pulse
            Compiled pulse for the RX gate.
        """
        # Placeholder implementation
        pass

    def _add_to_tlist(self, gate_times):
        """
        Add gate times to the global time list.

        Parameters
        ----------
        gate_times : list of float
            Time points to add to the global time list.
        """
        self._tlist.extend(gate_times)

    @abstractmethod
    def compile_gate(self, gate):
        """
        Abstract method to compile a single quantum gate into a Pulse.

        Parameters
        ----------
        gate : Gate
            Quantum gate to be compiled.

        Returns
        -------
        Pulse
            Compiled pulse for the gate.
        """
        pass

    def compile(self, circuit):
        """
        Compile a quantum circuit into a control pulse sequence.

        Parameters
        ----------
        circuit : QubitCircuit or list of Gate
            Quantum circuit to compile.
        schedule_mode : str, optional
            Scheduling mode for the pulses.

        Returns
        -------
        dict
            Dictionary containing compiled time sequence and pulse coefficients.
        """
        gates = circuit.gates if isinstance(circuit, QubitCircuit) else circuit

        for gate in gates:
            if gate.name not in self.gate_compilers:
                raise ValueError(f"Unsupported gate {gate.name}")
            pulse = self.gate_compilers[gate.name](gate)
            self._add_to_tlist(pulse.tlist)
            self._pulses[gate.targets].append(pulse)

        # Placeholder for scheduling and concatenation
        # self._schedule() and concatenate_pulses() need to be implemented
        # compiled_tlist, compiled_coeffs = ...

        scheduled_pulses = self.scheduler(self.pulses)

        # Concatenate pulses
        compiled_tlist, compiled_coeffs = self.concatenate_pulses(scheduled_pulses)

        return compiled_tlist, compiled_coeffs
