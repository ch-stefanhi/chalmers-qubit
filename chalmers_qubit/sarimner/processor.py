from qutip_qip.device import Processor

from chalmers_qubit.sarimner.model import SarimnerModel
from chalmers_qubit.sarimner.compiler import SarimnerCompiler


class SarimnerProcessor(Processor):
    """
    Description goes here

    Parameters
    ----------
    num_qubits: int
        The number of qubits in the system.
    dims: list, optional
        The dimension of each component system.
        Default value is a qubit system of ``dim=[2,2,2,...,2]``.
    """

    def __init__(self, 
                 num_qubits,
                 qubit_frequencies: list, 
                 anharmonicities: list, 
                 rotating_frame_frequencies: list = None, 
                 dims: list = None,
                 t1: list = None, 
                 t2: list = None):

        self.model = SarimnerModel(
                        num_qubits=num_qubits,
                        qubit_frequencies=qubit_frequencies, 
                        anharmonicities=anharmonicities, 
                        rotating_frame_frequencies=rotating_frame_frequencies, 
                        dims=dims,
                        t1=t1, 
                        t2=t2,
                    )
        
        super(SarimnerProcessor, self).__init__(model=self.model)
        self.native_gates = None
        self._default_compiler = SarimnerCompiler
        self.pulse_mode = "discrete"

    def load_circuit(self, qc, schedule_mode="ASAP", compiler=None):
        """
        The default routine of compilation.
        It first calls the :meth:`.transpile` to convert the circuit to
        a suitable format for the hardware model.
        Then it calls the compiler and save the compiled pulses.

        Parameters
        ----------
        qc : :class:`.QubitCircuit`
            Takes the quantum circuit to be implemented.

        schedule_mode: string
            "ASAP" or "ALAP" or None.

        compiler: subclass of :class:`.GateCompiler`
            The used compiler.

        Returns
        -------
        tlist, coeffs: dict of 1D NumPy array
            A dictionary of pulse label and the time sequence and
            compiled pulse coefficients.
        """
        # Choose a compiler and compile the circuit
        if compiler is None and self._default_compiler is not None:
            compiler = self._default_compiler(self.num_qubits, self.params)
        if compiler is not None:
            tlist, coeffs = compiler.compile(
                qc.gates, schedule_mode=schedule_mode
            )
        else:
            raise ValueError("No compiler defined.")
        # Save compiler pulses
        self.set_coeffs(coeffs)
        self.set_tlist(tlist)
        return tlist, coeffs