from typing import Optional
from qutip import propagator
from qutip_qip.noise import Noise
from qutip_qip.device import Processor, Model
from qutip_qip.compiler import GateCompiler

from chalmers_qubit.sarimner.compiler import SarimnerCompiler

class SarimnerProcessor(Processor):
    """
    Initialize a new SarimnerProcessor instance with a quantum model, an optional compiler, and noise models.

    Parameters
    ----------
    model : Model
        The quantum model that defines the physical properties and capabilities of the processor.
    compiler : GateCompiler, optional
        The compiler used to translate quantum gates into executable operations. If not provided,
        a default compiler specific to the model (SarimnerCompiler) is instantiated and used.
    noise : list, optional
        A list of noise models to be added to the processor. Each element in the list should be compatible
        with the processor's noise handling methods.

    Attributes
    ----------
    model : Model
        The model of the quantum processor, storing physical properties.
    _default_compiler : GateCompiler
        Holds the compiler instance being used, either the provided one or a default SarimnerCompiler.
    native_gates : NoneType
        Initially set to None, to be configured with the gate set natively supported by the processor.
    pulse_mode : str
        Mode of pulse operation, set to 'discrete' by default, indicating how pulses are treated in the processor.
    global_phase : float
        The global phase of the quantum state managed by the processor, initialized to 0.
    """

    def __init__(self,
                 model:Model,
                 compiler:Optional[GateCompiler] = None,
                 noise:Optional[list] = None):

        self.model = model

        if compiler is None:
            self._default_compiler = SarimnerCompiler(model=model)
        else:
            self._default_compiler = compiler

        if noise is not None:
            for elem in noise:
                self.add_noise(elem)

        super(SarimnerProcessor, self).__init__(model=self.model)
        self.native_gates = None
        self.pulse_mode = "continuous"
        self.spline_kind = "cubic"
        self.global_phase = 0

    def add_noise(self, noise):
        """
        Add a noise object to the processor.

        Parameters
        ----------
        noise : :class:`.Noise`
            The noise object defined outside the processor.
        """
        if isinstance(noise, Noise):
            self.model._add_noise(noise)
        else:
            raise TypeError("Input is not a Noise object.")

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
            compiler = self._default_compiler
        if compiler is not None:
            tlist, coeffs = compiler.compile(
                qc.gates, schedule_mode=schedule_mode
            )
        else:
            raise ValueError("No compiler defined.")

        # Update global phase
        self.global_phase = compiler.global_phase

        # Save compiler pulses
        self.set_coeffs(coeffs)
        self.set_tlist(tlist)
        return tlist, coeffs

    def run_state(
        self,
        init_state=None,
        analytical=False,
        states=None,
        noisy=True,
        solver="mesolve",
        qc=None,
        **kwargs):
        if qc is not None:
            self.load_circuit(qc)
        return super().run_state(init_state,analytical,states,noisy,solver,**kwargs)

    def run_propagator(self, qc=None, noisy=False, **kwargs):
        """
        NOT WORKING AFTER QUTIP UPDATE TO 5.0
        Parameters
        ----------
        qc: :class:`qutip.qip.QubitCircuit`, optional
            A quantum circuit. If given, it first calls the ``load_circuit``
            and then calculate the evolution.
        states: :class:`qutip.Qobj`, optional
         Old API, same as init_state.
        **kwargs
           Keyword arguments for the qutip solver.
        Returns
        -------
        evo_result: :class:`qutip.Result`
            If ``analytical`` is False,  an instance of the class
            :class:`qutip.Result` will be returned.
            If ``analytical`` is True, a list of matrices representation
            is returned.
        """
        if qc is not None:
            self.load_circuit(qc)
        # construct qobjevo for unitary evolution
        noisy_qobjevo, c_ops = self.get_qobjevo(noisy=noisy)

        # time steps
        tlist = self.get_full_tlist()
        H = noisy_qobjevo.to_list()

        # Compute drift Hamiltonians
        H_drift = 0
        drift = self._get_drift_obj()
        for drift_ham in drift.drift_hamiltonians:
            H_drift += drift_ham.get_qobj(self.dims)
        H[0] = H_drift

        # compute the propagator
        evo_result = propagator(H=H, t=tlist, c_ops=c_ops, **kwargs)
        return evo_result
