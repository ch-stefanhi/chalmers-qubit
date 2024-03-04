import numpy as np
from qutip_qip.compiler import GateCompiler, Instruction

__all__ = ["SarimnerCompiler"]

class SarimnerCompiler(GateCompiler):
    """
    Compiler for :class:`.ChalmersQubits`.
    Compiled pulse strength is in the unit of GHz.

    Supported native gates: "RX", "RZ", "CZ", "CCZS".

    Parameters
    ----------
    num_qubits: int
        The number of qubits in the system.

    params: dict
        A Python dictionary contains the name and the value of the parameters.
        See :meth:`.ChalmersQubitsModel` for the definition.

    Attributes
    ----------
    num_qubits: int
        The number of the component systems.

    params: dict
        A Python dictionary contains the name and the value of the parameters,
        such as laser frequency, detuning etc.

    gate_compiler: dict
        The Python dictionary in the form of {gate_name: decompose_function}.
        It saves the decomposition scheme for each gate.
    """

    def __init__(self, num_qubits:int, params:dict):
        super().__init__(num_qubits, params=params)
        self.params = params

        self.gate_compiler.update(
            {
                "RZ": self.rz_compiler,
                "RX": self.rx_compiler,
                "RY": self.ry_compiler,
                "CZ": self.cz_compiler,
                "CCZS": self.cczs_compiler,
                "GLOBALPHASE": self.globalphase_compiler,
            })
        # initialise the global phase for all the qubits to 0.
        self.phase = [0] * num_qubits

    def coupling(self, t, y: bool, args: dict) -> np.ndarray:
        omega_drive = args['drivefreq']
        delta = args['detuning']
        phi = args['phase']
        g = args['coupling_strength']
        if y:
            coeff = g * np.sin(omega_drive * t + phi) * \
                np.exp(1j * delta * t)
        else:
            coeff = g * np.sin(omega_drive * t + phi) * \
                np.exp(- 1j * delta * t)
        return coeff

    def drive_coeff(self, t: np.ndarray, y: bool, args: dict) -> np.ndarray:
        # Amplitude, needs to be optimized to get perfect pi-pulse or pi/2-pulse
        amp = args['amp']
        # DRAG-parameter, needs to be optimized to get no phase errors for a pi/2-pulse
        qscale = args['qscale']
        drive_freq = args['freq']
        phi = args['phase']
        L = args['gatetime']
        sigma = args['sigma']
        Omega_x = amp * np.exp(-pow(t - 0.5 * L, 2)
                               / (2 * pow(sigma, 2)))
        Omega_y = qscale * (t - 0.5 * L) / pow(sigma, 2) * Omega_x

        if y:
            coeff = (Omega_x * np.cos(drive_freq * t - phi)
                     + Omega_y * np.sin(drive_freq * t - phi))
        else:
            coeff = (Omega_x * np.sin(drive_freq * t + phi)
                     - Omega_y * np.cos(drive_freq * t + phi))
        return coeff

    def rz_compiler(self, gate, args):
        """
        Compiler for the Virtual-RZ gate

        Parameters
        ----------
        gate : :obj:`.Gate`:
            The quantum gate to be compiled.
        args : dict
            The compilation configuration defined in the attributes
            :obj:`.GateCompiler.args` or given as a parameter in
            :obj:`.GateCompiler.compile`.

        """
        q = gate.targets[0]  # target qubit
        self.phase[q] += gate.arg_value

    def ry_compiler(self, gate, args):
        """
        Compiler for the RX gate

        Parameters
        ----------
        gate : :obj:`.Gate`:
            The quantum gate to be compiled.
        args : dict
            The compilation configuration defined in the attributes
            :obj:`.GateCompiler.args` or given as a parameter in
            :obj:`.GateCompiler.compile`.

        Returns
        -------
        A list of :obj:`.Instruction`, including the compiled pulse
        information for this gate.
        """

        # Gaussian width in ns
        sigma = 10

        # Gate time in ns
        t_total = 50
        tlist = np.linspace(0, t_total, 5000)

        # target qubit
        q = gate.targets[0]

        # amplitude
        amp = 0.06345380720748712 * gate.arg_value / np.pi

        # parameters
        phase = self.phase[q]
        alpha = self.params['alpha'][q]
        omega_qubit = self.params['wq'][q]
        rotating_freq = self.params['wr'][q]
        omega_drive = rotating_freq - omega_qubit
        args = {'amp': amp,
                'qscale': -.5 / alpha,
                'phase': phase,
                'freq': omega_drive,
                'gatetime': t_total,
                'sigma': sigma}

        # set start and end to zero
        coeff_x = self.drive_coeff(tlist, False, args)
        coeff_y = self.drive_coeff(tlist, True, args)

        pulse_info = [
            # (control label, coeff)
            ("X" + str(q), coeff_x),
            ("Y" + str(q), coeff_y)
        ]

        return [Instruction(gate, tlist, pulse_info, t_total)]

    def rx_compiler(self, gate, args):
        """
        Compiler for the RX gate

        Parameters
        ----------
        gate : :obj:`.Gate`:
            The quantum gate to be compiled.
        args : dict
            The compilation configuration defined in the attributes
            :obj:`.GateCompiler.args` or given as a parameter in
            :obj:`.GateCompiler.compile`.

        Returns
        -------
        A list of :obj:`.Instruction`, including the compiled pulse
        information for this gate.
        """

        # Gaussian width in ns
        sigma = 10

        # Gate time in ns
        t_total = 50
        tlist = np.linspace(0, t_total, 5000)

        # target qubit
        q = gate.targets[0]

        # amplitude
        amp = 0.06345380720748712 * gate.arg_value / np.pi

        # parameters
        phase = self.phase[q]
        alpha = self.params['alpha'][q]
        omega_qubit = self.params['wq'][q]
        rotating_freq = self.params['wr'][q]
        omega_drive = rotating_freq - omega_qubit
        args = {'amp': amp,
                'qscale': -.5 / alpha,
                'phase': phase,
                'freq': omega_drive,
                'gatetime': t_total,
                'sigma': sigma}

        # set start and end to zero
        coeff_x = self.drive_coeff(tlist, True, args)
        coeff_y = self.drive_coeff(tlist, False, args)

        pulse_info = [
            # (control label, coeff)
            ("X" + str(q), coeff_x),
            ("Y" + str(q), coeff_y)
        ]

        return [Instruction(gate, tlist, pulse_info, t_total)]

    def cz_compiler(self, gate, args):
        """
        Compiler for CZ gate.

        Parameters
        ----------
        gate : :obj:`.Gate`:
            The quantum gate to be compiled.
        args : dict
            The compilation configuration defined in the attributes
            :obj:`.GateCompiler.args` or given as a parameter in
            :obj:`.GateCompiler.compile`.

        Returns
        -------
        A list of :obj:`.Instruction`, including the compiled pulse
        information for this gate.
        """
        q1 = gate.controls[0]
        q2 = gate.targets[0]

        omega1 = self.params["wq"][q1]
        omega1_rot = self.params["wr"][q1]
        omega2 = self.params["wq"][q2]
        omega2_rot = self.params["wr"][q2]
        alpha1 = self.params["alpha"][q1]
        coupling_strength = self.params["cpl_matrix"][q1,q2]

        t_total = np.sqrt(2) * np.pi / coupling_strength
        tlist = np.linspace(0, t_total, 5000)

        # Check if qubit q1 and q2 are coupled together
        if coupling_strength is 0:
            raise ValueError(f"Qubit {q1} and {q2} are not coupled together.")

        args = {'drivefreq': abs(omega1 + alpha1 - omega2),
                'detuning': (omega1_rot - omega2_rot),
                'phase': 0,
                'coupling_strength': coupling_strength}

        pulse_info = [("ab" + str(q1) + str(q2), self.coupling(tlist, True, args)),
                      ("ba" + str(q1) + str(q2), self.coupling(tlist, False, args))]

        return [Instruction(gate, tlist, pulse_info, t_total)]

    def cczs_compiler(self, gate, args):
        """
        Compiler for CCZS gate.

        Parameters
        ----------
        gate : :obj:`.Gate`:
            The quantum gate to be compiled.
        args : dict
            The compilation configuration defined in the attributes
            :obj:`.GateCompiler.args` or given as a parameter in
            :obj:`.GateCompiler.compile`.

        Returns
        -------
        A list of :obj:`.Instruction`, including the compiled pulse
        information for this gate.
        """
        q1 = gate.targets[0]  # this is the control qubit
        q2 = gate.targets[1]
        q3 = gate.targets[2]
        theta, phi, gamma = gate.arg_value

        # gate time
        t_total = np.pi / self.g
        tlist = np.linspace(0, t_total, 5000)

        omega1 = self.params["wq"][q1]
        omega2 = self.params["wq"][q2]
        omega3 = self.params["wq"][q3]
        omega1_rot = self.params["wr"][q1]
        omega2_rot = self.params["wr"][q2]
        omega3_rot = self.params["wr"][q3]
        alpha1 = self.params["alpha"][q1]

        omega1_drive = abs(omega1 + alpha1 - omega2)
        omega2_drive = abs(omega1 + alpha1 - omega3)

        # We do simultaneous drive
        args1 = {"drivefreq": omega1_drive,
                 "detuning": (omega1_rot - omega2_rot),
                 "phase": 0}
        args2 = {"drivefreq": omega2_drive,
                 "detuning": (omega1_rot - omega3_rot),
                 "phase": phi - np.pi}

        pulse_info = [
            ("ab" + str(q1) + str(q2), self.coupling(tlist, True, args1)),
            ("ba" + str(q1) + str(q2), self.coupling(tlist, False, args1)),
            ("ab" + str(q1) + str(q3), self.coupling(tlist, True, args2)),
            ("ba" + str(q1) + str(q3), self.coupling(tlist, False, args2))
        ]

        return [Instruction(gate, tlist, pulse_info)]

    def globalphase_compiler(self, gate, args):
        """
        Compiler for the GLOBALPHASE gate
        """
        self.global_phase += gate.arg_value