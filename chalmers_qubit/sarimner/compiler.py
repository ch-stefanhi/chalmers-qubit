import numpy as np

from qutip_qip.compiler import Instruction
from qutip_qip.operations import Gate

from ..base.compiler import Compiler


class SarimnerCompiler(Compiler):
    """
    Compiler for the Chalmers processor.

    Supported native gates: "RX", "RZ", "CZ", "CCZS".
    """

    def __init__(self, model, g: float = None):
        num_qubits = model.num_qubits
        if g is None:
            g = 2 * np.pi * 1 * 1e-3

        super().__init__(num_qubits)

        self.g = g
        self.compilers = {
            "RX": self.rx_compiler,
        }

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
        # sigma = 10

        # Gate time in ns
        # t_total = 50
        # tlist = np.linspace(0, t_total, 5000)

        # target qubit
        # q = gate.targets[0]
        # set start and end to zero
        # coeff_x = np.cos(4.2 * tlist)
        # coeff_y = np.cos(4.2 * tlist) * 0

        # pulse_info = [
        #     # (control label, coeff)
        #     ("sx" + str(q), coeff_x),
        #     ("sy" + str(q), coeff_y),
        # ]

        # return [Instruction(gate, tlist, pulse_info, t_total)]
