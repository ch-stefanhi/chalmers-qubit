from abc import ABC
from abc import abstractmethod

import numpy as np

from qutip_qip.circuit import QubitCircuit
from qutip_qip.pulse import Pulse
from scipy import interpolate


class Processor(ABC):
    """
    The abstract base class that defines the basic methods and their call
    signature for any processor that derives from this class.

    .. note::

        This is an abstract base class to outline the general API to define
        new processors but cannot be used by itself. It provides a series of
        methods and attributes expected while developing a new processor

    """

    def __init__(
        self,
        model=None,
        compiler=None,
        scheduler=None,
        transpiler=None,
        noise=None,
    ):
        self.model = model
        self.compiler = compiler
        self.scheduler = scheduler
        self.transpiler = transpiler
        self.noise = noise

        # The pulses will be created after a quantum circuit has been loaded
        # with self.load_circuit(qc)
        self.pulses = {}
        self.qc = None

    @abstractmethod
    def load_circuit(self, qc: QubitCircuit):
        """Load a quantum circuit to the processor and generate the pulses.

        Args:
            qc (QubitCircuit): _description_

        Returns:
            _type_: _description_
        """
        qc = self.transpile(qc)

        # Choose a compiler and compile the circuit
        if self.compiler is not None:
            tlist, coeffs = self.compiler.compile(qc.gates)
        else:
            raise ValueError("No compiler defined.")
        
        # Save compiler pulses
        return tlist, coeffs
    
    @abstractmethod
    def transpile(self, qc: QubitCircuit)-> (list, list):
        """
        Converts the circuit to one that can be executed on given hardware if
        the gates in the circuit are not in the supported gates.

        Parameters
        ----------
        qc: :class:`.QubitCircuit`
            The input quantum circuit.

        Returns
        -------
        qc: :class:`QubitCircuit`
            The transpiled quantum circuit.
        """
        # ...
        return qc

    def run_state(self, init_state=None, solver="mesolve", **kwargs):
        # Generate time-dependent Hamiltonian using the pulses

        # if solver == "mesolve":
        #     evo_result = mesolve(
        #         H=noisy_qobjevo, rho0=init_state, tlist=tlist, **kwargs
        #     )

        # elif solver == "mcsolve":
        #     evo_result = mcsolve(
        #         noisy_qobjevo, init_state, tlist=tlist, **kwargs
        #     )

        # return evo_result
        pass

    def plot_pulses(
        self,
        title=None,
        figsize=(12, 6),
        dpi=None,
        show_axis=False,
        rescale_pulse_coeffs=True,
        num_steps=1000,
    ):
        """
        Plot the ideal pulse coefficients.

        Parameters
        ----------
        title: str, optional
            Title for the plot.

        figsize: tuple, optional
            The size of the figure.

        dpi: int, optional
            The dpi of the figure.

        show_axis: bool, optional
            If the axis are shown.

        rescale_pulse_coeffs: bool, optional
            Rescale the hight of each pulses.

        num_steps: int, optional
            Number of time steps in the plot.

        Returns
        -------
        fig: matplotlib.figure.Figure
            The `Figure` object for the plot.

        axis: list of ``matplotlib.axes._subplots.AxesSubplot``
            The axes for the plot.

        Notes
        -----
        :meth:.Processor.plot_pulses` only works for array_like coefficients.
        """
        import matplotlib.pyplot as plt
        import matplotlib.gridspec as gridspec

        color_list = plt.rcParams["axes.prop_cycle"].by_key()["color"]

        # choose labels
        pulse_labels = [{pulse.label: pulse.label for pulse in self.pulses}]

        # If it is a nested list instead of a list of dict
        if isinstance(pulse_labels[0], list):
            for ind, pulse_group in enumerate(pulse_labels):
                pulse_labels[ind] = {i: latex for i, latex in enumerate(pulse_group)}

        # create a axis for each pulse
        fig = plt.figure(figsize=figsize, dpi=dpi)
        grids = gridspec.GridSpec(sum([len(d) for d in pulse_labels]), 1)
        grids.update(wspace=0.0, hspace=0.0)

        full_tlist = [pulse.tlist for pulse in self.pulses if pulse.tlist is not None]
        full_tlist = np.unique(np.sort(np.hstack(tlist)))

        tlist = np.linspace(0.0, full_tlist()[-1], num_steps)
        dt = tlist[1] - tlist[0]

        # make sure coeffs start and end with zero, for ax.fill
        tlist = np.hstack(([-dt * 1.0e-20], tlist, [tlist[-1] + dt * 1.0e-20]))

        coeffs = []
        for pulse in self.pulses:
            coeffs.append(pulse_interpolate(pulse, tlist))

        pulse_ind = 0
        axis = []
        for i, label_group in enumerate(pulse_labels):
            for j, (label, latex_str) in enumerate(label_group.items()):
                try:
                    pulse = self.find_pulse(label)
                    coeff = pulse_interpolate(pulse, tlist)
                except KeyError:
                    coeff = np.zeros(tlist.shape)
                grid = grids[pulse_ind]
                ax = plt.subplot(grid)
                axis.append(ax)
                ax.fill(tlist, coeff, color_list[i], alpha=0.7)
                ax.plot(tlist, coeff, color_list[i])
                if rescale_pulse_coeffs:
                    ymax = np.max(np.abs(coeff)) * 1.1
                else:
                    ymax = np.max(np.abs(coeffs)) * 1.1
                if ymax != 0.0:
                    ax.set_ylim((-ymax, ymax))

                # disable frame and ticks
                if not show_axis:
                    ax.set_xticks([])
                    ax.spines["bottom"].set_visible(False)
                ax.spines["top"].set_visible(False)
                ax.spines["right"].set_visible(False)
                ax.spines["left"].set_visible(False)
                ax.set_yticks([])
                ax.set_ylabel(latex_str, rotation=0)
                pulse_ind += 1
                if i == 0 and j == 0 and title is not None:
                    ax.set_title(title)
        fig.tight_layout()
        return fig, axis


def pulse_interpolate(pulse, tlist):
    """
    A function that calls Scipy interpolation routine. Used for plotting.
    """
    if pulse.tlist is None and pulse.coeff is None:
        coeff = np.zeros(len(tlist))
        return coeff
    if isinstance(pulse.coeff, bool):
        if pulse.coeff:
            coeff = np.ones(len(tlist))
        else:
            coeff = np.zeros(len(tlist))
        return coeff
    coeff = pulse.coeff
    if len(coeff) == len(pulse.tlist) - 1:  # for discrete pulse
        coeff = np.concatenate([coeff, [0]])

    if pulse.spline_kind == "step_func":
        kind = "previous"
    else:
        kind = "cubic"
    inter = interpolate.interp1d(
        pulse.tlist, coeff, kind=kind, bounds_error=False, fill_value=0.0
    )
    return inter(tlist)
