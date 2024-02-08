# chalmers-qubit

A simulation framework for Chalmers devices that can be used to simulate the
running of quantum algorithms with realistic noise. We follow [qutip-qip](https://qutip-qip.readthedocs.io/en/stable/)
to build a processor that can take in a quantum circuit (e.g., a QASM cicruit)
and performs a master equation simulation adding noise such as T1 and T2. It is
also possible to perform a Monte-Carlo trajectory simulation and customize the
processor to add various types of noise such as [ZZCrossTalk](https://qutip-qip.readthedocs.io/en/latest/apidoc/qutip_qip.noise.html#qutip_qip.noise.ZZCrossTalk).

The package is under development and testing. 

# Installation

The main requirement to use this package is [qutip-qip](https://qutip-qip.readthedocs.io/en/stable/)
based on [qutip](https://qutip-qip.readthedocs.io/en/stable/): The
Quantum Toolbox in Python. The requirements are already specified in the 
`setup.py` file and you can install the package `chalmers_qubit` simply by
downloading this folder or cloning this repository and running:

```
pip install .
```

It might be beneficial to install an editable version. In the editable version,
changes to the code are reflected system-wide without requiring a reinstallation.

```
pip install -e .
```

If you do not care about making changes to the source code and just want to
try out the package (e.g., from Google Colab), you can do a git+ install with

```
pip install git+https://github.com/aqp-mc2-chalmers/chalmers-qubit.git
```

# Usage

The usage of the package follows [qutip-qip](https://qutip-qip.readthedocs.io/en/stable/)
where first, a quantum circuit is defined using [`qutip-qip`](https://qutip-qip.readthedocs.io/en/stable/qip-simulator.html)
and then run on one of the custom Chalmers processors, e.g., the processor
called sarimner. The custom processor is defined 
in `chalmers_qubit.sarimner.processor` and can be initialized with a model,
compiler, noise, scheduler and a transpiler. 

Note that only gates with compilation instructions in `chalmers_qubit/sarimner/compiler.py`
will work for this particular processor.

Notebooks exploring the usage of the simulator is available in `doc/notebooks`. 

```
import numpy as np
from qutip import basis, tensor
from qutip_qip.circuit import QubitCircuit
from chalmers_qubit.sarimner.processor import SarimnerProcessor

# Define a circuit
qc = QubitCircuit(2)
qc.add_gate("X", targets=1)
qc.add_gate("SNOT", targets=0)

# Initial state to run a simulation
# The default assumptions is that each transmon is a qudit with 3 levels
init_state = tensor(basis(3, 0), basis(3, 0))

# Define a Model with model parameters
model = SarimnerModel(num_qubits = 2,
    wq=[1.0, 2.0],
    wr=[2.0, 3.0],
    alpha=[100, 200],
    t1=[1, 2],
    t2=[20, 30],
    zz_crosstalk_static=None,
)

# Define a compiler with gate parameters and times
compiler = SarimnerCompiler(model, g=2.0)

# Initialize the processor
processor = SarimnerProcessor(model=model, compiler=compiler)

# Load the circuit that generates the pulses and run the simulation
processor.load_circuit(qc)

tlist = np.linspace(0, 20, 300)
result = processor.run_state(init_state, tlist=tlist)
print("Final state", result.states[-1])

# Run the same circuit but with mcsolve (the initial state is a pure state)
init_state = tensor(basis(3, 0), basis(3, 0))
result = processor.run_state(psi0, tlist=tlist, solver="mcsolve")
print("Final state", result.states[-1])
```

It is also possible to import QASM circuits.

# Development

In order to add new custom pulses or modify the device, edit the processor, 
compiler or scheduler following the tutorials and detailed instructions in
[qutip-qip](https://qutip-qip.readthedocs.io/en/stable/).

The [tutorials](https://qutip.org/qutip-tutorials/) show examples of how to
customize the processor. If you have installed the package in the develop mode,
any changes to the processor, e.g., adding a new gate will be reflected
immediately system-wide without requiring a reinstallation of the package.

# Support

This package was built from contributions by Pontus Vikstål, Kamanasish Debnath
and Shahnawaz Ahmed.

Contact shahnawaz.ahmed95@gmail.com or anton.frisk.kockum@chalmers.se 
for help and support.
