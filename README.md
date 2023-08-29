# chalmers-qubit

A model for the Chalmers qubit to run realistic simulations of quantum circuits
using [qutip-qip](https://qutip-qip.readthedocs.io/en/stable/). 

# Installation

The package requires [qutip-qip](https://qutip-qip.readthedocs.io/en/stable/)
that is based on [qutip](https://qutip-qip.readthedocs.io/en/stable/): The
Quantum Toolbox in Python installed with

```
pip install qutip-qip
```

Then you can installing `chalmers-qubit` by downloading this package or cloning
it and running the following command that installs an editable version of the
package. The editable version can be modified in-place such that the changes
apply immediately to the code without requiring a re-install.

```
pip install -e .
```



# Usage

The usage of the package follows [qutip-qip](https://qutip-qip.readthedocs.io/en/stable/)
where first a quantum circuit is defined using [`qutip-qip`](https://qutip-qip.readthedocs.io/en/stable/qip-simulator.html)
and then run on the custom `ChalmersProcessor`. The custom processor is defined
in `chalmers_qubit.processor` and can be initialized with noise parameters,
e.g., T1 and T2 values for each qubit.

```
import numpy as np
from qutip import basis
from qutip_qip.circuit import QubitCircuit
from chalmers_qubit.device.processor import ChalmersProcessor

# Define a circuit
qc = QubitCircuit(3)
qc.add_gate("X", targets=2)
qc.add_gate("SNOT", targets=0)

# Run gate-level simulation
init_state = basis([2,2,2], [0,0,0])
ideal_result = qc.run(init_state)

# Run pulse-level simulation
processor = ChalmersProcessor(t2=30)
processor.load_circuit(qc)
tlist = np.linspace(0, 20, 300)
result = processor.run_state(init_state, tlist=tlist)
print("Final state", result.states[-1])
```

It is also possible to import QASM circuits. However, note that only gates with
compilation instructions in `chalmers_qubit/compiler.py` will work.

# Development

In order to add new custom pulses or modify the device, edit the processor, 
compiler or scheduler following the tutorials and detailed instruction in [qutip-qip](https://qutip-qip.readthedocs.io/en/stable/).

The [tutorials](https://qutip.org/qutip-tutorials/) show examples of how to
customize the `ChalmersProcessor`. Once changes are made, any script will automatically
run a quantum circuit with the new customized `ChalmersProcessor`

# Support

This package was built from contributions by Pontus Vikstål, Kamanasish Debnath
and Shahnawaz Ahmed.

Contact shahnawaz.ahmed95@gmail.com or anton.frisk.kockum@chalmers.se 
for help and support.
