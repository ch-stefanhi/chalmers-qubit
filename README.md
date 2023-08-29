# chalmers-qubit

A realistic model of the Chalmers device that can be used to run noisy simulations of quantum circuits
following [qutip-qip](https://qutip-qip.readthedocs.io/en/stable/). 

# Installation

The package requires [qutip-qip](https://qutip-qip.readthedocs.io/en/stable/) which is based on [qutip](https://qutip-qip.readthedocs.io/en/stable/): The
Quantum Toolbox in Python. The requirements are already specified in the `setup.py` file and you can install `chalmers-qubit` simply by downloading this package or cloning
it and running the following command:
```
python setup.py install
```

It might be beneficial to install an editable version which can be further customized using the following command:

```
python setup.py develop
```

# Usage

The usage of the package follows [qutip-qip](https://qutip-qip.readthedocs.io/en/stable/)
where first, a quantum circuit is defined using [`qutip-qip`](https://qutip-qip.readthedocs.io/en/stable/qip-simulator.html)
and then run on the custom `ChalmersProcessor`. The custom processor is defined
in `chalmers_qubit.processor` and can be initialized with noise parameters,
e.g., T1 and T2 values for each qubit.

Note that only gates with compilation instructions in `chalmers_qubit/compiler.py` will work.

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

It is also possible to import QASM circuits.

# Development

In order to add new custom pulses or modify the device, edit the processor, 
compiler or scheduler following the tutorials and detailed instructions in [qutip-qip](https://qutip-qip.readthedocs.io/en/stable/).

The [tutorials](https://qutip.org/qutip-tutorials/) show examples of how to
customize the processor. If you have installed the package in the develop mode, any changes to the processor, e.g., adding a new gate will be reflected immediately system-wide without requiring a reinstallation of the package.

# Support

This package was built from contributions by Pontus Vikstål, Kamanasish Debnath
and Shahnawaz Ahmed.

Contact shahnawaz.ahmed95@gmail.com or anton.frisk.kockum@chalmers.se 
for help and support.