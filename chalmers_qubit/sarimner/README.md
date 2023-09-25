# sarimner

A five qubit processor family with custom definitions for the following that 
builds on top of `qutip-qip`. The following modules contains code that defines
the custom simulator.

    - model.py: A [Model]() to represent the Hamiltonian and controls
    - processor.py: The [Processor]() that runs circuits and simulations
    - compiler.py: A custom [GateCompiler]() to go from circuit -> pulses
    - scheduler.py: We use the default scheduling algorithm in qutip-qip
    - noise.py: Various types of noises in defined as [Noise]() objects
    - transpiler.py: Any code for transpiling circuits
