Installation
============

The main requirement to use this package is `qutip-qip <https://qutip-qip.readthedocs.io/en/stable>`_. The requirements are already specified in the 
setup.py file and you can install the package chalmers_qubit simply by
downloading this folder or cloning this repository and running:

.. code-block:: bash

    pip install .

to get the minimal installation. You can instead use `'.[full]'` to install the package with all optional dependencies, such as matplotlib. Moreover, it might be beneficial to install an editable version. In the editable version, changes to the code are reflected system-wide without requiring a reinstallation.

.. code-block:: bash

    pip install -e '.[full]'

If you do not care about making changes to the source code and just want to try out the package (e.g., from Google Colab), you can do a git+ install with

.. code-block:: bash
    
    pip install git+https://github.com/aqp-mc2-chalmers/chalmers-qubit.git