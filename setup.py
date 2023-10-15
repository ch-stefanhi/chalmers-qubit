import setuptools

setuptools.setup(
    name="chalmers_qubit",
    version="0.1",
    author="Pontus Vikstål, Kamanasish Debnath, Shahnawaz Ahmed",
    description="A model of the Chalmers device to be used with qutip-qip",
    packages=setuptools.find_packages(),
    install_requires=[
        "qutip-qip",
    ],
    python_requires=">=3.9",
)
