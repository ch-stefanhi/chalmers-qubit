import numpy as np
from qutip import Qobj

# Ideal gate
def cczs(args):
    theta, phi, gamma = args
    U = np.array(
        [
            [1, 0, 0, 0, 0, 0, 0, 0],
            [0, 1, 0, 0, 0, 0, 0, 0],
            [0, 0, 1, 0, 0, 0, 0, 0],
            [0, 0, 0, 1, 0, 0, 0, 0],
            [0, 0, 0, 0, 1, 0, 0, 0],
            [
                0,
                0,
                0,
                0,
                0,
                -np.exp(-1j * gamma) * np.sin(theta / 2) ** 2 + np.cos(theta / 2) ** 2,
                (1 / 2) * (1 + np.exp(-1j * gamma)) * np.exp(-1j * phi) * np.sin(theta),
                0,
            ],
            [
                0,
                0,
                0,
                0,
                0,
                (1 / 2) * (1 + np.exp(-1j * gamma)) * np.exp(1j * phi) * np.sin(theta),
                -np.exp(-1j * gamma) * np.cos(theta / 2) ** 2 + np.sin(theta / 2) ** 2,
                0,
            ],
            [0, 0, 0, 0, 0, 0, 0, -np.exp(1j * gamma)],
        ],
        dtype="complex",
    )
    return Qobj(U, dims=[[2] * 3, [2] * 3])
