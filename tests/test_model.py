import pytest
from qutip import destroy, tensor, qeye
from chalmers_qubit.sarimner.model import SarimnerModel


@pytest.fixture
def sarimner_model():
    num_qubits = 3
    qubit_frequencies = [4.2e9] * num_qubits
    anharmonicities = [200e6] * num_qubits
    resonance_frequencies = [4.1e9] * num_qubits
    zz_crosstalk_static = {(0, 1): 0.5, (1, 2): 0.1}

    model = SarimnerModel(
        num_qubits,
        qubit_frequencies,
        anharmonicities,
        resonance_frequencies,
        zz_crosstalk_static=zz_crosstalk_static,
    )
    return model


def test_zz_coupling(sarimner_model):
    model = sarimner_model

    for (qubit1, qubit2), coupling_strength in model.zz_crosstalk.items():
        sigma_z1 = tensor(
            [
                qeye(model.dims[qubit1])
                if i != qubit1
                else destroy(model.dims[qubit1]).dag() * destroy(model.dims[qubit1])
                for i in range(model.num_qubits)
            ]
        )
        sigma_z2 = tensor(
            [
                qeye(model.dims[qubit2])
                if i != qubit2
                else destroy(model.dims[qubit2]).dag() * destroy(model.dims[qubit2])
                for i in range(model.num_qubits)
            ]
        )
        expected_zz_term = coupling_strength * sigma_z1 * sigma_z2

        zz_term, _ = model._drift[f"ZZ{qubit1}{qubit2}"]

        assert zz_term == expected_zz_term


def test_drift_hamiltonian(sarimner_model):
    num_qubits = 3
    for qubit in range(num_qubits):
        destroy_op = destroy(3)  # Assuming 3-level system for each qubit
        expected_freq_diff_term = (
            (
                sarimner_model.qubit_frequencies[qubit]
                - sarimner_model.resonance_frequencies[qubit]
            )
            * destroy_op.dag()
            * destroy_op
        )
        expected_anharmonicity_term = (
            (sarimner_model.anharmonicities[qubit] / 2.0)
            * destroy_op.dag() ** 2
            * destroy_op**2
        )

        freq_diff_term, _ = sarimner_model._drift[r"\Delta \omega" + str(qubit)]
        anharmonicity_term, _ = sarimner_model._drift[r"\alpha" + str(qubit)]

        assert freq_diff_term == expected_freq_diff_term
        assert anharmonicity_term == expected_anharmonicity_term
