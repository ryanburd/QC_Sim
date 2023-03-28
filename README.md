# Quantum Computer Simulator

This repository is for a quantum computer simulator that can run quantum circuits created by the user.

The user can:
- Create quantum circuits consisting of a chosen number of qubits and classical bits (for storing the measurement results of qubits).
- Apply Pauli-X, -Y, -Z, Hadamard, phase, pi/8, Rx, Ry, Rz, U, controlled-X, -Y, -Z, -Rx, -Ry, -Rz, -U, and SWAP gates to the circuit.
- Create a diagram of the circuit.
- Measure the final state of the quantum circuit after measurement in the computational basis.
- Create a histogram of the result of many shots.

## Tutorial

A tutorial of the QPU simulator can be found in tutorial.py, along with the circuit diagram, tutorial_diagram.png, and the results histogram, tutorial_hist.png. The script is also provided below:

```py
# Import the QPU simulator which creates the quantum circuit and contains all gate operations. Import numpy for using pi in rotation gates.
import QPUsimulator as QPU
import numpy as np

# Create a quantum circuit with a chosen number of qubits. The same number of classica bits will be created for qubit measurement output. The output of each qubit is stored in its similarly indexed classical bit.
numQubits = 3
circuit = QPU.Circuit(numQubits)

# Single qubit gates require a target qubit to act on.
circuit.X(0)
circuit.Y(1)
circuit.Z(2)
circuit.H(0)
circuit.S(1)
circuit.T(2)

# Rotation gates also require angles: theta for Rx, Ry, and Rz; theta, phi, lambda for U
circuit.RX(0, np.pi)
circuit.RY(1, np.pi)
circuit.RZ(2, np.pi)
circuit.U(1, 0, np.pi, 0)

# Add a barrier to separate the circuit into segments for visual purposes.
circuit.barrier()

# Controlled gates require ([controls], target) qubits to act on. Note that the controls must be passed as a list, even if there is only one control qubit.
circuit.CX([0], 1)
circuit.CY([1], 2)
circuit.CZ([0, 1], 2)
circuit.CRX([1], 0, np.pi)
circuit.CRY([2], 1, np.pi)
circuit.CRZ([2, 1], 0, np.pi)
circuit.CU([0], 1, 0, np.pi, 0)

# SWAP gates require two qubits to swap states.
circuit.SWAP(0, 2)

# Measurements (in the computational basis) require a qubit to measure. The output is stored in the same indexed classical bit.
circuit.measure(0)
circuit.measure(1)
circuit.measure(2)

# Choose the number of shots for the circuit. Run the circuit 'shots' times and obtain the results as a list. Set hist=True to create a histogram of the results.
shots = 1024
results = circuit.run(shots, hist=True)

# Display the quantum circuit.
circuit.display_circuit()
```

## Notes
Current limitations being addressed in ongoing updates:
- Only pure states can be considered since the circuit math is completed with state vectors. This will be updated to using density operators to represent the qubits' states so that mixed states can be considered as well.

Future updates will include common algorithms preprogrammed for easy use.
