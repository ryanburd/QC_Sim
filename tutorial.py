# Import the QPU simulator which creates the quantum circuit and contains all gate operations. Import matplotlib for making a histogram of the results.
import QPUsimulator as QPU
import numpy as np
import matplotlib.pyplot as plt

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

# Controlled gates require (control, target) qubits to act on.
circuit.CX(0, 1)
circuit.CY(1, 2)
circuit.CZ(0, 1)
circuit.CRX(1, 0, np.pi)
circuit.CRY(2, 1, np.pi)
circuit.CRZ(2, 0, np.pi)
circuit.CU(0, 1, 0, np.pi, 0)

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