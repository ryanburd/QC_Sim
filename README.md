# Quantum Computer Simulator

This repository is for a quantum computer simulator that can run quantum circuits created by the user.

The user can:
- Create quantum circuits consisting of a chosen number of qubits and classical bits (for storing the measurement results of qubits).
- Apply Pauli-X, -Y, -Z, Hadamard, phase, Rx, Ry, Rz, U, controlled-X, -Y, -Z, =phase, -Rx, -Ry, -Rz, -U, and SWAP gates to the circuit.
- Create a diagram of the circuit.
- Measure the final state of the quantum circuit after measurement in the computational basis.
- Create a histogram of the result of many shots.

Preprogrammed algorithms available to the user (tutorials are provided for each):
- Deutsch-Jozsa
- Quantum Fourier Transform (QFT)
- Inverse QFT (IQFT)
- Quantum phase estimation (QPE)

## Tutorial

A tutorial of the QPU simulator can be found in General_tutorial.py, along with the circuit diagram, General_tutorial_diagram.png, and the results histogram, General_tutorial_hist.png. The script is also provided below:

```py
# This tutorial shows how to use the QPU simulator to make your own circuits using quantum logic gates, display a diagram of your circuit, and run your circuit for many shots to obtain the results.
#
# In addition to the quantum logic gates shown here, some common algorithms have been preprogrammed for your convenience. See the tutorials on those algorithms to learn how to use them.

################################################################################

# Import the QPU simulator which creates the quantum circuit and contains all gate operations. Import numpy for using pi in rotation gates.
import Simulator as QPU
import numpy as np

# Create a quantum circuit with a chosen number of qubits. The same number of classical bits will be created for qubit measurement output. The output of each qubit is stored in its similarly indexed classical bit.
numQubits = 3
circuit = QPU.Circuit(numQubits)

# Single qubit gates require a target qubit to act on. For convenience, if you would like to apply the same single qubit gate to multiple qubits at the same time, all the targets can be passed in a list. Commented out gates are provided to show you other similar gates that can be applied.
circuit.X([0, 1, 2])
# circuit.Y(1)
# circuit.Z(2)
circuit.H(0)

# The general phase gate 'P' also requires an angle theta. For convenience, 'S' and 'T' gates have been preprogrammed, which correspond the P gates with theta = pi/2 and pi/4, respectively. The 'S' and 'T' gates do not take a theta argument.
circuit.P(1, theta=np.pi)
circuit.S(2)
# circuit.T(2)

# Rotation gates also require angles: theta for Rx, Ry, and Rz; theta, phi, lambda for U. Commented out gates are provided to show you other similar gates that can be applied.
circuit.RX(0, theta=np.pi)
# circuit.RY(1, theta=np.pi)
# circuit.RZ(2, theta=np.pi)
circuit.U(1, theta=0, phi=np.pi, lambd=0)

# Add a barrier to separate the circuit into segments for visual purposes.
circuit.barrier()

# Controlled gates require ([controls], target) qubits to act on. Note that the controls must be passed as a list, even if there is only one control qubit. Commented out gates are provided to show you other similar gates that can be applied.
circuit.CX([0], 1)
circuit.CY([0, 1], 2)
# circuit.CZ([0, 1], 2)
circuit.CP([1], 0, theta=np.pi)
circuit.CRX([2, 0], 1, theta=np.pi)
# circuit.CRY([0, 1], 2, theta=np.pi)
# circuit.CRZ([2, 1], 0, theta=np.pi)
circuit.CU([0], 1, theta=0, phi=np.pi, lambd=0)

# SWAP gates require two qubits to swap states.
circuit.SWAP(0, 2)

# Measurements (in the computational basis) require a qubit (or list of qubits) to measure. The output is stored in the same indexed classical bit.
circuit.measure([0, 1, 2])

# Display the quantum circuit. Each qubit is labeled with a Q and an index and has a horizontal line extending to the right, representing the qubit's wire in the circuit. Classical bits are labeled with a C and index and have a double horizontal line extending to the right representing its wire in the circuit. Quantum logic gates are presented by boxes with a letter indicating the gate type. Control qubits have a filled circle with a vertical line extending to the target qubit's gate. SWAP gates are shown as X's for the qubits to be swapped with a vertical line connecting them. The output of a qubit measurement is shown with an arrow pointing to an open circle on the wire for the classical bit that stores the output.
circuit.display_circuit()

# Choose the number of shots for the circuit. Run the circuit 'shots' times and obtain the results as a list. Set hist=True to create a histogram of the results.
shots = 1024
results = circuit.run(shots, hist=True)
```

## Notes
Current limitations being addressed in ongoing updates:
- Only pure states can be considered since the circuit math is completed with state vectors. This will be updated to using density operators to represent the qubits' states so that mixed states can be considered as well.
- Preprogrammed algorithms must use the first n qubits indexed 0 to n-1. For the Deutsch-Jozsa and QPE algorithms, the extra qubit must be after the input qubits. Future updates will allow the user to identify which qubits within the circuit they would like to apply an algorithm to.

Future updates will include common algorithms preprogrammed for easy use.