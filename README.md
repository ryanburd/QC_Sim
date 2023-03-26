# Quantum Computer Simulator

This repository is for a quantum computer simulator that can run quantum circuits created by the user.

The user can:
- Create quantum circuits consisting of a chosen number of qubits and classical bits (for storing the measurement results of qubits).
- Apply Pauli-X, -Y, -Z, Hadamard, phase, pi/8, xontrolled-X, -Y, -Z, SWAP, and Toffoli gates to the circuit.
- Create a diagram of the circuit.
- Measure the final state of the quantum circuit after measurement in the computational basis.

## Example

An example use of the QPU simulator can be found in example.py. The script is also provided below:

```
# Import the QPU simulator which creates the quantum circuit and contains all gate operations. Import matplotlib for making a histogram of the results.
import QPUsimulator as QPU
import matplotlib.pyplot as plt

# Create a quantum circuit with a chosen number of qubits and classical bits (for qubit measurement output).
numQubits = 3
numCbits = numQubits
circuit = QPU.Circuit(numQubits, numCbits)

# Apply Hadamard gates to each qubit. The passed number is the qubit to apply the gate to.
circuit.H(0)
circuit.H(1)
circuit.H(2)

# Add a barrier to end this section (for display purposes only).
circuit.barrier()

# Apply controlled-Z gates to the first qubit controlled by the other 2 qubits. The first passed number is the control qubit. The second number is the target qubit.
circuit.CZ(2, 0)
circuit.CZ(1, 0)

# Add a barrier to end this section.
circuit.barrier()

# Apply a Hadamard and Pauli X gate to each qubit. Then control a Z gate on the first qubit with both remaining qubits acting as controls for the one Z gate. Apply a Pauli X and Hadamard gate to each qubit.
circuit.H(0)
circuit.H(1)
circuit.H(2)
circuit.X(0)
circuit.X(1)
circuit.X(2)
# The first 2 passed qubits are the control qubits. The third passed number is the target qubit.
circuit.Toff(1, 2, 0)
circuit.X(0)
circuit.X(1)
circuit.X(2)
circuit.H(0)
circuit.H(1)
circuit.H(2)

# Add a barrier.
circuit.barrier()

# Choose the number of shots for the circuit. Create a list to store the result of each shot.
shots = 1
results = []

# For each shot, measure the circuit, obtain the result from the classical bits, and store in the results list.
for shot in range(shots):
    # The first passed number is the qubit to measure. The second passed number is the classicla bit to store the result in.
    circuit.measure(0, 0)
    circuit.measure(1, 1)
    circuit.measure(2, 2)

    result0 = circuit.cbits[0].state
    result1 = circuit.cbits[1].state
    result2 = circuit.cbits[2].state

    results.append(str(result0)+str(result1)+str(result2))

# Display the quantum circuit.
circuit.display_circuit()

# Create a histogram of the results.
plt.hist(results)
plt.title('Histogram of results of %i shots of the quantum circuit'%shots)
plt.xlabel('Quantum circuit state')
plt.show()
```

## Notes
Current limitations being addressed in ongoing updates:
- Only pure states can be considered since the circuit math is completed with state vectors. This will be updated to using density operators to represent the qubits' states so that mixed states can be considered as well.
- The circuit can only be run once (single shot). This is due to the circuit's state being updated when a gate is applied to the circuit. This will be fixed by only creating a list of gates when gates are applied by the user, and a new run() function will apply the math of all the gates to the circuit at the end, which can be easily looped for multiple shots of the circuit.
- New gates will be added soon: rotation gates, general U gate, and controlled versions of any single qubit gate with multiple control qubits.
