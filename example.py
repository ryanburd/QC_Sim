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