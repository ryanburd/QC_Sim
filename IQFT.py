# Import the QPU simulator which creates the quantum circuit and contains all gate operations. Import numpy for using pi in rotation gates.
import QPUsimulator as QPU
import numpy as np

# Create a quantum circuit with a chosen number of qubits. The same number of classica bits will be created for qubit measurement output. The output of each qubit is stored in its similarly indexed classical bit.
numQubits = 3
circuit = QPU.Circuit(numQubits)

circuit.H([0, 1, 2])
circuit.P(0, theta=5*np.pi/4)
circuit.P(1, theta=5*np.pi/2)
circuit.P(2, theta=5*np.pi)

# for qubit in range(numQubits):
#     if qubit == numQubits-qubit-1:
#         break
#     if qubit >= 0.5*numQubits:
#         break
#     circuit.SWAP(qubit, numQubits-qubit-1)

# for qubit in range(numQubits):
#     for control in range(qubit-1, -1, -1):
#         circuit.CP([control], qubit, theta=-np.pi/2**(qubit-control))
#     circuit.H(qubit)

circuit.IQFT()

for qubit in range(numQubits):
    circuit.measure(qubit)

# Display the quantum circuit.
circuit.display_circuit()

# Choose the number of shots for the circuit. Run the circuit 'shots' times and obtain the results as a list. Set hist=True to create a histogram of the results.
shots = 1024
results = circuit.run(shots, hist=True)