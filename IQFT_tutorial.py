# This tutorial shows how to use the inverse Quantum Fourier Transform (IQFT) algorithm with the QPU simulator.

################################################################################

# Import the QPU simulator which creates the quantum circuit and contains all gate operations.
# Import IQFT from Algorithms if you want to access the algorithm directly, rather than through QPUsimulator. Computationally, using the algorithm through Simulator.py or Algorithms.py is the same. When displaying the circuit, Simulator.py simplifies the diagram to show a general "IQFT" block over all the qubits involved to represent the algorithm; Algorithms.py shows all the individual gates completed in the algorithm.
# Import numpy for using np.pi.
import Simulator as QPU
from Algorithms import IQFT
import numpy as np

# Create a quantum circuit with a chosen number of qubits. The same number of classical bits will be created for qubit measurement output. The output of each qubit is stored in its similarly indexed classical bit.
numQubits = 3
circuit = QPU.Circuit(numQubits)

# Initialize the circuit to a state of your choosing in the Fourier basis. For example, the state |5> = |101> in the Fourier basis is prepared by applying an H to all qubits and then the phases 5pi/4, 5pi/2, and 5pi, respectively.
circuit.H(range(numQubits))
circuit.P(0, theta=5*np.pi/4)
circuit.P(1, theta=5*np.pi/2)
circuit.P(2, theta=5*np.pi)

# To apply the algorithm to your circuit, you can apply it just as you would a qubit gate. You may provide the qubit indices to perform the IQFT on as a list using algQubits. Note that the qubits involved must be sequential and ordered from least significant (lowest index) to most significant (highest index). To perform IQFT on all qubits within the circuit, you may leave this argument as the default by passing nothing.
circuit.IQFT(algQubits=range(numQubits))

# Measure each input qubit. 
circuit.measure(range(numQubits))

# Display the quantum circuit.
circuit.display_circuit()

# Choose the number of shots for the circuit. Run the circuit 'shots' times and obtain the results as a list. Set hist=True to create a histogram of the results. You should see that the results give |101> = |5> with 100% probability.
shots = 1024
results = circuit.run(shots, hist=True)

# If you would like to see all the gates performed by the algorithm, apply the IQFT directly through Algorithms.py. Note: the function IQFT was imported from Algorithms.py above, so we can use this function directly now.
# Computatinally, this method and the method above are exactly the same.
numQubits = 3
circuit = QPU.Circuit(numQubits)

circuit.H(range(numQubits))
circuit.P(0, theta=5*np.pi/4)
circuit.P(1, theta=5*np.pi/2)
circuit.P(2, theta=5*np.pi)

IQFT(circuit, range(numQubits))

for qubit in range(numQubits):
    circuit.measure(qubit)

circuit.display_circuit()

shots = 1024
results = circuit.run(shots, hist=True)