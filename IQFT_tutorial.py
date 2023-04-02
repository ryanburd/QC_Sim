# This tutorial shows how to use the inverse Quantum Fourier Transform (IQFT) algorithm with the QPU simulator.

################################################################################

# Import the QPU simulator which creates the quantum circuit and contains all gate operations. Import IQFT from Algorithms if you want to access the algorithm directly, rather than through QPUsimulator. Note this is only a syntax difference when calling the algorithm; the circuit operation will be exactly the same for either method. Import numpy for using np.pi.
import Simulator as QPU
from Algorithms import IQFT
import numpy as np

# Create a quantum circuit with a chosen number of qubits. The same number of classica bits will be created for qubit measurement output. The output of each qubit is stored in its similarly indexed classical bit.
numQubits = 3
circuit = QPU.Circuit(numQubits)

# Initialize the circuit to a state of your choosing in the Fourier basis. For example, the state |5> = |101> in the Fourier basis is prepared by applying an H to all qubits and then the phases 5pi/4, 5pi/2, and 5pi, respectively.
circuit.H([0, 1, 2])
circuit.P(0, theta=5*np.pi/4)
circuit.P(1, theta=5*np.pi/2)
circuit.P(2, theta=5*np.pi)

# To apply the algorithm to your circuit, you can apply it just as you would a qubit gate.
circuit.IQFT()

# Measure each input qubit. 
circuit.measure(range(numQubits))

# Display the quantum circuit.
circuit.display_circuit()

# Choose the number of shots for the circuit. Run the circuit 'shots' times and obtain the results as a list. Set hist=True to create a histogram of the results. You should see that the results give |101> = |5> with 100% probability.
shots = 1024
results = circuit.run(shots, hist=True)

# You can apply the algorithm using any oracle by accessing Algorithms.py directly if it is imported. In this case, pass the circuit as the only argument to the IQFT function. For circuit operations, this is identical to applying the algorithm using the gate syntax used above. Only the syntax for calling the algorithm changes. Use whichever syntax is more intuitive to you.
numQubits = 3
circuit = QPU.Circuit(numQubits)

circuit.H([0, 1, 2])
circuit.P(0, theta=5*np.pi/4)
circuit.P(1, theta=5*np.pi/2)
circuit.P(2, theta=5*np.pi)

IQFT(circuit)

for qubit in range(numQubits):
    circuit.measure(qubit)

circuit.display_circuit()

shots = 1024
results = circuit.run(shots, hist=True)