# This tutorial shows how to use Grover's algorithm with the QPU simulator.

################################################################################

# Import the QPU simulator which creates the quantum circuit and contains all gate operations.
# Import Grover from Algorithms if you want to access the algorithm directly, rather than through QPUsimulator. Computationally, using the algorithm through Simulator.py or Algorithms.py is the same. When displaying the circuit, Simulator.py simplifies the diagram to show a general "Grover" block over all the qubits involved to represent the algorithm; Algorithms.py shows all the individual gates completed in the algorithm.
import Simulator as QPU
from Algorithms import Grover

# Create a quantum circuit with a chosen number of qubits. The same number of classical bits will be created for qubit measurement output. The output of each qubit is stored in its similarly indexed classical bit.
numQubits = 3
circuit = QPU.Circuit(numQubits)

# Define the oracle as a python function. Use the oracle to mark states for amplification.
def tutorialOracle(circuit):

    # For this 3 qubit oracle, apply controlled-Z gates to the last qubit controlled by the other 2 qubits. This marks the |101> and |110> states.
    circuit.CZ([0], 2)
    circuit.CZ([1], 2)

    return

# Apply the algorithm using the tutorial oracle, which marks the |101> and |110> states for amplification.
circuit.Grover(oracle=tutorialOracle)

# Measure the qubits.
circuit.measure([0, 1, 2])

# Display the quantum circuit.
circuit.display_circuit()

# Choose the number of shots for the circuit. Run the circuit 'shots' times and obtain the results. Set hist=True to create a histogram of the results.
shots = 1024
results = circuit.run(shots, hist=True)


# There is an example oracle pre-programmed in Algorithms.py that marks the state |101> and |110> in a 3 qubit circuit.
numQubits = 3
circuit = QPU.Circuit(numQubits)

# Apply the algorithm using the example 3 qubit oracle, which marks the |101> and |110> states for amplification.
circuit.Grover(oracle='example')

circuit.measure([0, 1, 2])

circuit.display_circuit()

shots = 1024
results = circuit.run(shots, hist=True)


# If you would like to see all the gates performed by the algorithm, apply the Grover algorithm directly through Algorithms.py. Note: the function Grover was imported from Algorithms.py above, so we can use this function directly now.
# Computatinally, this method and the method above are exactly the same.
numQubits = 3
circuit = QPU.Circuit(numQubits)

Grover(circuit, oracle='example')

circuit.measure([0, 1, 2])

circuit.display_circuit()

shots = 1024
results = circuit.run(shots, hist=True)