# This tutorial shows how to use the Deutsch-Jozsa algorithm with the QPU simulator.

################################################################################

# Import the QPU simulator which creates the quantum circuit and contains all gate operations.
# Import DeutschJozsa from Algorithms if you want to access the algorithm directly, rather than through QPUsimulator. Computationally, using the algorithm through Simulator.py or Algorithms.py is the same. When displaying the circuit, Simulator.py simplifies the diagram to show a general "DJ" block over all the qubits involved to represent the algorithm; Algorithms.py shows all the individual gates completed in the algorithm.
import Simulator as QPU
from Algorithms import DeutschJozsa

# Create a quantum circuit with a chosen number of input qubits for the oracle. An additional qubit will be added as the single output qubit. The same number of classical bits will be created for qubit measurement output. The output of each qubit is stored in its similarly indexed classical bit.
numInputQubits = 3
circuit = QPU.Circuit(numInputQubits+1)

# To use your own oracle, define the oracle as a python function taking only the circuit as an argument. Be sure your oracle is either constant or balanced! For this tutorial, a balanced oracle will be used. (See Algorithms.py for an example of a constant oracle.)
def tutorialOracle(circuit, algQubits=None):

    # If no qubits are provided in algQubits, use all the qubits in the circuit. Set numQubits as the length of the qubits involved.
    if algQubits == None:
        algQubits = list(range(circuit.numQubits))
    numQubits = len(algQubits)

    # Select input qubits that will be flipped before applying controlled-X gates.
    inputFlips = [0, 1]
    circuit.X(inputFlips)

    # Apply controlled-X gates from each input qubit to the output qubit as the target.
    for control in algQubits[:-1]:
        circuit.CX([control], algQubits[-1])

    # Undo the input qubit flips.
    circuit.X(inputFlips)

    return

# To apply the algorithm to your circuit, you can apply it just as you would a qubit gate. Pass the function containing your oracle as the 'oracle' argument.
circuit.DeutschJozsa(oracle=tutorialOracle, oracleType='balanced')

# Measure each input qubit. For constant oracles, all inputs should measure |0> vs all |1>'s for balanced oracles.
for qubit in range(numInputQubits):
    circuit.measure(qubit)

# Display the quantum circuit.
circuit.display_circuit()

# Choose the number of shots for the circuit. Run the circuit 'shots' times and obtain the results as a list. Set hist=True to create a histogram of the results.
# The Deutsch-Jozsa only needs a single shot to solve this problem of whether the oracle is constant or balanced. Here, many shots will be run to show that the same result is achieved every time.
shots = 1024
results = circuit.run(shots, hist=True)

# Algorithms.py also inlcudes example constant and balanced oracles if you do not want to define your own.
#
# To use the example constant oracle, pass the string 'constant' as the 'oracle' argument. Choose the constant for the output qubit with 'constant OracleOutput'. The default is 0.
numInputQubits = 3
circuit = QPU.Circuit(numInputQubits+1)

circuit.DeutschJozsa(oracle='constant', constantOracleOutput=1)

for qubit in range(numInputQubits):
    circuit.measure(qubit)

circuit.display_circuit()

shots = 1024
results = circuit.run(shots, hist=True)

# The example balanced oracle is actually the exact same oracle defined above as 'tutorialOracle'. To use this oracle without defining it yourself, pass the string 'balanced' as the 'oracle' argument. Choose which input qubits to flip in the oracle by passing their indicies as a list in 'balancedInputFlips'. The default is no input qubits.
numInputQubits = 3
circuit = QPU.Circuit(numInputQubits+1)

circuit.DeutschJozsa(oracle='balanced', balancedInputFlips=[0, 1])

for qubit in range(numInputQubits):
    circuit.measure(qubit)

circuit.display_circuit()

shots = 1024
results = circuit.run(shots, hist=True)

# If you would like to see all the gates performed by the algorithm, apply the DeutschJozsa directly through Algorithms.py. Note: the function DeutschJozsa was imported from Algorithms.py above, so we can use this function directly now.
# Computatinally, this method and the method above are exactly the same.
numInputQubits = 3
circuit = QPU.Circuit(numInputQubits+1)

DeutschJozsa(circuit, oracle=tutorialOracle)

for qubit in range(numInputQubits):
    circuit.measure(qubit)

circuit.display_circuit()

shots = 1024
results = circuit.run(shots, hist=True)