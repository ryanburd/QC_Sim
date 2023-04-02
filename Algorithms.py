# This file contains all the algorithms preprogrammed into the QPUsimulator. See the individual tutorials for an example on how to use each algorithm in your circuits.

import Simulator
import numpy as np

# An example constant oracle for use in the Deutsch-Jozsa algorithm. The constant output of 0 or 1 can be specified with the 'output' argument, which has 0 as default. All inputs should be measured in the |0> state.
def constantOracle(circuit, numQubits, output):

    # If the desired constant output from the oracle is 1, apply an X gate to the output qubit.
    if output == 1:
        circuit.X(numQubits-1)

    return

# An example balanced oracle for use in the Deutsch-Jozsa algorithm. The set of inputs yielding 0 vs 1 can be changed by flipping the states of any of the input qubits using the 'inputFlips' argument. All inputs should be measured in the |1> state.
def balancedOracle(circuit, numQubits, inputFlips=[]):

    # Apply X gates to select input qubits to change which set of inputs result in an output=0 (and vice versa for output=1)
    circuit.X(inputFlips)

    circuit.barrier()

    # Apply controlled-X gates for each input qubit acting as a single control for the target (output qubit)
    for control in range(numQubits-1):
        circuit.CX([control], numQubits-1)

    circuit.barrier()

    # Reapply the X gates on the select input qubits from above
    circuit.X(inputFlips)

    return

# Deutsch-Jozsa algorithm: this algorithm takes a provided oracle in the form of f({x1, x2, ..., xn}) = {0, 1} where each input x as 0 or 1. The oracle must be constant, i.e. all possible inputs yield the same output of 0 or 1, or balanced, i.e. half of all possible inputs yield 0 and the other half all yield 1. This algorithm can identify if the oracle is constant or balanced with 100% accuracy in just a single shot of the circuit. For constant oracles, all input qubits will be measured in the |0> state at the end of the algorithm, vs all |1>'s for balanced oracles.
# 
# While not an interesting algorithm from a practical perspective, this algorithm does show that certain problems can be completed much more efficiently on a quanutm computer than a classical computer. On a classical computer, it would take 2^(n-1)+1 shots to identify the oracle as constant or balanced with 100% accuracy, with n being the number of inputs (x's) to the oracle.
#
# Example constant and balanced oracles are provided above. These can be selected by setting the 'oracle' argument to the string 'constant' or 'balanced' respecively. The user may also create their own oracle and pass it (as a python function) to the 'oracle' argument.
#
# For the example constant oracle, the constant output of 0 or 1 can be specified with the 'constantOracleOutput' argument, which has 0 as default.
#
# For the example balanced oracle, the set of inputs yielding 0 vs 1 can be changed by flipping the states of any of the input qubits using the 'balancedInputFlips' argument.
def DeutschJozsa(circuit, oracle, constantOracleOutput=0, balancedInputFlips=[]):

    numQubits = circuit.numQubits

    circuit.barrier()

    # Apply an X gate to the output qubit. Then apply H gates to all qubits. This will initialize all input qubits in the |+> state and the output qubit in the |-> state.
    circuit.X(numQubits-1)
    for Qidx in range(numQubits):
        circuit.H(Qidx)

    # Apply the oracle specified when calling the algorithm function. Pass the appropriate arguments.
    if oracle == 'balanced':
        balancedOracle(circuit, numQubits, balancedInputFlips)
    elif oracle == 'constant':
        constantOracle(circuit, numQubits, constantOracleOutput)
    else:
        oracle(circuit)

    # Apply H gates to all input qubits to put the back in the computational basis.
    for Qidx in range(numQubits-1):
        circuit.H(Qidx)

    circuit.barrier()

    return

# Quantum Fourier Transform (QFT)
def QFT(circuit):

    numQubits = circuit.numQubits

    for qubit in range(numQubits-1, -1, -1):
        circuit.H(qubit)
        for control in range(qubit):
            circuit.CP([control], qubit, theta=np.pi/2**(qubit-control))

    for qubit in range(numQubits):
        if qubit == numQubits-qubit-1:
            break
        if qubit >= 0.5*numQubits:
            break
        circuit.SWAP(qubit, numQubits-qubit-1)
    
    return

# Inverse Quantum Fourier Transform (IQFT)
def IQFT(circuit):

    numQubits = circuit.numQubits

    for qubit in range(numQubits):
        if qubit == numQubits-qubit-1:
            break
        if qubit >= 0.5*numQubits:
            break
        circuit.SWAP(qubit, numQubits-qubit-1)

    for qubit in range(numQubits):
        for control in range(qubit-1, -1, -1):
            circuit.CP([control], qubit, theta=-np.pi/2**(qubit-control))
        circuit.H(qubit)
    
    return

def main():
    return

if __name__ == "__main__":
    main()