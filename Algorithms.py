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

# Quantum Fourier Transform (QFT): this algorithm converts qubits in the computational basis into the Fourier basis. This is commonly used as a sub-step within other algorithms.
#
# You may provide the number of qubits to perform the QFT on using numQubits. Note that the qubits involved must be sequential and ordered from least significant (lowest index) to most significant (highest index). To perform QFT on all qubits within the circuit, you may leave this argument as the default, and the function will get the number of qubits in the circuit.
def QFT(circuit, algQubits=None):

    # If no qubits are provided in algQubits, use all the qubits in the circuit. Set numQubits as the length of the qubits involved.
    if algQubits == None:
        algQubits = list(range(circuit.numQubits))
    numQubits = len(algQubits)

    # Start with the most significant qubit. Apply an H gate to convert the qubit into the Fourier basis. Turn the qubit the appropriate angle using controlled-P gates. Apply controlled-P gates with the current qubit as the target and all other qubits with lower index as control qubits. Theta starts at pi/2^n, where 'n' is the number of controlled-P gates to be applied to the target, and is doubled with each successive controlled-P gate applied to the target, ending at pi/2.
    for qubit in range(numQubits-1, algQubits[0]-1, -1):
        circuit.H(qubit)
        for control in range(qubit):
            circuit.CP([control], qubit, theta=np.pi/2**(qubit-control))

    # Swap the qubit order. Conversion between the computational and Fourier bases reverses the qubit order of significance.

    # First, create a list of the qubit indices reversed.
    reversedQubits = list(reversed(algQubits))

    # Loop over the qubits to match up each qubit with its SWAP partner.
    for Qidx, qubit in enumerate(algQubits):

        # If current qubit's index is greater than the qubit's index from the reversed list, all qubits have been swapped. If the current qubit's index equals the qubit's index in the reversed list, this qubit will not have its order changed. For both cases, all SWAPs have been completed, so end the loop.
        if qubit >= reversedQubits[Qidx]:
            break

        # Swap the current qubit with the last unswapped qubit (which will have the same list index as the algQubits list).
        circuit.SWAP(qubit, reversedQubits[Qidx])
    
    return

# Inverse Quantum Fourier Transform (IQFT): this algorithm converts qubits in the Fourier basis into the computational basis. This is commonly used as a sub-step within other algorithms.
#
# You may provide the qubit indices to perform the IQFT on using algQubits. Note that the qubits involved must be sequential and ordered from least significant (lowest index) to most significant (highest index). To perform IQFT on all qubits within the circuit, you may leave this argument as the default, and the function will use all the qubits in the circuit.
def IQFT(circuit, algQubits=None):

    # If no qubits are provided in algQubits, use all the qubits in the circuit. Set numQubits as the length of the qubits involved.
    if algQubits == None:
        algQubits = list(range(circuit.numQubits))
    numQubits = len(algQubits)

    # Swap the qubit order. Conversion between the computational and Fourier bases reverses the qubit order of significance.

    # First, create a list of the qubit indices reversed.
    reversedQubits = list(reversed(algQubits))

    # Loop over the qubits to match up each qubit with its SWAP partner.
    for Qidx, qubit in enumerate(algQubits):

        # If current qubit's index is greater than the qubit's index from the reversed list, all qubits have been swapped. If the current qubit's index equals the qubit's index in the reversed list, this qubit will not have its order changed. For both cases, all SWAPs have been completed, so end the loop.
        if qubit >= reversedQubits[Qidx]:
            break

        # Swap the current qubit with the last unswapped qubit (which will have the same list index as the algQubits list).
        circuit.SWAP(qubit, reversedQubits[Qidx])

    # Start with the least significant qubit. Turn the qubit the appropriate angle using controlled-P gates. Apply controlled-P gates with the current qubit as the target and all other qubits with lower index as control qubits. Theta starts at -pi/2, and is halved with each successive controlled-P gate applied to the target, ending at -pi/2^n, where 'n' is the number of controlled-P gates to be applied to the target. Apply an H gate to convert the qubit into the computational basis.
    for qubit in algQubits:
        for control in range(qubit-1, algQubits[0]-1, -1):
            circuit.CP([control], qubit, theta=-np.pi/2**(qubit-control))
        circuit.H(qubit)
    
    return

# Quantum phase estimation (QPE): this algorithm estimates the angle theta within the eigenvalue problem U|psi> = e^(2pi*i*theta)|psi>. The more qubits are included in the algorithm, the higher the precision of the algorithm (at the expense of higher computational cost). This is commonly used as a sub-step within other algorithms.
#
# The angle lambd = 2pi*theta must be passed as an argument.
#
# You may provide the number of qubits to perform the QPE on using numPrecisionQubits. To perform QPE on all qubits within the circuit (minus the final qubit which represents |psi>), you may leave this argument as the default, and the function will get the number of qubits in the circuit.
def QPE(circuit, lambd, algQubits=None):

    # If no qubits are provided in algQubits, use all the qubits in the circuit.
    if algQubits == None:
        algQubits = list(range(circuit.numQubits))

    # Use the first n-1 qubits in algQubits as the precision qubits for the algorithm. The nth qubit will be the qubit that stores psi.
    precisionQubits = algQubits[:-1]
    psiQubit = algQubits[-1]

    # Initialize the |psi> qubit in the |1> state.
    circuit.X(psiQubit)

    # Convert the precision qubits into the Fourier basis with H gates.
    circuit.H(precisionQubits)

    # Turn the |psi> qubit using each precision qubit as the control in controlled-P gates. The angle for each turn is lambd = 2pi*theta. For each precision qubit, apply 2^n controlled-P gates, where n is the index of the precision qubit.
    for control in precisionQubits:
        for repeat in range(2**control):
                circuit.CP([control], psiQubit, lambd)

    # Apply the inverse QFT to the precision qubits to convert them back into the computational basis.
    circuit.IQFT(algQubits=precisionQubits)

    return

def expGroverOracle(circuit):

    # For this 3 qubit oracle, apply controlled-Z gates to the last qubit controlled by the other 2 qubits. This marks the |101> and |110> states.
    circuit.CZ([0], 2)
    circuit.CZ([1], 2)

    return

def Grover(circuit, numQubits=0, oracle='example'):

    # If no number of qubits to apply the QFT are passed to the function, get the number of qubits in the circuit to use all of them in the algorithm.
    if numQubits == 0:
        numQubits = circuit.numQubits
    else:
        pass

    # Apply Hadamard gates to each qubit to create a superposition of all possible states.
    circuit.H(range(numQubits))

    circuit.barrier()

    # Apply the oracle which marks the selected states.

    if oracle == 'example':
        expGroverOracle(circuit)
    else:
        oracle(circuit)
    
    circuit.barrier()

    # Amplify the marked states. Apply a Hadamard and Pauli X gate to each qubit. Then control a Z gate on the last qubit with all remaining qubits acting as controls for the one Z gate. Apply a Pauli X and Hadamard gate to each qubit.
    circuit.H(range(numQubits))
    circuit.X(range(numQubits))
    circuit.CZ(range(numQubits-1), numQubits-1)
    circuit.X(range(numQubits))
    circuit.H(range(numQubits))

    circuit.barrier()

    return

def main():
    return

if __name__ == "__main__":
    main()