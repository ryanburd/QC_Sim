## This file allows the user to create a quantum circuit, apply common qubit gates, and collect the measurement result of the circuit with many shots. A diagram of the circuit can be created.

import Algorithms
import numpy as np
import matplotlib.pyplot as plt
from itertools import product as CartesianProduct

# Define common matrices used for gate operations. Since the rotation matrices need to receive angles, these matrices are packaged into a function instead of a dictionary, though this function essential acts as a dictionary.
def gateMatrix(gateType, angles=[0, 0, 0]):
    [theta, phi, lambd] = angles
    if gateType == 'I':
        return np.eye(2)
    if gateType == 'X':
        return np.array([[0, 1],
                         [1, 0]])
    if gateType == 'Y':
        return np.array([[ 0, -1j],
                         [1j,   0]])
    if gateType == 'Z':
        return np.array([[1,  0],
                         [0, -1]])
    if gateType == 'H':
        return 1/np.sqrt(2)*np.array([[1,  1],
                                      [1, -1]])
    if gateType == 'S':
        return np.array([[1,  0],
                         [0, 1j]])
    if gateType == 'T':
        return np.array([[1,                  0],
                         [0, np.exp(1j*np.pi/4)]])
    if gateType == 'P':
        return np.array([[1,                  0],
                         [0, np.exp(1j*theta)]])
    if gateType == 'RX':
        return np.array([[    np.cos(theta/2), -1j*np.sin(theta/2)],
                         [-1j*np.sin(theta/2),     np.cos(theta/2)]])
    if gateType == 'RY':
        return np.array([[np.cos(theta/2), -1*np.sin(theta/2)],
                         [np.sin(theta/2),    np.cos(theta/2)]])
    if gateType == 'RZ':
        return np.array([[np.exp(-1j*theta/2),                  0],
                         [                  0, np.exp(1j*theta/2)]])
    if gateType == 'U':
        return np.array([[np.cos(theta/2), -np.exp(1j*lambd)*np.sin(theta/2)],[np.exp(1j*phi)*np.sin(theta/2), np.exp(1j*(phi+lambd))*np.cos(theta/2)]])

    # Projection matrices into the computational basis. Used for constructing controlled-U gates and measurements.
    if gateType == 'P0':
        return np.array([[1, 0],
                         [0, 0]])
    if gateType == 'P1':
        return np.array([[0, 0],
                         [0, 1]])

# Creates a qubit object, which stores all gates applied to the qubit and all connections the qubit is a part of (e.g. as a control for another target qubit's gate).
class Qubit:

    def __init__(self):
        self.gates = []
        self.gatePos = []
        self.gateAngles = []
        self.connections = []
        self.connectTo = []
        self.connectPos = []
        self.earliestPos = 1

# Creates a classical bit object, which stores the bit's state (0 or 1) and all connections the bit is a part of (e.g. as storage for the result of measurement on a qubit).
class Cbit:

    def __init__(self, state):
        self.state = state
        self.connections = []
        self.connectTo = []
        self.connectPos = []
        self.earliestPos = 1

# Creates an instance of a quantum circuit with a provided number of quantum bits and classical bits and allows the user
# to apply qubit gates to the circuit
class Circuit:

    # Create the provided number of qubits and classical bits upon instance initialization
    def __init__(self, numQubits):

        # Create a list of qubits. Each instance of the class Qubit will store the gates applied to the qubit. This is useful for creating a diagram of the circuit.
        self.numQubits = numQubits
        self.qubits = [Qubit() for qubit in range(numQubits)]
        
        # Form the state of all the qubits in the circuit. Assume all qubits are initialized in the |0> state, [1, 0]
        self.state = np.array([1])
        for qubit in self.qubits:
            self.state = np.tensordot([1, 0], self.state, axes=0).reshape(len(self.state)*2, 1)
        
        # Create a list of classical bits, each initialized in the 0 state.
        self.numCbits = numQubits
        self.cbits = [Cbit(0) for cbit in range(numQubits)]

    ## Gate functions below add their respective gates to the ongoing list of gates defined for each qubit. When running the circuit with run(), the gate lists are collected and applied to the circuit's initial state vector.

    ## SINGLE QUBIT GATES ##

    # Pauli-X gate
    def X(self, targets):

        # Append the gate onto the running list of gates for the target qubits. If only one target is passed as an int, place it in a list to avoid an error when iterating over the for loop. Use the qubit's earliest position for the gate position, then increment the earliest position for potential future gates.
        if type(targets) == int:
            targets = [targets]
        for target in targets:
            self.qubits[target].gates.append('X')
            angles = [None, None, None]
            self.qubits[target].gateAngles.append(angles)
            self.qubits[target].gatePos.append(self.qubits[target].earliestPos)
            self.qubits[target].earliestPos += 1

        return self

    # Pauli-Y gate
    def Y(self, targets):

        # Append the gate onto the running list of gates for the target qubits. If only one target is passed as an int, place it in a list to avoid an error when iterating over the for loop. Use the qubit's earliest position for the gate position, then increment the earliest position for potential future gates.
        if type(targets) == int:
            targets = [targets]
        for target in targets:
            self.qubits[target].gates.append('Y')
            angles = [None, None, None]
            self.qubits[target].gateAngles.append(angles)
            self.qubits[target].gatePos.append(self.qubits[target].earliestPos)
            self.qubits[target].earliestPos += 1

        return self

    # Pauli-Z gate
    def Z(self, targets):

        # Append the gate onto the running list of gates for the target qubits. If only one target is passed as an int, place it in a list to avoid an error when iterating over the for loop. Use the qubit's earliest position for the gate position, then increment the earliest position for potential future gates.
        if type(targets) == int:
            targets = [targets]
        for target in targets:
            self.qubits[target].gates.append('Z')
            angles = [None, None, None]
            self.qubits[target].gateAngles.append(angles)
            self.qubits[target].gatePos.append(self.qubits[target].earliestPos)
            self.qubits[target].earliestPos += 1

        return self
    
    # Hadamard gate
    def H(self, targets):

        # Append the gate onto the running list of gates for the target qubits. If only one target is passed as an int, place it in a list to avoid an error when iterating over the for loop. Use the qubit's earliest position for the gate position, then increment the earliest position for potential future gates.
        if type(targets) == int:
            targets = [targets]
        for target in targets:
            self.qubits[target].gates.append('H')
            angles = [None, None, None]
            self.qubits[target].gateAngles.append(angles)
            self.qubits[target].gatePos.append(self.qubits[target].earliestPos)
            self.qubits[target].earliestPos += 1

        return self
    
    # Phase gate
    def S(self, targets):

        # Append the gate onto the running list of gates for the target qubits. If only one target is passed as an int, place it in a list to avoid an error when iterating over the for loop. Use the qubit's earliest position for the gate position, then increment the earliest position for potential future gates.
        if type(targets) == int:
            targets = [targets]
        for target in targets:
            self.qubits[target].gates.append('S')
            angles = [None, None, None]
            self.qubits[target].gateAngles.append(angles)
            self.qubits[target].gatePos.append(self.qubits[target].earliestPos)
            self.qubits[target].earliestPos += 1

        return self
    
    # pi/8 gate
    def T(self, targets):

        # Append the gate onto the running list of gates for the target qubits. If only one target is passed as an int, place it in a list to avoid an error when iterating over the for loop. Use the qubit's earliest position for the gate position, then increment the earliest position for potential future gates.
        if type(targets) == int:
            targets = [targets]
        for target in targets:
            self.qubits[target].gates.append('T')
            angles = [None, None, None]
            self.qubits[target].gateAngles.append(angles)
            self.qubits[target].gatePos.append(self.qubits[target].earliestPos)
            self.qubits[target].earliestPos += 1

        return self
    
    # phase gate
    def P(self, targets, theta):

        # Append the gate onto the running list of gates for the target qubits. If only one target is passed as an int, place it in a list to avoid an error when iterating over the for loop. Use the qubit's earliest position for the gate position, then increment the earliest position for potential future gates.
        if type(targets) == int:
            targets = [targets]
        for target in targets:
            self.qubits[target].gates.append('P')
            angles = [theta, None, None]
            self.qubits[target].gateAngles.append(angles)
            self.qubits[target].gatePos.append(self.qubits[target].earliestPos)
            self.qubits[target].earliestPos += 1

        return self

    # R_X gate
    def RX(self, targets, theta):

        # Append the gate onto the running list of gates for the target qubits. If only one target is passed as an int, place it in a list to avoid an error when iterating over the for loop. gateAngles stores lists of 3 angles, so place theta and 2 Nones in a list before appending it. Use the qubit's earliest position for the gate position, then increment the earliest position for potential future gates.
        if type(targets) == int:
            targets = [targets]
        for target in targets:
            self.qubits[target].gates.append('RX')
            angles = [theta, None, None]
            self.qubits[target].gateAngles.append(angles)
            self.qubits[target].gatePos.append(self.qubits[target].earliestPos)
            self.qubits[target].earliestPos += 1

        return self
    
    # R_Y gate
    def RY(self, targets, theta):

        # Append the gate onto the running list of gates for the target qubits. If only one target is passed as an int, place it in a list to avoid an error when iterating over the for loop. gateAngles stores lists of 3 angles, so place theta and 2 Nones in a list before appending it. Use the qubit's earliest position for the gate position, then increment the earliest position for potential future gates.
        if type(targets) == int:
            targets = [targets]
        for target in targets:
            self.qubits[target].gates.append('RY')
            angles = [theta, None, None]
            self.qubits[target].gateAngles.append(angles)
            self.qubits[target].gatePos.append(self.qubits[target].earliestPos)
            self.qubits[target].earliestPos += 1

    # R_Z gate
    def RZ(self, targets, theta):

        # Append the gate onto the running list of gates for the target qubits. If only one target is passed as an int, place it in a list to avoid an error when iterating over the for loop. gateAngles stores lists of 3 angles, so place theta and 2 Nones in a list before appending it. Use the qubit's earliest position for the gate position, then increment the earliest position for potential future gates.
        if type(targets) == int:
            targets = [targets]
        for target in targets:
            self.qubits[target].gates.append('RZ')
            angles = [theta, None, None]
            self.qubits[target].gateAngles.append(angles)
            self.qubits[target].gatePos.append(self.qubits[target].earliestPos)
            self.qubits[target].earliestPos += 1

    # U gate
    def U(self, targets, theta, phi, lambd):

        # Append the gate onto the running list of gates for the target qubits. If only one target is passed as an int, place it in a list to avoid an error when iterating over the for loop. gateAngles stores lists of 3 angles, so place theta and 2 Nones in a list before appending it. Use the qubit's earliest position for the gate position, then increment the earliest position for potential future gates.
        if type(targets) == int:
            targets = [targets]
        for target in targets:
            self.qubits[target].gates.append('U')
            angles = [theta, phi, lambd]
            self.qubits[target].gateAngles.append(angles)
            self.qubits[target].gatePos.append(self.qubits[target].earliestPos)
            self.qubits[target].earliestPos += 1

    ## TWO QUBIT GATES ##

    # Controlled-X gate
    def CX(self, controls, target):

        # Append the gate onto the running list of gates for the target qubit. Append the connection type onto the running list of connections for the control qubit and which qubit it is controlling (the target).
        self.qubits[target].gates.append('X')
        angles = [None, None, None]
        self.qubits[target].gateAngles.append(angles)
        for control in controls:
            self.qubits[control].connections.append('C')
            self.qubits[control].connectTo.append(target)

        # Get the earliest possible gate position within the circuit for each qubit between the target and control (inclusive). The max of this list will be used for the gate position for both the target and control. Then increase the earliest position for all qubits.
        earliestPositions = [self.qubits[idx].earliestPos for idx in range(min(min(controls), target), max(max(controls), target)+1)]
        position = max(earliestPositions)
        self.qubits[target].gatePos.append(position)
        for control in controls:
            self.qubits[control].connectPos.append(position)
        for Qidx in range(min(min(controls), target), max(max(controls), target)+1):
            self.qubits[Qidx].earliestPos = position + 1

        return self
    
    # Controlled-Y gate
    def CY(self, controls, target):

        # Append the gate onto the running list of gates for the target qubit. Append the connection type onto the running list of connections for the control qubit and which qubit it is controlling (the target).
        self.qubits[target].gates.append('Y')
        angles = [None, None, None]
        self.qubits[target].gateAngles.append(angles)
        for control in controls:
            self.qubits[control].connections.append('C')
            self.qubits[control].connectTo.append(target)

        # Get the earliest possible gate position within the circuit for each qubit between the target and control (inclusive). The max of this list will be used for the gate position for both the target and control. Then increase the earliest position for all qubits.
        earliestPositions = [self.qubits[idx].earliestPos for idx in range(min(min(controls), target), max(max(controls), target)+1)]
        position = max(earliestPositions)
        self.qubits[target].gatePos.append(position)
        for control in controls:
            self.qubits[control].connectPos.append(position)
        for Qidx in range(min(min(controls), target), max(max(controls), target)+1):
            self.qubits[Qidx].earliestPos = position + 1

        return self
    
    # Controlled-Z gate
    def CZ(self, controls, target):

        # Append the gate onto the running list of gates for the target qubit. Append the connection type onto the running list of connections for the control qubit and which qubit it is controlling (the target).
        self.qubits[target].gates.append('Z')
        angles = [None, None, None]
        self.qubits[target].gateAngles.append(angles)
        for control in controls:
            self.qubits[control].connections.append('C')
            self.qubits[control].connectTo.append(target)

        # Get the earliest possible gate position within the circuit for each qubit between the target and control (inclusive). The max of this list will be used for the gate position for both the target and control. Then increase the earliest position for all qubits.
        earliestPositions = [self.qubits[idx].earliestPos for idx in range(min(min(controls), target), max(max(controls), target)+1)]
        position = max(earliestPositions)
        self.qubits[target].gatePos.append(position)
        for control in controls:
            self.qubits[control].connectPos.append(position)
        for Qidx in range(min(min(controls), target), max(max(controls), target)+1):
            self.qubits[Qidx].earliestPos = position + 1

        return self
    
    # Controlled-P gate
    def CP(self, controls, target, theta):

        # Append the gate onto the running list of gates for the target qubit. Append the connection type onto the running list of connections for the control qubit and which qubit it is controlling (the target).
        self.qubits[target].gates.append('P')
        angles = [theta, None, None]
        self.qubits[target].gateAngles.append(angles)
        for control in controls:
            self.qubits[control].connections.append('C')
            self.qubits[control].connectTo.append(target)

        # Get the earliest possible gate position within the circuit for each qubit between the target and control (inclusive). The max of this list will be used for the gate position for both the target and control. Then increase the earliest position for all qubits.
        earliestPositions = [self.qubits[idx].earliestPos for idx in range(min(min(controls), target), max(max(controls), target)+1)]
        position = max(earliestPositions)
        self.qubits[target].gatePos.append(position)
        for control in controls:
            self.qubits[control].connectPos.append(position)
        for Qidx in range(min(min(controls), target), max(max(controls), target)+1):
            self.qubits[Qidx].earliestPos = position + 1

        return self

    # Controlled-RX gate
    def CRX(self, controls, target, theta):

        # Append the gate onto the running list of gates for the target qubit. Append the connection type onto the running list of connections for the control qubit and which qubit it is controlling (the target).
        self.qubits[target].gates.append('RX')
        angles = [theta, None, None]
        self.qubits[target].gateAngles.append(angles)
        for control in controls:
            self.qubits[control].connections.append('C')
            self.qubits[control].connectTo.append(target)

        # Get the earliest possible gate position within the circuit for each qubit between the target and control (inclusive). The max of this list will be used for the gate position for both the target and control. Then increase the earliest position for all qubits.
        earliestPositions = [self.qubits[idx].earliestPos for idx in range(min(min(controls), target), max(max(controls), target)+1)]
        position = max(earliestPositions)
        self.qubits[target].gatePos.append(position)
        for control in controls:
            self.qubits[control].connectPos.append(position)
        for Qidx in range(min(min(controls), target), max(max(controls), target)+1):
            self.qubits[Qidx].earliestPos = position + 1

        return self

        # Controlled-RY gate
    def CRY(self, controls, target, theta):

        # Append the gate onto the running list of gates for the target qubit. Append the connection type onto the running list of connections for the control qubit and which qubit it is controlling (the target).
        self.qubits[target].gates.append('RY')
        angles = [theta, None, None]
        self.qubits[target].gateAngles.append(angles)
        for control in controls:
            self.qubits[control].connections.append('C')
            self.qubits[control].connectTo.append(target)

        # Get the earliest possible gate position within the circuit for each qubit between the target and control (inclusive). The max of this list will be used for the gate position for both the target and control. Then increase the earliest position for all qubits.
        earliestPositions = [self.qubits[idx].earliestPos for idx in range(min(min(controls), target), max(max(controls), target)+1)]
        position = max(earliestPositions)
        self.qubits[target].gatePos.append(position)
        for control in controls:
            self.qubits[control].connectPos.append(position)
        for Qidx in range(min(min(controls), target), max(max(controls), target)+1):
            self.qubits[Qidx].earliestPos = position + 1

        return self

    # Controlled-RZ gate
    def CRZ(self, controls, target, theta):

        # Append the gate onto the running list of gates for the target qubit. Append the connection type onto the running list of connections for the control qubit and which qubit it is controlling (the target).
        self.qubits[target].gates.append('RZ')
        angles = [theta, None, None]
        self.qubits[target].gateAngles.append(angles)
        for control in controls:
            self.qubits[control].connections.append('C')
            self.qubits[control].connectTo.append(target)

        # Get the earliest possible gate position within the circuit for each qubit between the target and control (inclusive). The max of this list will be used for the gate position for both the target and control. Then increase the earliest position for all qubits.
        earliestPositions = [self.qubits[idx].earliestPos for idx in range(min(min(controls), target), max(max(controls), target)+1)]
        position = max(earliestPositions)
        self.qubits[target].gatePos.append(position)
        for control in controls:
            self.qubits[control].connectPos.append(position)
        for Qidx in range(min(min(controls), target), max(max(controls), target)+1):
            self.qubits[Qidx].earliestPos = position + 1

        return self

        # Controlled-U gate
    def CU(self, controls, target, theta, phi, lambd):

        # Append the gate onto the running list of gates for the target qubit. Append the connection type onto the running list of connections for the control qubit and which qubit it is controlling (the target).
        self.qubits[target].gates.append('U')
        angles = [theta, phi, lambd]
        self.qubits[target].gateAngles.append(angles)
        for control in controls:
            self.qubits[control].connections.append('C')
            self.qubits[control].connectTo.append(target)

        # Get the earliest possible gate position within the circuit for each qubit between the target and control (inclusive). The max of this list will be used for the gate position for both the target and control. Then increase the earliest position for all qubits.
        earliestPositions = [self.qubits[idx].earliestPos for idx in range(min(min(controls), target), max(max(controls), target)+1)]
        position = max(earliestPositions)
        self.qubits[target].gatePos.append(position)
        for control in controls:
            self.qubits[control].connectPos.append(position)
        for Qidx in range(min(min(controls), target), max(max(controls), target)+1):
            self.qubits[Qidx].earliestPos = position + 1

        return self

    # SWAP gate
    def SWAP(self, target1, target2):

        # Append the gate onto the running list of gates for the target qubit (which we'll use target2 as for consistency with controlled gates). Append the connection type onto the running list of connections for the control qubit (target1) and which qubit it is controlling (target2).
        self.qubits[target2].gates.append('SWAP')
        angles = [None, None, None]
        self.qubits[target2].gateAngles.append(angles)
        self.qubits[target1].connections.append('SWAP')
        self.qubits[target1].connectTo.append(target2)

        # Get the earliest possible gate position within the circuit for each qubit between the targets (inclusive). The max of this list will be used for the gate position for both targets. Then increase the earliest position for all qubits.
        earliestPositions = [self.qubits[idx].earliestPos for idx in range(min(target1, target2), max(target1, target2)+1)]
        position = max(earliestPositions)
        self.qubits[target2].gatePos.append(position)
        self.qubits[target1].connectPos.append(position)
        for Qidx in range(min(target1, target2), max(target1, target2)+1):
            self.qubits[Qidx].earliestPos = position + 1

        return self
    
    ## OTHER CIRCUIT FUNCTIONS ##

    # Add a barrier to the circuit. The state vector does not change. A barrier is purely for visual purposes when displaying the circuit to divide the circuit into segments.
    def barrier(self):

        # Append a barrier to the first qubit. Use the last qubit as the connector to extend the barrier across the entire circuit. The max earliest position for all qubits is the position of the barrier. All qubits' earliest position is then updated to the position after the barrier.
        self.qubits[0].gates.append('B')
        angles = [None, None, None]
        self.qubits[0].gateAngles.append(angles)
        self.qubits[-1].connections.append('B')
        self.qubits[-1].connectTo.append(0)
        earliestPosition = max([qubit.earliestPos for qubit in self.qubits])
        self.qubits[0].gatePos.append(earliestPosition)
        self.qubits[-1].connectPos.append(earliestPosition)
        for qubit in self.qubits:
            qubit.earliestPos = earliestPosition + 1

    # Measure a qubit and store the result in a classical bit. This is a measurement in the computational basis (projection into the 0 or 1 state).
    def measure(self, targets):

        # Append the gate onto the running list of gates for the target qubit. Append the connection type onto the running list of connections for the control qubit and which qubit it is controlling (the target).
        if type(targets) == int:
            targets = [targets]
        for target in targets:
            self.qubits[target].gates.append('M')
            angles = [None, None, None]
            self.qubits[target].gateAngles.append(angles)
            self.cbits[target].connections.append('O')
            self.cbits[target].connectTo.append(target)

            # Get the earliest possible gate position within the circuit for each qubit between the target and control (inclusive). The max of this list will be used for the gate position for both the target and control. Then increase the earliest position for all qubits between the target and control (inclusive).
            earliestPositions = [qubit.earliestPos for qubit in self.qubits[target:]]
            for cbit in self.cbits[:target]:
                earliestPositions.append(cbit.earliestPos)
            position = max(earliestPositions)
            self.qubits[target].gatePos.append(position)
            self.cbits[target].connectPos.append(position)
            for qubit in self.qubits[target:]:
                qubit.earliestPos = position + 1
            for cbit in self.cbits:
                cbit.earliestPos = position + 1

        return self
    
    ## PREPROGRAMMED ALGORITHMS ##

    # Deutsch-Jozsa algorithm: this algorithm takes a provided oracle in the form of f({x1, x2, ..., xn}) = {0, 1} where each input x as 0 or 1. The oracle must be constant, i.e. all possible inputs yield the same output of 0 or 1, or balanced, i.e. half of all possible inputs yield 0 and the other half all yield 1. This algorithm can identify if the oracle is constant or balanced with 100% accuracy in just a single shot of the circuit. For constant oracles, all input qubits will be measured in the |0> state at the end of the algorithm, vs all |1>'s for balanced oracles.
    # 
    # While not an interesting algorithm from a practical perspective, this algorithm does show that certain problems can be completed much more efficiently on a quanutm computer than a classical computer. On a classical computer, it would take 2^(n-1)+1 shots to identify the oracle as constant or balanced with 100% accuracy, with n being the number of inputs (x's) to the oracle.
    #
    # See Algorithms.py to learn how to use an example oracle, or create your own constant or balanced oracle and pass it to 'oracle' as a python function.
    def DeutschJozsa(self, oracle, constantOracleOutput=0, balancedInputFlips=[]):
        Algorithms.DeutschJozsa(self, oracle, constantOracleOutput, balancedInputFlips)
        return
    
    # Quantum Fourier Transform (QFT): this algorithm converts qubits in the computational basis into the Fourier basis. This is commonly used as a sub-step within other algorithms.
    #
    # You may provide the number of qubits to perform the QFT on using numQubits. Note that the qubits involved must be sequential and ordered from least significant (lowest index) to most significant (highest index). To perform QFT on all qubits within the circuit, you may leave this argument as the default, and the function will get the number of qubits in the circuit.
    def QFT(self, numQubits=0):
        Algorithms.QFT(self, numQubits)
        return
    
    # Inverse Quantum Fourier Transform (IQFT): this algorithm converts qubits in the Fourier basis into the computational basis. This is commonly used as a sub-step within other algorithms.
    #
    # You may provide the number of qubits to perform the IQFT on using numQubits. Note that the qubits involved must be sequential and ordered from least significant (lowest index) to most significant (highest index). To perform IQFT on all qubits within the circuit, you may leave this argument as the default, and the function will get the number of qubits in the circuit.
    def IQFT(self, numQubits=0):
        Algorithms.IQFT(self, numQubits)
        return
    
    # Quantum phase estimation (QPE): this algorithm estimates the angle theta within the eigenvalue problem U|psi> = e^(2pi*i*theta)|psi>. The more qubits are included in the algorithm, the higher the precision of the algorithm (at the expense of higher computational cost). This is commonly used as a sub-step within other algorithms.
    #
    # The angle lambd = 2pi*theta must be passed as an argument.
    #
    # You may provide the number of qubits to perform the QPE on using numPrecisionQubits. To perform QPE on all qubits within the circuit (minus the final qubit which represents |psi>), you may leave this argument as the default, and the function will get the number of qubits in the circuit.
    def QPE(self, lambd, numPrecisionQubits=0):
        Algorithms.QPE(self, lambd, numPrecisionQubits)
        return

    # Assign the gate label, box, and connection property to be used for displaying the circuit
    def format_gate(self, gate, angles=[0, 0, 0]):
        # User-defined phase gate
        if gate == 'P':
            [theta, phi, lambd] = angles
            gateLabel = 'P\n(%s)'%str(round(theta, 2))
            size = 10
            bbox=dict(boxstyle='square', facecolor='white')
            arrowprops=dict()
        # Rotation gates
        elif gate in {'RX', 'RY', 'RZ'}:
            [theta, phi, lambd] = angles
            letters = list(gate)
            gateLabel = '$%s_%s$\n(%s)'%(letters[0], letters[1], str(round(theta, 2)))
            size = 10
            bbox=dict(boxstyle='square', facecolor='white')
            arrowprops=dict()
        # U gates
        elif gate == 'U':
            [theta, phi, lambd] = angles
            gateLabel = 'U\n(%s,%s,%s)'%(str(round(theta, 2)), str(round(phi, 2)), str(round(lambd, 2)))
            size = 10
            bbox=dict(boxstyle='square', facecolor='white')
            arrowprops=dict()
        # SWAP gates
        elif gate == 'SWAP':
            gateLabel = 'x'
            size = 25
            bbox=dict(boxstyle='square', pad=0, facecolor='none', edgecolor='none')
            arrowprops=dict(arrowstyle="-", edgecolor='black', linewidth=2)
        # Controls in controlled-gates
        elif gate == 'C':
            gateLabel = ' '
            size = 5
            bbox=dict(boxstyle='circle', facecolor='black', edgecolor='black')
            arrowprops=dict(arrowstyle="-", edgecolor='black', linewidth=2)
        # Output in measurement gates
        elif gate == 'O':
            gateLabel = ' '
            size = 10
            bbox=dict(boxstyle='circle', facecolor='none', edgecolor='black', linewidth=2)
            arrowprops=dict(arrowstyle="<|-", edgecolor='black', linewidth=2)
        elif gate == 'B':
            gateLabel = ' '
            size = 10
            bbox=dict(boxstyle='square', facecolor='gray', edgecolor='none')
            arrowprops=dict(arrowstyle="-", edgecolor='grey', linewidth=9)
        else:
            gateLabel = gate
            size = 15
            bbox=dict(boxstyle='square', facecolor='white')
            arrowprops=dict()
        
        return [gateLabel, size, bbox, arrowprops]

    # Create a figure of the circuit.
    def display_circuit(self):

        # Create the figure
        fig = plt.figure()
        ax = fig.add_subplot(111)

        # Get the circuit length
        circuitLength = 0
        for qubit in self.qubits:
            if len(qubit.gates) > 0:
                lastPos = qubit.gatePos[-1]
                if lastPos > circuitLength:
                    circuitLength = lastPos

        # Set some style parameters
        bitLabelPosition = 0
        bitLabelFontSize = 15
        ax.set(xlim=(0, circuitLength+1), ylim=(-1*(self.numQubits+self.numCbits), 1))
        ax.set_axis_off()

        # Display each classical bit. These are displayed first for proper layer ordering when displaying connections between qubits and classical bits.
        offset = 0.05
        for Bidx, cbit in enumerate(self.cbits):
            # Display the classical bit labels and a double line to represent their wires. "offset" creates a small spacing between the two plotted lines for each bit to give the double line visual.
            ax.annotate('$C_%s$'%Bidx, xy=(bitLabelPosition, -1*(Bidx+self.numQubits)), size=bitLabelFontSize, va='center', ha='center', bbox=dict(boxstyle='square', facecolor='white', edgecolor='none'))
            ax.plot([bitLabelPosition, circuitLength], [-1*(Bidx+self.numQubits)+offset, -1*(Bidx+self.numQubits)+offset], color='black')
            ax.plot([bitLabelPosition, circuitLength], [-1*(Bidx+self.numQubits)-offset, -1*(Bidx+self.numQubits)-offset], color='black')

            # Display each connection using the properties from format_gate
            for Cidx, connection in enumerate(cbit.connections):
                [connectLabel, size, bbox, arrowprops] = self.format_gate(connection)
                ax.annotate(connectLabel, xy=(cbit.connectPos[Cidx], -1*cbit.connectTo[Cidx]), xytext=(cbit.connectPos[Cidx], -1*(Bidx+self.numQubits)), size=size, va='center', ha='center', bbox=bbox, arrowprops=arrowprops)

        # Display each qubit
        for Qidx, qubit in enumerate(self.qubits):
            # Display the qubit labels and a horizontal line to represent the wire for each qubit's circuit.
            ax.annotate('$Q_%s$'%Qidx, xy=(bitLabelPosition, -1*Qidx), size=bitLabelFontSize, va='center', ha='center', bbox=dict(boxstyle='square', facecolor='white', edgecolor='none'))
            ax.plot([bitLabelPosition, circuitLength], [-1*Qidx, -1*Qidx], color='black')

        # Display each qubit connection using the properties from format_gate. These are displayed next for proper layer ordering of connections between qubits.
        for Qidx, qubit in enumerate(self.qubits):
            for Cidx, connection in enumerate(qubit.connections):
                [connectLabel, size, bbox, arrowprops] = self.format_gate(connection)
                ax.annotate(connectLabel, xy=(qubit.connectPos[Cidx], -1*qubit.connectTo[Cidx]), xytext=(qubit.connectPos[Cidx], -1*Qidx), size=size, va='center', ha='center', bbox=bbox, arrowprops=arrowprops)

        # Display each qubit gate using the properties from format_gate.
        for Qidx, qubit in enumerate(self.qubits):
            for Gidx, gate in enumerate(qubit.gates):
                if gate in {'P', 'RX', 'RY', 'RZ', 'U'}:
                    angles = qubit.gateAngles[Gidx]
                    [gateLabel, size, bbox, arrowprops] = self.format_gate(gate, angles)
                else:
                    [gateLabel, size, bbox, arrowprops] = self.format_gate(gate)
                ax.annotate(gateLabel, xy=(qubit.gatePos[Gidx], -1*Qidx), size=size, va='center', ha='center', bbox=bbox)

        plt.show()

        return
    
    def run(self, shots, hist=False):

        # Get the circuit length
        circuitLength = 0
        for qubit in self.qubits:
            if len(qubit.gates) > 0:
                lastPos = qubit.gatePos[-1]
                if lastPos > circuitLength:
                    circuitLength = lastPos

        # For the gates, create a matrix of I's with circuitLength rows and numQubits columns. Update the matrix with the gates and connections applied to each qubit, placing the gate in the corresponding qubit's column and the row of the circuit position where the gate/connection is applied. Do the same for the angles of each gate (only non-empty lists for rotation matrices)
        qubit_gates = [['I' for qubit in self.qubits] for pos in range(circuitLength)]
        gate_angles = [[[] for qubit in self.qubits] for pos in range(circuitLength)]
        for Qidx, qubit in enumerate(self.qubits):
            for Gidx, gate in enumerate(qubit.gates):
                # Subtract 1 since plotted gate positions start at 1 but Python indexing starts at 0
                gatePos = qubit.gatePos[Gidx] - 1
                qubit_gates[gatePos][Qidx] = gate
                gate_angles[gatePos][Qidx] = qubit.gateAngles[Gidx]
            for Cidx, connection in enumerate(qubit.connections):
                # Subtract 1 since plotted gate positions start at 1 but Python indexing starts at 0
                connectPos = qubit.connectPos[Cidx] - 1
                qubit_gates[connectPos][Qidx] = connection
        
        # Create an empty list of results with size equal to the total number of shots. Repeat the circuit's gates applications and measurements for each shot, storing the result from each shot in the list 'results'
        results = ['']*shots
        originalState = self.state
        for shot in range(shots):

            self.state = originalState

            for pos in range(circuitLength):

                # Get the next circuit position's list of gates
                gates = qubit_gates[pos]

                # Skip over barriers since they do not change the circuit's state
                if gates[0] == 'B':
                    continue

                # For controlled gates:
                if 'C' in gates:
                    # 
                    numControls = gates.count('C')
                    numOutcomes = 2**numControls
                    Coutcomes = [np.array([1]) for outcome in range(numOutcomes)]
                    outcomeCombos = CartesianProduct([0, 1], repeat=numControls)
                    outcomeCombos = [list(outcome) for outcome in outcomeCombos]

                    controlNum = 0
                    for Qidx, qubit in enumerate(self.qubits):
                        gateType = gates[Qidx]
                        angles = gate_angles[pos][Qidx]

                        if gateType == 'C':
                            idx = 0
                            for combo in outcomeCombos:
                                outcome = combo[controlNum]
                                if outcome == 0:
                                    Coutcomes[idx] = np.kron(gateMatrix('P0'), Coutcomes[idx])
                                else: # outcome == 1
                                    Coutcomes[idx] = np.kron(gateMatrix('P1'), Coutcomes[idx])
                                idx += 1
                            controlNum += 1

                        elif gateType != 'I':
                            idx = 0
                            for combo in outcomeCombos:
                                if all(combo):
                                    Coutcomes[idx] = np.kron(gateMatrix(gateType, angles), Coutcomes[idx])
                                else:
                                    Coutcomes[idx] = np.kron(gateMatrix('I'), Coutcomes[idx])
                                idx += 1

                        else: # gateType == 'I'
                            idx = 0
                            for combo in outcomeCombos:
                                Coutcomes[idx] = np.kron(gateMatrix('I'), Coutcomes[idx])
                                idx += 1

                    kronMatrix = np.sum(Coutcomes, axis=0)
                
                # For SWAP gates:
                elif 'SWAP' in gates:
                    # Create the Kronecker product matrix defining the gate operations. First create 3 separate matrices for the two target qubits to receive X, Y, and Z gates together. All other qubits get an identity. Then add the 3 matrices along with an identity and divide the sum by 2.
                    kronXMatrix = np.array([1])
                    kronYMatrix = np.array([1])
                    kronZMatrix = np.array([1])
                    for Qidx, qubit in enumerate(self.qubits):
                        gateType = gates[Qidx]
                        if gateType == 'SWAP':
                            kronXMatrix = np.kron(gateMatrix('X'), kronXMatrix)
                            kronYMatrix = np.kron(gateMatrix('Y'), kronYMatrix)
                            kronZMatrix = np.kron(gateMatrix('Z'), kronZMatrix)
                        else: # gateType == 'I'
                            kronXMatrix = np.kron(gateMatrix('I'), kronXMatrix)
                            kronYMatrix = np.kron(gateMatrix('I'), kronYMatrix)
                            kronZMatrix = np.kron(gateMatrix('I'), kronZMatrix)
                    kronMatrix = 0.5*(np.eye(np.size(kronXMatrix, 0)) + kronXMatrix + kronYMatrix + kronZMatrix)
                
                # For measurements (in computational basis):
                elif 'M' in gates:
                    # Create 2 Kronecker product matrices for the projection of the target qubit into the 0 or 1 state. All other qubits get an identity.
                    kron0Matrix = np.array([1])
                    kron1Matrix = np.array([1])
                    measuredQubit = 0
                    for Qidx, qubit in enumerate(self.qubits):
                        gateType = gates[Qidx]
                        if gateType == 'M':
                            kron0Matrix = np.kron(gateMatrix('P0'), kron0Matrix)
                            kron1Matrix = np.kron(gateMatrix('P1'), kron1Matrix)
                            measuredQubit = Qidx
                        else: # gateType == 'I'
                            kron0Matrix = np.kron(gateMatrix('I'), kron0Matrix)
                            kron1Matrix = np.kron(gateMatrix('I'), kron1Matrix)
                    
                    # For each of the 0 and 1 state projection matrices, apply the projection to the circuit's current state to get the resulting state. Transpose the circuit's current state (without the projection) and apply it to the resulting state (with the projection) to get the probability of the measurement outcome.
                    state0 = np.dot(kron0Matrix, self.state)
                    prob0 = np.vdot(state0, state0)
                    state1 = np.dot(kron1Matrix, self.state)
                    prob1 = np.vdot(state1, state1)
                    
                    # Generate a random number between 0 and 1. If it is less than the probability of the target qubit being in the 0 state, set the classical bit to 0 and update the circuit's state with the projection-into-0 state from above (normalized with the square root of the probability of measuring 0). Otherwise, set the classical bit to 1 and update the circuit's state with the projection-into-1 state (normalized).
                    if np.random.rand(1) < prob0:
                        self.cbits[measuredQubit].state = 0
                        self.state = state0 / np.sqrt(prob0)
                    else:
                        self.cbits[measuredQubit].state = 1
                        self.state = state1 / np.sqrt(prob1)

                # For single qubit gates:
                else:
                    #  Create the Kronecker product matrix defining the gate operations.
                    kronMatrix = np.array([1])
                    for Qidx, qubit in enumerate(self.qubits):
                        gateType = gates[Qidx]
                        if gateType in {'P', 'RX', 'RY', 'RZ', 'U'}:
                            angles = gate_angles[pos][Qidx]
                            kronMatrix = np.kron(gateMatrix(gateType, angles), kronMatrix)
                        else:
                            kronMatrix = np.kron(gateMatrix(gateType), kronMatrix)
                
                # The measurement operation above changes the circuit's state with the elif statement. If the current circuit position does not contain a measurement, apply the gates to the circuit's state, updating the state.
                if 'M' not in gates:
                    self.state = np.dot(kronMatrix, self.state)

            result = ''.join(reversed([str(cbit.state) for cbit in self.cbits]))
            # Reverse the bit order so bit 0 is on the far right
            result = '|' + result + '>'
            results[shot] = result

        if hist:
            labels, counts = np.unique(results, return_counts=True)
            plt.bar(labels, counts/shots, align='center')
            plt.title('Histogram of results of %i shots of the quantum circuit'%shots)
            plt.xlabel('Quantum circuit state')
            plt.xticks(rotation=45)
            plt.gca().set_xticks(labels)
            plt.ylabel('Probability')
            plt.show()

        return results

def main():
    return

if __name__ == "__main__":
    main()
