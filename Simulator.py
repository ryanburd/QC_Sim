## This file allows the user to create a quantum circuit, apply common qubit gates, and collect the measurement result of the circuit with many shots. A diagram of the circuit can be created.

import Algorithms
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Ellipse
from itertools import product as CartesianProduct
import tkinter

# Define common matrices used for gate operations. Since the rotation matrices need to receive angles, these matrices are packaged into a function instead of a dictionary, though this function essentially acts as a dictionary.
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

# Creates a qubit object, which stores all gates applied to the qubit, all connections the qubit is a part of (e.g. as a control for another target qubit's gate), and all algorithms that the qubit initiates. Note: while multiple qubits will be involved in an algorithm, only the highest index qubit inolved will receive the algorithm (and additional properties) appended to its lists. For circuit display purposes, it is only necessary to use one qubit to track algorithms, and the display code is written such that using the highest index qubit as the tracker is easiest.
class Qubit:

    def __init__(self):

        # gates = type of gate; gatePos = position of the gate along the circuit wire; gateAngles = theta, phi, lambd angles for phase and rotation gates
        self.gates = []
        self.gatePos = []
        self.gateAngles = []

        # connections = type of connection; connectTo = qubit index that the current qubit will connect to (such as as a control); connectPos = position of the connection along the circuit wire
        self.connections = []
        self.connectTo = []
        self.connectPos = []

        # algorithms = type of algorithm; algQubits = qubit indices involved in the algorithm; algNumQubits = number of qubits involved in the algorithm; algStart = circuit position where the algorithm begins; algEnd = circuit position where the algorithm ends
        self.algorithms = []
        self.algQubits = []
        self.algNumQubits = []
        self.algStart = []
        self.algEnd = []

        # The next available position along the qubit's circuit wire where a new gate can go. This is updated as more gates and algorithms are applied to the whole circuit and is used to determine gatePos, connectPos, and algStart above.
        self.earliestPos = 1

# Creates a classical bit object, which stores the bit's state (0 or 1) and all connections the bit is a part of (e.g. as storage for the result of measurement on a qubit).
class Cbit:

    def __init__(self, state):

        # state = state of the bit, i.e. 0 or 1; connections = type of connection; connectTo = qubit index that the current bit will connect to (such as as a measurement output storage); connectPos = position of the connection along the circuit wire
        self.state = state
        self.connections = []
        self.connectTo = []
        self.connectPos = []

        # The next available position along the bit's circuit wire where a new gate can go. This is updated as more gates and algorithms are applied to the whole circuit and is used to determine connectPos above.
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

        # If only one target is provided, place the index in a list. This is to remain consistent with situations where lists of targets are provided and avoids an error in the code below.
        if type(targets) == int:
            targets = [targets]

        # For each target, append the gate onto the running list of gates for the target qubit. No angles are needed, so a list of None's are appended as a placeholder. Use the qubit's earliest position for the gate position, then increment the earliest position for future gates.
        for target in targets:
            self.qubits[target].gates.append('X')
            angles = [None, None, None]
            self.qubits[target].gateAngles.append(angles)
            self.qubits[target].gatePos.append(self.qubits[target].earliestPos)
            self.qubits[target].earliestPos += 1

        return self

    # Pauli-Y gate
    def Y(self, targets):

        # If only one target is provided, place the index in a list. This is to remain consistent with situations where lists of targets are provided and avoids an error in the code below.
        if type(targets) == int:
            targets = [targets]

        # For each target, append the gate onto the running list of gates for the target qubit. No angles are needed, so a list of None's are appended as a placeholder. Use the qubit's earliest position for the gate position, then increment the earliest position for future gates.
        for target in targets:
            self.qubits[target].gates.append('Y')
            angles = [None, None, None]
            self.qubits[target].gateAngles.append(angles)
            self.qubits[target].gatePos.append(self.qubits[target].earliestPos)
            self.qubits[target].earliestPos += 1

        return self

    # Pauli-Z gate
    def Z(self, targets):

        # If only one target is provided, place the index in a list. This is to remain consistent with situations where lists of targets are provided and avoids an error in the code below.
        if type(targets) == int:
            targets = [targets]

        # For each target, append the gate onto the running list of gates for the target qubit. No angles are needed, so a list of None's are appended as a placeholder. Use the qubit's earliest position for the gate position, then increment the earliest position for future gates.
        for target in targets:
            self.qubits[target].gates.append('Z')
            angles = [None, None, None]
            self.qubits[target].gateAngles.append(angles)
            self.qubits[target].gatePos.append(self.qubits[target].earliestPos)
            self.qubits[target].earliestPos += 1

        return self
    
    # Hadamard gate
    def H(self, targets):

        # If only one target is provided, place the index in a list. This is to remain consistent with situations where lists of targets are provided and avoids an error in the code below.
        if type(targets) == int:
            targets = [targets]

        # For each target, append the gate onto the running list of gates for the target qubit. No angles are needed, so a list of None's are appended as a placeholder. Use the qubit's earliest position for the gate position, then increment the earliest position for future gates.
        for target in targets:
            self.qubits[target].gates.append('H')
            angles = [None, None, None]
            self.qubits[target].gateAngles.append(angles)
            self.qubits[target].gatePos.append(self.qubits[target].earliestPos)
            self.qubits[target].earliestPos += 1

        return self
    
    # Phase gate
    def S(self, targets):

        # If only one target is provided, place the index in a list. This is to remain consistent with situations where lists of targets are provided and avoids an error in the code below.
        if type(targets) == int:
            targets = [targets]

        # For each target, append the gate onto the running list of gates for the target qubit. No angles are needed, so a list of None's are appended as a placeholder. Use the qubit's earliest position for the gate position, then increment the earliest position for future gates.
        for target in targets:
            self.qubits[target].gates.append('S')
            angles = [None, None, None]
            self.qubits[target].gateAngles.append(angles)
            self.qubits[target].gatePos.append(self.qubits[target].earliestPos)
            self.qubits[target].earliestPos += 1

        return self
    
    # pi/8 gate
    def T(self, targets):

        # If only one target is provided, place the index in a list. This is to remain consistent with situations where lists of targets are provided and avoids an error in the code below.
        if type(targets) == int:
            targets = [targets]

        # For each target, append the gate onto the running list of gates for the target qubit. No angles are needed, so a list of None's are appended as a placeholder. Use the qubit's earliest position for the gate position, then increment the earliest position for future gates.
        for target in targets:
            self.qubits[target].gates.append('T')
            angles = [None, None, None]
            self.qubits[target].gateAngles.append(angles)
            self.qubits[target].gatePos.append(self.qubits[target].earliestPos)
            self.qubits[target].earliestPos += 1

        return self
    
    # phase gate
    def P(self, targets, theta):

        # If only one target is provided, place the index in a list. This is to remain consistent with situations where lists of targets are provided and avoids an error in the code below.
        if type(targets) == int:
            targets = [targets]

        # For each target, append the gate onto the running list of gates for the target qubit. Append theta and None's for phi and lambd. Use the qubit's earliest position for the gate position, then increment the earliest position for future gates.
        for target in targets:
            self.qubits[target].gates.append('P')
            angles = [theta, None, None]
            self.qubits[target].gateAngles.append(angles)
            self.qubits[target].gatePos.append(self.qubits[target].earliestPos)
            self.qubits[target].earliestPos += 1

        return self

    # R_X gate
    def RX(self, targets, theta):

        # If only one target is provided, place the index in a list. This is to remain consistent with situations where lists of targets are provided and avoids an error in the code below.
        if type(targets) == int:
            targets = [targets]

        # For each target, append the gate onto the running list of gates for the target qubit. Append theta and None's for phi and lambd. Use the qubit's earliest position for the gate position, then increment the earliest position for future gates.
        for target in targets:
            self.qubits[target].gates.append('RX')
            angles = [theta, None, None]
            self.qubits[target].gateAngles.append(angles)
            self.qubits[target].gatePos.append(self.qubits[target].earliestPos)
            self.qubits[target].earliestPos += 1

        return self
    
    # R_Y gate
    def RY(self, targets, theta):

        # If only one target is provided, place the index in a list. This is to remain consistent with situations where lists of targets are provided and avoids an error in the code below.
        if type(targets) == int:
            targets = [targets]

        # For each target, append the gate onto the running list of gates for the target qubit. Append theta and None's for phi and lambd. Use the qubit's earliest position for the gate position, then increment the earliest position for future gates.
        for target in targets:
            self.qubits[target].gates.append('RY')
            angles = [theta, None, None]
            self.qubits[target].gateAngles.append(angles)
            self.qubits[target].gatePos.append(self.qubits[target].earliestPos)
            self.qubits[target].earliestPos += 1

    # R_Z gate
    def RZ(self, targets, theta):

        # If only one target is provided, place the index in a list. This is to remain consistent with situations where lists of targets are provided and avoids an error in the code below.
        if type(targets) == int:
            targets = [targets]

        # For each target, append the gate onto the running list of gates for the target qubit. Append theta and None's for phi and lambd. Use the qubit's earliest position for the gate position, then increment the earliest position for future gates.
        for target in targets:
            self.qubits[target].gates.append('RZ')
            angles = [theta, None, None]
            self.qubits[target].gateAngles.append(angles)
            self.qubits[target].gatePos.append(self.qubits[target].earliestPos)
            self.qubits[target].earliestPos += 1

    # U gate
    def U(self, targets, theta, phi, lambd):

        # If only one target is provided, place the index in a list. This is to remain consistent with situations where lists of targets are provided and avoids an error in the code below.
        if type(targets) == int:
            targets = [targets]

        # For each target, append the gate onto the running list of gates for the target qubit. Append theta, phi, and lambd. Use the qubit's earliest position for the gate position, then increment the earliest position for future gates.
        for target in targets:
            self.qubits[target].gates.append('U')
            angles = [theta, phi, lambd]
            self.qubits[target].gateAngles.append(angles)
            self.qubits[target].gatePos.append(self.qubits[target].earliestPos)
            self.qubits[target].earliestPos += 1

    ## TWO QUBIT GATES ##

    # Controlled-X gate
    def CX(self, controls, target):

        # Append the gate onto the running list of gates for the target qubit. No angles are needed, so a list of None's are appended as a placeholder. Append a control, 'C', to the list of connections and the index of the target to the list of connectTo for the control qubits.
        self.qubits[target].gates.append('X')
        angles = [None, None, None]
        self.qubits[target].gateAngles.append(angles)
        for control in controls:
            self.qubits[control].connections.append('C')
            self.qubits[control].connectTo.append(target)

        # Get the earliest possible gate position within the circuit for each qubit between the target and control (inclusive). The max of this list will be used for the gate position for both the target and control. Then increment the earliest position for all qubits.
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

        # Append the gate onto the running list of gates for the target qubit. No angles are needed, so a list of None's are appended as a placeholder. Append a control, 'C', to the list of connections and the index of the target to the list of connectTo for the control qubits.
        self.qubits[target].gates.append('Y')
        angles = [None, None, None]
        self.qubits[target].gateAngles.append(angles)
        for control in controls:
            self.qubits[control].connections.append('C')
            self.qubits[control].connectTo.append(target)

        # Get the earliest possible gate position within the circuit for each qubit between the target and control (inclusive). The max of this list will be used for the gate position for both the target and control. Then increment the earliest position for all qubits.
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

        # Append the gate onto the running list of gates for the target qubit. No angles are needed, so a list of None's are appended as a placeholder. Append a control, 'C', to the list of connections and the index of the target to the list of connectTo for the control qubits.
        self.qubits[target].gates.append('Z')
        angles = [None, None, None]
        self.qubits[target].gateAngles.append(angles)
        for control in controls:
            self.qubits[control].connections.append('C')
            self.qubits[control].connectTo.append(target)

        # Get the earliest possible gate position within the circuit for each qubit between the target and control (inclusive). The max of this list will be used for the gate position for both the target and control. Then increment the earliest position for all qubits.
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

        # Append the gate onto the running list of gates for the target qubit. Append theta and None's for phi and lambd. Append a control, 'C', to the list of connections and the index of the target to the list of connectTo for the control qubits.
        self.qubits[target].gates.append('P')
        angles = [theta, None, None]
        self.qubits[target].gateAngles.append(angles)
        for control in controls:
            self.qubits[control].connections.append('C')
            self.qubits[control].connectTo.append(target)

        # Get the earliest possible gate position within the circuit for each qubit between the target and control (inclusive). The max of this list will be used for the gate position for both the target and control. Then increment the earliest position for all qubits.
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

        # Append the gate onto the running list of gates for the target qubit. Append theta and None's for phi and lambd. Append a control, 'C', to the list of connections and the index of the target to the list of connectTo for the control qubits.
        self.qubits[target].gates.append('RX')
        angles = [theta, None, None]
        self.qubits[target].gateAngles.append(angles)
        for control in controls:
            self.qubits[control].connections.append('C')
            self.qubits[control].connectTo.append(target)

        # Get the earliest possible gate position within the circuit for each qubit between the target and control (inclusive). The max of this list will be used for the gate position for both the target and control. Then increment the earliest position for all qubits.
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

        # Append the gate onto the running list of gates for the target qubit. Append theta and None's for phi and lambd. Append a control, 'C', to the list of connections and the index of the target to the list of connectTo for the control qubits.
        self.qubits[target].gates.append('RY')
        angles = [theta, None, None]
        self.qubits[target].gateAngles.append(angles)
        for control in controls:
            self.qubits[control].connections.append('C')
            self.qubits[control].connectTo.append(target)

        # Get the earliest possible gate position within the circuit for each qubit between the target and control (inclusive). The max of this list will be used for the gate position for both the target and control. Then increment the earliest position for all qubits.
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

        # Append the gate onto the running list of gates for the target qubit. Append theta and None's for phi and lambd. Append a control, 'C', to the list of connections and the index of the target to the list of connectTo for the control qubits.
        self.qubits[target].gates.append('RZ')
        angles = [theta, None, None]
        self.qubits[target].gateAngles.append(angles)
        for control in controls:
            self.qubits[control].connections.append('C')
            self.qubits[control].connectTo.append(target)

        # Get the earliest possible gate position within the circuit for each qubit between the target and control (inclusive). The max of this list will be used for the gate position for both the target and control. Then increment the earliest position for all qubits.
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

        # Append the gate onto the running list of gates for the target qubit. Append theta, phi, and lambd. Append a control, 'C', to the list of connections and the index of the target to the list of connectTo for the control qubits.
        self.qubits[target].gates.append('U')
        angles = [theta, phi, lambd]
        self.qubits[target].gateAngles.append(angles)
        for control in controls:
            self.qubits[control].connections.append('C')
            self.qubits[control].connectTo.append(target)

        # Get the earliest possible gate position within the circuit for each qubit between the target and control (inclusive). The max of this list will be used for the gate position for both the target and control. Then increment the earliest position for all qubits.
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

        # Append the gate onto the running list of gates for the target qubit (which we'll use target2 as for consistency with controlled gates). No angles are needed, so a list of None's are appended as a placeholder. Append the connection type onto the running list of connections for the control qubit (target1) and which qubit it is controlling (target2).
        self.qubits[target2].gates.append('SWAP')
        angles = [None, None, None]
        self.qubits[target2].gateAngles.append(angles)
        self.qubits[target1].connections.append('SWAP')
        self.qubits[target1].connectTo.append(target2)

        # Get the earliest possible gate position within the circuit for each qubit between the targets (inclusive). The max of this list will be used for the gate position for both targets. Then increment the earliest position for all qubits.
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

        # Append a barrier to the highest index qubit. Similar to tracking algorithms, only one qubit needs to act as the tracker, and the cirucit display code is written such that using the last qubit is easiest. The max earliest position for all qubits is the position of the barrier. All qubits' earliest position is then updated to the position after the barrier.
        self.qubits[-1].gates.append('B')
        angles = [None, None, None]
        self.qubits[-1].gateAngles.append(angles)
        earliestPosition = max([qubit.earliestPos for qubit in self.qubits])
        self.qubits[-1].gatePos.append(earliestPosition)
        for qubit in self.qubits:
            qubit.earliestPos = earliestPosition + 1

    # Measure a qubit and store the result in a classical bit. This is a measurement in the computational basis (projection into the 0 or 1 state).
    def measure(self, targets):

        # If only one target is provided, place the index in a list. This is to remain consistent with situations where lists of targets are provided and avoids an error in the code below.
        if type(targets) == int:
            targets = [targets]

        # For each target, append the gate onto the running list of gates for the target qubit. No angles are needed, so a list of None's are appended as a placeholder. Append an output, 'O', to the list of connections and the index of the target to the list of connectTo for the classical bit that will store the measurement outcome. For simplicity, the classical bit with the same index as the target qubit will be used.
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
    # You may provide the qubit indices to perform the IQFT on using algQubits. Note that the qubits involved must be sequential and ordered from least significant (lowest index) to most significant (highest index). To perform IQFT on all qubits within the circuit, you may leave this argument as the default, and the function will create a list of all qubit indices in the circuit.
    def IQFT(self, algQubits=None):

        # If no qubits are provided in algQubits, use all the qubits in the circuit. Set numQubits as the length of the qubits involved.
        if algQubits == None:
            algQubits = list(range(self.numQubits))
        numQubits = len(algQubits)

        # Use the highest index qubit as the algorithm tracker. Append the algorithm type, qubit indices involved, and number of qubits involved to the lists for the tracker.
        algTracker = max(algQubits)
        self.qubits[algTracker].algorithms.append('IQFT')
        self.qubits[algTracker].algQubits.append(algQubits)
        self.qubits[algTracker].algNumQubits.append(numQubits)

        # Get the earliest position for all qubits in the circuit (not just the algorithm qubits). The max will be used as the starting point for the algorithm. Append this value to algStart for the tracker. Increase the earliest position for all qubits in the circuit to the algStart + 1.
        earliestPosition = max([qubit.earliestPos for qubit in self.qubits])
        self.qubits[algTracker].algStart.append(earliestPosition)
        for qubit in self.qubits:
            qubit.earliestPos = earliestPosition + 1

        # Apply the algorithm. See Algorithms.py.
        Algorithms.IQFT(self, algQubits)

        # Get the earliest position for all qubits in the circuit (not just the algorithm qubits). The max will be used as the end point for the algorithm. Append this value to algEnd for the tracker. Increase the earliest position for all qubits in the circuit to the algEnd + 1.
        earliestPosition = max([qubit.earliestPos for qubit in self.qubits])
        self.qubits[algTracker].algEnd.append(earliestPosition)
        for qubit in self.qubits:
            qubit.earliestPos = earliestPosition + 1

        return
    
    # Quantum phase estimation (QPE): this algorithm estimates the angle theta within the eigenvalue problem U|psi> = e^(2pi*i*theta)|psi>. The more qubits are included in the algorithm, the higher the precision of the algorithm (at the expense of higher computational cost). This is commonly used as a sub-step within other algorithms.
    #
    # The angle lambd = 2pi*theta must be passed as an argument.
    #
    # You may provide the number of qubits to perform the QPE on using numPrecisionQubits. To perform QPE on all qubits within the circuit (minus the final qubit which represents |psi>), you may leave this argument as the default, and the function will get the number of qubits in the circuit.
    def QPE(self, lambd, numPrecisionQubits=0):
        Algorithms.QPE(self, lambd, numPrecisionQubits)
        return
    
    def Grover(self, numQubits=0, oracle='example'):
        Algorithms.Grover(self, numQubits, oracle)
        return

    # Calculate the dimensions in data coordinates of the boxes used to represent gates when displaying the circuit. xy = center of the box, in data coordinates; sizeX (sizeY) = total pixel size of the box's text in the X (Y) axis, accounting for number of letters in the text (X) and number of lines of text (Y); ax = figure axis.
    def gate_size(self, xy, sizeX, sizeY, ax):

        # Separate the center x and y coordinates. Convert them to display coordinates.
        xData, yData = xy
        xDisplay, yDisplay = ax.transData.transform([xData,yData])

        # The width (height) in display coordinates will be 2 times the text width (height).
        wDisplay = sizeX*2
        hDisplay = sizeY*2

        # Calculate the minimum x and y coordinates of the box in display coordinates (lower left corner). This is the center is half of the width (height) for the x (y) coordinate. Then convert these to data coordinates.
        xMinDisplay, yMinDisplay = xDisplay-wDisplay/2, yDisplay-hDisplay/2
        xMinData, yMinData = ax.transData.inverted().transform((xMinDisplay,yMinDisplay))

        # Calculate the maximum x and y coordinates by adding the width (height) to the minimum x (y) coordinate in display coordinates and converting to data coordinates. Then subtract the max and min for x (y) to get the width (height) in data coordinates.
        xMaxData, yMaxData = ax.transData.inverted().transform((xMinDisplay+wDisplay,yMinDisplay+hDisplay))
        wData = xMaxData - xMinData
        hData = yMaxData - yMinData

        # Return the lower left box coordinates, width, and height, all in data coordinates.
        return [(xMinData, yMinData), wData, hData]

    # Assign the gate label, box parameters, and connection parameters to be used for displaying the circuit. gate = gate type; xy = center of the box, in data coordinates; ax = figure axis; zorder = layer order for rendering the boxes in the circuit diagram; angles = theta, phi, lambd, when needed for phase and rotation gates.
    def format_gate(self, gate, xy, ax, zorder, angles=[0, 0, 0]):

        # User-defined phase gate
        if gate == 'P':

            # Gate label: 'P' with the theta value below it, with specified text size.
            [theta, phi, lambd] = angles
            thetaStr = str(round(theta, 2))
            gateLabel = 'P\n(%s)'%thetaStr
            textSize = 10

            # Get the character length of the gate label, using theta since it will be >= character length of 'P.' Ignore the parentheses since they do not add much width. If theta is a single digit whole number, add 1 to textLen. If theta is a decimal, subtract 1. These changes were made heuristically so that an appropriate padding is added around the text when rendering the box for the gate.
            textLen = len(thetaStr)
            if textLen == 1: textLen += 1
            else: textLen -= 1

            # Set the width and height of the gate label text in display coordinates, i.e. textSize times the number of characters in the X and Y directions, respectively. For the text width, the factor of 0.8 was determined heuristically for appropriate padding when rendering the box.
            textWidth = textSize*textLen*0.8
            textHeight = textSize*2

            # Get the gate coordinates, width, and height in data coordinates using the function gate_size.
            gateXY, gateWidth, gateHeight = self.gate_size(xy, textWidth, textHeight, ax)

            # Create the box object with appropriate style parameters. Since this is not a connection, leave the arrow properties blank.
            gateBox = Rectangle(gateXY, gateWidth, gateHeight, fc='white', ec='black', zorder=zorder)
            arrowprops=dict()

        # Rotation gates
        elif gate in {'RX', 'RY', 'RZ'}:

            # Gate label: 'R' with the axis of rotation as a subscript and the theta value below it, with specified text size. Split the rotation label into its two letters, 'R' and the axis of rotation, to make the axis a subscript.
            [theta, phi, lambd] = angles
            thetaStr = str(round(theta, 2))
            letters = list(gate)
            gateLabel = '$%s_%s$\n(%s)'%(letters[0], letters[1], thetaStr)
            textSize = 10

            # Get the character length of the gate label, using theta since it will usually be >= character length of gate type. Ignore the parentheses since they do not add much width. If theta is a single digit whole number, add 1 to textLen. If theta is a decimal, subtract 1. These changes were made heuristically so that an appropriate padding is added around the text when rendering the box for the gate.
            textLen = len(thetaStr)
            if textLen == 1: textLen += 1
            else: textLen -= 1

            # Set the width and height of the gate label text in display coordinates, i.e. textSize times the number of characters in the X and Y directions, respectively. For the text width, the factor of 0.8 was determined heuristically for appropriate padding when rendering the box.
            textWidth = textSize*textLen*0.8
            textHeight = textSize*2

            # Get the gate coordinates, width, and height in data coordinates using the function gate_size.
            gateXY, gateWidth, gateHeight = self.gate_size(xy, textWidth, textHeight, ax)

            # Create the box object with appropriate style parameters. Since this is not a connection, leave the arrow properties blank.
            gateBox = Rectangle(gateXY, gateWidth, gateHeight, fc='white', ec='black', zorder=zorder)
            arrowprops=dict()

        # U gates
        elif gate == 'U':

            # Gate label: 'U' with the angle values below it, with specified text size.
            [theta, phi, lambd] = angles
            thetaStr = str(round(theta, 2))
            phiStr = str(round(phi, 2))
            lambdStr = str(round(lambd, 2))
            gateLabel = 'U\n(%s,%s,%s)'%(thetaStr, phiStr, lambdStr)
            textSize = 10

            # Get the character length of the gate label, using the sum of the angles' text lengths since they will be > character length of 'U.' Ignore the parentheses since they do not add much width. Remove the decimal points. Add 2 to the text length sum. These changes were made heuristically so that an appropriate padding is added around the text when rendering the box for the gate.
            if '.' in thetaStr:
                thetaStr = thetaStr.replace('.','')
            if '.' in phiStr:
                phiStr = phiStr.replace('.','')
            if '.' in lambdStr:
                lambdStr = lambdStr.replace('.','')
            textLen = len(thetaStr+phiStr+lambdStr)+2

            # Set the width and height of the gate label text in display coordinates, i.e. textSize times the number of characters in the X and Y directions, respectively. For the text width, the factor of 0.55 was determined heuristically for appropriate padding when rendering the box.
            textWidth = textSize*textLen*0.55
            textHeight = textSize*2

            # Get the gate coordinates, width, and height in data coordinates using the function gate_size.
            gateXY, gateWidth, gateHeight = self.gate_size(xy, textWidth, textHeight, ax)

            # Create the box object with appropriate style parameters. Since this is not a connection, leave the arrow properties blank.
            gateBox = Rectangle(gateXY, gateWidth, gateHeight, fc='white', ec='black', zorder=zorder)
            arrowprops=dict()

        # SWAP gates
        elif gate == 'SWAP':

            # Gate label: 'x' with the specified text size. This size can be used for the text width and height since it is a single character and a single line.
            gateLabel = 'x'
            textSize = 25
            textWidth = textSize
            textHeight = textSize

            # Get the gate coordinates, width, and height in data coordinates using the function gate_size.
            gateXY, gateWidth, gateHeight = self.gate_size(xy, textWidth, textHeight, ax)

            # Create the box object with appropriate style parameters. Use a solid black line for the connection to the swapped qubit.
            gateBox = Rectangle(gateXY, gateWidth, gateHeight, fc='none', ec='none', zorder=zorder)
            arrowprops=dict(arrowstyle="-", edgecolor='black', linewidth=2)

        # Controls in controlled-gates
        elif gate == 'C':

            # Gate label: blank. A filled circle with the specified size will be used as the operation symbol. The symbol width and height are the same since the symbol is a circle.
            gateLabel = ' '
            textSize = 8
            symWidth = textSize
            symHeight = textSize

            # Get the gate coordinates, width, and height in data coordinates using the function gate_size.
            gateXY, gateWidth, gateHeight = self.gate_size(xy, symWidth, symHeight, ax)

            # Create the circle object with appropriate style parameters. Use a solid black line for the connection to the target qubit.
            gateBox = Ellipse(xy, gateWidth, gateHeight, fc='black', ec='black', zorder=zorder)
            arrowprops=dict(arrowstyle="-", edgecolor='black', linewidth=2)

        # Output in measurement gates
        elif gate == 'O':

            # Gate label: blank. An empty circle with the specified size will be used as the operation symbol. The symbol width and height are the same since the symbol is a circle.
            gateLabel = ' '
            textSize = 12
            symWidth = textSize
            symHeight = textSize

            # Get the gate coordinates, width, and height in data coordinates using the function gate_size.
            gateXY, gateWidth, gateHeight = self.gate_size(xy, symWidth, symHeight, ax)

            # Create the circle object with appropriate style parameters. Use a solid black line with an arrow that points to the circle to represent the measurement output from the qubit being stored in the classical bit.
            gateBox = Ellipse(xy, gateWidth, gateHeight, fc='none', ec='black', lw=2, zorder=zorder)
            arrowprops=dict(arrowstyle="<|-", edgecolor='black', linewidth=2)

        # Barriers
        elif gate == 'B':

            # Gate label: blank. A grey box spanning all the qubits will be used for the barrier symbol. The barrier width and height are the same, chosen heuristically.
            gateLabel = ' '
            textSize = 10
            barrierWidth = textSize
            barrierHeight = textSize

            # Get the gate coordinates, width, and height in data coordinates using the function gate_size.
            gateXY, gateWidth, gateHeight = self.gate_size(xy, barrierWidth, barrierHeight, ax)

            # Create the box object with appropriate style parameters. Since gateHeight would assume only one qubit has the barrier, add the total number of qubits in the circuit (-1 since gateHeight already includes one qubit) to stretch the barrier over the entire circuit. Since this is not a connection, leave the arrow properties blank.
            gateBox = Rectangle(gateXY, gateWidth, gateHeight+self.numQubits-1, fc='gray', ec='none', zorder=zorder)
            arrowprops=dict()

        # Gates with no rotation angles.
        else:

            # Gate label: the gate type, with specified text size. Since single letters are used for these gate types, text width and height equal the text size.
            gateLabel = gate
            textSize = 15
            textWidth = textSize
            textHeight = textSize

            # Get the gate coordinates, width, and height in data coordinates using the function gate_size.
            gateXY, gateWidth, gateHeight = self.gate_size(xy, textWidth, textHeight, ax)

            # Create the box object with appropriate style parameters. Since this is not a connection, leave the arrow properties blank.
            gateBox = Rectangle(gateXY, gateWidth, gateHeight, fc='white', ec='black', zorder=zorder)
            arrowprops=dict()
        
        return [gateLabel, textSize, arrowprops, gateBox]
    
    # Assign the algorithm label and box parameters to be used for displaying the circuit. algorithm = algorithm type; numQubits = number of qubits involved in the algorithm; xy = center of the box, in data coordinates; ax = figure axis; zorder = layer order for rendering the boxes in the circuit diagram.
    def format_algorithm(self, algorithm, numQubits, xy, ax, zorder):

        # Algorithm label: algorithm type. The box width is the text size times the character length of the algorithm type, with a factor of 0.5 determined heuristically for appropriate padding. The box height is just the text size since only a single line is used.
        algLabel = algorithm
        textSize = 15
        boxWidth = textSize*len(algorithm)*0.5
        boxHeight = textSize

        # Get the algorithm box coordinates, width, and height in data coordinates using the function gate_size.
        algXY, algWidth, algHeight = self.gate_size(xy, boxWidth, boxHeight, ax)

        # Create the box object with appropriate style parameters. Since algHeight would assume only one qubit has the barrier, add the total number of qubits involved in the algorithm (-1 since algHeight already includes one qubit) to stretch the box over all qubits involved. Since this is not a connection, leave the arrow properties blank.
        algBox = Rectangle(algXY, algWidth, algHeight+numQubits-1, fc='white', ec='black', zorder=zorder)
        arrowprops=dict()

        return [algLabel, textSize, arrowprops, algBox]

    # Create a figure showing a diagram of the circuit.
    def display_circuit(self):

        # Get the screen size and dpi to scale the figure window.
        win = tkinter.Tk()
        screenWidth = win.winfo_screenwidth()
        screenHeight = win.winfo_screenheight()
        dpi = win.winfo_fpixels('1i')
        win.withdraw() # don't show the tkinter window

        # Create the figure with size 0.75 times the screen width and height.
        fig = plt.figure(figsize=(screenWidth/dpi*0.75, screenHeight/dpi*0.75))
        ax = fig.add_subplot(111)

        # Get the circuit length. For each qubit, get the last gate position applied to the qubit. Update circuitLength to be the highest gate position in the circuit.
        circuitLength = 0
        for qubit in self.qubits:
            if len(qubit.gates) > 0:
                lastPos = qubit.gatePos[-1]
                if lastPos > circuitLength:
                    circuitLength = lastPos
        
        # Algorithms will be displayed as simple boxes without showing the individual gates comprising the algorithm (for simplicity; use Algorithms.py directly to display the individial gates). Calculate a circuit length offset to account for gates applied after algorithms being shifted to earlier positions in the diagram. For each qubit, get the number of positions taken up by each algorithm that the qubit is a tracker for. Add the number to the running total.
        circuitLengthOffset = 0
        for qubit in self.qubits:
            for Aidx, algorithm in enumerate(qubit.algorithms):
                circuitLengthOffset += qubit.algEnd[Aidx]-qubit.algStart[Aidx]

        # Set some style parameters. Start the bit labels at x=0. The diagram will display the x axis from 0 to the circuit length, with the offset subtracted out. Along the y axis, qubits will be placed at the negative of their index so that the lowest index qubit will be at the top of the diagram. Classical bits will be below the qubits. Set the y axis to be from the negative of the total number of qubits and classical bits to 1, which prevents a buffer with the top qubit.
        bitLabelPosition = 0
        bitLabelFontSize = 15
        ax.set(xlim=(0, circuitLength-circuitLengthOffset+1), ylim=(-1*(self.numQubits+self.numCbits), 1))
        ax.set_axis_off()

        # Begin the circuit element rendering order at 3. This will be increased when necessary to ensure proper display ordering of the circuit elements.
        zorder = 3

        ## Display all the gate operations in the circuit.

        # algorithmOn is a boolean that tracks whether the current circuit position is part of an algorithm. algTracker is the qubit that is the tracker for the current algorithm when algorithmOn is True. posOffset keeps track of how many positions algorithms take up so that gates after algorithms can be shifted to earlier positions in the diagram.
        algorithmOn = False
        algTracker = None
        posOffset = 0

        # Loop over the circuit position until maxGatePos (highest position of a gate among all qubits) is exceeded.
        position = 1
        maxGatePos = max([qubit.gatePos[-1] for qubit in self.qubits])
        while position <= maxGatePos:

            # If algorithmOn is True, the current position is part of an algorithm.
            if algorithmOn:
                # If the current position is the algorithm end, set algorithmOn to False.
                if position in self.qubits[algTracker].algEnd:
                    algorithmOn = False
                # Increment the posOffset since the current position will not be displayed in the diagram.
                posOffset += 1
            
            # If algorithmOn is False, loop over all the qubits and display their gate, connection, or algorithm being tracked at the current circuit position.
            else:
                for Qidx, qubit in reversed(list(enumerate(self.qubits))):

                    # Coordinates to place the qubit's gate at; x = current position minus the current posOffset; y = negative of the current qubit index.
                    xy = (position-posOffset, -1*Qidx)

                    # If the current position is an algorithm, display the algorithm box over all qubits involved.
                    if position in qubit.algStart:

                        # Set algorithmOn to True and the algTracker to the current qubit. Get the index of this algorithm among all the current qubit's algorithms tracked (Aidx). Get the algorithm type (alg) and number of qubits involved (algNumQubits).
                        algorithmOn = True
                        algTracker = Qidx
                        Aidx = qubit.algStart.index(position)
                        alg = qubit.algorithms[Aidx]
                        algNumQubits = qubit.algNumQubits[Aidx]

                        # Format the algorithm box using format_algorithm. Add gateBox as a patch.
                        [gateLabel, textSize, arrowprops, gateBox] = self.format_algorithm(alg, algNumQubits, xy, ax, zorder)
                        ax.add_patch(gateBox)

                        # Add the gate label as an annotation over the box. To center the label vertically, find the mean of the first and last qubit involved in the algorithm and negate the result.
                        y = -1*np.mean([list(qubit.algQubits[Aidx])[0], list(qubit.algQubits[Aidx])[-1]])
                        ax.annotate(gateLabel, xy=(xy[0],y), size=textSize, va='center', ha='center', zorder=zorder)

                        # No algorithm gates should be displayed. Break out of the qubit for loop since algorithmOn is now True. The for loop will not be reentered until the entire algorithm gate sequence has been skipped over (in the display only).
                        break

                    # If the current position is in the qubit's connectPos, display the connection.
                    if position in qubit.connectPos:

                        # Reduce the zorder to render the connection below the target gate.
                        zorder -= 1

                        # Get the index of this connection among all the current qubit's connections (Cidx). Get the connection type (connection) and the target qubit to connect to (connectTo).
                        Cidx = qubit.connectPos.index(position)
                        connection = qubit.connections[Cidx]
                        connectTo = qubit.connectTo[Cidx]

                        # Format the connection symbol using format_gate. Add connectSym as a patch.
                        [connectLabel, textSize, arrowprops, connectSym] = self.format_gate(connection, xy, ax, zorder)
                        ax.add_patch(connectSym)

                        # Add the connection label and connector as an annotation. xy = position of the target; xytext = position of the connection symbol.
                        ax.annotate(connectLabel, xy=(position, -1*connectTo), xytext=xy, size=textSize, va='center', ha='center', arrowprops=arrowprops, zorder=zorder)

                        # Increase the zorder back to the gate layer.
                        zorder += 1
                    
                    # If the current position is within the qubit's gatePos, display the gate.
                    if position in qubit.gatePos:

                        # Get the index of this gate among all the current qubit's gates (Gidx). Get the gate type (gate) and angles for phase or rotation gates (angles).
                        Gidx = qubit.gatePos.index(position)
                        gate = qubit.gates[Gidx]
                        angles = qubit.gateAngles[Gidx]

                        # Format the gate box using format_gate. Add gateBox as a patch. If the gate is a phase or rotation gate, provide the angles to the function as well.
                        if gate in {'P', 'RX', 'RY', 'RZ', 'U'}:
                            angles = qubit.gateAngles[Gidx]
                            [gateLabel, textSize, arrowprops, gateBox] = self.format_gate(gate, xy, ax, zorder, angles)
                        else:
                            [gateLabel, textSize, arrowprops, gateBox] = self.format_gate(gate, xy, ax, zorder)
                        ax.add_patch(gateBox)

                        # Add the gate label as an annotation.
                        ax.annotate(gateLabel, xy=xy, size=textSize, va='center', ha='center', zorder=zorder)

                # Display each classical bit connection using the properties from format_gate
                for Bidx, cbit in reversed(list(enumerate(self.cbits))):

                    # Coordinates to place the bit's operation at; x = current position minus the current posOffset; y = negative of the current bit index + total number of qubits in the circuit.
                    xy = (position-posOffset, -1*(Bidx+self.numQubits))

                    # If the current position is in the bit's connectPos, display the connection.
                    if position in cbit.connectPos:

                        # Reduce the zorder to render the connection below the target gate.
                        zorder -= 1

                        # Get the index of this connection among all the current bit's connections (Cidx). Get the connection type (connection) and the target qubit to connect to (connectTo).
                        Cidx = cbit.connectPos.index(position)
                        connection = cbit.connections[Cidx]
                        connectTo = cbit.connectTo[Cidx]

                        # Format the connection symbol using format_gate. Add connectSym as a patch.
                        [connectLabel, textSize, arrowprops, connectSym] = self.format_gate(connection, xy, ax, zorder)
                        ax.add_patch(connectSym)

                        # Add the connection label and connector as an annotation. xy = position of the qubit being connected to; xytext = position of the connection symbol.
                        ax.annotate(connectLabel, xy=(xy[0], -1*connectTo), xytext=xy, size=textSize, va='center', ha='center', arrowprops=arrowprops, zorder=zorder)

                        # Increase the zorder back to the gate layer.
                        zorder += 1

            # Increment the position for the next loop iteration.
            position += 1

        # Display each classical bit. These are displayed first for proper layer ordering when displaying connections between qubits and classical bits.
        offset = 0.05
        for Bidx, cbit in enumerate(self.cbits):
            # Display the classical bit labels and a double line to represent their wires. "offset" creates a small spacing between the two plotted lines for each bit to give the double line visual.
            ax.plot([bitLabelPosition, circuitLength-posOffset], [-1*(Bidx+self.numQubits)+offset, -1*(Bidx+self.numQubits)+offset], color='black', zorder=1)
            ax.plot([bitLabelPosition, circuitLength-posOffset], [-1*(Bidx+self.numQubits)-offset, -1*(Bidx+self.numQubits)-offset], color='black', zorder=1)
            ax.annotate('$C_%s$'%Bidx, xy=(bitLabelPosition, -1*(Bidx+self.numQubits)), size=bitLabelFontSize, va='center', ha='center', bbox=dict(boxstyle='square', facecolor='white', edgecolor='none'), zorder=2)

        # Display each qubit label and wire.
        for Qidx, qubit in enumerate(self.qubits):
            # Display the qubit labels and a horizontal line to represent the wire for each qubit's circuit.
            ax.plot([bitLabelPosition, circuitLength-posOffset], [-1*Qidx, -1*Qidx], color='black', zorder=1)
            ax.annotate('$Q_%s$'%Qidx, xy=(bitLabelPosition, -1*Qidx), size=bitLabelFontSize, va='center', ha='center', bbox=dict(boxstyle='square', facecolor='white', edgecolor='none'), zorder=2)

        plt.show()

        return
    
    # Run the circuit to calculate the final state of the qubits.
    def run(self, shots, hist=False):

        # Get the circuit length. For each qubit, get the last gate position applied to the qubit. Update circuitLength to be the highest gate position in the circuit.
        circuitLength = 0
        for qubit in self.qubits:
            if len(qubit.gates) > 0:
                lastPos = qubit.gatePos[-1]
                if lastPos > circuitLength:
                    circuitLength = lastPos

        # For the gates, create a matrix of I's with circuitLength rows and numQubits columns. Update the matrix with the gates and connections applied to each qubit, placing the gate in the corresponding qubit's column and the row of the circuit position where the gate/connection is applied. Do the same for the angles of each gate.
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

        # Store for initial state of the circuit (usually |0> for each qubit) to reset the circuit state at the start of each shot.
        initialState = self.state

        # Create an empty list of results with size equal to the total number of shots. Repeat the circuit's gate applications and measurements for each shot, storing the result from each shot in the list 'results'.
        results = ['']*shots
        for shot in range(shots):

            # Reset the circuit state
            self.state = initialState

            # Loop over each circuit position, creating the Kronecker matrix using all the gates at the position and applying it to the current circuit state (self.state).
            for pos in range(circuitLength):

                # Get the current position's list of gates.
                gates = qubit_gates[pos]

                # Skip over barriers since they do not change the circuit's state. Go to the next circuit position.
                if gates[-1] == 'B':
                    continue

                # Skip over algorithm indicators since they are only used for displaying the circuit. The actual gates within the algorithm start at the next circuit position. Go to the next position.
                if 'IQFT' in gates[-1]:
                    continue

                # If there is a controlled gate within the current position:
                if 'C' in gates:

                    # Get the number of control qubits for the controlled gate.
                    numControls = gates.count('C')

                    # To calculate the matrix that represents a controlled gate, we have to calculate a Kronecker matrix for each possible combination of controlled qubit outcomes and then sum the matrices. These Kronecker matrices are calculated from the projection matrices into either the |0> or |1> state for each control qubit and the target's intended gate if all control qubits are projected into the |1> state (identity matrix otherwise).

                    # If there are numControls control qubits that can each be measured in either the |0> or |1> state, then there are 2**numControls different outcomes for the controlled gate. Create a list that will contain the Kronecker matrix for each possible outcome.
                    numOutcomes = 2**numControls
                    outcomeKroneckers = [np.array([1]) for outcome in range(numOutcomes)]

                    # Each Kronecker matrix within outcomeKroneckers is associated with a unique combination of control qubit measurement outcomes. Create a list of each combination of control qubit outcomes using the Cartesian product. E.g. if there are 2 control qubits, the list would contain: '00', '01', '10', '11'
                    controlOutcomeCombos = CartesianProduct([0, 1], repeat=numControls)
                    controlOutcomeCombos = [list(outcome) for outcome in controlOutcomeCombos]

                    # Within each combo in the list controlOutcomeCombos, the first number represents the outcome for the first control qubit, the second number for the second control qubit, etc. To keep track of which index to use for each control qubit, the variable controlNum will track how many control qubits we have already encountered for the current controlled gate. The variable starts at 0 so that the first control qubit will use index 0, and the variable will be incremented every time a control qubit is encountered.
                    controlNum = 0

                    # Loop over each qubit to get its gate for the current circuit position.
                    for Qidx, qubit in enumerate(self.qubits):

                        # Get the gate type and the angles for it (for phase and rotation gates).
                        gateType = gates[Qidx]
                        angles = gate_angles[pos][Qidx]

                        # If the qubit is a control:
                        if gateType == 'C':

                            # Loop over the different control qubit outcome combinations. The variable idx will track which matrix within outcomeKroneckers to apply the current qubit's projection matrix to.
                            idx = 0
                            for combo in controlOutcomeCombos:

                                # Get the outcome within the current combo for the current control qubit. Use controlNum as the index (see explanation of controlNum above).
                                outcome = combo[controlNum]

                                # For an outcome of 0 for the current control qubit, use the projection into |0> for the Kronecker matrix.
                                if outcome == 0:
                                    outcomeKroneckers[idx] = np.kron(gateMatrix('P0'), outcomeKroneckers[idx])
                                
                                # For an outcome of 1 for the current control qubit, use the projection into |1> for the Kronecker matrix.
                                else: # outcome == 1
                                    outcomeKroneckers[idx] = np.kron(gateMatrix('P1'), outcomeKroneckers[idx])

                                # Increment idx so that the projection matrix for the next combo in controlOutcomeCombo will be applied to the next matrix within outcomeKroneckers.
                                idx += 1

                            # The current control qubit is done, so increment controlNum so that the next control qubit encountered will use the next index within each combo in controlOutcomeCombo.
                            controlNum += 1

                        # If the gate type is not a control or an identity, this is the target qubit.
                        elif gateType != 'I':

                            # Loop over the possible outcome combinations of the control qubits.
                            idx = 0
                            for combo in controlOutcomeCombos:

                                # If all control qubits measure |1> within the current combo, apply the intended target gate to the corresponding matrix within outcomeKroneckers.
                                if all(combo):
                                    outcomeKroneckers[idx] = np.kron(gateMatrix(gateType, angles), outcomeKroneckers[idx])
                                
                                # Otherwise, when at least one control qubit measures to 0, apply the identity matrix.
                                else:
                                    outcomeKroneckers[idx] = np.kron(gateMatrix('I'), outcomeKroneckers[idx])

                                # Go to the next matrix within outcomeKroneckers.
                                idx += 1

                        # For all other qubits in the circuit that are not involved within the controlled gate, apply an identity matrix to each matrix within outcomeKroneckers.
                        else:
                            idx = 0
                            for combo in controlOutcomeCombos:
                                outcomeKroneckers[idx] = np.kron(gateMatrix('I'), outcomeKroneckers[idx])
                                idx += 1

                    # With the Kronecker matrix for each combination of control qubit outcomes calculated, sum the matrices to get the final matrix that represents the operation of the controlled gate.
                    kronMatrix = np.sum(outcomeKroneckers, axis=0)
                
                # For SWAP gates:
                elif 'SWAP' in gates:
                    
                    # SWAP gates can be decomposed into 1/2 the sum of Kronecker matrices that apply an identity to each qubit, an X gate to each qubit, a Y gate to each qubit, and a Z gate to each qubit.

                    # Create 3 separate Kronecker matrices for the two target qubits to receive X, Y, and Z gates together. All other qubits will get an identity.
                    kronXMatrix = np.array([1])
                    kronYMatrix = np.array([1])
                    kronZMatrix = np.array([1])

                    for Qidx, qubit in enumerate(self.qubits):
                        
                        # For each qubit, get the gate type.
                        gateType = gates[Qidx]

                        # For SWAP gates, apply an X, Y, and Z gate to the corresponding Kronecker matrix.
                        if gateType == 'SWAP':
                            kronXMatrix = np.kron(gateMatrix('X'), kronXMatrix)
                            kronYMatrix = np.kron(gateMatrix('Y'), kronYMatrix)
                            kronZMatrix = np.kron(gateMatrix('Z'), kronZMatrix)
                        
                        # Otherwise, the gate type will be an identity. Apply the identity to all 3 Kronecker matrices.
                        else:
                            kronXMatrix = np.kron(gateMatrix('I'), kronXMatrix)
                            kronYMatrix = np.kron(gateMatrix('I'), kronYMatrix)
                            kronZMatrix = np.kron(gateMatrix('I'), kronZMatrix)
                    
                    # Add the 3 Kronecker matrices to an identity matrix of the same size and take 1/2 the sum. This is the final matrix that represents the SWAP gate operation.
                    kronMatrix = 0.5*(np.eye(np.size(kronXMatrix, 0)) + kronXMatrix + kronYMatrix + kronZMatrix)
                
                # For measurements (in computational basis):
                elif 'M' in gates:

                    # Create 2 Kronecker product matrices for the projection of the target qubit into the |0> or |1> state. All other qubits will get an identity.
                    kron0Matrix = np.array([1])
                    kron1Matrix = np.array([1])

                    measuredQubit = 0
                    for Qidx, qubit in enumerate(self.qubits):

                        # For each qubit, get the gate type.
                        gateType = gates[Qidx]

                        # For measurements, apply the projection matrix into |0> and |1> to the corresponding Kronecker matrix.
                        if gateType == 'M':
                            kron0Matrix = np.kron(gateMatrix('P0'), kron0Matrix)
                            kron1Matrix = np.kron(gateMatrix('P1'), kron1Matrix)

                            # Store the index of the qubit being measured in measuredQubit.
                            measuredQubit = Qidx
                        
                        # Otherwise, the gate type is an identity. Apply an identity to each Kronecker matrix.
                        else:
                            kron0Matrix = np.kron(gateMatrix('I'), kron0Matrix)
                            kron1Matrix = np.kron(gateMatrix('I'), kron1Matrix)
                    
                    # For each of the |0> and |1> state projection matrices, apply the projection to the circuit's current state to get the resulting state. Transpose the circuit's current state (without the projection) and apply it to the resulting state (with the projection) to get the probability of the measurement outcome.
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

                        # Get the gate type.
                        gateType = gates[Qidx]

                        # For phase and rotation gates, get the angles needed to define the gate and apply the gate to the Kronecker matrix.
                        if gateType in {'P', 'RX', 'RY', 'RZ', 'U'}:
                            angles = gate_angles[pos][Qidx]
                            kronMatrix = np.kron(gateMatrix(gateType, angles), kronMatrix)
                        
                        # Other gates don't need angles and can be applied to the Kronecker matrix.
                        else:
                            kronMatrix = np.kron(gateMatrix(gateType), kronMatrix)
                
                # The measurement operation above changes the circuit's state within the elif statement. If the current circuit position does not contain a measurement, apply the gates to the circuit's state, updating the state.
                if 'M' not in gates:
                    self.state = np.dot(kronMatrix, self.state)

            # Create a string containing the classical bit states at the end of the circuit. Reverse the bit order so bit 0 is on the far right. Style the list as a ket since thise is the state of the qubits, despite being stored in the classical bits.
            result = ''.join(reversed([str(cbit.state) for cbit in self.cbits]))
            result = '|' + result + '>'

            # Store the current result in the 'shot' index in the list of all results.
            results[shot] = result

        # If you want to create a histogram of your results:
        if hist:

            # Create a list 'labels' that contains each unique final circuit state within the list of all results. Get the number of times each unique state was obtained and store in the list 'counts'.
            labels, counts = np.unique(results, return_counts=True)

            # Create a bar graph of the unique states with the number of counts of each state normalized to the total number of shots taken for the circuit. The bar graph thus gives the percent chance of obtaining each unique state.
            plt.bar(labels, counts/shots, align='center')

            # Apply labels to the graph and axes.
            plt.title('Histogram of results of %i shots of the quantum circuit'%shots)
            plt.xlabel('Quantum circuit state')
            plt.ylabel('Probability')

            # Use the circuit states stored in 'labels' as the x-axis labels for the bars. Rotate the labels to prevent overlap.
            plt.xticks(rotation=45)
            plt.gca().set_xticks(labels)
            
            # Display the histogram.
            plt.show()

        # Return the list of results for further computation, if needed.
        return results

def main():
    return

if __name__ == "__main__":
    main()
