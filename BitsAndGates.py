## This file allows the user to create a quantum circuit with a provided number of qubits and classical bits.
## Common qubit gates to operate on the quantum circuit are defined.

import numpy as np
import matplotlib.pyplot as plt

# Define common matrices used for gate operations.
Xmatrix = np.array([[0, 1],
                    [1, 0]])
Ymatrix = np.array([[ 0, -1j],
                    [1j,   0]])
Zmatrix = np.array([[1,  0],
                    [0, -1]])
Hmatrix = 1/np.sqrt(2)*np.array([[1,  1],
                                 [1, -1]])
Smatrix = np.array([[1,  0],
                    [0, 1j]])
Tmatrix = np.array([[1,                  0],
                    [0, np.exp(1j*np.pi/4)]])

# Projection matrices into the computational basis. Used for constructing controlled-U gates and measurements.
proj0matrix = np.array([[1, 0],
                        [0, 0]])
proj1matrix = np.array([[0, 0],
                        [0, 1]])

class Qubit:

    def __init__(self):
        self.gates = []
        self.gatePos = []
        self.connections = []
        self.connectTo = []
        self.connectPos = []
        self.earliestPos = 1

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
    def __init__(self, numQubits, numCbits):

        # Create a list of qubits. Each instance of the class Qubit will store the gates applied to the qubit. This is useful for creating a diagram of the circuit.
        self.qubits = [Qubit() for qubit in range(numQubits)]
        
        # Form the state of all the qubits in the circuit. Assume all qubits are initialized in the |0> state, [1, 0]
        self.state = np.array([1])
        for qubit in self.qubits:
            self.state = np.tensordot(self.state, [1, 0], axes=0).reshape(len(self.state)*2, 1)
        
        # Create a list of classical bits, each initialized in the 0 state.
        self.cbits = [Cbit(0) for cbit in range(numCbits)]

    ## SINGLE QUBIT GATES ##

    # Pauli-X gate
    def X(self,target):

        # Create the Kronecker product matrix defining the gate operation. The target qubit gets an X gate, all other qubits get identity.
        kronMatrix = np.array([1])
        for qIndex, qubit in enumerate(self.qubits):
            if qIndex == target:
                kronMatrix = np.kron(kronMatrix, Xmatrix)
            else:
                kronMatrix = np.kron(kronMatrix, np.eye(2))
        
        # Apply the gate to the circuit's state, updating the state.
        self.state = np.dot(kronMatrix, self.state)

        # Append the gate onto the running list of gates for the target qubit. Use the qubit's earliest position for the gate position, then increment the earliest position for potential future gates.
        self.qubits[target].gates.append('X')
        self.qubits[target].gatePos.append(self.qubits[target].earliestPos)
        self.qubits[target].earliestPos += 1

        return self

    # Pauli-Y gate
    def Y(self,target):

        # Create the Kronecker product matrix defining the gate operation. The target qubit gets a Y gate, all other qubits get identity.
        kronMatrix = np.array([1])
        for qIndex, qubit in enumerate(self.qubits):
            if qIndex == target:
                kronMatrix = np.kron(kronMatrix, Ymatrix)
            else:
                kronMatrix = np.kron(kronMatrix, np.eye(2))
        
        # Apply the gate to the circuit's state, updating the state.
        self.state = np.dot(kronMatrix, self.state)

        # Append the gate onto the running list of gates for the target qubit. Use the qubit's earliest position for the gate position, then increment the earliest position for potential future gates.
        self.qubits[target].gates.append('Y')
        self.qubits[target].gatePos.append(self.qubits[target].earliestPos)
        self.qubits[target].earliestPos += 1

        return self

    # Pauli-Z gate
    def Z(self,target):

        # Create the Kronecker product matrix defining the gate operation. The target qubit gets a Z gate, all other qubits get identity.
        kronMatrix = np.array([1])
        for qIndex, qubit in enumerate(self.qubits):
            if qIndex == target:
                kronMatrix = np.kron(kronMatrix, Zmatrix)
            else:
                kronMatrix = np.kron(kronMatrix, np.eye(2))
        
        # Apply the gate to the circuit's state, updating the state.
        self.state = np.dot(kronMatrix, self.state)

        # Append the gate onto the running list of gates for the target qubit. Use the qubit's earliest position for the gate position, then increment the earliest position for potential future gates.
        self.qubits[target].gates.append('Z')
        self.qubits[target].gatePos.append(self.qubits[target].earliestPos)
        self.qubits[target].earliestPos += 1

        return self
    
    # Hadamard gate
    def H(self,target):

        # Create the Kronecker product matrix defining the gate operation. The target qubit gets an H gate, all other qubits get identity.
        kronMatrix = np.array([1])
        for qIndex, qubit in enumerate(self.qubits):
            if qIndex == target:
                kronMatrix = np.kron(kronMatrix, Hmatrix)
            else:
                kronMatrix = np.kron(kronMatrix, np.eye(2))
        
        # Apply the gate to the circuit's state, updating the state.
        self.state = np.dot(kronMatrix, self.state)

        # Append the gate onto the running list of gates for the target qubit. Use the qubit's earliest position for the gate position, then increment the earliest position for potential future gates.
        self.qubits[target].gates.append('H')
        self.qubits[target].gatePos.append(self.qubits[target].earliestPos)
        self.qubits[target].earliestPos += 1

        return self
    
    # Phase gate
    def S(self,target):

        # Create the Kronecker product matrix defining the gate operation. The target qubit gets an S gate, all other qubits get identity.
        kronMatrix = np.array([1])
        for qIndex, qubit in enumerate(self.qubits):
            if qIndex == target:
                kronMatrix = np.kron(kronMatrix, Smatrix)
            else:
                kronMatrix = np.kron(kronMatrix, np.eye(2))
        
        # Apply the gate to the circuit's state, updating the state.
        self.state = np.dot(kronMatrix, self.state)

        # Append the gate onto the running list of gates for the target qubit. Use the qubit's earliest position for the gate position, then increment the earliest position for potential future gates.
        self.qubits[target].gates.append('S')
        self.qubits[target].gatePos.append(self.qubits[target].earliestPos)
        self.qubits[target].earliestPos += 1

        return self
    
    # pi/8 gate
    def T(self,target):

        # Create the Kronecker product matrix defining the gate operation. The target qubit gets a T gate, all other qubits get identity.
        kronMatrix = np.array([1])
        for qIndex, qubit in enumerate(self.qubits):
            if qIndex == target:
                kronMatrix = np.kron(kronMatrix, Tmatrix)
            else:
                kronMatrix = np.kron(kronMatrix, np.eye(2))
        
        # Apply the gate to the circuit's state, updating the state.
        self.state = np.dot(kronMatrix, self.state)

        # Append the gate onto the running list of gates for the target qubit. Use the qubit's earliest position for the gate position, then increment the earliest position for potential future gates.
        self.qubits[target].gates.append('T')
        self.qubits[target].gatePos.append(self.qubits[target].earliestPos)
        self.qubits[target].earliestPos += 1

        return self
    
    ## TWO QUBIT GATES ##

    # Controlled-X gate
    def CX(self, control, target):

        # Create the Kronecker product matrix defining the gate operation. First create separate matrices for projecting the control qubit into the 0 or 1 state. For the projection into the 1 state, the target qubit gets an X gate. All other qubits, including the target in the projection into the 0 state for the control, get an identity. Then add the two matrices together.
        kron0Matrix = np.array([1])
        kron1Matrix = np.array([1])
        for qIndex, qubit in enumerate(self.qubits):
            if qIndex == control:
                kron0Matrix = np.kron(kron0Matrix, proj0matrix)
                kron1Matrix = np.kron(kron1Matrix, proj1matrix)
            elif qIndex == target:
                kron0Matrix = np.kron(kron0Matrix, np.eye(2))
                kron1Matrix = np.kron(kron1Matrix, Xmatrix)
            else:
                kron0Matrix = np.kron(kron0Matrix, np.eye(2))
                kron1Matrix = np.kron(kron1Matrix, np.eye(2))
        kronMatrix = kron0Matrix + kron1Matrix

        # Apply the gate to the circuit's state, updating the state.
        self.state = np.dot(kronMatrix, self.state)

        # Append the gate onto the running list of gates for the target qubit. Append the connection type onto the running list of connections for the control qubit and which qubit it is controlling (the target).
        self.qubits[target].gates.append('X')
        self.qubits[control].connections.append('C')
        self.qubits[control].connectTo.append(target)

        # Get the earliest possible gate position within the circuit for each qubit between the target and control (inclusive). The max of this list will be used for the gate position for both the target and control. Then increase the earliest position for all qubits between the target and control (inclusive).
        earliestPositions = [self.qubits[idx].earliestPos for idx in range(min(control, target), max(control, target)+1)]
        position = max(earliestPositions)
        self.qubits[target].gatePos.append(position)
        self.qubits[control].connectPos.append(position)
        for idx in range(min(control, target), max(control, target)+1):
            self.qubits[idx].earliestPos = position + 1

        return self
    
    # Controlled-Y gate
    def CY(self, control, target):

        # Create the Kronecker product matrix defining the gate operation. First create separate matrices for projecting the control qubit into the 0 or 1 state. For the projection into the 1 state, the target qubit gets a Y gate. All other qubits, including the target in the projection into the 0 state for the control, get an identity. Then add the two matrices together.
        kron0Matrix = np.array([1])
        kron1Matrix = np.array([1])
        for qIndex, qubit in enumerate(self.qubits):
            if qIndex == control:
                kron0Matrix = np.kron(kron0Matrix, proj0matrix)
                kron1Matrix = np.kron(kron1Matrix, proj1matrix)
            elif qIndex == target:
                kron0Matrix = np.kron(kron0Matrix, np.eye(2))
                kron1Matrix = np.kron(kron1Matrix, Ymatrix)
            else:
                kron0Matrix = np.kron(kron0Matrix, np.eye(2))
                kron1Matrix = np.kron(kron1Matrix, np.eye(2))
        kronMatrix = kron0Matrix + kron1Matrix

        # Apply the gate to the circuit's state, updating the state.
        self.state = np.dot(kronMatrix, self.state)

        # Append the gate onto the running list of gates for the target qubit. Append the connection type onto the running list of connections for the control qubit and which qubit it is controlling (the target).
        self.qubits[target].gates.append('Y')
        self.qubits[control].connections.append('C')
        self.qubits[control].connectTo.append(target)

        # Get the earliest possible gate position within the circuit for each qubit between the target and control (inclusive). The max of this list will be used for the gate position for both the target and control. Then increase the earliest position for all qubits between the target and control (inclusive).
        earliestPositions = [self.qubits[idx].earliestPos for idx in range(min(control, target), max(control, target)+1)]
        position = max(earliestPositions)
        self.qubits[target].gatePos.append(position)
        self.qubits[control].connectPos.append(position)
        for idx in range(min(control, target), max(control, target)+1):
            self.qubits[idx].earliestPos = position + 1

        return self
    
    # Controlled-Z gate
    def CZ(self, control, target):

        # Create the Kronecker product matrix defining the gate operation. First create separate matrices for projecting the control qubit into the 0 or 1 state. For the projection into the 1 state, the target qubit gets a Z gate. All other qubits, including the target in the projection into the 0 state for the control, get an identity. Then add the two matrices together.
        kron0Matrix = np.array([1])
        kron1Matrix = np.array([1])
        for qIndex, qubit in enumerate(self.qubits):
            if qIndex == control:
                kron0Matrix = np.kron(kron0Matrix, proj0matrix)
                kron1Matrix = np.kron(kron1Matrix, proj1matrix)
            elif qIndex == target:
                kron0Matrix = np.kron(kron0Matrix, np.eye(2))
                kron1Matrix = np.kron(kron1Matrix, Zmatrix)
            else:
                kron0Matrix = np.kron(kron0Matrix, np.eye(2))
                kron1Matrix = np.kron(kron1Matrix, np.eye(2))
        kronMatrix = kron0Matrix + kron1Matrix

        # Apply the gate to the circuit's state, updating the state.
        self.state = np.dot(kronMatrix, self.state)

        # Append the gate onto the running list of gates for the target qubit. Append the connection type onto the running list of connections for the control qubit and which qubit it is controlling (the target).
        self.qubits[target].gates.append('Z')
        self.qubits[control].connections.append('C')
        self.qubits[control].connectTo.append(target)

        # Get the earliest possible gate position within the circuit for each qubit between the target and control (inclusive). The max of this list will be used for the gate position for both the target and control. Then increase the earliest position for all qubits between the target and control (inclusive).
        earliestPositions = [self.qubits[idx].earliestPos for idx in range(min(control, target), max(control, target)+1)]
        position = max(earliestPositions)
        self.qubits[target].gatePos.append(position)
        self.qubits[control].connectPos.append(position)
        for idx in range(min(control, target), max(control, target)+1):
            self.qubits[idx].earliestPos = position + 1

        return self
    
    # SWAP gate
    def SWAP(self, target1, target2):

        # Create the Kronecker product matrix defining the gate operation. First create 3 separate matrices for the two target qubits to receive X, Y, and Z gates together. All other qubits get an identity. Then add the 3 matrices along with an identity and divide the sum by 2.
        kronXMatrix = np.array([1])
        kronYMatrix = np.array([1])
        kronZMatrix = np.array([1])
        for qIndex, qubit in enumerate(self.qubits):
            if qIndex == target1 or qIndex == target2:
                kronXMatrix = np.kron(kronXMatrix, Xmatrix)
                kronYMatrix = np.kron(kronYMatrix, Ymatrix)
                kronZMatrix = np.kron(kronZMatrix, Zmatrix)
            else:
                kronXMatrix = np.kron(kronXMatrix, np.eye(2))
                kronYMatrix = np.kron(kronYMatrix, np.eye(2))
                kronZMatrix = np.kron(kronZMatrix, np.eye(2))
        kronMatrix = 0.5*(np.eye(np.size(kronXMatrix, 0)) + kronXMatrix + kronYMatrix + kronZMatrix)

        # Apply the gate to the circuit's state, updating the state.
        self.state = np.dot(kronMatrix, self.state)

        # Append the gate onto the running list of gates for the target qubit (which we'll use target2 as for consistency with controlled gates). Append the connection type onto the running list of connections for the control qubit (target1) and which qubit it is controlling (target2).
        self.qubits[target2].gates.append('SWAP')
        self.qubits[target1].connections.append('SWAP')
        self.qubits[target1].connectTo.append(target2)

        # Get the earliest possible gate position within the circuit for each qubit between the targets (inclusive). The max of this list will be used for the gate position for both targets. Then increase the earliest position for all qubits between the targets (inclusive).
        earliestPositions = [self.qubits[idx].earliestPos for idx in range(min(target1, target2), max(target1, target2)+1)]
        position = max(earliestPositions)
        self.qubits[target2].gatePos.append(position)
        self.qubits[target1].connectPos.append(position)
        for idx in range(min(target1, target2), max(target1, target2)+1):
            self.qubits[idx].earliestPos = position + 1

        return self
    
    ## THREE QUBIT GATES ##

    # Toffoli gate
    def Toff(self, control1, control2, target):

        # Create the Kronecker product matrix defining the gate operation. First create 4 separate matrices for the 4 possible projection combos of the 2 control qubits. Apply the appropriate projection matrix for each control qubit for each of the 4 matrices. The target only gets an X gate when both controls are projected into the 1 state; the target gets an identity otherwise. All other qubits get an identity. Then add the 4 matrices together.
        kron00Matrix = np.array([1])
        kron01Matrix = np.array([1])
        kron10Matrix = np.array([1])
        kron11Matrix = np.array([1])
        for qIndex, qubit in enumerate(self.qubits):
            if qIndex == control1:
                kron00Matrix = np.kron(kron00Matrix, proj0matrix)
                kron01Matrix = np.kron(kron01Matrix, proj0matrix)
                kron10Matrix = np.kron(kron10Matrix, proj1matrix)
                kron11Matrix = np.kron(kron11Matrix, proj1matrix)
            elif qIndex == control2:
                kron00Matrix = np.kron(kron00Matrix, proj0matrix)
                kron01Matrix = np.kron(kron01Matrix, proj1matrix)
                kron10Matrix = np.kron(kron10Matrix, proj0matrix)
                kron11Matrix = np.kron(kron11Matrix, proj1matrix)
            elif qIndex == target:
                kron00Matrix = np.kron(kron00Matrix, np.eye(2))
                kron01Matrix = np.kron(kron01Matrix, np.eye(2))
                kron10Matrix = np.kron(kron10Matrix, np.eye(2))
                kron11Matrix = np.kron(kron11Matrix, Xmatrix)
            else:
                kron00Matrix = np.kron(kron00Matrix, np.eye(2))
                kron01Matrix = np.kron(kron01Matrix, np.eye(2))
                kron10Matrix = np.kron(kron10Matrix, np.eye(2))
                kron11Matrix = np.kron(kron11Matrix, np.eye(2))
        kronMatrix = kron00Matrix + kron01Matrix + kron10Matrix + kron11Matrix

        # Apply the gate to the circuit's state, updating the state.
        self.state = np.dot(kronMatrix, self.state)

        # Append the gate onto the running list of gates for the target qubit. Append the connection type onto the running lists of connections for the control qubits and which qubit they are controlling (the target).
        self.qubits[target].gates.append('X')
        self.qubits[control1].connections.append('C')
        self.qubits[control1].connectTo.append(target)
        self.qubits[control2].connections.append('C')
        self.qubits[control2].connectTo.append(target)

        # Get the earliest possible gate position within the circuit for each qubit between the target and controls (inclusive). The max of this list will be used for the gate position for both the target and controls. Then increase the earliest position for all qubits between the target and controls (inclusive).
        earliestPositions = [self.qubits[idx].earliestPos for idx in range(min(control1, control2, target), max(control1, control2, target)+1)]
        position = max(earliestPositions)
        self.qubits[target].gatePos.append(position)
        self.qubits[control1].connectPos.append(position)
        self.qubits[control2].connectPos.append(position)
        for idx in range(min(control1, control2, target), max(control1, control2, target)+1):
            self.qubits[idx].earliestPos = position + 1

        return self
    
    ## OTHER FUNCTIONS ##

    # Add a barrier to the circuit. The state vector does not change. A barrier is purely for visual purposes when displaying the circuit to divide the circuit into segments.
    def barrier(self):

        # Append a barrier to the first qubit. Use the last qubit as the connector to extend the barrier across the entire circuit. The max earliest position for all qubits is the position of the barrier. All qubits' earliest position is then updated to the position after the barrier.
        self.qubits[0].gates.append('B')
        self.qubits[-1].connections.append('B')
        self.qubits[-1].connectTo.append(0)
        earliestPosition = max([qubit.earliestPos for qubit in self.qubits])
        self.qubits[0].gatePos.append(earliestPosition)
        self.qubits[-1].connectPos.append(earliestPosition)
        for qubit in self.qubits:
            qubit.earliestPos = earliestPosition + 1

    # Measure a qubit and store the result in a classical bit. This is a measurement in the computational basis (projection into the 0 or 1 state).
    def measure(self, target, output):

        # Create 2 Kronecker product matrices for the projection of the target qubit into the 0 or 1 state. All other qubits get an identity.
        kron0Matrix = np.array([1])
        kron1Matrix = np.array([1])
        for qIndex, qubit in enumerate(self.qubits):
            if qIndex == target:
                kron0Matrix = np.kron(kron0Matrix, proj0matrix)
                kron1Matrix = np.kron(kron1Matrix, proj1matrix)
            else:
                kron0Matrix = np.kron(kron0Matrix, np.eye(2))
                kron1Matrix = np.kron(kron1Matrix, np.eye(2))
        
        # For each of the 0 and 1 state projection matrices, apply the projection to the circuit's current state to get the resulting state. Transpose the circuit's current state (without the projection) and apply it to the resulting state (with the projection) to get the probability of the measurement outcome.
        state0 = np.dot(kron0Matrix, self.state)
        prob0 = np.dot(self.state.T, state0)[0][0]
        state1 = np.dot(kron1Matrix, self.state)
        prob1 = np.dot(self.state.T, state1)[0][0]
        
        # Generate a random number between 0 and 1. If it is less than the probability of the target qubit being in the 0 state, set the classical bit to 0 and update the circuit's state with the projection-into-0 state from above (normalized with the square root of the probability of measuring 0). Otherwise, set the classical bit to 1 and update the circuit's state with the projection-into-1 state (normalized).
        if np.random.rand(1) < prob0:
            self.cbits[output].state = 0
            self.state = state0 / np.sqrt(prob0)
        else:
            self.cbits[output].state = 1
            self.state = state1 / np.sqrt(prob1)

        # Append the gate onto the running list of gates for the target qubit. Append the connection type onto the running list of connections for the control qubit and which qubit it is controlling (the target).
        self.qubits[target].gates.append('M')
        self.cbits[output].connections.append('O')
        self.cbits[output].connectTo.append(target)

        # Get the earliest possible gate position within the circuit for each qubit between the target and control (inclusive). The max of this list will be used for the gate position for both the target and control. Then increase the earliest position for all qubits between the target and control (inclusive).
        earliestPositions = [qubit.earliestPos for qubit in self.qubits[target:]]
        for cbit in self.cbits[:output]:
            earliestPositions.append(cbit.earliestPos)
        position = max(earliestPositions)
        self.qubits[target].gatePos.append(position)
        self.cbits[output].connectPos.append(position)
        for qubit in self.qubits[target:]:
            qubit.earliestPos = position + 1
        for cbit in self.cbits[:output]:
            cbit.earliestPos = position + 1

        return self
    
    # Assign the gate label, box, and connection property to be used for displaying the circuit
    def format_gate(self, gate):
        # SWAP gates
        if gate == 'SWAP':
            gateLabel = 'x'
            size = 30
            bbox=dict(boxstyle='square', pad=0, facecolor='none', edgecolor='none')
            arrowprops=dict(arrowstyle="-", edgecolor='black', linewidth=2)
        # Controls in controlled-gates
        elif gate == 'C':
            gateLabel = ' '
            size = 7
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
            size = 20
            bbox=dict(boxstyle='square', facecolor='gray', edgecolor='none')
            arrowprops=dict(arrowstyle="-", edgecolor='grey', linewidth=18)
        else:
            gateLabel = gate
            size = 20
            bbox=dict(boxstyle='square', facecolor='white')
            arrowprops=dict()
        
        return [gateLabel, size, bbox, arrowprops]

    # Create a figure of the circuit.
    def display_circuit(self):

        # Create the figure and set some style parameters.
        fig = plt.figure()
        ax = fig.add_subplot(111)
        circuitLength = max([qubit.gatePos[-1] for qubit in self.qubits])
        numQubits = len(self.qubits)
        numCbits = len(self.cbits)
        bitLabelPosition = 0.5
        ax.set(xlim=(0, circuitLength+1), ylim=(-1*(numQubits+numCbits), 1))
        ax.set_axis_off()

        # Display each classical bit. These are displayed first for proper layer ordering when displaying connections between qubits and classical bits.
        offset = 0.05
        for Bidx, cbit in enumerate(self.cbits):
            # Display the classical bit labels and a double line to represent their wires. "offset" creates a small spacing between the two plotted lines for each bit to give the double line visual.
            ax.annotate('$C_%s$'%Bidx, xy=(bitLabelPosition, -1*(Bidx+numQubits)), size=20, va='center', ha='center', bbox=dict(boxstyle='square', facecolor='white', edgecolor='none'))
            ax.plot([bitLabelPosition, circuitLength], [-1*(Bidx+numQubits)+offset, -1*(Bidx+numQubits)+offset], color='black')
            ax.plot([bitLabelPosition, circuitLength], [-1*(Bidx+numQubits)-offset, -1*(Bidx+numQubits)-offset], color='black')

            # Display each connection using the properties from format_gate
            for Cidx, connection in enumerate(cbit.connections):
                [connectLabel, size, bbox, arrowprops] = self.format_gate(connection)
                ax.annotate(connectLabel, xy=(cbit.connectPos[Cidx], -1*cbit.connectTo[Cidx]), xytext=(cbit.connectPos[Cidx], -1*(Bidx+numQubits)), size=size, va='center', ha='center', bbox=bbox, arrowprops=arrowprops)

        # Display each qubit
        for Qidx, qubit in enumerate(self.qubits):
            # Display the qubit labels and a horizontal line to represent the wire for each qubit's circuit.
            ax.annotate('$Q_%s$'%Qidx, xy=(bitLabelPosition, -1*Qidx), size=20, va='center', ha='center', bbox=dict(boxstyle='square', facecolor='white', edgecolor='none'))
            ax.plot([bitLabelPosition, circuitLength], [-1*Qidx, -1*Qidx], color='black')

        # Display each qubit connection using the properties from format_gate. These are displayed next for proper layer ordering of connections between qubits.
        for Qidx, qubit in enumerate(self.qubits):
            for Cidx, connection in enumerate(qubit.connections):
                [connectLabel, size, bbox, arrowprops] = self.format_gate(connection)
                ax.annotate(connectLabel, xy=(qubit.connectPos[Cidx], -1*qubit.connectTo[Cidx]), xytext=(qubit.connectPos[Cidx], -1*Qidx), size=size, va='center', ha='center', bbox=bbox, arrowprops=arrowprops)

        # Display each qubit gate using the properties from format_gate.
        for Qidx, qubit in enumerate(self.qubits):
            for Gidx, gate in enumerate(qubit.gates):
                [gateLabel, size, bbox, arrowprops] = self.format_gate(gate)
                ax.annotate(gateLabel, xy=(qubit.gatePos[Gidx], -1*Qidx), size=size, va='center', ha='center', bbox=bbox)

        plt.show()

        return

def main():
    numQubits = 3
    test = Circuit(numQubits, numQubits)
    test.H(0)
    test.barrier()
    test.X(1)
    test.CZ(0, 2)
    test.S(1)
    test.T(0)
    test.barrier()
    test.SWAP(1, 2)
    test.Y(2)
    test.barrier()
    test.barrier()
    test.Z(1)
    test.measure(0, 0)
    test.measure(1, 1)
    test.measure(2, 2)
    test.display_circuit()

if __name__ == "__main__":
    main()