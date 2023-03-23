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

# Creates an instance of a quantum circuit with a provided number of quantum bits and classical bits and allows the user
# to apply qubit gates to the circuit
class Circuit:

    # Create the provided number of qubits and classical bits upon instance initialization
    def __init__(self, numQubits, numCbits):

        # Create a list of qubits, each intialized in the |0> state, [1, 0].
        qubit_list = np.zeros([numQubits, 2], dtype=np.complex_)
        qubit_list[:, 0] = 1
        self.qubits = qubit_list
        
        # Form the state of all the qubits in the circuit.
        self.state = np.array([1])
        for qubit in self.qubits:
            self.state = np.tensordot(self.state, qubit, axes=0).reshape(len(self.state)*2, 1)
        
        # Create a list of classical bits, each initialized in the 0 state.
        self.cbits = np.zeros([numCbits, 1], dtype=np.int8)

        # Create a blank list that will store all gates when they are applied to the circuit.
        self.gates = []

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

        # Append the gate onto the running list of gates in the circuit.
        self.gates.append(['X', target])

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

        # Append the gate onto the running list of gates in the circuit.
        self.gates.append(['Y', target])

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

        # Append the gate onto the running list of gates in the circuit.
        self.gates.append(['Z', target])

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

        # Append the gate onto the running list of gates in the circuit.
        self.gates.append(['H', target])

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

        # Append the gate onto the running list of gates in the circuit.
        self.gates.append(['S', target])

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

        # Append the gate onto the running list of gates in the circuit.
        self.gates.append(['T', target])

        return self
    
    ## TWO QUBIT GATES ##

    # Controlled-X gate
    def CX(self, target, control):

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

        # Append the gate onto the running list of gates in the circuit.
        self.gates.append(['X', target, 'C', control])

        return self
    
    # Controlled-Y gate
    def CY(self, target, control):

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

        # Append the gate onto the running list of gates in the circuit.
        self.gates.append(['Y', target, 'C', control])

        return self
    
    # Controlled-Z gate
    def CZ(self, target, control):

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

        # Append the gate onto the running list of gates in the circuit.
        self.gates.append(['Z', target, 'C', control])

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

        # Append the gate onto the running list of gates in the circuit.
        self.gates.append(['SWAP', target1, 'SWAP', target2])

        return self
    
    ## THREE QUBIT GATES ##

    # Toffoli gate
    def Toff(self, target, control1, control2):

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

        # Append the gate onto the running list of gates in the circuit.
        self.gates.append(['X', target, 'C', control1, 'C', control2])

        return self
    
    ## OTHER FUNCTIONS ##

    # Measure a qubit and store the result in a classical bit. This is a measurement in the computational basis (projection into the 0 or 1 state).
    def measure(self, target, cbit):

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
            self.cbits[cbit] = 0
            self.state = state0 / np.sqrt(prob0)
        else:
            self.cbits[cbit] = 1
            self.state = state1 / np.sqrt(prob1)

        # Append the gate onto the running list of gates in the circuit.
        self.gates.append(['M', target, cbit])

        return self
    
    # Create a figure of the circuit.
    def display_circuit(self):

        # Create the figure and set some style parameters.
        fig = plt.figure()
        ax = fig.add_subplot(111)
        plt.rcParams.update({'font.size': 20})
        numQubits = len(self.qubits)
        numCbits = len(self.cbits)
        ax.set(xlim=(0, 5), ylim=(-1*(numQubits+numCbits), 1))
        ax.set_axis_off()

        # Display the qubit labels and a horizontal line to represent the wire for each qubit's circuit.
        numGates = len(self.gates)
        bitLabelPosition = 0.5
        for qIndex, qubit in enumerate(self.qubits):
            ax.annotate('$Q_%s$'%qIndex, xy=(bitLabelPosition, -1*qIndex), va='center', ha='center', bbox=dict(boxstyle='square', facecolor='white', edgecolor='none'))
            ax.plot([bitLabelPosition, numGates], [-1*qIndex, -1*qIndex], color='black')
        
        # Display the classical bit labels and a double line to represent their wires. "offset" creates a small spacing between the two plotted lines for each bit to give the double line visual.
        offset = 0.05
        for cIndex, cbit in enumerate(self.cbits):
            ax.annotate('$C_%s$'%cIndex, xy=(bitLabelPosition, -1*(cIndex+numQubits)), va='center', ha='center', bbox=dict(boxstyle='square', facecolor='white', edgecolor='none'))
            ax.plot([bitLabelPosition, numGates], [-1*(cIndex+numQubits)+offset, -1*(cIndex+numQubits)+offset], color='black')
            ax.plot([bitLabelPosition, numGates], [-1*(cIndex+numQubits)-offset, -1*(cIndex+numQubits)-offset], color='black')

        # Display each gate applied to the circuit.
        for gIndex, gate in enumerate(self.gates):
            # For plotting, the negative of the qubit number will be used for the gate's y value so that qubit 0 can be at the top with additional qubits descending.
            target = -1*gate[1]
            gateType = gate[0]
            gateNumber = gIndex+1

            # For multi-qubit gates, first plot the control gates. This puts the line connecting the control to the target gate behind the target gate (only needed for visual appeal).
            if len(gate) > 2:
                for index in range(int(len(gate) / 2)):
                    # Each gate function returns the target first, so skip the first gate element.
                    if index == 0:
                        continue
                    # Style control gates with a filled dot and connect it with a line to the target gate position.
                    if gate[index*2] == 'C':
                        control = -1*gate[index*2+1]
                        ax.annotate('•', xy=(gateNumber, target), xytext=(gateNumber, control), va='center', ha='center', bbox=dict(boxstyle='Circle, pad=0', facecolor='black'), arrowprops=dict(arrowstyle="-"))
                    # Style SWAP gates with an X and a line to the swapped qubit's position.
                    elif gate[index*2] == 'SWAP':
                        control = -1*gate[index*2+1]
                        ax.annotate('X', xy=(gateNumber, target), xytext=(gateNumber, control), va='center', ha='center', arrowprops=dict(arrowstyle="-"))
            
            # If the target gate is a SWAP, style it with an X.
            if gateType == 'SWAP':
                ax.annotate('X', xy=(gateNumber, target), va='center', ha='center')
            # All other target gates are styled with their gate label in a box.
            else:
                ax.annotate(gateType, xy=(gateNumber, target), va='center', ha='center', bbox=dict(boxstyle='square', fc='white'))

        plt.show()

        return

def main():
    numQubits = 3
    test = Circuit(numQubits, numQubits)
    test.X(0)
    test.H(1)
    test.SWAP(0, 1)
    test.Toff(1, 0, 2)
    test.display_circuit()

if __name__ == "__main__":
    main()