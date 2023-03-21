## This file allows the user to create a quantum circuit with a provided number of qubits and classical bits.
## Common qubit gates to operate on the quantum circuit are defined.

import numpy as np

# Creates an instance of a quantum circuit with a provided number of quantum bits and classical bits. The user can display the state of the bits and/or the circuit of all bits and gates.
class Circuit:
    # Create the provided number of qubits and classical bits upon instance initialization
    def __init__(self,numQubits,numCbits):
        self.numQubits = numQubits
        self.qubits = []
        for qubit in range(numQubits):
            self.qubits.append(Qubit())

        self.numCbits = numCbits
        self.cbits = []
        for cbit in range(numCbits):
            self.cbits.append(Cbit())
    
    # Print the state of all qubits in the command line.
    def display_states(self):
        for qubit in range(self.numQubits):
            print('Q%i:'%qubit, self.qubits[qubit].state)

        for cbit in range(self.numCbits):
            print('C%i:'%cbit, self.cbits[cbit].state)
        return

    # Print the circuit of all applied gates in the command line.
    def display_circuit(self):
        for qubit in self.qubits:
            for gate in qubit.gates:
                print('┌───┐')
                print('│ %s │'%gate)
                print('└───┘')
        return

# Creates an instance of a classical bit. The state is initialized to 0.
class Cbit:
    def __init__(self):
        self.state = 0

# Creates an instance of a quantum bit, or qubit. The state is initialized to |0>, [1, 0]. The user can apply gates to the qubit, which will be tracked.
class Qubit:

    def __init__(self):
        self.state = [1, 0]
        self.gates = []
        return

    ## SINGLE QUBIT GATES ##

    # Pauli-X gate
    def X(self):
        self.state = np.dot(np.array([[0,1],
                                      [1,0]]),self.state)
        self.gates.append('X')
        return self

    # Pauli-Y gate
    def Y(self):
        self.state = np.dot(np.array([[0,-1j],
                                      [1j,0]]),self.state)
        self.gates.append('Y')
        return self

    # Pauli-Z gate
    def Z(self):
        self.state = np.dot(np.array([[1,0],
                                      [0,-1]]),self.state)
        self.gates.append('Z')
        return self
    
    # Hadamard gate
    def H(self):
        self.state = np.dot(np.array([[1,1],
                                      [1,-1]])/np.sqrt(2),self.state)
        self.gates.append('H')
        return self
    
    # Phase gate
    def S(self):
        self.state = np.dot(np.array([[1,0],
                                      [0,1j]]),self.state)
        self.gates.append('S')
        return self
    
    # pi/8 gate
    def T(self):
        self.state = np.dot(np.array([[1,0],
                                      [0,np.exp(1j*np.pi/4)]]),self.state)
        self.gates.append('T')
        return self
    
    # ## TWO QUBIT GATES ##

    # # Controlled-X gate
    # def CX(self,c_qubit,t_qubit):
    #     tensor = np.tensordot(self.qubits[:,c_qubit],self.qubits[:,t_qubit],axes=0).reshape(4,1)
    #     CX_tensor = np.dot(np.array([[1,0,0,0],
    #                                  [0,1,0,0],
    #                                  [0,0,0,1],
    #                                  [0,0,1,0]]),tensor)
    #     if self.qubits[0,c_qubit] != 0: self.qubits[:,t_qubit] = CX_tensor[0:2,0]/self.qubits[0,c_qubit]
    #     else: self.qubits[:,t_qubit] = CX_tensor[2:4,0]/self.qubits[1,c_qubit]
    #     self.gates.append(['CX',c_qubit,t_qubit])
    #     return self
    
    # # Controlled-Y gate
    # def CY(self,c_qubit,t_qubit):
    #     tensor = np.tensordot(self.qubits[:,c_qubit],self.qubits[:,t_qubit],axes=0).reshape(4,1)
    #     CY_tensor = np.dot(np.array([[1,0,0,0],
    #                                  [0,1,0,0],
    #                                  [0,0,0,-1j],
    #                                  [0,0,1j,0]]),tensor)
    #     if self.qubits[0,c_qubit] != 0: self.qubits[:,t_qubit] = CY_tensor[0:2,0]/self.qubits[0,c_qubit]
    #     else: self.qubits[:,t_qubit] = CY_tensor[2:4,0]/self.qubits[1,c_qubit]
    #     self.gates.append(['CY',c_qubit,t_qubit])
    #     return self
    
    # # Controlled-Z gate
    # def CZ(self,c_qubit,t_qubit):
    #     tensor = np.tensordot(self.qubits[:,c_qubit],self.qubits[:,t_qubit],axes=0).reshape(4,1)
    #     CZ_tensor = np.dot(np.array([[1,0,0,0],
    #                                  [0,1,0,0],
    #                                  [0,0,1,0],
    #                                  [0,0,0,-1]]),tensor)
    #     if self.qubits[0,c_qubit] != 0: self.qubits[:,t_qubit] = CZ_tensor[0:2,0]/self.qubits[0,c_qubit]
    #     else: self.qubits[:,t_qubit] = CZ_tensor[2:4,0]/self.qubits[1,c_qubit]
    #     self.gates.append(['CZ',c_qubit,t_qubit])
    #     return self
    
    # # SWAP gate
    # def SWAP(self,qubit1,qubit2):
    #     coeff = [self.qubits[0,qubit1],self.qubits[1,qubit1],self.qubits[0,qubit2],self.qubits[1,qubit2]]
    #     self.qubits[0,qubit1], self.qubits[1,qubit1], self.qubits[0,qubit2], self.qubits[1,qubit2] = coeff[2], coeff[3], coeff[0], coeff[1]
    #     self.gates.append(['SWAP',qubit1,qubit2])
    #     return self
    
    # ## THREE QUBIT GATES ##

    # # Toffoli gate
    # def Toff(self,c_qubit1,c_qubit2,t_qubit):
    #     tensor = np.tensordot(self.qubits[:,c_qubit1],self.qubits[:,c_qubit2],axes=0).reshape(4,1)
    #     tensor = np.tensordot(tensor,self.qubits[:,t_qubit],axes=0).reshape(8,1)
    #     Toff_tensor = (np.dot(np.array([[1,0,0,0,0,0,0,0],
    #                                     [0,1,0,0,0,0,0,0],
    #                                     [0,0,1,0,0,0,0,0],
    #                                     [0,0,0,1,0,0,0,0],
    #                                     [0,0,0,0,1,0,0,0],
    #                                     [0,0,0,0,0,1,0,0],
    #                                     [0,0,0,0,0,0,0,1],
    #                                     [0,0,0,0,0,0,1,0]]),tensor))
    #     if self.qubits[0,c_qubit1] != 0:
    #         new_tensor = Toff_tensor[0:4,0]/self.qubits[0,c_qubit1]
    #     else:
    #         new_tensor = Toff_tensor[4:8,0]/self.qubits[1,c_qubit1]
    #     if self.qubits[0,c_qubit2] != 0:
    #         self.qubits[:,t_qubit] = new_tensor[0:2]/self.qubits[0,c_qubit2]
    #     else:
    #         self.qubits[:,t_qubit] = new_tensor[2:4]/self.qubits[1,c_qubit2]
    #     self.gates.append(['Toff',c_qubit1,c_qubit2,t_qubit])
    #     return self
    
    ## OTHER FUNCTIONS ##

    # Measure a qubit and store the result in a classical bit.
    def measure(self, cbit):
        qubit_probabilities = self.state*np.conjugate(self.state)
        if np.random.rand(1) < qubit_probabilities[0]:
            cbit = 0
            self.state = [1, 0]
        else:
            cbit = 1
            self.state = [0, 1]
        self.gates.append(['M'])
        return self

def main():
    test = Circuit(3,3)
    test.qubits[0].X()
    test.qubits[0].X()
    test.display_circuit()
    # test.measure(0,0)

if __name__ == "__main__":
    main()