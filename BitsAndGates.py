## This file allows the user to create a quantum circuit with a provided number of qubits and classical bits.
## Common qubit gates to operate on the quantum circuit are defined.

import numpy as np

# Creates an instance of a quantum circuit with a provided number of quantum bits and classical bits and allows the user
# to apply qubit gates to the circuit
class Circuit:

    # Create the provided number of qubits and classical bits upon instance initialization
    def __init__(self,num_qubits,num_cbits):
        qubit_list = np.zeros([2,num_qubits],dtype=np.complex_)
        qubit_list[0,:] = 1
        self.qubits = qubit_list
        self.cbits = np.zeros([1,num_cbits],dtype=np.int8)
        self.gates = []

    ## SINGLE QUBIT GATES ##

    # Pauli-X gate
    def X(self,qubit_index):
        self.qubits[:,qubit_index] = np.dot(np.array([[0,1],
                                                      [1,0]]),self.qubits[:,qubit_index])
        self.gates.append(['X',qubit_index])
        return self

    # Pauli-Y gate
    def Y(self,qubit_index):
        self.qubits[:,qubit_index] = np.dot(np.array([[0,-1j],
                                                      [1j,0]]),self.qubits[:,qubit_index])
        self.gates.append(['Y',qubit_index])
        return self

    # Pauli-Z gate
    def Z(self,qubit_index):
        self.qubits[:,qubit_index] = np.dot(np.array([[1,0],
                                                      [0,-1]]),self.qubits[:,qubit_index])
        self.gates.append(['Z',qubit_index])
        return self
    
    # Hadamard gate
    def H(self,qubit_index):
        self.qubits[:,qubit_index] = np.dot(np.array([[1,1],
                                                      [1,-1]])/np.sqrt(2),self.qubits[:,qubit_index])
        self.gates.append(['H',qubit_index])
        return self
    
    # Phase gate
    def S(self,qubit_index):
        self.qubits[:,qubit_index] = np.dot(np.array([[1,0],
                                                      [0,1j]]),self.qubits[:,qubit_index])
        self.gates.append(['S',qubit_index])
        return self
    
    # pi/8 gate
    def T(self,qubit_index):
        self.qubits[:,qubit_index] = np.dot(np.array([[1,0],
                                                      [0,np.exp(1j*np.pi/4)]]),self.qubits[:,qubit_index])
        self.gates.append(['T',qubit_index])
        return self
    
    ## TWO QUBIT GATES ##

    # Controlled-X gate
    def CX(self,c_qubit,t_qubit):
        tensor = np.tensordot(self.qubits[:,c_qubit],self.qubits[:,t_qubit],axes=0).reshape(4,1)
        CX_tensor = np.dot(np.array([[1,0,0,0],
                                     [0,1,0,0],
                                     [0,0,0,1],
                                     [0,0,1,0]]),tensor)
        if self.qubits[0,c_qubit] != 0: self.qubits[:,t_qubit] = CX_tensor[0:2,0]/self.qubits[0,c_qubit]
        else: self.qubits[:,t_qubit] = CX_tensor[2:4,0]/self.qubits[1,c_qubit]
        self.gates.append(['CX',c_qubit,t_qubit])
        return self
    
    # Controlled-Y gate
    def CY(self,c_qubit,t_qubit):
        tensor = np.tensordot(self.qubits[:,c_qubit],self.qubits[:,t_qubit],axes=0).reshape(4,1)
        CY_tensor = np.dot(np.array([[1,0,0,0],
                                     [0,1,0,0],
                                     [0,0,0,-1j],
                                     [0,0,1j,0]]),tensor)
        if self.qubits[0,c_qubit] != 0: self.qubits[:,t_qubit] = CY_tensor[0:2,0]/self.qubits[0,c_qubit]
        else: self.qubits[:,t_qubit] = CY_tensor[2:4,0]/self.qubits[1,c_qubit]
        self.gates.append(['CY',c_qubit,t_qubit])
        return self
    
    # Controlled-Z gate
    def CZ(self,c_qubit,t_qubit):
        tensor = np.tensordot(self.qubits[:,c_qubit],self.qubits[:,t_qubit],axes=0).reshape(4,1)
        CZ_tensor = np.dot(np.array([[1,0,0,0],
                                     [0,1,0,0],
                                     [0,0,1,0],
                                     [0,0,0,-1]]),tensor)
        if self.qubits[0,c_qubit] != 0: self.qubits[:,t_qubit] = CZ_tensor[0:2,0]/self.qubits[0,c_qubit]
        else: self.qubits[:,t_qubit] = CZ_tensor[2:4,0]/self.qubits[1,c_qubit]
        self.gates.append(['CZ',c_qubit,t_qubit])
        return self
    
    # SWAP gate
    def SWAP(self,qubit1,qubit2):
        coeff = [self.qubits[0,qubit1],self.qubits[1,qubit1],self.qubits[0,qubit2],self.qubits[1,qubit2]]
        self.qubits[0,qubit1], self.qubits[1,qubit1], self.qubits[0,qubit2], self.qubits[1,qubit2] = coeff[2], coeff[3], coeff[0], coeff[1]
        self.gates.append(['SWAP',qubit1,qubit2])
        return self
    
    ## THREE QUBIT GATES ##

    # Toffoli gate
    def Toff(self,c_qubit1,c_qubit2,t_qubit):
        tensor = np.tensordot(self.qubits[:,c_qubit1],self.qubits[:,c_qubit2],axes=0).reshape(4,1)
        tensor = np.tensordot(tensor,self.qubits[:,t_qubit],axes=0).reshape(8,1)
        Toff_tensor = (np.dot(np.array([[1,0,0,0,0,0,0,0],
                                        [0,1,0,0,0,0,0,0],
                                        [0,0,1,0,0,0,0,0],
                                        [0,0,0,1,0,0,0,0],
                                        [0,0,0,0,1,0,0,0],
                                        [0,0,0,0,0,1,0,0],
                                        [0,0,0,0,0,0,0,1],
                                        [0,0,0,0,0,0,1,0]]),tensor))
        if self.qubits[0,c_qubit1] != 0:
            new_tensor = Toff_tensor[0:4,0]/self.qubits[0,c_qubit1]
        else:
            new_tensor = Toff_tensor[4:8,0]/self.qubits[1,c_qubit1]
        if self.qubits[0,c_qubit2] != 0:
            self.qubits[:,t_qubit] = new_tensor[0:2]/self.qubits[0,c_qubit2]
        else:
            self.qubits[:,t_qubit] = new_tensor[2:4]/self.qubits[1,c_qubit2]
        self.gates.append(['Toff',c_qubit1,c_qubit2,t_qubit])
        return self
    
    ## OTHER FUNCTIONS ##

    # Measure a qubit and store the result in a classical bit.
    def measure(self, qubit, cbit):
        qubit_probabilities = self.qubits[:, qubit]*np.conjugate(self.qubits[:, qubit])
        if np.random.rand(1) < qubit_probabilities[0]:
            self.cbits[:, cbit] = 0
            self.qubits[:, qubit] = [1, 0]
        else:
            self.cbits[:, cbit] = 1
            self.qubits[:, qubit] = [0, 1]
        self.gates.append(['measure',qubit])
        return self

def main():
    test = Circuit(3,3)
    test.X(2)
    test.X(0)
    test.H(0)
    print(test.qubits)
    test.measure(0,0)
    print(test.qubits)
    print(test.cbits)
    print(test.gates)

if __name__ == "__main__":
    main()