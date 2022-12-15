## This file allows the user to create a quantum circuit with a provided number of qubits and classical bits.
## Common qubit gates to operate on the quantum circuit are defined.

import numpy as np

# Creates an instance of a quantum circuit with a provided number of quantum bits and classical bits and allows the user
# to apply qubit gates to the circuit
class Circuit:

    # Create the provided number of qubits and classical bits upon instance initialization
    def __init__(self,num_qubits,num_cbits):
        qubit_list = np.zeros([num_qubits,2],dtype=np.complex_)
        qubit_list[:,0] = 1
        self.qubits = qubit_list
        self.cbits = np.zeros([num_cbits,1],dtype=np.int8)

    # Pauli-X gate
    def X(self,qubit_index):
        self.qubits[qubit_index,:] = np.dot(np.array([[0,1],[1,0]]),self.qubits[qubit_index])
        return self

    # Pauli-Y gate
    def Y(self,qubit_index):
        self.qubits[qubit_index,:] = np.dot(np.array([[0,-1j],[1j,0]]),self.qubits[qubit_index])
        return self

    # Pauli-Z gate
    def Z(self,qubit_index):
        self.qubits[qubit_index,:] = np.dot(np.array([[1,0],[0,-1]]),self.qubits[qubit_index])
        return self
    
    def H(self,qubit_index):
        self.qubits[qubit_index,:] = np.dot(np.array([[1,1],[1,1]])/np.sqrt(2),self.qubits[qubit_index])
        return self