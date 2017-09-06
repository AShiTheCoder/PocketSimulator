# PocketSimulator
PocketSimulator is a space-efficient simulator for arbitrary quantum circuits. Many existing quantum computation simulators operate on the state-vector representation of quantum states, using time and space exponential to the # of qubits involved. In contrast, PocketSimulator adapts a path-integral-esque approach to simulation by iterating over different computation paths successively, achieving *linear* space usage in the # of qubits + the # of gates. 
## About
Shi has [shown](https://arxiv.org/abs/quant-ph/0205115) that Hadamard and Toffoli gates together are sufficient for universal quantum computation (ie. any quantum computation can be computed with a circuit composed of only of Hadamard/Toffoli gates).
Because Toffoli gates are universal for classical computation and therefore...
## Getting Started
PocketSimulator is contained entirely in PathIntegral.cpp; experimental implementations of other simulation methods are located in LabTests.cpp. Any compiler or environment capable of c++ program execution is suitable for running PocketSimulator.
## Running Circuits
PocketSimulator takes three main arguments for custom simulation:
- **n**: # of qubits
- **length**: # of "splitting" gates in the circuit
- **gates.txt**: text file encoding the computation to be simulated
n and length can be modified in the main() method of PathIntegral.cpp; gates.txt can be edited directly to input a desired algorithm or sequence. 

Currently, the following gates are supported in PocketSimulator:

Gate | Abbreviation | Arguments | Example
---|---|---|---
Hadamard | h | target qubit (h \_) | h 0
Toffoli | t | 2 controls and a target qubit (t c1 c2 \_) | t 0 1 2
U (2^n phase) | U | ???? | ????
(\*) Pauli X | x | target qubit (x \_) | x 0
(\*) Pauli Y | y | target qubit (y \_) | y 0
(\*) Pauli Z | z | target qubit (z \_) | z 0
