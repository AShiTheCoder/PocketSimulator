# PocketSimulator
A space-efficient simulator for arbitrary quantum circuits
## About
Many existing quantum computation simulators operate on the state-vector representation of quantum states, using time and space exponential to the # of qubits involved. In contrast, PocketSimulator approaches simulation with an algorithm reminiscent of Feynman's path integral formulation of quantum mechanics. It successively iterates over every computation path, achieving *linear* space usage in the # of qubits + the # of gates.
## Getting Started
PocketSimulator is contained entirely in PathIntegral.cpp; experimental implementations of other simulation methods are located in LabTests.cpp. Any compiler or environment capable of c++ program execution is suitable for running PocketSimulator.
## Running Circuits
PocketSimulator takes three main arguments for custom simulation:
- **n**: # of qubits
- **length**: # of "changing" gates in the circuit (all gates excluding those which purely add a relative phase)
- **gates.txt**: text file encoding the computation to be simulated

n and length can be modified in the main() method of PathIntegral.cpp; gates.txt can be edited directly to input a desired algorithm or sequence.

Currently, the following gates are/will be supported in PocketSimulator:

Gate | Abbreviation | Arguments | Example
---|---|---|---
Hadamard | h | target qubit (h \_) | h 0
Toffoli | t | 2 controls and a target qubit (t c1 c2 \_) | t 0 1 2
U (2^n phase) | U | ???? | ????
(\*) Pauli X | x | target qubit (x \_) | x 0
(\*) Pauli Y | y | target qubit (y \_) | y 0
(\*) Pauli Z | z | target qubit (z \_) | z 0
(\*) General Phase | R | numerator denominator target (R n d \_) | R 3 7 0 [phase of e^(i\*pi\*6/7)]

(\*) coming soon
