# PocketSimulator
A space-efficient simulator for arbitrary quantum circuits
## About
Many existing quantum computation simulators operate on the state-vector representation of quantum states, using time and space exponential to the # of qubits involved. In contrast, PocketSimulator approaches simulation with an algorithm reminiscent of Feynman's path integral formulation of quantum mechanics. It successively iterates over every computation path, achieving *linear* space usage in the # of qubits + the # of gates.
## Getting Started
The main algorithm for which PocketSimulator was developed lies in pathIntegral.cpp; experimental implementations of other simulation methods are located in their respective files. Various settings for testing these algorithms on pre-made or user-inputted circuits are modified in main.cpp. Any compiler or environment capable of c++ program execution is suitable for running PocketSimulator.
## Running Circuits
### Parameters
PocketSimulator takes several arguments for custom simulation:
- **N**: # of qubits
- **length**: # of "changing" gates in the circuit (all gates excluding those which purely add a relative phase)
- **gates.txt**: text file encoding the computation to be simulated
- **startState** and **endState**

gates.txt can be edited directly to input a desired algorithm (sequence of quantum gates), while the remaining parameters are inputted directly into the main() method of PathIntegral.cpp. DID THE PROJECT MOVE WORK???
### Execution
Upon execution of algorithm U, PocketSimulator will return a probability amplitude <endState|U|startState>, as well as memory and time usage details of the execution as returned by the system method `getrusage()` ([documentation here](http://pubs.opengroup.org/onlinepubs/009695399/functions/getrusage.html)).
### Gates
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

