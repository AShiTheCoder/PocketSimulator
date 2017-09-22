# PocketSimulator
A space-efficient simulator for arbitrary quantum circuits
## About
Many existing quantum computation simulators operate on the state-vector representation of quantum states, using time and space exponential to the # of qubits involved. In contrast, PocketSimulator approaches simulation with an algorithm reminiscent of Feynman's path integral formulation of quantum mechanics. It successively iterates over every computation path, achieving *linear* space usage in the # of qubits + the # of gates.
## Getting Started
The main algorithm for which PocketSimulator was developed lies in pathIntegral.cpp; experimental implementations of other simulation methods are located in their respective files. Various settings for testing these algorithms on pre-made or user-inputted circuits are modified in ```main.cpp```. Any compiler or environment capable of C++ program execution is suitable for running PocketSimulator.
## Running Circuits
### Pre-built circuits
PocketSimulator comes pre-loaded with features for writing and testing certain types of circuits using the various algorithms in the project, accessed using the ```int``` variables ```circuitSetting``` and ```algorithmSetting```.

```circuitSetting``` details:

Value | Effect
---|---
0 | Execute a user-inputted circuit from gates.txt
1 | Write and execute a Hadamard-Toffoli layered circuit, consisting of two n-Hadamard layers surrounding a randomly generated collection of n^2 toffoli gates, for a total of n^2 + 2n gates.
2 | Write and execute a QFT-layered circuit, consisting of two quantum Fourier transforms surrounding a randomly generated collection of n^2 toffoli gates, for a total of n^2 + 2n gates.
3 | Write and execute a Draper adder circuit.

```algorithmSetting``` details:

Value | Effect
---|---
0 | Simulate using the recursive path-summing algorithm (```pathIntegral.cpp```)
1 | Simulate using the state vector algorithm (```stateVector.cpp```)
2 | Simulate using the recursive Aaronson method (```savitch.cpp```)

### Parameters
PocketSimulator takes several arguments for simulation:
- **N**: # of qubits
- **nonPhaseGates**: # of "changing" gates in the circuit (all gates excluding those which purely add a relative phase). Used only for custom user-inputted circuits (```circuitSetting = 0```)
- **gates.txt**: text file encoding the computation to be simulated
- **startState** and **endState**

```gates.txt``` can be edited directly to input a desired algorithm (sequence of quantum gates), while the remaining parameters are edited directly in ```main.cpp```.

### Gates
Currently, the following gates are/will be supported in PocketSimulator:

Gate | Abbreviation | Arguments | Example
---|---|---|---
Hadamard | h | target qubit (h \_) | h 0
Toffoli | t | 2 controls and a target qubit (t c1 c2 \_) | t 0 1 2
U (1/(2^a) phase) | U | phase power a, target qubit (U a \_) | U 3 0 [phase of e^(i\*pi/4]
u (-1/(2^a) phase) | u | phase power a, target qubit (U a \_) | U 4 0 [phase of e^(-i\*pi/8]
(\*) Pauli X | x | target qubit (x \_) | x 0
(\*) Pauli Y | y | target qubit (y \_) | y 0
(\*) Pauli Z | z | target qubit (z \_) | z 0
(\*) General Phase | R | numerator denominator target (R n d \_) | R 3 7 0 [phase of e^(i\*pi\*6/7)]

(\*) coming soon

Gates are inputted in chronological order in ```gates.txt```. A control bit with value of either 0 or 1 must precede each gate inputted; if it equals 1, the gate is modified as a controlled operation and will need to specify an additional argument for the control qubit number.

### Execution
Upon execution of algorithm U, PocketSimulator will return a complex probability amplitude <endState|U|startState>, as well as time used (in seconds) by the execution as returned by the system method `getrusage()` ([documentation here](http://pubs.opengroup.org/onlinepubs/009695399/functions/getrusage.html)).

## About
This project was created for the 2017 Siemens Competition.
