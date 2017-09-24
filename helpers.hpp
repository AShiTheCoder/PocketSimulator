//
//  helpers.hpp
//  PocketSimulator
//
//  Created on 9/20/17.
//  Copyright © 2017. All rights reserved.
//

#ifndef helpers_hpp
#define helpers_hpp

#include <stdio.h>
#include <string>
using namespace std;

string binString(int x, int N); //helper method for int x to a binary string

int reverseBit(int x, int len); //returns x reversed in binary (bitstring length l)

string randToff(int *indices, int N); //generates a random Toffoli gate in indices[3]

/* randControlToff: generate a random Toffoli with control bits in the first a bits and target bit in the rest of the N-bit register (N-a choices) */
string randControlToff(int *indices, int a, int N);

/* paradigmCircuit: writes a "HSP standard method" circuit for N qubits, using the first a qubits as the control register (the one we will eventually measure) writing N random toffoli gates.
 
 HSP Standard Method:
 1) Put a-bit register into superposition with Hadamards
 2) Compute the function f(a --> b), storing the result in a b-bit register (a + b = N)
 (in this case, our function is a collection of random toffoli gates.
 3) Perform the appropriate actions on the a-bit register to solve the problem (in this case, we use the QFT that Shor's algorithm uses).
 
 All Toffoli gates are randomly generated within the control restrictions. */
string paradigmCircuit(int a, int N);

/* writeCircuit: writes layered/dispersed circuits for n qubits as a string (with Hadamard + Toffoli).
 All Toffoli gates are randomly generated. */
string writeCircuit(int length, bool QFT, int N);

string writeHlayer(int N); //writes a Hadamard layer on n qubits

string writeQFT(int N); //writes a QFT circuit on n qubits

string writeAdder(int N); //Writes a draper adder circuit

int bitDiff(int a, int b); //Returns bit difference between a and b

void resetCounter(bool *reached, int n); //helper method for layer separation in Aaronson's Savitch implementation

#endif /* helpers_hpp */
