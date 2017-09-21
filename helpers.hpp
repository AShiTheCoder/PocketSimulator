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

/* writeCircuit: writes layered/dispersed circuits for n qubits as a string (with Hadamard + Toffoli).
 All Toffoli gates are randomly generated. */
string writeCircuit(int length, bool layered, int N);

string writeQFT(int N); //writes a QFT circuit on n qubits

string writeAdder(int N); //Writes a draper adder circuit

int bitDiff(int a, int b); //Returns bit difference between a and b

void resetCounter(bool *reached, int n); //helper method for layer separation in Aaronson's Savitch implementation

#endif /* helpers_hpp */
