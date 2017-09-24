//
//  helpers.cpp
//  PocketSimulator
//
//  Created on 9/20/17.
//  Copyright © 2017. All rights reserved.
//
#include "helpers.hpp"
#include <string>
#include <iostream>
#include <math.h>

using namespace std;

//----------------------------------AUXILIARY METHODS--------------------------------------

string binString(int x, int N){ //helper method for int x to a binary string
    int a = x;
    string result = "";
    while (a > 0){
        result = to_string(a % 2) + result;
        a /= 2;
    }
    while (result.length() < N){
        result = "0" + result;
    }
    return result;
}

int reverseBit(int x, int len){ //returns x reversed in binary (bitstring length l)
    int a = x, result = 0, counter = len - 1;
    while (a > 0){
        result += (a % 2) * pow(2, counter);
        counter--;
        a /= 2;
    }
    return result;
}

string randToff(int *indices, int N){ //generates a random Toffoli gate in indices[3]
    if (N < 3) {
        cout << "not enough qubits\n";
        return "fail";
    }
    
    int i = 0, temp;
    string result = "0 t";
    bool success;
    while (i < 3){
        success = true;
        temp = rand() % N;
        for (int j = 0; j < i; j++){
            if (indices[j] == temp) success = false;
        }
        if (success){
            result = result + " " + to_string(temp);
            indices[i] = temp;
            i++;
        }
    }
    return result;
}

/* writeCircuit: writes layered/dispersed circuits for n qubits as a string (with Hadamard + Toffoli).
 All Toffoli gates are randomly generated. */
string writeCircuit(int length, bool QFT, int N){
    string out = "";
    int toff[3];
    if (QFT) out += writeQFT(N);
    else {
        for (int i = 0; i < N; i++){
            out = out + "0 h " + to_string(i) + "\n";
        }
    }
    for (int i = 0; i < length; i++){
        out = out + randToff(toff, N) + "\n";
    }
    if (QFT) out += writeQFT(N);
    else {
        for (int i = 0; i < N; i++){
            out = out + "0 h " + to_string(i) + "\n";
        }
    }
    return out;
}

string writeQFT(int N){ //writes a QFT circuit on n qubits
    string out = "";
    for (int i = 0; i < N; i++){
        out = out + "0 h " + to_string(i) + "\n"; //hadamard
        for (int j = 2; j <= N - i; j++){
            out = out+"1 U "+to_string(j)+" "+to_string(i + j - 1)+" "+to_string(i)+"\n";
        }
    }
    return out;
}

string writeAdder(int N){ //Writes a draper adder circuit
    string out = "";
    for (int i = N/2; i < N; i++){
        out = out + "0 h " + to_string(i) + "\n"; //hadamard
        for (int j = 2; j <= N - i; j++){
            out = out+"1 U "+to_string(j)+" "+to_string(i + j - 1)+" "+to_string(i)+"\n";
        }
    }
    for (int i = 0; i < N/2; i++){
        for (int j = 0; j < N/2 - i; j++){
            out = out+"1 U "+to_string(i + 1)+" "+to_string(N/2 + j)+" "+to_string(j + i)+"\n";
        }
    }
    for (int i = N - 1; i >= N/2; i--){
        for (int j = N - i; j >= 2; j--){
            out = out+"1 u "+to_string(j)+" "+to_string(i + j - 1)+" "+to_string(i)+"\n";
        }
        out = out + "0 h " + to_string(i) + "\n"; //hadamard
    }
    return out;
}

int bitDiff(int a, int b){
    return __builtin_popcount(a ^ b);
}

void resetCounter(bool *reached, int n){
    for (int i = 0; i < n; i++){
        reached[i] = false;
    }
}
