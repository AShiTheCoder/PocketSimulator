//
//  stateVector.cpp
//  PocketSimulator
//
//  Created on 9/20/17.
//  Copyright © 2017. All rights reserved.
//
#include <iostream>
#include <complex>
#include <fstream>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <math.h>
#define _USE_MATH_DEFINES

#include "helpers.hpp"
#include "stateVector.hpp"

using namespace std;
//STATE VECTOR VARIABLES
#define MAX_LAYERS 1000
complex<double> amps[INT_MAX/2048]; //used for amplitude storage

//---------------------------------STATE VECTOR EVOLUTION----------------------------------

/* First simulation algorithm: tracking entire state vector
 Takes time T*exp(O(n)) and space exp(O(n)) [T = total # of gates]
 
 PARAMETERS:
 in: file input stream to read gates from
 n: number of qubits
 startState: starting state of qubit register
 verbose: set to true to print intermediate amplitude values between each gate,
 false to only print the end amplitudes
 
 MODIFIED VERBOSE: TRUE = PRINT ALL END AMPLITUDES, FALSE = ONLY PRINTS "DONE"
 (because of very large state spaces yielding massive console outputs, verbose was adjusted from the previous definition.) */

void stateVector(string gatePath, int N, int startState, int endState, bool verbose, bool showRuntime){
    cout << "Comparison algorithm: [stateVector]\nSimulation in progress........\n";
    ifstream in = ifstream(gatePath);
    int spaceSize = (int)pow(2,N);
    
    for (int i = 0; i < spaceSize; i++){ //initialize amps array
        amps[i] = 0;
    }
    amps[startState] = 1; //amplitude of the starting state is one
    in.clear();
    in.seekg(0, ios::beg); //go to beginning of file
    
    char gate;
    int control, c1, c2, target, ctrl;
    int index, toff;
    int Hplus, H, zeroI, oneI;
    double temp;
    complex<double> ctemp;
    complex<double> phase;
    int phasePow;
    
    in >> control >> gate;
    while (!in.eof()){
        switch (gate) {
            case 'h': //hadamard gate
            {
                if (verbose) cout << "hadamard detected\n";
                in >> target;
                complex<double> zero, one;
                Hplus = pow(2, N - target), H = Hplus/2;
                
                //iterate over state vectors where target qubit == 0, updating both vectors (where the target = 0 and = 1) together
                for (int i = 0; i < spaceSize; i += Hplus){
                    for (int j = 0; j < H; j++){
                        zeroI = i + j;
                        oneI = zeroI | (1 << (N - target - 1));
                        zero = 1/sqrt(2) * amps[zeroI];
                        one = 1/sqrt(2) * amps[oneI];
                        amps[zeroI] = zero + one;
                        amps[oneI] = zero - one;
                    }
                }
                break;
            }
            case 't':
            {
                if (verbose) cout << "toffoli detected\n";
                in >> c1 >> c2 >> target;
                if (c1 > c2){
                    temp = c1;
                    c1 = c2;
                    c2 = temp;
                }
                
                /* Iterate over state vectors where both control qubits are 1. The previous approach looped over all states, but this one saves time by factor of 4 by only iterating over those where toffoli actually does something. */
                int inci = pow(2, N - c1), incj = pow(2, N - c2), C1 = inci/2, C2 = incj/2;
                for (int i = 0; i < spaceSize; i += inci){ //increments qubits before c1
                    for (int j = 0; j < C1; j += incj){ //increments qubits between c1 and c2
                        for (int k = 0; k < C2; k++){ //increments qubits after c2
                            index = k + C1 + j + C2 + i; //index where c1, c2 == 1
                            if (((index >> (N - target - 1)) & 1) == 0){
                                toff = index ^ (1 << (N - target - 1));
                                ctemp = amps[toff];
                                amps[toff] = amps[index];
                                amps[index] = ctemp;
                            }
                        }
                    }
                }
                break;
            }
            case 'U':
            {
                in >> phasePow;
                phase = polar(1.0, 1/pow(2, phasePow) * 2 * M_PI);
                if (control){ // controlled gate case
                    in >> ctrl >> target;
                    if (target > ctrl){
                        temp = ctrl;
                        ctrl = target;
                        target = temp;
                    }
                    int inci = pow(2, N - target), incj = pow(2, N - ctrl), C = inci/2, targBit = incj/2;
                    for (int i = 0; i < spaceSize; i += inci){
                        for (int j = 0; j < C; j += incj){
                            for (int k = 0; k < targBit; k++) amps[i + C + j + targBit + k] *= phase;
                        }
                    }
                } else { // non-controlled gate case
                    in >> target;
                    int inc = pow(2, N - target), targBit = inc/2;
                    for (int i = 0; i < spaceSize; i += inc){
                        for (int j = 0; j < targBit; j++) amps[i + targBit + j] *= phase;
                    }
                }
                break;
            }
            case 'u':
            {
                in >> phasePow;
                phase = polar(1.0, -1/pow(2, phasePow) * 2 * M_PI);
                if (control){ // controlled gate case
                    in >> ctrl >> target;
                    if (target > ctrl){
                        temp = ctrl;
                        ctrl = target;
                        target = temp;
                    }
                    int inci = pow(2, N - target), incj = pow(2, N - ctrl), C = inci/2, targBit = incj/2;
                    for (int i = 0; i < spaceSize; i += inci){
                        for (int j = 0; j < C; j += incj){
                            for (int k = 0; k < targBit; k++) amps[i + C + j + targBit + k] *= phase;
                        }
                    }
                } else { // non-controlled gate case
                    in >> target;
                    int inc = pow(2, N - target), targBit = inc/2;
                    for (int i = 0; i < spaceSize; i += inc){
                        for (int j = 0; j < targBit; j++) amps[i + targBit + j] *= phase;
                    }
                }
                break;
            }
            default:
            {
                cout << "Incompatible gate type: " << gate << "\n";
                break;
            }
        }
        in >> control >> gate;
    }
    if (verbose){
        for (int i = 0; i < spaceSize; i++){
            cout << binString(i, N) << ": " << amps[i] << "\n";
        }
    }
    cout << "Finished computation on " << N << " qubits\n";
    cout << "<" << binString(endState, N) << "|Circuit|" << binString(startState, N)
    << ">(vector) = " << amps[endState] << "\n";
    
    if (showRuntime){ //Print time usage
        cout.precision(7);
        struct rusage usage;
        getrusage(RUSAGE_SELF, &usage);
        long totaluTime = (usage.ru_stime.tv_sec + usage.ru_utime.tv_sec) * 1000000 + usage.ru_stime.tv_usec + usage.ru_utime.tv_usec;
        double totalTime = totaluTime/ (double) 1000000;
        cout << "Runtime: " << totalTime << " seconds\n";
        //        cout << "Memory usage: " << usage.ru_maxrss / (double) memConst << " qunits [1 qunit ≈ 1 mb]\n\n";
        //Memory usage details removed due to unclear units
    }
    cout << "\n";
}
