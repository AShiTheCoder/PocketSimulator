/*  Main method kernel for PocketSimulator and other simulation methods for comparison
    
    Note: In a binary representation, qubits 0 through n-1 are represented from leftmost
    digit to rightmost digit (i.e. 6 = 110 = 0th qubit:1, 1st qubit:1, 2nd qubit:0)
 
    PocketSimulator
    Copyright © 2017. All rights reserved.
*/
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
#include "savitch.hpp"
#include "pathIntegral.hpp"

using namespace std;

//----------------------------------CONTROL PANEL--------------------------------------
//Input simulation parameters
int N = 5;
int startState = rand()%(int)pow(2,N), endState = rand()%(int)pow(2,N);
bool showRuntime = true; //controls whether runtime details are printed on console
string gatePath = "/Users/AShi/Documents/Repos/PocketSimulator/PocketSimulator/gates.txt"; //Directory path to gate file
ifstream in = ifstream(gatePath);

//SIMULATION SETTING VARIABLES

/* Changing the 'circuitSetting' variable allows you to choose between executing+writing different circuits.
 
 0 = Execute a user-inputted circuit from gates.txt
 
 1 = write and execute a Hadamard-Toffoli layered circuit, consisting of two n-Hadamard layers surrounding a randomly generated collection of n^2 toffoli gates, for a total of n^2 + 2n gates.
 
 2 = write and execute a QFT-layered circuit, consisting of two quantum Fourier transforms surrounding a randomly generated collection of n^2 toffoli gates, for a total of n^2 + 2n gates.
 
 3 = write and execute a Draper adder circuit (used in SEQCSim)
 
 The algorithmSetting variable controls whether to run the PocketSimulator recursive algorithm (= 0), the classic state vector implementation (= 1), or Aaronson's simulation algorithm (= 2). */

int circuitSetting = 2; //Circuit setting control
int algorithmSetting = 2; //Algorithm setting control

//VARIABLES FOR SETTING 0 ONLY: user-inputted circuit
int nonPhaseGates = 0; //Number of gates in circuit EXCLUDING PHASE GATES
bool cmplx = false; //Change this to true if there are complex-valued gates in the circuit

//------------------------------------MAIN METHOD----------------------------------------------

int main(int argc, const char * argv[]){
    srand((int)time(0)); //Set seed for generating random Toffoli gates
    cout << fixed;

    switch (circuitSetting){
        case 0: //Execute user-inputted circuit from gates.txt
        {
            switch(algorithmSetting){
                case 0: pathIntegral(gatePath, N, startState, endState, nonPhaseGates, cmplx, showRuntime); break;
                case 1: stateVector(gatePath, N, startState, endState, false, showRuntime); break;
                case 2: savitch(gatePath, N, startState, endState, false, showRuntime); break;
                default: break;
            }
            break;
        }
        case 1: //Write and execute layered-Hadamard circuit
        {
            ofstream out (gatePath);
            
            /* The circuit consists of two n-Hadamard layers surrounding a randomly generated collection of [length] toffoli gates, for a total of [length] + 2n gates. */
            int length = (int)pow(N,2);
            string circuit = writeCircuit(length, false, N);
            out << circuit;
            out.close();
            
            switch(algorithmSetting){
                case 0: pathIntegral(gatePath, N, startState, endState, pow(N,2) + 2*N, false, showRuntime); break;
                case 1: stateVector(gatePath, N, startState, endState, false, showRuntime); break;
                case 2: savitch(gatePath, N, startState, endState, false, showRuntime); break;
                default: break;
            }
            
            break;
        }
        case 2: //Write and execute layered-QFT circuit
        {
            ofstream out (gatePath);
            
            /* The circuit consists of 2 QFT (Quantum Fourier Transform) layers surrounding [length] random toffoli gates, for a total of [length] + 2n gates. */
            int length = (int)pow(N,2);
            string circuit = writeCircuit(length, true, N);
            out << circuit;
            out.close();
            
            switch(algorithmSetting){
                case 0: pathIntegral(gatePath, N, startState, endState, pow(N,2) + 2*N, true, showRuntime); break;
                case 1: stateVector(gatePath, N, startState, endState, false, showRuntime); break;
                case 2: savitch(gatePath, N, startState, endState, false, showRuntime); break;
                default: break;
            }
            break;
        }
        case 3: //Write and execute draper adder
        {
            ofstream out (gatePath);
            string circuit = writeAdder(N);
            out << circuit;
            out.close();
            
            int a = rand()%(int)pow(2,N/2), b = rand()%(int)pow(2,N/2), sum = (a + b)%(int)pow(2,N/2);
            startState = a*pow(2, N/2) + b, endState = startState - b + sum;
            
            switch(algorithmSetting){
                case 0: pathIntegral(gatePath, N, startState, endState, N, true, showRuntime); break;
                case 1: stateVector(gatePath, N, startState, endState, false, showRuntime); break;
                case 2: savitch(gatePath, N, startState, endState, false, showRuntime); break;
                default: break;
            }
            
            cout << "Confirmed addition of " << a << " + " << b << " = " << sum << " (modulo " << (int)pow(2, N/2) << ")\n";
            break;
        }
        default: break;
    }
    
    if (showRuntime){
        cout << "Press enter once memory/time data has been collected.\n";
        string line;
        getline(cin, line);
    }
    return 0;
}
