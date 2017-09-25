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
/* SIMULATION SETTING VARIABLES
 N: number of qubits to simulate
 startState, endState: states to simulate and compute <endState|C|startState>
 showRuntime: toggle algorithm time details at end of simulation
 gatePath: directory path to gate input file
 
 Changing the 'circuitSetting' variable allows you to choose between executing+writing different circuits.
 
 0 = Execute a user-inputted circuit from gates.txt
 
 1 = write and execute a Hadamard-Toffoli layered circuit, consisting of two n-Hadamard layers surrounding a randomly generated collection of n toffoli gates.
 
 2 = write and execute a QFT layered circuit, consisting of two QFT circuits surrounding a randomly generated collection of n toffoli gates.
 
 3 = write and execute an "HSP standard method" circuit.
 
 4 = write and execute a QFT circuit.
 
 5 = write and execute a Draper adder circuit (used in SEQCSim)
 
 The algorithmSetting variable controls whether to run the PocketSimulator recursive algorithm (= 0), the classic state vector implementation (= 1), or Aaronson's simulation algorithm (= 2). */

int N = 18;
int startState = rand()%(int)pow(2,N), endState = rand()%(int)pow(2,N);
bool showRuntime = true; //controls whether runtime details are printed on console
string gatePath = "/Users/AShi/Documents/Repos/PocketSimulator/PocketSimulator/gates.txt"; //Directory path to gate file
ifstream in = ifstream(gatePath);

int circuitSetting = 3; //Circuit setting control
int algorithmSetting = 2; //Algorithm setting control

//VARIABLE FOR SETTING 0 ONLY: user-inputted circuit
int nonPhaseGates = 0; //Number of gates in circuit EXCLUDING PHASE GATES
//------------------------------------------------------------------------------------------------

//------------------------------------MAIN METHOD-------------------------------------------------

int main(int argc, const char * argv[]){
    srand((int)time(0)); //Set seed for generating random Toffoli gates
    cout << fixed;
    ofstream out (gatePath);
    string circuit;
    
    switch (circuitSetting){
        case 0: //Execute user-inputted circuit from gates.txt
        {
            break;
        }
        case 1: //Write and execute layered-Hadamard circuit
        {
            /* The circuit consists of two n-Hadamard layers surrounding a randomly generated collection of n toffoli gates, for a total of 3n gates. */
            nonPhaseGates = 3*N;
            circuit = writeCircuit(N, false, N);
            out << circuit;
            out.close();
            cout << "Circuit type: [layered Hadamard]\n";
            break;
        }
        case 2: //Write and execute layered-QFT circuit
        {
            /* The circuit consists of two n-Hadamard layers surrounding a randomly generated collection of n toffoli gates, for a total of 3n gates. */
            nonPhaseGates = 3*N;
            circuit = writeCircuit(N, true, N);
            out << circuit;
            out.close();
            cout << "Circuit type: [layered QFT]\n";
            break;
        }
        case 3: //Write and execute an "HSP standard method" circuit
        {
            nonPhaseGates = (int)(2*N/3)*2 + N;
            startState = 0;
            circuit = paradigmCircuit(2*N/3, N);
            out << circuit;
            out.close();
            cout << "Circuit type: [HSP standard method]\n";
            break;
        }
        case 4: //Write and execute QFT circuit
        {
            /* The circuit consists of a QFT (Quantum Fourier Transform) circuit with N branching gates. */
            nonPhaseGates = N;
            circuit = writeQFT(N);
            out << circuit;
            out.close();
            cout << "Circuit type: [QFT]\n";
            break;
        }
        case 5: //Write and execute draper adder
        {
            nonPhaseGates = N;
            circuit = writeAdder(N);
            out << circuit;
            out.close();
            
            int a = rand()%(int)pow(2,N/2), b = rand()%(int)pow(2,N/2), sum = (a + b)%(int)pow(2,N/2);
            startState = a*pow(2, N/2) + b, endState = startState - b + sum;
            
            cout << "Circuit type: [Draper adder]\n";
            cout << "Confirming addition of " << a << " + " << b << " = " << sum << " (modulo " << (int)pow(2, N/2) << ")\n";
            break;
        }
        default: break;
    }
    
    switch(algorithmSetting){
        case 0: pathIntegral(gatePath, N, startState, endState, nonPhaseGates, showRuntime); break;
        case 1: stateVector(gatePath, N, startState, endState, false, showRuntime); break;
        case 2: savitch(gatePath, N, startState, endState, false, showRuntime); break;
        default: break;
    }
    
    if (showRuntime){
        cout << "Press enter once memory/time data has been collected.\n";
        string line;
        getline(cin, line);
    }
    return 0;
}
