//
//  savitch.cpp
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
#define MAX_LAYERS 1000 //Set a cap on the maximum simulation layer count

#include "helpers.hpp"
#include "savitch.hpp"

using namespace std;

//AARONSON VARIABLES
bool bitReached[32];
string layerGates[MAX_LAYERS];
int layers[MAX_LAYERS];
extern ifstream in;

//----------------------------------AARONSON RECURSION-------------------------------------

/* This is the algorithm described in Aaronson/Chen's paper (arXiv:1612.05903 [quant-ph]) based off of Savitch's Theorem, section 4. It simulates a quantum circuit with a recursive procedure in time O(n*(2d)^(n+1)) and space O(nlog(d)). d is the circuit depth, or the number of gate groups applied chronologically where each gate group only acts on a qubit 0 or 1 times. (effectively we assume d ~ T, the total # of gates).
 UPDATE: takes less time now due to some adjustments
 
 There is a modified tradeoff algorithm (also from Aaronson/Chen) exchanging time and space resources where for a chosen integer k, the algorithm takes time O(n*2^(n+1)*d^(k+1)) and space O(2^(n-k)logd). However, note that algorithmThree is approximately a special case where k = n, and when k = 0 both space and time are of complexity exp(n).
 
 The tradeoff algorithm is not implemented here because despite the extra flexibility, time is the limiting factor preventing scaling to more qubits in all three algorithms and is still lower bounded by O(nexp(n)) in the tradeoff, a larger complexity than even algorithm one.
 
 V1: first version
 V2: added zero-term checking
 V3: added out-of-reach path pruning */

complex<double> savitchRecur(int N, int beginD, int endD, int startS, int endS, int *layers, bool verbose){ //Recursive subalgorithm for algorithm three
    complex<double> result = complex<double>(0);
    if (verbose) cout << beginD << "(" << startS << ") to " << endD << "(" << endS << ")\n";
    if (beginD == endD) { //base case
        result = 1;
        in.clear();
        in.seekg(ios::beg);
        stringstream gates;
        gates << layerGates[beginD];
        
        char gate;
        int c1, c2, target, endBit, control, c;
        int qubits = startS;
        int phasePow;
        
        gates >> control >> gate;
        while (!gates.eof()){ //act on gates inside this layer
            switch (gate){
                case 'h': //hadamard
                {
                    gates >> target;
                    endBit = (endS >> (N - target - 1)) & 1;
                    if ((((qubits >> (N - target - 1)) & 1) == 1) && endBit == 1){
                        //if the hadamarded bit is 1 in both the start and end states, the end amplitude will be negative
                        result *= -1;
                    }
                    //set the hadamarded bit in the register (qubits) to whatever it equals in the end state
                    qubits = qubits & ~(1 << (N - target - 1));
                    qubits += (endBit << (N - target - 1));
                    result /= sqrt(2); //hadamard-adjusted amplitude
                    break;
                }
                case 't': //toffoli
                {
                    gates >> c1 >> c2 >> target;
                    int add = ((qubits >> (N - c1 - 1)) & 1) * ((qubits >> (N - c2 - 1)) & 1);
                    qubits = qubits ^ (add << (N - target - 1)); //update register
                    break;
                }
                case 'U':
                {
                    gates >> phasePow;
                    complex<double> phase = polar(1.0, 1/pow(2, phasePow) * 2 * M_PI);
                    if (control){ // controlled gate case
                        gates >> c >> target;
                        if (((qubits >> (N - c - 1)) & 1) && ((qubits >> (N - target - 1)) & 1)) result *= phase;
                    } else { // non-controlled gate case
                        gates >> target;
                        if ((qubits >> (N - target - 1)) & 1) result *= phase;
                    }
                    break;
                }
                case 'u':
                {
                    gates >> phasePow;
                    complex<double> phase = polar(1.0, -1/pow(2, phasePow) * 2 * M_PI);
                    if (control){ // controlled gate case
                        gates >> c >> target;
                        if (((qubits >> (N - c - 1)) & 1) && ((qubits >> (N - target - 1)) & 1)) result *= phase;
                    } else { // non-controlled gate case
                        gates >> target;
                        if ((qubits >> (N - target - 1)) & 1) result *= phase;
                    }
                    break;
                }
                default:
                {
                    cout << "Incompatible gate type: " << gate << "\n";
                    break;
                }
            }
            gates >> control >> gate;
        }
        // <endS|qubits> == 0 if endS ≠ qubits [|qubits> = Layer|startS>]
        if (qubits != endS) result = 0;
        return result;
    } else { //recursive case
        /* Compute <y|C|x> by summing all <y|C_1|i><i|C_2|x> for i = {0,1}^n.
         Compute the two sub terms recursively. */
        for (int i = 0; i < pow(2, N); i++){
            if (bitDiff(startS, i) <= (layers[(beginD + endD)/2 + 1] - layers[beginD]) &&
                bitDiff(i, endS) <= (layers[endD + 1] - layers[(beginD + endD)/2 + 1])){
                complex<double> termOne = savitchRecur(N, beginD, (beginD + endD)/2, startS, i, layers, verbose);
                if (termOne != complex<double>(0)){ //only compute second term if first is nonzero
                    result += (termOne * savitchRecur(N, (beginD + endD)/2 + 1, endD, i, endS, layers, verbose));
                }
            }
        }
    }
    return result;
}

void savitch(string gatePath, int N, int startState, int endState, bool verbose, bool showRuntime){
    cout << "Comparison algorithm: [Aaronson's Savitch]\n" << N << " qubit simulation in progress........\n";
    char gate;
    int c1, c2, target, depth = 0, gCount = 0, ctrl, c;
    string gateBuffer = "";
    int phasePow;
    
    in = ifstream(gatePath);
    in.clear();
    in.seekg(ios::beg);
    in >> ctrl >> gate;
    resetCounter(bitReached, N);
    layers[0] = 0;
    while (!in.eof()){ //separate circuit into layers, record layer dividers in layers[], record gate sequences in layerGates[]
        switch (gate){
            case 'h':
            {
                in >> target;
                if (bitReached[target]){
                    layerGates[depth] = gateBuffer;
                    depth++;
                    layers[depth] = gCount;
                    resetCounter(bitReached, N);
                    gateBuffer = "";
                }
                bitReached[target] = true;
                gateBuffer = gateBuffer + "0 h " + to_string(target) + " ";
                break;
            }
            case 't':
            {
                in >> c1 >> c2 >> target;
                if (bitReached[c1] || bitReached[c2] || bitReached[target]){
                    layerGates[depth] = gateBuffer;
                    depth++;
                    layers[depth] = gCount;
                    resetCounter(bitReached, N);
                    gateBuffer = "";
                }
                bitReached[c1] = true;
                bitReached[c2] = true;
                bitReached[target] = true;
                gateBuffer = gateBuffer + "0 t " + to_string(c1) + " " + to_string(c2) + " " +to_string(target) + " ";
                break;
            }
            case 'U':
            {
                in >> phasePow;
                if (ctrl){ // controlled gate case
                    in >> c >> target;
                    if (bitReached[c] || bitReached[target]){
                        layerGates[depth] = gateBuffer;
                        depth++;
                        layers[depth] = gCount;
                        resetCounter(bitReached, N);
                        gateBuffer = "";
                    }
                    bitReached[c] = true;
                    bitReached[target] = true;
                    gateBuffer = gateBuffer + "1 U " + to_string(phasePow) + " " + to_string(c) + " " +to_string(target) + " ";
                } else { // non-controlled gate case
                    in >> target;
                    if (bitReached[target]){
                        layerGates[depth] = gateBuffer;
                        depth++;
                        layers[depth] = gCount;
                        resetCounter(bitReached, N);
                        gateBuffer = "";
                    }
                    bitReached[target] = true;
                    gateBuffer = gateBuffer + "0 U " + to_string(phasePow) + " " +to_string(target) + " ";
                }
                break;
            }
            case 'u':
            {
                in >> phasePow;
                if (ctrl){ // controlled gate case
                    in >> c >> target;
                    if (bitReached[c] || bitReached[target]){
                        layerGates[depth] = gateBuffer;
                        depth++;
                        layers[depth] = gCount;
                        resetCounter(bitReached, N);
                        gateBuffer = "";
                    }
                    bitReached[c] = true;
                    bitReached[target] = true;
                    gateBuffer = gateBuffer + "1 u " + to_string(phasePow) + " " + to_string(c) + " " +to_string(target) + " ";
                } else { // non-controlled gate case
                    in >> target;
                    if (bitReached[target]){
                        layerGates[depth] = gateBuffer;
                        depth++;
                        layers[depth] = gCount;
                        resetCounter(bitReached, N);
                        gateBuffer = "";
                    }
                    bitReached[target] = true;
                    gateBuffer = gateBuffer + "0 u " + to_string(phasePow) + " " +to_string(target) + " ";
                }
                break;
            }
            default:
            {
                cout << "Incompatible gate type: " << gate << "\n";
                break;
            }
        }
        in >> ctrl >> gate;
        gCount++;
    }
    layerGates[depth] = gateBuffer;
    depth++;
    layers[depth] = gCount;
    cout << "Divided into " << depth << " layers\n";
    
    complex<double> result = savitchRecur(N, 0, depth - 1, startState, endState, layers, verbose); //call recursive algorithm
    cout << "<" << binString(endState, N) << "|Circuit|" << binString(startState, N) << ">: " << result << "\n";
    
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
