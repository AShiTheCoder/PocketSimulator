//
//  pathIntegral.cpp
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

#include "pathIntegral.hpp"
#include "helpers.hpp"
using namespace std;

// Global memory storage for recursive calls
extern int N;
extern int startState, endState;
extern ifstream in;
int currState;
double amplitudes[50];
complex<double> cmplxAmps[50];
//---------------------------------PATH INTEGRAL SUMMING-----------------------------------

/* A recursive path-summing simulation algorithm
 Takes time O(t2^h) and space O(h) + O(n) [t = # of "non-branching" gates, h = # of "branching" gates]
 
 V1: first version
 V2: changed to DFS procedure
 V3: added out-of-reach path pruning
 V4: added QFT: controlled-U gates, complex numbers, phase accumulation
 V5: globalized variables to minimize space usage, rearranged parameters */

void pathStep(streampos pos, int changesLeft, double currPhase, int currDepth){
    char gate;
    int control, c1, c2, target, oneFactor = 1;
    
    in.clear();
    in.seekg(pos); //move to correct position
    in >> control >> gate; //read control boolean and gate
    
    streampos currPos;
    while (!in.eof()){
        switch (gate){ //check the type of gate read
            case 'h': //Hadamard gate
            {
                changesLeft--;
                in >> target;
                currPos = in.tellg();
                //|0><+| case: amp is always +1
                //|1><-| case: if the target qubit is a 1, amp turns negative; stays positive otherwise
                if (((currState >> (N - target - 1)) & 1) == 1) oneFactor = -1;
                
                if (bitDiff(currState, endState) <= (changesLeft + 1)){ //is the end state reachable?
                    currState &= ~(1 << (N - target - 1)); //change state register to 0 branch
                    pathStep(currPos, changesLeft, 1/sqrt(2) * currPhase, currDepth + 1); //travel down the 0 branch
                    amplitudes[currDepth] = amplitudes[currDepth + 1]; //record 0 branch's intermediate amplitude
                    currState |= (1 << (N - target - 1)); //change state register to 1 branch
                    pathStep(currPos, changesLeft, oneFactor/sqrt(2) * currPhase, currDepth + 1); //travel down the 1 branch
                    amplitudes[currDepth] += amplitudes[currDepth + 1]; //add 1 branch's intermediate amplitude
                    if (oneFactor == 1) currState &= ~(1 << (N - target - 1));
                    else currState |= (1 << (N - target - 1));
                    //reset the state to its pre-branching value (necessary "cleanup" for other calls)
                } else amplitudes[currDepth] = 0; //otherwise, terminate computation prematurely
                return;
                break;
            }
            case 't': //Toffoli gate
            {
                changesLeft--;
                in >> c1 >> c2 >> target;
                currPos = in.tellg();
                int add = ((currState >> (N - c1 - 1)) & 1) * ((currState >> (N - c2 - 1)) & 1);
                if (bitDiff(currState, endState) <= (changesLeft + 1)){ //is the end state reachable?
                    currState ^= (add << (N - target - 1)); //Toffoli state
                    pathStep(currPos, changesLeft, currPhase, currDepth); //Step forwards and compute
                    currState ^= (add << (N - target - 1)); //Un-toffoli state
                } else amplitudes[currDepth] = 0;
                return;
            }
            default: break;
        }
        in >> control >> gate;
    }
    
    if (currState == endState){ //inner product <a|C|b> is 0 unless end state |a> matches start state |b>
        amplitudes[currDepth] = currPhase;
    } else amplitudes[currDepth] = 0;
}

void complexPathStep(streampos pos, int changesLeft, complex<double> currPhase, int currDepth){
    char gate;
    int control, c = -1, c1, c2, phasePow, target, oneFactor = 1;
    
    in.clear();
    in.seekg(pos); //move to correct position
    in >> control >> gate; //read control boolean and gate
    
    streampos currPos;
    while (!in.eof()){
        switch (gate){ //check the type of gate read
            case 'h': //Hadamard gate
            {
                changesLeft--;
                in >> target;
                currPos = in.tellg();
                //|0><+| case: amp is always +1
                //|1><-| case: if the target qubit is a 1, amp turns negative; stays positive otherwise
                if (((currState >> (N - target - 1)) & 1) == 1) oneFactor = -1;
                
                if (bitDiff(currState, endState) <= (changesLeft + 1)){ //is the end state reachable?
                    currState &= ~(1 << (N - target - 1)); //travel down the 0 branch
                    complexPathStep(currPos, changesLeft, 1/sqrt(2) * currPhase, currDepth + 1);
                    cmplxAmps[currDepth] = cmplxAmps[currDepth + 1];
                    currState |= (1 << (N - target - 1)); //travel down the 1 branch
                    complexPathStep(currPos, changesLeft, oneFactor/sqrt(2) * currPhase, currDepth + 1);
                    cmplxAmps[currDepth] += cmplxAmps[currDepth + 1];
                    if (oneFactor == 1) currState &= ~(1 << (N - target - 1));
                    else currState |= (1 << (N - target - 1)); //reset the state
                } else cmplxAmps[currDepth] = 0; //otherwise, terminate computation prematurely
                return;
                break;
            }
            case 't': //Toffoli gate
            {
                changesLeft--;
                in >> c1 >> c2 >> target;
                currPos = in.tellg();
                int add = ((currState >> (N - c1 - 1)) & 1) * ((currState >> (N - c2 - 1)) & 1);
                if (bitDiff(currState, endState) <= (changesLeft + 1)){ //is the end state reachable?
                    currState ^= (add << (N - target - 1)); //Toffoli state
                    complexPathStep(currPos, changesLeft, currPhase, currDepth); //Step forwards and compute
                    currState ^= (add << (N - target - 1)); //Un-toffoli state
                } else cmplxAmps[currDepth] = 0;
                return;
            }
            case 'U':
            {
                in >> phasePow;
                complex<double> phase = polar(1.0, 1/pow(2, phasePow) * 2 * M_PI);
                if (control){ // controlled gate case
                    in >> c >> target;
                    if (((currState >> (N - c - 1)) & 1) && ((currState >> (N - target - 1)) & 1)) currPhase *= phase;
                } else { // non-controlled gate case
                    in >> target;
                    if ((currState >> (N - target - 1)) & 1) currPhase *= phase;
                }
                break;
            }
            case 'u':
            {
                in >> phasePow;
                complex<double> phase = polar(1.0, -1/pow(2, phasePow) * 2 * M_PI);
                if (control){ // controlled gate case
                    in >> c >> target;
                    if (((currState >> (N - c - 1)) & 1) && ((currState >> (N - target - 1)) & 1)) currPhase *= phase;
                } else { // non-controlled gate case
                    in >> target;
                    if ((currState >> (N - target - 1)) & 1) currPhase *= phase;
                }
                break;
            }
            default: break;
        }
        in >> control >> gate;
    }
    
    if (currState == endState){ //inner product <a|C|b> is 0 unless end state |a> matches start state |b>
        cmplxAmps[currDepth] = currPhase;
    } else cmplxAmps[currDepth] = 0;
}

void pathIntegral(string gatePath, int n, int startS, int endS, int numChanges, bool cmplx, bool printMem){
    cout << "Main Method: [PocketSimulator]\nSimulation in progress........\n\n\n";
    currState = startS, startState = startS, endState = endS;
    N = n;
    in = ifstream(gatePath);
    
    //initial recursive call (the "root" of the path tree)
    if (cmplx){
        complexPathStep(ios::beg, numChanges, 1, 0);
        cout << "<" << binString(endS, N) << "|Circuit|" << binString(startS, N) << "> = " << cmplxAmps[0] << "\n";
    } else {
        pathStep(ios::beg, numChanges, 1, 0);
        cout << "<" << binString(endS, N) << "|Circuit|" << binString(startS, N) << "> = " << amplitudes[0] << "\n";
    }
    
    if (printMem){ //Print time usage
        cout.precision(7);
        struct rusage usage;
        getrusage(RUSAGE_SELF, &usage);
        long totaluTime = (usage.ru_stime.tv_sec + usage.ru_utime.tv_sec) * 1000000 + usage.ru_stime.tv_usec + usage.ru_utime.tv_usec;
        double totalTime = totaluTime/ (double) 1000000;
        cout << "Runtime: " << totalTime << " seconds\n";
        //        cout << "Memory usage: " << usage.ru_maxrss / (double) memConst << " qunits [1 qunit ≈ 1 mb]\n\n";
        //Memory usage details removed due to unclear units
    }
}
