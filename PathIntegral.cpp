/*  A path integral algorithm for simulating general quantum circuits
    
    Note: In a binary representation, qubits 0 through n-1 are represented from leftmost
    digit to rightmost digit (i.e. 6 = 110 = 0th qubit:1, 1st qubit:1, 2nd qubit:0)
 
    QuantumLab
    Created by Andrew Shi on 7/8/17.
    Copyright © 2017 Andrew Shi. All rights reserved.
*/
#include <iostream>
#include <complex>
#include <fstream>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <math.h>
#define _USE_MATH_DEFINES
#define MAX_LAYERS 1000

using namespace std;
double amps[INT_MAX]; //used for amplitude storage
long memConst = 1024 * 1024; //constant for interpreting rusage's memory output
int special = 0, cool = 0;

//----------------------------------AUXILIARY METHODS--------------------------------------

string binString(int x, int len){ //helper method for int x to binary string of length len
    int a = x;
    string result = "";
    while (a > 0){
        result = to_string(a % 2) + result;
        a /= 2;
    }
    while (result.length() < len){
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

string randToff(int *indices, int n){ //generates a random toffoli gate for n qubit circuit
    if (n < 3) {
        cout << "not enough qubits\n";
        return "fail";
    }
    
    int i = 0, temp;
    string result = "0 t";
    bool success;
    while (i < 3){
        success = true;
        temp = rand() % n;
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

void algorithmOne(ifstream &in, int n, int startState, bool verbose){
    int spaceSize = (int)pow(2,n);
    
    for (int i = 0; i < spaceSize; i++){ //initialize amps array
        amps[i] = 0;
    }
    amps[startState] = 1; //amplitude of the starting state is one
    in.clear();
    in.seekg(0, ios::beg); //go to beginning of file
    
    char gate;
    int c1, c2, target;
    in >> gate;
    while (!in.eof()){
        switch (gate) {
            case 'h': //hadamard gate
            {
                if (verbose) cout << "hadamard detected\n";
                in >> target;
                double zero, one;
                int Hplus = pow(2, n - target), H = Hplus/2, zeroI, oneI;
                
                //iterate over state vectors where target qubit == 0, updating both vectors (where the target = 0 and = 1) together
                for (int i = 0; i < spaceSize; i += Hplus){
                    for (int j = 0; j < H; j++){
                        zeroI = i + j;
                        oneI = zeroI | (1 << (n - target - 1));
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
                double temp;
                int index, toff;
                in >> c1 >> c2 >> target;
                if (c1 > c2){
                    temp = c1;
                    c1 = c2;
                    c2 = temp;
                }
                
                /* Iterate over state vectors where both control qubits are 1. The previous approach looped over all states, but this one saves time by factor of 4 by only iterating over those where toffoli actually does something. */
                int inci = pow(2, n - c1), incj = pow(2, n - c2), C1 = inci/2, C2 = incj/2;
                for (int i = 0; i < spaceSize; i += inci){ //increments qubits before c1
                    for (int j = 0; j < C1; j += incj){ //increments qubits between c1 and c2
                        for (int k = 0; k < C2; k++){ //increments qubits after c2
                            index = k + C1 + j + C2 + i; //index where c1, c2 == 1
                            if (((index >> (n - target - 1)) & 1) == 0){
                                toff = index ^ (1 << (n - target - 1));
                                temp = amps[toff];
                                amps[toff] = amps[index];
                                amps[index] = temp;
                            }
                        }
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
        in >> gate;
    }
    if (verbose){
        for (int i = 0; i < spaceSize; i++){
            cout << binString(i, n) << ": " << amps[i] << "\n";
        }
    }
    cout << "Finished computation on " << n << " qubits\n";
    cout << "<" << binString(startState, n) << "|Circuit|" << binString(startState, n)
    << "> = " << amps[startState] << "\n";
    
    //Print time/memory usage
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    long totaluTime = (usage.ru_stime.tv_sec + usage.ru_utime.tv_sec) * 1000000 + usage.ru_stime.tv_usec + usage.ru_utime.tv_usec;
    double totalTime = totaluTime/ (double) 1000000;
    cout << "Runtime: " << totalTime << " seconds\n";
    cout << "Memory usage: " << usage.ru_maxrss / (double) memConst << " qunits [1 qunit ≈ 1 mb]\n\n";
}

/* writeCircuit: returns a random circuit for n qubits as a string (in Hadamard + toffoli).
 If superRandom = true: all gates will be random (1/2 chance for either toffoli or hadamard), circuit will have (length) gates.
 If superRandom = false, a row of hadamards (one for each qubit) will start and end the circuit, with (length) toffoli gates in between. 
 All toffoli gates are randomly generated. */
string writeCircuit(int length, int n, bool superRandom){
    string out = "";
    int toff[3];
    if (superRandom){
        for (int i = 0; i < length; i++){
            if (rand() % 2){
                out = out + "0 h " + to_string(rand() % n);
            } else {
                out += randToff(toff, n);
            }
            out = out + "\n";
        }
    } else {
        for (int i = 0; i < n; i++){
            out = out + "0 h " + to_string(i) + "\n";
        }
        for (int i = 0; i < length; i++){
            out = out + randToff(toff, n) + "\n";
        }
        for (int i = 0; i < n; i++){
            out = out + "0 h " + to_string(i) + "\n";
        }
    }
    return out;
}

string writeQFT(int n){ //writes a QFT circuit on n qubits
    string out = "";
    for (int i = 0; i < n; i++){
        out = out + "0 h " + to_string(i) + "\n"; //hadamard
        for (int j = 2; j <= n - i; j++){
            out = out+"1 U "+to_string(j)+" "+to_string(i + j - 1)+" "+to_string(i)+"\n";
        }
    }
    return out;
}

int bitDiff(int a, int b){
    int x = a ^ b;
    std::bitset<sizeof(int) * CHAR_BIT> y(x);
    return (int)y.count();
}

//---------------------------------PATH INTEGRAL SUMMING-----------------------------------

/* Second simulation algorithm: summing hadamard possibilites
 Takes time T*2^h and space O(h) + O(n) [T = total # of gates, h = # of hadamards]
 UPDATE: takes less time now due to recursive restructuring (~2^h)
 Verbose will print the resulting state for *every* hadamard string (and the amplitude).
 
 V1: first version
 V2: changed to DFS procedure 
 V3: added out-of-reach path pruning 
 V4: added QFT: controlled-U gates, complex numbers, phase accumulation */

double pathStep(ifstream &in, streampos pos, int n, int startState, int currState, double currPhase, int endState, int changesLeft){
    char gate;
    int control, c1, c2, target, qubits = currState, oneFactor = 1, changeCounter = changesLeft;
    in.clear();
    in.seekg(pos); //move to correct position
    in >> control >> gate;
    
    streampos endPos;
    while (!in.eof()){
        switch (gate){
            case 'h':
            {
                changeCounter--;
                in >> target;
                endPos = in.tellg();
                int branchZero = qubits, branchOne = qubits;
                //|0><+| case: amp is always +1
                //|1><-| case: if the target qubit is a 1, amp turns negative; stays positive otherwise
                if (((branchOne >> (n - target - 1)) & 1) == 1) oneFactor = -1;
                branchZero = branchZero & ~(1 << (n - target - 1));
                branchOne = branchOne | (1 << (n - target - 1));
                
                if (bitDiff(qubits, endState) <= (changeCounter + 1)){
                    return 1/sqrt(2) * pathStep(in, endPos, n, startState, branchZero, currPhase, endState, changeCounter) + oneFactor/sqrt(2) * pathStep(in, endPos, n, startState, branchOne, currPhase, endState, changeCounter);
                } else return 0;
                break;
            }
            case 't':
            {
                changeCounter--;
                in >> c1 >> c2 >> target;
                int add = ((qubits >> (n - c1 - 1)) & 1) * ((qubits >> (n - c2 - 1)) & 1);
                qubits = qubits ^ (add << (n - target - 1));
                //if (bitDiff(qubits, endState) > changeCounter) return 0;
                break;
            }
            default: break;
        }
        in >> control >> gate;
    }
    //    cout << "RETURNED\n";
    if (qubits == endState){ //inner product <a|C|b> is 0 unless end state |a> matches start state |b>
        return currPhase;
    } else return 0;
}

complex<double> complexPathStep(ifstream &in, streampos pos, int n, int startState, int currState, complex<double> currPhase, int endState, int changesLeft){
    char gate;
    int control, c = -1, c1, c2, phasePow, target, qubits = currState, oneFactor = 1, changeCounter = changesLeft;
    in.clear();
    in.seekg(pos); //move to correct position
    in >> control >> gate;
    
    streampos endPos;
    while (!in.eof()){
        switch (gate){
            case 'h':
            {
                changeCounter--;
                in >> target;
                endPos = in.tellg();
                int branchZero = qubits, branchOne = qubits;
                //|0><+| case: amp is always +1
                //|1><-| case: if the target qubit is a 1, amp turns negative; stays positive otherwise
                if (((branchOne >> (n - target - 1)) & 1) == 1) oneFactor = -1;
                branchZero = branchZero & ~(1 << (n - target - 1));
                branchOne = branchOne | (1 << (n - target - 1));
                
                if (bitDiff(qubits, endState) <= (changeCounter + 1)){
                    return 1/sqrt(2) * complexPathStep(in, endPos, n, startState, branchZero, currPhase, endState, changeCounter) + oneFactor/sqrt(2) * complexPathStep(in, endPos, n, startState, branchOne, currPhase, endState, changeCounter);
                } else return 0;
                break;
            }
            case 't':
            {
                changeCounter--;
                in >> c1 >> c2 >> target;
                int add = ((qubits >> (n - c1 - 1)) & 1) * ((qubits >> (n - c2 - 1)) & 1);
                qubits = qubits ^ (add << (n - target - 1));
                //if (bitDiff(qubits, endState) > changeCounter) return 0;
                break;
            }
            case 'U':
            {
                in >> phasePow;
                complex<double> phase = polar(1.0, 1/pow(2, phasePow) * 2 * M_PI);
                if (control){
                    in >> c >> target;
                    if (((qubits >> (n - c - 1)) & 1) && ((qubits >> (n - target - 1)) & 1)) currPhase *= phase;
                } else {
                    in >> target;
                    if ((qubits >> (n - target - 1)) & 1) currPhase *= phase;
                }
                break;
            }
            default: break;
        }
        in >> control >> gate;
    }
    //    cout << "RETURNED\n";
    if (qubits == endState){ //inner product <a|C|b> is 0 unless end state |a> matches start state |b>
        return currPhase;
    } else return 0;
}

void pathIntegral(ifstream &in, int n, int startState, int endState, int numChanges, bool verbose, bool printMem, bool cmplx){
    //initial recursive call (the "root" of the path tree)
    if (cmplx){
        complex<double> amplitude = complexPathStep(in, ios::beg, n, startState, startState, 1, endState, numChanges);
        if (verbose) cout << "<" << binString(endState, n) << "|Circuit|" << binString(startState, n) << "> = " << amplitude << "\n";
        else cout << amplitude;
    } else {
        double amplitude = pathStep(in, ios::beg, n, startState, startState, 1, endState, numChanges);
        if (verbose) cout << "<" << binString(endState, n) << "|Circuit|" << binString(startState, n) << "> = " << amplitude << "\n";
        else cout << amplitude;
    }
    
    if (printMem){ //Print time/memory usage
        cout.precision(7);
        struct rusage usage;
        getrusage(RUSAGE_SELF, &usage);
        long totaluTime = (usage.ru_stime.tv_sec + usage.ru_utime.tv_sec) * 1000000 + usage.ru_stime.tv_usec + usage.ru_utime.tv_usec;
        double totalTime = totaluTime/ (double) 1000000;
        cout << "Runtime: " << totalTime << " seconds\n";
        cout << "Memory usage: " << usage.ru_maxrss / (double) memConst << " qunits [1 qunit ≈ 1 mb]\n\n";
    }
}

//--------------------------------------EXPERIMENT----------------------------------------

//-------------------------------------CONTROL PANEL---------------------------------------

int main(int argc, const char * argv[]){
    ifstream in ("/Users/AShi/Documents/Repos/QuantumLab/QuantumLab/gates.txt");
    srand((int)time(0)); //set rng seed
    int setting = 0, n = 3; //control panel setting and number of qubits
    /*0 = execute circuit
     1 = write QFT
     2 = write random circuit */
    
//    string Uset[] = {"t 0 1 2\nh 0\nh 1\nh 2\nt 0 1 2\nh 0\nh 1\nh 2\n","h 2\n","h 2\nt 1 2 0\nh 2\n","h 2\nt 0 2 1\nh 2\n","t 1 2 0\n","h 2\nt 1 2 0\nh 2\n","t 0 2 1\n"};
    int x;
    int woo = 0;
    complex<double> y;
    for (int i = 0; i < woo; i++){
        y = complex<double>(i, i) * complex<double>(i, i);
    }
    cout << "done2\n";
    for (int i = 0; i < woo; i++){
        x = i * i;
    }
    cout << "done1\n";
    
    switch (setting){
        case 0: //read and execute circuit from gates.txt
        {
//            cout.precision(10);
            cout << fixed;
            int start = 0, end = start;
//          start = rand()%(int)pow(2,n), end = start;
//            for (int i = 0; i < 8; i++){
//                for (int j = 0; j < 8; j++){
//                    pathIntegral(in, n, j, reverseBit(i, n), n, false, false);
//                    cout << " ";
//                }
//                cout << "\n\n";
//            }
            pathIntegral(in, n, start, end, (int)pow(n,2) + 2*n, true, true, false);
            //pathIntegral(in, n, start, reverseBit(end, n), n, true, true); //QFT
            break;

        }
        case 1: //write Quantum Fourier Transform on n qubits
        {
            ofstream out ("/Users/AShi/Documents/Repos/QuantumLab/QuantumLab/gates.txt");
            string circuit = writeQFT(n);
            out << circuit;
            out.close();
            break;
        }
        case 2: //write circuit into gates.txt
        {
            ofstream out ("/Users/AShi/Documents/Repos/QuantumLab/QuantumLab/gates.txt");
            int length = (int)pow(n,2);
            string circuit = writeCircuit(length, n, false);
            out << circuit;
            out.close();
            break;
        }
        default: break;
    }
    return 0;
}
