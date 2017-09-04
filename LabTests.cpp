/*  Algorithms for simulating general quantum circuits
    First simulation algorithm: tracking amplitudes (algorithmOne)
    Second simulation algorithm: summing over hadamard possibilities (algorithmTwo)
    Third simulation algorithm: recursively tracing through layers (algorithmThree) (due to Aaronson/Chen)

    Note: In a binary representation, qubits 0 through n-1 are represented from leftmost
    digit to rightmost digit (i.e. 6 = 110 = 0th qubit:1, 1st qubit:1, 2nd qubit:0)

    QuantumLab
    Created by Andrew Shi on 7/8/17.
    Copyright © 2017 Andrew Shi. All rights reserved.
*/
#include <iostream>
#include <sstream>
#include <complex>
#include <fstream>
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <sys/time.h>
#include <sys/resource.h>
#include <bitset>
#include <climits>
#define MAX_LAYERS 1000

using namespace std;
double amps[INT_MAX]; //used for amplitude storage
long memConst = 1024 * 1024; //constant for interpreting rusage's memory output
bool bitReached[32];
string layerGates[MAX_LAYERS];
int layerMarks[MAX_LAYERS];
int layers[MAX_LAYERS];
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

string randToff(int *indices, int n){ //generates a random toffoli gate for n qubit circuit
    if (n < 3) {
        cout << "not enough qubits\n";
        return "fail";
    }

    int i = 0, temp;
    string result = "t";
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
                out = out + "h " + to_string(rand() % n);
            } else {
                out += randToff(toff, n);
            }
            out = out + "\n";
        }
    } else {
        for (int i = 0; i < n; i++){
            out = out + "h " + to_string(i) + "\n";
        }
        for (int i = 0; i < length; i++){
            out = out + randToff(toff, n) + "\n";
        }
        for (int i = 0; i < n; i++){
            out = out + "h " + to_string(i) + "\n";
        }
    }
    return out;
}

void resetCounter(bool *reached, int n){
    for (int i = 0; i < n; i++){
        reached[i] = false;
    }
}

void printLayers(ifstream &in, int *layers, int depth){
    in.clear();
    in.seekg(ios::beg);
    int gateCount = 0;
    string gate;
    for (int i = 0; i < depth; i++){
        cout << "Layer " << i << ":\n";
        while (gateCount < layers[i]){
            getline(in, gate);
            cout << gate << "\n";
            gateCount++;
        }
        cout << "\n";
    }
}

int bitDiff(int a, int b){
    int x = a ^ b;
    std::bitset<sizeof(int) * CHAR_BIT> y(x);
    return (int)y.count();
}

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

//---------------------------------PATH INTEGRAL SUMMING-----------------------------------

/* Second simulation algorithm: summing hadamard possibilites
 Takes time T*exp(O(h)) and space O(h) + O(n) [T = total # of gates, h = # of hadamards]
 UPDATE: takes less time now due to recursive restructuring

 This algorithm has the same parameter definitions as before; chosenState is the state for which an amplitude will be recorded (aka startState).
 Verbose will print the resulting state for *every* hadamard string (and the amplitude).

 V1: first version
 V2: changed to DFS procedure
 V3: added out-of-reach path pruning */

double algTwoRecur(ifstream &in, streampos pos, int n, int startState, int currState, int gatesLeft, bool verbose){
    char gate;
    int c1, c2, target, qubits = currState, oneFactor = 1, gateCounter = gatesLeft;
    in.clear();
    in.seekg(pos); //move to correct position
    in >> gate;

    streampos endPos;
    while (!in.eof()){
        gateCounter--;
        switch (gate){
            case 'h':
            {
                in >> target;
                endPos = in.tellg();
                int branchZero = qubits, branchOne = qubits;
                //|0><+| case: amp is always +1
                //|1><-| case: if the target qubit is a 1, amp turns negative; stays positive otherwise
                if (((branchOne >> (n - target - 1)) & 1) == 1) oneFactor = -1;
                branchZero = branchZero & ~(1 << (n - target - 1));
                branchOne = branchOne | (1 << (n - target - 1));

                if (bitDiff(qubits, startState) <= (gateCounter + 1)){
                    return 1/sqrt(2) * algTwoRecur(in, endPos, n, startState, branchZero, gateCounter, verbose) + oneFactor/sqrt(2) * algTwoRecur(in, endPos, n, startState, branchOne, gateCounter, verbose);
                } else return 0;
                break;
            }
            case 't':
            {
                in >> c1 >> c2 >> target;
                int add = ((qubits >> (n - c1 - 1)) & 1) * ((qubits >> (n - c2 - 1)) & 1);
                qubits = qubits ^ (add << (n - target - 1));
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
    if (qubits == startState){ //inner product <a|C|b> is 0 unless end state |a> matches start state |b>
        return 1;
    } else return 0;
}

void algorithmTwo(ifstream &in, int n, int startState, int numGates, bool verbose){
    double amplitude = algTwoRecur(in, ios::beg, n, startState, startState, numGates, verbose);
    cout << "<" << binString(startState, n) << "|Circuit|" << binString(startState, n)
    << "> = " << amplitude << "\n";

    //Print time/memory usage
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    long totaluTime = (usage.ru_stime.tv_sec + usage.ru_utime.tv_sec) * 1000000 + usage.ru_stime.tv_usec + usage.ru_utime.tv_usec;
    double totalTime = totaluTime/ (double) 1000000;
    cout << "Runtime: " << totalTime << " seconds\n";
    cout << "Memory usage: " << usage.ru_maxrss / (double) memConst << " qunits [1 qunit ≈ 1 mb]\n\n";
}

void oldalgorithmTwo(ifstream &in, int n, int startState, bool verbose){
    char gate;
    int c1, c2, target, qubits, hCount = 0, h;
    double amplitude = 0, temp;
    in.clear();
    in.seekg(0, ios::beg); //move to beginning of file
    in >> gate;

    int counter = 0;

    while (!in.eof()){ //count the number of total hadamards in hCount
        switch (gate){
            case 'h':
            {
                in >> target;
                hCount++;
            }
            case 't':
            {
                in >> c1 >> c2 >> target;
            }
            default:
            {
            cout << "Incompatible gate type: " << gate << "\n";
            break;
            }
        }
        in >> gate;
    }
    in.clear();
    in.seekg(0, ios::beg);

    for (int i = 0; i < (int)pow(2, hCount); i++){ //enumerate over all hadamard branches
        qubits = startState, temp = 1, h = 0;
        in >> gate;

        while (!in.eof()){
            switch (gate){
                case 'h': //hadamard gate
                {
                    in >> target;
                    temp *= 1/sqrt(2);
                    if (((i >> (hCount - h - 1)) & 1) == 1){ //hadamard bit at position h is 1
                        //if the target qubit is a 1, amp turns negative; stays positive otherwise
                        if (((qubits >> (n - target - 1)) & 1) == 1) temp *= -1;
                        //map the target bit to 1
                        qubits = qubits | (1 << (n - target - 1));
                    } else { //hadamard bit at pos is 0
                        //map target bit to 0
                        qubits = qubits & ~(1 << (n - target - 1));
                    }
                    h++; //increase h index
                    break;
                }
                case 't': //toffoli gate
                {
                    in >> c1 >> c2 >> target;
                    int add = ((qubits >> (n - c1 - 1)) & 1) * ((qubits >> (n - c2 - 1)) & 1);
                    qubits = qubits ^ (add << (n - target - 1));
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
        if (qubits == startState){ //inner product <a|C|b> is 0 unless end state |a> matches start state |b>
            amplitude += temp; //accumulate amplitude of desired state
            counter++;
        }
        //reset file reader
        in.clear();
        in.seekg(0, ios::beg);
        if (verbose){
            cout << "Hadamard string " << binString(i, hCount) << " produced " << binString(qubits, n) << " with amplitude " << temp << "\n";

        }
    }
    cout << "<" << binString(startState, n) << "|Circuit|" << binString(startState, n)
    << "> = " << amplitude << "\n";

    //    cout << "Count: " << counter << "\n";

    //Print time/memory usage
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    long totaluTime = (usage.ru_stime.tv_sec + usage.ru_utime.tv_sec) * 1000000 + usage.ru_stime.tv_usec + usage.ru_utime.tv_usec;
    double totalTime = totaluTime/ (double) 1000000;
    cout << "Runtime: " << totalTime << " seconds\n";
    cout << "Memory usage: " << usage.ru_maxrss / (double) memConst << " qunits [1 qunit ≈ 1 mb]\n\n";
}

//----------------------------------AARONSON RECURSION-------------------------------------

/* This is the algorithm described in Aaronson/Chen's paper (arXiv:1612.05903 [quant-ph]), section 4. It simulates a quantum circuit with a recursive procedure in time O(n*(2d)^(n+1)) and space O(nlog(d)). d is the circuit depth, or the number of gate groups applied chronologically where each gate group only acts on a qubit 0 or 1 times. (effectively we assume d ~ T, the total # of gates).
 UPDATE: takes less time now due to some adjustments

 There is a modified tradeoff algorithm (also from Aaronson/Chen) exchanging time and space resources where for a chosen integer k, the algorithm takes time O(n*2^(n+1)*d^(k+1)) and space O(2^(n-k)logd). However, note that algorithmThree is approximately a special case where k = n, and when k = 0 both space and time are of complexity exp(n).

 The tradeoff algorithm is not implemented here because despite the extra flexibility, time is the limiting factor preventing scaling to more qubits in all three algorithms and is still lower bounded by O(nexp(n)) in the tradeoff, a larger complexity than even algorithm one.

 V1: first version
 V2: added zero-term checking
 V3: added out-of-reach path pruning */

double algThreeRecur(ifstream &in, int n, int beginD, int endD, int startS, int endS, int *layers, bool verbose){ //Recursive subalgorithm for algorithm three
    double result = 0;
    if (verbose) cout << beginD << "(" << startS << ") to " << endD << "(" << endS << ")\n";
    if (beginD == endD) { //base case
        result = 1;
        in.clear();
        in.seekg(ios::beg);
        stringstream gates;
        gates << layerGates[beginD];

        char gate;
        int c1, c2, target, endBit;
        int qubits = startS;
        gates >> gate;
        while (!gates.eof()){ //act on gates inside this layer
            switch (gate){
                case 'h': //hadamard
                {
                    gates >> target;
                    endBit = (endS >> (n - target - 1)) & 1;
                    if ((((startS >> (n - target - 1)) & 1) == 1) && endBit == 1){
                        //if the hadamarded bit is 1 in both the start and end states, the end amplitude will be negative
                        result *= -1;
                    }
                    //set the hadamarded bit in the register (qubits) to whatever it equals in the end state
                    qubits = qubits & ~(1 << (n - target - 1));
                    qubits += (endBit << (n - target - 1));
                    result *= 1/sqrt(2); //hadamard-adjusted amplitude
                    break;
                }
                case 't': //toffoli
                {
                    gates >> c1 >> c2 >> target;
                    int add = ((startS >> (n - c1 - 1)) & 1) * ((startS >> (n - c2 - 1)) & 1);
                    qubits = qubits ^ (add << (n - target - 1)); //update register
                    break;
                }
                default:
                {
                    cout << "Incompatible gate type: " << gate << "\n";
                    break;
                }
            }
            gates >> gate;
        }
        // <endS|qubits> == 0 if endS ≠ qubits [|qubits> = Layer|startS>]
        if (qubits != endS) result = 0;
        return result;
    } else { //recursive case
        /* Compute <y|C|x> by summing all <y|C_1|i><i|C_2|x> for i = {0,1}^n.
         Compute the two sub terms recursively. */
        for (int i = 0; i < pow(2, n); i++){
            if (bitDiff(startS, i) <= (layers[(beginD + endD)/2 + 1] - layers[beginD]) &&
                bitDiff(i, endS) <= (layers[endD + 1] - layers[(beginD + endD)/2 + 1])){
                double termOne = algThreeRecur(in, n, beginD, (beginD + endD)/2, startS, i, layers, verbose);
                if (termOne != 0){ //only compute second term if first is nonzero
                result += (termOne * algThreeRecur(in, n, (beginD + endD)/2 + 1, endD, i, endS, layers, verbose));
                }
            }
        }
    }
    return result;
}

void algorithmThree(ifstream &in, int n, int startState, int endState, bool verbose){
    char gate;
    int c1, c2, target, depth = 0, gCount = 0;
    string gateBuffer = "";

    in.clear();
    in.seekg(ios::beg);
    in >> gate;
    resetCounter(bitReached, n);
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
                    resetCounter(bitReached, n);
                    gateBuffer = "";
                }
                bitReached[target] = true;
                gateBuffer = gateBuffer + "h " + to_string(target) + " ";
                break;
            }
            case 't':
            {
                in >> c1 >> c2 >> target;
                if (bitReached[c1] || bitReached[c2] || bitReached[target]){
                    layerGates[depth] = gateBuffer;
                    depth++;
                    layers[depth] = gCount;
                    resetCounter(bitReached, n);
                    gateBuffer = "";
                }
                bitReached[c1] = true;
                bitReached[c2] = true;
                bitReached[target] = true;
                gateBuffer = gateBuffer + "t " + to_string(c1) + " " + to_string(c2) + " " +to_string(target) + " ";
                break;
            }
            default:
            {
                cout << "Incompatible gate type: " << gate << "\n";
                break;
            }
        }
        in >> gate;
        gCount++;
    }
    layerGates[depth] = gateBuffer;
    depth++;
    layers[depth] = gCount;
    cout << "Divided into " << depth << " layers\n";

    double result = algThreeRecur(in, n, 0, depth - 1, startState, endState, layers, verbose); //call recursive algorithm
    cout << "<" << binString(endState, n) << "|Circuit|" << binString(startState, n) << ">: " << result << "\n";

    //Print time/memory usage
    struct rusage usage;
    getrusage(RUSAGE_SELF, &usage);
    long totaluTime = (usage.ru_stime.tv_sec + usage.ru_utime.tv_sec) * 1000000 + usage.ru_stime.tv_usec + usage.ru_utime.tv_usec;
    double totalTime = totaluTime/ (double) 1000000;
    cout << "Runtime: " << totalTime << " seconds\n";
    cout << "Memory usage: " << usage.ru_maxrss / (double) memConst << " qunits [1 qunit ≈ 1 mb]\n\n";
}

//-------------------------------------CONTROL PANEL---------------------------------------

int main(int argc, const char * argv[]){
    ifstream in ("/Users/AShi/Documents/Repos/QuantumLab/QuantumLab/gates.txt");
    srand((int)time(0)); //set rng seed
    int setting = 0;
    /*0 = execute circuit
     1 = write circuit
     2 = special mode (resource testing!)
     3 = special mode no. 2 (Uset!! from Uncle's paper [used in proof of Thm. 3.2]) */

    int n = 10, length = pow(n,2); //number of qubits and length of circuit
    int r = 0; //starting state
    r = rand() % (int)pow(2,n); //random start state for testing, can be commented out
    string Uset[] = {"t 0 1 2\nh 0\nh 1\nh 2\nt 0 1 2\nh 0\nh 1\nh 2\n",
        "h 2\n","h 2\nt 1 2 0\nh 2\n","h 2\nt 0 2 1\nh 2\n","t 1 2 0\n",
        "h 2\nt 1 2 0\nh 2\n","t 0 2 1\n"};

    switch (setting){
        case 0: //read and execute circuit from gates.txt
        {
            cout.precision(10);
            algorithmOne(in, n, r, false);
            algorithmTwo(in, n, r, 16, false);
//            oldalgorithmTwo(in, n, r, false);
//            algorithmThree(in, n, r, r, false);
            break;

        }
        case 1: //write circuit into gates.txt
        {
            ofstream out ("/Users/AShi/Documents/Repos/QuantumLab/QuantumLab/gates.txt");
            string circuit = writeCircuit(length, n, false);
            out << circuit;
            out.close();
            break;
        }
        case 2:
        {
            ofstream out ("/Users/AShi/Documents/Repos/QuantumLab/QuantumLab/gates.txt");
            int maxCircuit = 25;
            for (int i = 3; i <= maxCircuit; i++){
                string circuit = writeCircuit(length, i, false);
                out << circuit;
                out.close();
//                r = rand() % (int)pow(2,i);
                algorithmOne(in, i, r, false);
//                algorithmTwo(in, i, r, false);
                if (i != maxCircuit) out.open("/Users/AShi/Documents/Repos/QuantumLab/QuantumLab/gates.txt", ofstream::out | ofstream::trunc);
            }
            break;
        }
        case 3:
        {
            //write a circuit of [length] number of U-gates, from the set defined above
            ofstream out ("/Users/AShi/Documents/Repos/QuantumLab/QuantumLab/gates.txt");
            for (int i = 0; i < length; i++){
                out << Uset[rand() % 7];
            }
            out.close();

            n = 3;
            algorithmOne(in, n, r, true);
//            algorithmTwo(in, n, r, false);
            break;
        }
        default: break;
    }
    return 0;
}
