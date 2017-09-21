//
//  savitch.hpp
//  PocketSimulator
//
//  Created on 9/20/17.
//  Copyright © 2017. All rights reserved.
//

#ifndef savitch_hpp
#define savitch_hpp

#include <stdio.h>

double savitchRecur(int N, int beginD, int endD, int startS, int endS, int *layers, bool verbose); //Recursive subalgorithm for algorithm three

void savitch(string gatePath, int N, int startState, int endState, bool verbose);

#endif /* savitch_hpp */
