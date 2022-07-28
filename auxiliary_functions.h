/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   auxiliary_functions.h
 * Author: Marshall
 *
 * Created on June 7, 2022, 5:04 PM
 */

#ifndef AUXILIARY_FUNCTIONS_H
#define AUXILIARY_FUNCTIONS_H

#include <cstdlib>
#include <cstring>
#include <iostream>
#include <cmath>
#include "structures.h"
#include "constants.h"

#define SIGNUM(x) ((x > 0) - (x < 0))
#define ODDEVEN(n) (((n & 1) == 1) ? -1 : 1)

void fileOpenError(char* errorFileName, int lineNumber, const char* fileName);
int returnSize(char fileName[]);
char* concat(const char* s1, const char* s2);

// Outlier detection functions
template<typename T>
int partition5(T* input, int left, int right);
template<typename T>
int pivot(T* input, int left, int right);
template<typename T>
int partition(T* input, int left, int right, int pivotIndex, int n);
template<typename T>
int select(T* input, int left, int right, int n);
template<typename T>
float findTrueMedian(T* input, int inputLength);

#include "outlier_detection.h"

//bool substepBoundaryConditionsSorter_Ascending (SubstepBoundaryConditions lhs, 
//        SubstepBoundaryConditions rhs);
//
//bool substepBoundaryConditionsSorter_Descending (SubstepBoundaryConditions lhs, 
//        SubstepBoundaryConditions rhs);

int compareFloat(const void* a, const void* b);

bool compareOutlierInfo(outlierInfo a, outlierInfo b);

inline int quadrantRotationFunction(int n){
    return (SIGNUM(1 - ODDEVEN(n)) * ODDEVEN((n - 1) >> 1));
};
std::complex<double> qrf_complex(int n);

int gcd3(int a, int b, int c);
int gcd(int a, int b);
void cart2sph(double& r, double& theta, double& phi, double dx, double dy, double dz);
void cart2sph_exact(double& r, double& theta, double& phi, double dx, double dy, double dz);

#endif /* AUXILIARY_FUNCTIONS_H */

