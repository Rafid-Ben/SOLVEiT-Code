/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "auxiliary_functions.h"

/**
 * Prints error if file cannot be opened
 * 
 * @param errorFileName Name of file
 * @param lineNumber Line number
 * @param fileName File where error occurred
 */
void fileOpenError(char* errorFileName, int lineNumber, const char* fileName) {
    printf("Can't open %s file on line %d in file %s.\n", errorFileName, lineNumber, fileName);
    exit(1);
}

/**
 * Concatenate two const char* strings and return the result
 * @param s1 First string
 * @param s2 Second string
 * @return Result of concatenation
 */
char* concat(const char* s1, const char* s2) {
    char* result = (char*)malloc(strlen(s1) + strlen(s2) + 1);
    strcpy(result, s1);
    strcat(result, s2);
    return result;
}

/**
 * Returns the number of scattered points for building interpolation scheme
 * @param fileName Name of text file containing scattered points
 * @return Number of nodes
 */
int returnSize(char fileName[]) {
    FILE *fp = fopen(fileName, "r");
    if (!fp) {
        printf("Can't open file\n");
        return 0;
    }

    char buf[1024]; //allocating bits to read

    int node_num2 = 0; //start line counter to 0
    while (fgets(buf, 1024, fp)) {
        node_num2++;
    }
    return node_num2;
}

/**
 * Compares two floats to determine which is larger (carry-over from C function)
 * @param a - first float value
 * @param b - second float value
 * @return 1 if a > b, -1 if a < b, 0 otherwise
 */
int compareFloat (const void * a, const void * b) {
    if (*(float*)a > *(float*)b) return 1;
    else if (*(float*)a < *(float*)b) return -1;
    else return 0;
}

/**
 * Compares the signedMorton value for two outlierInfo structs
 * Used to sort array of outlierInfo structs
 * @param a - first outlierInfo value
 * @param b - second outlierInfo value
 * @return True if a < b, False otherwise
 */
bool compareOutlierInfo(outlierInfo a, outlierInfo b) {
    return (a.signedMorton < b.signedMorton);
}

/**
 * Computes a complex number indicating rotation given an input quadrant
 * @param n - quadrant
 * @return Complex number to use to rotate to/from desired quadrant
 */
std::complex<double> qrf_complex(int n){
    int c = quadrantRotationFunction(n + 1);
    int s = quadrantRotationFunction(n);
    std::complex<double> r(c, s);
    return r;
}

/**
 * Returns the Greatest Common Divisor of three integer values
 * @param a - first integer value
 * @param b - second integer value
 * @param c - third integer value
 * @return Greatest Common Divisor of three inputs
 */
int gcd3(int a, int b, int c){
    return (gcd(a, gcd(b, c)));
}

/**
 * Returns the Greatest Common Divisor of two integer values
 * @param a - first integer value
 * @param b - second integer value
 * @return Greatest Common Divisor of two inputs
 */
int gcd(int a, int b){
    int temp;
    a = abs(a);
    b = abs(b);
    while (b != 0){
        temp = b;
        b = a % b;
        a = temp;
    }
    return a;
}

/**
 * Convert Cartesian to spherical coordinates, using the softening parameter when
 * calculating the distance
 * @param r - output r coordinate
 * @param theta - output theta coordinate
 * @param phi - output phi coordinate
 * @param dx - input x coordinate
 * @param dy - input y coordinate
 * @param dz - input z coordinate
 */
void cart2sph(double& r, double& theta, double& phi, double dx, double dy, double dz) {
    r = sqrtf(dx * dx + dy * dy + dz * dz + PhysicalConstant::SOFTENING);
//    r = sqrtf(dx * dx + dy * dy + dz * dz) + eps;
    theta = acosf(dz / r);
    if (theta == 0 && dz != 0){
        // First order approximation of arctan(sqrt(dx * dx + dy * dy) / dz);
        theta = sqrtf(dx * dx + dy * dy) / dz;
    }
    if (fabs(dx) + fabs(dy) < PhysicalConstant::SOFTENING) {
        phi = 0;
    }
    else if (fabs(dx) < PhysicalConstant::SOFTENING) {
        phi = dy / fabs(dy) * PhysicalConstant::PI * 0.5;
    }
    else if (dx > 0) {
        phi = atanf(dy / dx);
    }
    else {
        phi = atanf(dy / dx) + PhysicalConstant::PI;
    }
}

/**
 * Convert Cartesian to spherical coordinates, using the exact relations
 * @param r - output r coordinate
 * @param theta - output theta coordinate
 * @param phi - output phi coordinate
 * @param dx - input x coordinate
 * @param dy - input y coordinate
 * @param dz - input z coordinate
 */
void cart2sph_exact(double& r, double& theta, double& phi, double dx, double dy, double dz){
    r = sqrtf(dx * dx + dy * dy + dz * dz);
    if (r != 0){
        theta = acosf(dz / r);
        if (theta == 0 && dz != 0){
            // First order approximation of arctan(sqrt(dx * dx + dy * dy) / dz);
            theta = sqrtf(dx * dx + dy * dy) / dz;
        }
    }
    else{
        //theta = pi / 2;
        // Note 100% sure about this case, but it results in the correct rotation matrix, which is the only place it's used
        theta = 0;
    }
    if ((dx == 0) && (dy == 0)){
        phi = 0;
    }
    else if (dx == 0){
        phi = dy / fabs(dy) * PhysicalConstant::PI * 0.5;
    }
    else if (dx > 0){
        phi = atanf(dy / dx);
    }
    else{
        phi = atanf(dy / dx) + PhysicalConstant::PI;
    }
}
