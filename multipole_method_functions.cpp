/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

#include "multipole_method_functions.h"

// Generate Morton index for a box center to use in M2L translation
void MultipoleMethodFunction::morton1(vec3<int> boxIndex3D, int& boxIndex, int numLevel) {
    int i, nx, ny, nz;
    boxIndex = 0;
    for (i = 0; i < numLevel; i++) {
        nx = boxIndex3D.x % 2;
        boxIndex3D.x >>= 1;
        boxIndex += nx * (1 << (3 * i + 1));

        ny = boxIndex3D.y % 2;
        boxIndex3D.y >>= 1;
        boxIndex += ny * (1 << (3 * i));

        nz = boxIndex3D.z % 2;
        boxIndex3D.z >>= 1;
        boxIndex += nz * (1 << (3 * i + 2));
    }
}

// Returns 3D box index from Morton index
void MultipoleMethodFunction::unmorton(int boxIndex, vec3<int>& boxIndex3D) {
    int i, j, k, n, mortonIndex3D[3];

    for (i = 0; i < 3; i++) mortonIndex3D[i] = 0;
    n = boxIndex;
    k = 0;
    i = 0;
    while (n != 0) {
        j = 2 - k;
        mortonIndex3D[j] += (n % 2) * (1 << i);
        n >>= 1;
        k = (k + 1) % 3;
        if (k == 0) i++;
    }
    boxIndex3D.x = mortonIndex3D[1];
    boxIndex3D.y = mortonIndex3D[2];
    boxIndex3D.z = mortonIndex3D[0];
}

// Generate Morton3 index for a box center
void MultipoleMethodFunction::morton3(vec3<int> boxIndex3D, int& boxIndex, int numLevel) {
    int i, nx, ny, nz;
    boxIndex = 0;
    for (i = 0; i < numLevel; i++) {
        nx = boxIndex3D.x % 3;
        boxIndex3D.x /= 3;
        boxIndex += nx * pow(3, 3 * i + 1);

        ny = boxIndex3D.y % 3;
        boxIndex3D.y /= 3;
        boxIndex += ny * pow(3, 3 * i);

        nz = boxIndex3D.z % 3;
        boxIndex3D.z /= 3;
        boxIndex += nz * pow(3, 3 * i + 2);
    }
}

// Returns 3D box index from Morton3 index
void MultipoleMethodFunction::unmorton3(int boxIndex, vec3<int>& boxIndex3D) {
    int i, j, k, n, powThree, mortonIndex3D[3];

    for (i = 0; i < 3; i++) mortonIndex3D[i] = 0;
    n = boxIndex;
    k = 0;
    i = 0;
    powThree = 1;
    while (n != 0) {
        j = 2 - k;
        mortonIndex3D[j] += (n % 3) * powThree;
        n /= 3;
        k = (k + 1) % 3;
        if (k == 0) powThree *= 3;
    }
    boxIndex3D.x = mortonIndex3D[1];
    boxIndex3D.y = mortonIndex3D[2];
    boxIndex3D.z = mortonIndex3D[0];
}

// Generate Morton11 index for a box center
void MultipoleMethodFunction::morton11(vec3<int> boxIndex3D, int& boxIndex, int numLevel) {
    int i, nx, ny, nz;
    boxIndex = 0;
    if (numLevel == 1) {
        boxIndex += boxIndex3D.x * 11;
        boxIndex += boxIndex3D.y;
        boxIndex += boxIndex3D.z * 121;
    } else {
        for (i = 0; i < numLevel; i++) {
            nx = boxIndex3D.x % 11;
            boxIndex3D.x /= 11;
            boxIndex += nx * pow(11, 3 * i + 1);

            ny = boxIndex3D.y % 11;
            boxIndex3D.y /= 11;
            boxIndex += ny * pow(11, 3 * i);

            nz = boxIndex3D.z % 11;
            boxIndex3D.z /= 11;
            boxIndex += nz * pow(11, 3 * i + 2);
        }
    }
}

// Returns 3D box index from Morton11 index
void MultipoleMethodFunction::unmorton11(int boxIndex, vec3<int>& boxIndex3D) {
    int i, j, k, n, powEleven, mortonIndex3D[3];

    for (i = 0; i < 3; i++) mortonIndex3D[i] = 0;
    n = boxIndex;
    k = 0;
    i = 0;
    powEleven = 1;
    while (n != 0) {
        j = 2 - k;
        if (i == 0) mortonIndex3D[j] += (n % 11);
        else mortonIndex3D[j] += (n % 11) * powEleven;
        n /= 11;
        k = (k + 1) % 3;
        if (k == 0) powEleven *= 11;
    }
    boxIndex3D.x = mortonIndex3D[1];
    boxIndex3D.y = mortonIndex3D[2];
    boxIndex3D.z = mortonIndex3D[0];
}

// Generate Morton index for a box center for a general box size
void MultipoleMethodFunction::mortonGeneral(vec3<int> boxIndex3D, int& boxIndex, int numLevel, int boxSize) {
    int i, nx, ny, nz;
    boxIndex = 0;
    for (i = 0; i < numLevel; i++) {
        nx = boxIndex3D.x % boxSize;
        boxIndex3D.x /= boxSize;
        boxIndex += nx * pow(boxSize, 3 * i + 1);

        ny = boxIndex3D.y % boxSize;
        boxIndex3D.y /= boxSize;
        boxIndex += ny * pow(boxSize, 3 * i);

        nz = boxIndex3D.z % boxSize;
        boxIndex3D.z /= boxSize;
        boxIndex += nz * pow(boxSize, 3 * i + 2);
    }
}

// Returns 3D box index from Morton index for a general box size
void MultipoleMethodFunction::unmortonGeneral(int boxIndex, vec3<int>& boxIndex3D, int boxSize) {
    int i, j, k, n, powBoxSize, mortonIndex3D[3];

    for (i = 0; i < 3; i++) mortonIndex3D[i] = 0;
    n = boxIndex;
    k = 0;
    i = 0;
    powBoxSize = 1;
    while (n != 0) {
        j = 2 - k;
        if (i == 0) mortonIndex3D[j] += (n % boxSize);
        else mortonIndex3D[j] += (n % boxSize) * powBoxSize;
        n /= boxSize;
        k = (k + 1) % 3;
        if (k == 0) powBoxSize *= boxSize;
    }
    boxIndex3D.x = mortonIndex3D[1];
    boxIndex3D.y = mortonIndex3D[2];
    boxIndex3D.z = mortonIndex3D[0];
}

void MultipoleMethodFunction::signedUnmorton3(unsigned int boxIndex, vec3<int>& boxIndex3D) {
    int i, j, k, n, powThree, mortonIndex3D[3];
    int xflip, yflip, zflip;
    xflip = 0;
    yflip = 0;
    zflip = 0;
    const unsigned int xCheck = 1 << 30;
    const unsigned int yCheck = 1 << 29;
    const unsigned int zCheck = 1 << 31;

    for (i = 0; i < 3; i++) mortonIndex3D[i] = 0;
    n = boxIndex;

    if ((n & xCheck) == xCheck) {
        xflip = 1;
        n &= ~xCheck;
    }
    if ((n & yCheck) == yCheck) {
        yflip = 1;
        n &= ~yCheck;
    }
    if ((n & zCheck) == zCheck) {
        zflip = 1;
        n &= ~zCheck;
    }

    k = 0;
    i = 0;
    powThree = 1;
    while (n != 0) {
        j = 2 - k;
        mortonIndex3D[j] += (n % 3) * powThree;
        n /= 3;
        k = (k + 1) % 3;
        if (k == 0) powThree *= 3;
    }

    if (xflip == 1) {
        boxIndex3D.x = -mortonIndex3D[1] - 1;
    } else {
        boxIndex3D.x = mortonIndex3D[1];
    }

    if (yflip == 1) {
        boxIndex3D.y = -mortonIndex3D[2] - 1;
    } else {
        boxIndex3D.y = mortonIndex3D[2];
    }

    if (zflip == 1) {
        boxIndex3D.z = -mortonIndex3D[0] - 1;
    } else {
        boxIndex3D.z = mortonIndex3D[0];
    }
}

void MultipoleMethodFunction::signedMorton3(vec3<int> boxIndex3D, unsigned int& boxIndex) {
    int i, nx, ny, nz;
    boxIndex = 0;
    i = 0;
    const unsigned int xCheck = 1 << 30;
    const unsigned int yCheck = 1 << 29;
    const unsigned int zCheck = 1 << 31;

    if (boxIndex3D.x < 0) {
        boxIndex3D.x = -boxIndex3D.x - 1;
        boxIndex |= xCheck;
    }
    if (boxIndex3D.y < 0) {
        boxIndex3D.y = -boxIndex3D.y - 1;
        boxIndex |= yCheck;
    }
    if (boxIndex3D.z < 0) {
        boxIndex3D.z = -boxIndex3D.z - 1;
        boxIndex |= zCheck;
    }

    while (boxIndex3D.x != 0 || boxIndex3D.y != 0 || boxIndex3D.z != 0) {
        nx = boxIndex3D.x % 3;
        boxIndex3D.x /= 3;
        boxIndex += nx * pow(3, 3 * i + 1);

        ny = boxIndex3D.y % 3;
        boxIndex3D.y /= 3;
        boxIndex += ny * pow(3, 3 * i);

        nz = boxIndex3D.z % 3;
        boxIndex3D.z /= 3;
        boxIndex += nz * pow(3, 3 * i + 2);

        i++;
        if (i > 6) {
            printf("Error: particle too far away.");
        }
    }
}

void MultipoleMethodFunction::balancedTernary(vec3<int> boxIndex3D, int& boxIndex){
    int nx, ny, nz, i;
    i = 0;
    boxIndex = 0;
    int xflip, yflip, zflip;
    
    if (boxIndex3D.x < 0){
        xflip = 1;
        boxIndex3D.x = -boxIndex3D.x;
    } else {
        xflip = 0;
    }
    
    if (boxIndex3D.y < 0){
        yflip = 1;
        boxIndex3D.y = -boxIndex3D.y;
    }
    else {
        yflip = 0;
    }
    
    if (boxIndex3D.z < 0){
        zflip = 1;
        boxIndex3D.z = -boxIndex3D.z;
    }
    else {
        zflip = 0;
    }
    
    while (i < 6  && (boxIndex3D.x != 0 || boxIndex3D.y != 0 || boxIndex3D.z != 0)) {
        
        nx = boxIndex3D.x % 3;
        nx += xflip * ((nx == 1) - (nx == 2));
        boxIndex3D.x += (1 - xflip) * (nx == 2) + xflip * (nx == 1);
        boxIndex3D.x /= 3;
        boxIndex += nx * pow(3, 3 * i + 1);

        ny = boxIndex3D.y % 3;
        ny += yflip * ((ny == 1) - (ny == 2));
        boxIndex3D.y += (1 - yflip) * (ny == 2) + yflip * (ny == 1);
        boxIndex3D.y /= 3;
        boxIndex += ny * pow(3, 3 * i);

        nz = boxIndex3D.z % 3;
        nz += zflip * ((nz == 1) - (nz == 2));
        boxIndex3D.z += (1 - zflip) * (nz == 2) + zflip * (nz == 1);
        boxIndex3D.z /= 3;
        boxIndex += nz * pow(3, 3 * i + 2);
        
        i++;
    }
}

void MultipoleMethodFunction::undoBalancedTernary(int boxIndex, vec3<int>& boxIndex3D){
    
    int i, j, k, n, powBoxSize, mortonIndex3D[3];
    int remainder;

    for (i = 0; i < 3; i++) mortonIndex3D[i] = 0;
    n = boxIndex;
    k = 0;
    i = 0;
    powBoxSize = 1;
    while (n != 0) {
        j = 2 - k;
        remainder = (n % 3);
        if (remainder == 2){
            mortonIndex3D[j] -= powBoxSize;
        }
        else if (remainder == 1){
            mortonIndex3D[j] += powBoxSize;
        }
        n /= 3;
        k = (k + 1) % 3;
        if (k == 0) powBoxSize *= 3;
    }
    
    boxIndex3D.x = mortonIndex3D[1];
    boxIndex3D.y = mortonIndex3D[2];
    boxIndex3D.z = mortonIndex3D[0];
}