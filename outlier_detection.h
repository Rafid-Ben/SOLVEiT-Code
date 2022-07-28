/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   outlier_detection.h
 * Author: Marshall
 *
 * Created on January 21, 2022, 4:51 PM
 */

#ifndef OUTLIER_DETECTION_H
#define OUTLIER_DETECTION_H

/*Note that this needs to be a separate file, due to using templates.
 * If try to write function definitions in .cpp file, then the compiler won't
 * like the templates and will say that there are undefined reference errors.
 */

template<typename T>
int partition5(T* input, int left, int right) {
    int i, j;
    T temp;
    i = left + 1;
    while (i <= right) {
        j = i;
        while ((j > left) && (input[j - 1] > input[j])) {
            temp = input[j];
            input[j] = input[j - 1];
            input[j - 1] = temp;
            j--;
        }
        i++;
    }
    return ((left + right) / 2);
}

template<typename T>
int pivot(T* input, int left, int right) {
    if (right - left < 5) {
        return partition5(input, left, right);
    }
    int i, subRight, median5;
    T temp;
    for (i = left; i < right; i += 5) {
        subRight = i + 4;
        if (subRight > right) {
            subRight = right;
        }
        median5 = partition5(input, i, subRight);
        temp = input[left + (i - left) / 5];
        input[left + (i - left) / 5] = input[median5];
        input[median5] = temp;
    }
    int mid = (right - left) / 10 + left + 1;
    return select(input, left, left + (right - left) / 5, mid);
}

template<typename T>
int partition(T* input, int left, int right, int pivotIndex, int n) {
    T pivotValue = input[pivotIndex];
    T temp = input[right];
    input[right] = input[pivotIndex];
    input[pivotIndex] = temp;
    int storeIndex = left;
    int i;
    for (i = left; i < right; i++) {
        if (input[i] < pivotValue) {
            temp = input[i];
            input[i] = input[storeIndex];
            input[storeIndex] = temp;
            storeIndex++;
        }
    }
    int storeIndexEq = storeIndex;
    for (i = storeIndex; i < right; i++) {
        if (input[i] == pivotValue) {
            temp = input[i];
            input[i] = input[storeIndexEq];
            input[storeIndexEq] = input[i];
            storeIndexEq++;
        }
    }
    temp = input[right];
    input[right] = input[storeIndexEq];
    input[storeIndexEq] = temp;
    if (n < storeIndex) {
        return storeIndex;
    }
    if (n <= storeIndexEq) {
        return n;
    }
    return storeIndexEq;
}

template<typename T>
int select(T* input, int left, int right, int n) {
    int pivotIndex;
    do {
        if (left == right) {
            return left;
        }
        pivotIndex = pivot(input, left, right);
        pivotIndex = partition(input, left, right, pivotIndex, n);
        if (n == pivotIndex) {
            return n;
        } else if (n < pivotIndex) {
            right = pivotIndex - 1;
        } else {
            left = pivotIndex + 1;
        }
    } while (true);
}

template<typename T>
float findTrueMedian(T* input, int inputLength) {
    if ((inputLength & 1) == 1) {
        int medIndex = select(input, 0, inputLength - 1, inputLength / 2);
        float median = input[medIndex];
        return median;
    } else {
        int lowerIndex = select(input, 0, inputLength - 1, (inputLength - 1) / 2);
        T lower = input[lowerIndex];
        // Alternatively, should just be able to find the minimum of the upper half of the table
        int upperIndex = select(input, inputLength / 2, inputLength - 1, inputLength / 2);
        T upper = input[upperIndex];


        return (((float) lower + (float) upper) / 2);
    }
}


#endif /* OUTLIER_DETECTION_H */

