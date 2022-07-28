/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   snapshot.h
 * Author: Marshall
 *
 * Created on January 23, 2022, 10:01 AM
 */

#ifndef SNAPSHOT_H
#define SNAPSHOT_H

#include <cassert>
#include <cstdlib>
#include <cmath>
#include <cstring>

template<typename T>
class Array {
public:
    T* vals;
    size_t maxElements;
    size_t numElements;
    Array(size_t numElements, size_t maxElements) : numElements{numElements},
            maxElements{maxElements} {
        vals = new T[maxElements];
        this->setZero(0, maxElements);
    }
            
    Array(size_t n = 0) : numElements{n}, maxElements{n} {
        vals = new T[maxElements];
        this->setZero(0, maxElements);
    }
    ~Array(){
        delete[] vals;
    };
    
    void setNull(size_t firstElement, size_t lastElement){
        for (int i = firstElement; i < lastElement; i++){
//            this->vals[i] = NULL;
        }
    }
    void setZero(size_t firstElement, size_t lastElement){
        for (int i = firstElement; i < lastElement; i++){
            this->vals[i] = 0;
        }
    }
    
    void changeNumElements(size_t n){
        if (n <= this->maxElements){
            this->setZero(this->numElements, n);
            this->numElements = n;
        }
        else{
            size_t newSize = (n > 2 * this->maxElements) ? n : 2 * this->maxElements;
            T* temp = new T[newSize];
            for (int i = this->numElements; i < newSize; i++){
                temp[i] = 0;
            }
            std::memcpy(temp, this->vals, this->numElements * sizeof(T));
            this->numElements = n;
            this->maxElements = newSize;
            delete[] this->vals;
            this->vals = temp;
        }
    }
    
    void changeMaxElements(size_t n){
        T* temp = new T[n];
        for (int i = this->numElements; i < n; i++){
            temp[i] = 0;
        }
        this->numElements = (this->numElements > n) ? n : this->numElements;
        this->maxElements = n;
        std::memcpy(temp, this->vals, this->numElements * sizeof(T));
        delete[] this->vals;
        this->vals = temp;
    }
    
//    T& operator[] (int i) {
//        assert(i >= 0 && i < this->numElements);
//        return this->vals[i];
//    }
    
    T& operator[] (int i) const{
        assert(i >= 0 && i < this->numElements);
        return this->vals[i];
    }
    
    template<typename S>
    Array<T>& operator=(S s){
        for (int i = 0; i < this->numElements; i++){
            this->vals[i] = s;
        }
        return *this;
    }
    
    Array<T>& operator=(const Array<T>& other){
        if (this != &other){
//            assert(this->numElements == other.numElements);
//            delete[] this->vals;
            this->changeMaxElements(other.maxElements);
            this->changeNumElements(other.numElements);
            for (int i = 0; i < this->numElements; i++){
                this->vals[i] = other[i];
            }
        }
        return *this;
    }
    
    Array<T>(const Array<T> &other){
        if (this != &other){
            this->numElements = other.numElements;
            this->maxElements = other.maxElements;
//            delete[] this->vals;
            this->vals = new T[this->maxElements];
            for (int i = 0; i < this->numElements; i++){
                this->vals[i] = other[i];
            }
        }
    }
    
    Array<T>& operator+=(Array<T> o){
        assert(this->numElements == o.numElements);
        for (int i = 0; i < this->numElements; i++){
            this->vals[i] += o[i];
        }
        return *this;
    }
    Array<T>& operator-=(Array<T> o){
        assert(this->numElements == o.numElements);
        for (int i = 0; i < this->numElements; i++){
            this->vals[i] -= o[i];
        }
        return *this;
    }
    Array<T>& operator*=(Array<T> o){
        assert(this->numElements == o.numElements);
        for (int i = 0; i < this->numElements; i++){
            this->vals[i] *= o[i];
        }
        return *this;
    }
    Array<T>& operator/=(Array<T> o){
        assert(this->numElements == o.numElements);
        for (int i = 0; i < this->numElements; i++){
            this->vals[i] /= o[i];
        }
        return *this;
    }
    
    template<typename S>
    Array<T>& operator+=(S s){
        for (int i = 0; i < this->numElements; i++){
            this->vals[i] += s;
        }
        return *this;
    }
    template<typename S>
    Array<T>& operator-=(S s){
        for (int i = 0; i < this->numElements; i++){
            this->vals[i] -= s;
        }
        return *this;
    }
    template<typename S>
    Array<T>& operator*=(S s){
        for (int i = 0; i < this->numElements; i++){
            this->vals[i] *= s;
        }
        return *this;
    }
    template<typename S>
    Array<T>& operator/=(S s){
        for (int i = 0; i < this->numElements; i++){
            this->vals[i] /= s;
        }
        return *this;
    }
};

template<typename T>
class Array<T*> {
public:
    T** vals;
    size_t maxElements;
    size_t numElements;
    Array(size_t numElements, size_t maxElements) : numElements{numElements},
            maxElements{maxElements} {
        vals = new T*[maxElements];
        this->setNull(0, maxElements);
    }
            
    Array(size_t n = 0) : numElements{n}, maxElements{n} {
        vals = new T*[maxElements];
        this->setNull(0, maxElements);
    }
    ~Array(){
        delete[] vals;
    };
    
    void setNull(size_t firstElement, size_t lastElement){
        for (int i = firstElement; i < lastElement; i++){
            this->vals[i] = NULL;
        }
    }
    void setZero(size_t firstElement, size_t lastElement){
        for (int i = firstElement; i < lastElement; i++){
//            this->vals[i] = 0;
        }
    }
    
    void changeNumElements(size_t n){
        if (n <= this->maxElements){
            this->setNull(this->numElements, n);
            this->numElements = n;
        }
        else{
            size_t newSize = (n > 2 * this->maxElements) ? n : 2 * this->maxElements;
            T** temp = new T*[newSize];
            for (int i = this->numElements; i < newSize; i++){
                temp[i] = NULL;
            }
            std::memcpy(temp, this->vals, this->numElements * sizeof(T));
            this->numElements = n;
            this->maxElements = newSize;
            delete[] this->vals;
            this->vals = temp;
        }
    }
    
    void changeMaxElements(size_t n){
        T** temp = new T*[n];
        for (int i = this->numElements; i < n; i++){
            temp[i] = NULL;
        }
        this->numElements = (this->numElements > n) ? n : this->numElements;
        this->maxElements = n;
        std::memcpy(temp, this->vals, this->numElements * sizeof(T));
        delete[] this->vals;
        this->vals = temp;
    }
    
    T*& operator[] (int i) {
        assert(i >= 0 && i < this->numElements);
        return this->vals[i];
    }
    
    T*& operator[] (int i) const{
        assert(i >= 0 && i < this->numElements);
        return this->vals[i];
    }
    
    template<typename S>
    Array<T*>& operator=(S s){
        for (int i = 0; i < this->numElements; i++){
            this->vals[i] = s;
        }
        return *this;
    }
    
    Array<T*>& operator=(const Array<T*>& other){
        if (this != &other){
//            assert(this->numElements == other.numElements);
//            delete[] this->vals;
            this->changeMaxElements(other.maxElements);
            this->changeNumElements(other.numElements);
            for (int i = 0; i < this->numElements; i++){
                this->vals[i] = other[i];
            }
        }
        return *this;
    }
    
    Array<T*>(const Array<T*> &other){
        if (this != &other){
            this->numElements = other.numElements;
            this->maxElements = other.maxElements;
//            delete[] this->vals;
            this->vals = new T[this->maxElements];
            for (int i = 0; i < this->numElements; i++){
                this->vals[i] = other[i];
            }
        }
    }
};

template<typename T>
Array<T> operator+(const Array<T>& a, const Array<T>& b){
    assert(a.numElements == b.numElements);
    Array<T> temp(a.numElements, a.maxElements);
    for (int i = 0; i < a.numElements; i++){
        temp[i] = a[i] + b[i];
    }
    return temp;
}
template<typename T>
Array<T> operator-(const Array<T>& a, const Array<T>& b){
    assert(a.numElements == b.numElements);
    Array<T> temp(a.numElements, a.maxElements);
    for (int i = 0; i < a.numElements; i++){
        temp[i] = a[i] - b[i];
    }
    return temp;
}
template<typename T>
Array<T> operator*(const Array<T>& a, const Array<T>& b){
    assert(a.numElements == b.numElements);
    Array<T> temp(a.numElements, a.maxElements);
    for (int i = 0; i < a.numElements; i++){
        temp[i] = a[i] * b[i];
    }
    return temp;
}
template<typename T>
Array<T> operator/(const Array<T>& a, const Array<T>& b){
    assert(a.numElements == b.numElements);
    Array<T> temp(a.numElements, a.maxElements);
    for (int i = 0; i < a.numElements; i++){
        temp[i] = a[i] / b[i];
    }
    return temp;
}

template<typename T>
Array<T> operator+(const Array<T>& a, T s){
    Array<T> temp(a.numElements, a.maxElements);
    for (int i = 0; i < a.numElements; i++){
        temp[i] = a[i] + s;
    }
    return temp;
}
template<typename T>
Array<T> operator+(T s, const Array<T>& a){
    Array<T> temp(a.numElements, a.maxElements);
    for (int i = 0; i < a.numElements; i++){
        temp[i] = a[i] + s;
    }
    return temp;
}
template<typename T>
Array<T> operator-(const Array<T>& a, T s){
    Array<T> temp(a.numElements, a.maxElements);
    for (int i = 0; i < a.numElements; i++){
        temp[i] = a[i] - s;
    }
    return temp;
}
template<typename T>
Array<T> operator-(T s, const Array<T>& a){
    Array<T> temp(a.numElements, a.maxElements);
    for (int i = 0; i < a.numElements; i++){
        temp[i] = s - a[i];
    }
    return temp;
}
template<typename T>
Array<T> operator*(const Array<T>& a, T s){
    Array<T> temp(a.numElements, a.maxElements);
    for (int i = 0; i < a.numElements; i++){
        temp[i] = a[i] * s;
    }
    return temp;
}
template<typename T>
Array<T> operator*(T s, const Array<T>& a){
    Array<T> temp(a.numElements, a.maxElements);
    for (int i = 0; i < a.numElements; i++){
        temp[i] = s * a[i];
    }
    return temp;
}
template<typename T>
Array<T> operator/(const Array<T>& a, T s){
    Array<T> temp(a.numElements, a.maxElements);
    for (int i = 0; i < a.numElements; i++){
        temp[i] = a[i] / s;
    }
    return temp;
}
template<typename T>
Array<T> operator/(T s, const Array<T>& a){
    Array<T> temp(a.numElements, a.maxElements);
    for (int i = 0; i < a.numElements; i++){
        temp[i] = s / a[i];
    }
    return temp;
}

template<typename T, typename S>
Array<T> operator+(const Array<T>& a, S s){
    Array<T> temp(a.numElements, a.maxElements);
    for (int i = 0; i < a.numElements; i++){
        temp[i] = a[i] + s;
    }
    return temp;
}
template<typename T, typename S>
Array<T> operator+(S s, const Array<T>& a){
    Array<T> temp(a.numElements, a.maxElements);
    for (int i = 0; i < a.numElements; i++){
        temp[i] = a[i] + s;
    }
    return temp;
}
template<typename T, typename S>
Array<T> operator-(const Array<T>& a, S s){
    Array<T> temp(a.numElements, a.maxElements);
    for (int i = 0; i < a.numElements; i++){
        temp[i] = a[i] - s;
    }
    return temp;
}
template<typename T, typename S>
Array<T> operator-(S s, const Array<T>& a){
    Array<T> temp(a.numElements, a.maxElements);
    for (int i = 0; i < a.numElements; i++){
        temp[i] = s - a[i];
    }
    return temp;
}
template<typename T, typename S>
Array<T> operator*(const Array<T>& a, S s){
    Array<T> temp(a.numElements, a.maxElements);
    for (int i = 0; i < a.numElements; i++){
        temp[i] = a[i] * s;
    }
    return temp;
}
template<typename T, typename S>
Array<T> operator*(S s, const Array<T>& a){
    Array<T> temp(a.numElements, a.maxElements);
    for (int i = 0; i < a.numElements; i++){
        temp[i] = s * a[i];
    }
    return temp;
}
template<typename T, typename S>
Array<T> operator/(const Array<T>& a, S s){
    Array<T> temp(a.numElements, a.maxElements);
    for (int i = 0; i < a.numElements; i++){
        temp[i] = a[i] / s;
    }
    return temp;
}
template<typename T, typename S>
Array<T> operator/(S s, const Array<T>& a){
    Array<T> temp(a.numElements, a.maxElements);
    for (int i = 0; i < a.numElements; i++){
        temp[i] = s / a[i];
    }
    return temp;
}

template<typename T>
Array<T> averageArray(const Array<T>& a, const Array<T>& b){
    assert(a.numElements == b.numElements);
    Array<T> temp(a.numElements, a.maxElements);
    for (int i = 0; i < a.numElements; i++){
        temp[i] = (a[i] + b[i]) / 2.0;
    }
    return temp;
}

float LSFf_temp(float x, float alpha);

float LSFf_inv_temp(float x, float alpha);

Array<float> LSF_Array(const Array<float>& a, float alpha);
Array<float> LSF_Inv_Array(const Array<float>& a, float alpha);

template<typename T>
T maxArray(const Array<T>& a){
    T temp = a[0];
    for (int i = 0; i < a.numElements; i++){
        if (a[i] > temp){
            temp = a[i];
        }
    }
    return temp;
}
template<typename T>
T minArray(const Array<T>& a){
    T temp = a[0];
    for (int i = 0; i < a.numElements; i++){
        if (a[i] < temp){
            temp = a[i];
        }
    }
    return temp;
}
template<typename T>
unsigned int minBitArray(const Array<T>& a){
    unsigned int tempMin, temp;
    std::memcpy(&tempMin, &a[0], sizeof(T));
    for (int i = 0; i < a.numElements; i++){
        std::memcpy(&temp, &a[i], sizeof(T));
        if (temp < tempMin){
            tempMin = temp;
        }
    }
    return tempMin;
}
template<typename T>
unsigned int maxBitArray(const Array<T>& a){
    unsigned int tempMax, temp;
    std::memcpy(&tempMax, &a[0], sizeof(T));
    for (int i = 0; i < a.numElements; i++){
        std::memcpy(&temp, &a[i], sizeof(T));
        if (temp > tempMax){
            tempMax = temp;
        }
    }
    return tempMax;
}
template<typename T>
T absMinArray(const Array<T>& a){
    T temp = fabs(a[0]);
    for (int i = 0; i < a.numElements; i++){
        if (fabs(a[i]) < temp){
            temp = fabs(a[i]);
        }
    }
    return temp;
}
template<typename T>
T nonzeroAbsMinArray(const Array<T>& a){
    T temp = fabs(a[0]);
    for (int i = 0; i < a.numElements; i++){
        if (fabs(a[i]) < temp && fabs(a[i]) != 0){
            temp = fabs(a[i]);
        }
    }
    return temp;
}



class Snapshot {
public:
    int posDims = 3;
    int velDims = 3;
    int numContinuous;
    int numDiscrete;
    Array<float>* position;
    Array<float>* velocity;
    Array<float> indices;
    Array<float>* continuousParams;
    Array<float>* discreteParams;
    Snapshot(){
        position = NULL;
        velocity = NULL;
        continuousParams = NULL;
        discreteParams = NULL;
        numDiscrete = 0;
        numContinuous = 0;
    }
    Snapshot(size_t nBodies, int numContinuous, int numDiscrete, int posDims = 3, int velDims = 3) : 
        numContinuous{numContinuous}, numDiscrete{numDiscrete}, posDims{posDims},
        velDims{velDims}, indices{nBodies} {
        int i;
        position = new Array<float>[posDims];
        for (i = 0; i < posDims; i++){
            position[i] = Array<float>(nBodies);
        }
        velocity = new Array<float>[velDims];
        for (i = 0; i < velDims; i++){
            velocity[i] = Array<float>(nBodies);
        }
        continuousParams = new Array<float>[numContinuous];
        for (i = 0; i < numContinuous; i++){
            continuousParams[i] = Array<float>(nBodies);
        }
        discreteParams = new Array<float>[numDiscrete];
        for (i = 0; i < numDiscrete; i++){
            discreteParams[i] = Array<float>(nBodies);
        }
    }
    ~Snapshot(){
        delete[] position;
        delete[] velocity;
        delete[] continuousParams;
        delete[] discreteParams;
    }
    
    void setNumParticles(size_t n){
        int i;
        indices.changeNumElements(n);
        for (i = 0; i < posDims; i++){
            position[i].changeNumElements(n);
        }
        for (i = 0; i < velDims; i++){
            velocity[i].changeNumElements(n);
        }
        for (i = 0; i < numContinuous; i++){
            continuousParams[i].changeNumElements(n);
        }
        for (i = 0; i < numDiscrete; i++){
            discreteParams[i].changeNumElements(n);
        }
    }
    
    Snapshot& operator=(const Snapshot& other){
        int i;
        if (this != &other){
            this->posDims = other.posDims;
            delete[] position;
            position = new Array<float>[this->posDims];
            this->velDims = other.velDims;
            delete[] velocity;
            velocity = new Array<float>[this->velDims];
            this->numContinuous = other.numContinuous;
            delete[] continuousParams;
            continuousParams = new Array<float>[this->numContinuous];
            this->numDiscrete = other.numDiscrete;
            delete[] discreteParams;
            discreteParams = new Array<float>[this->numDiscrete];
            
            this->indices = other.indices;
            if (other.position != NULL){
                for (i = 0; i < posDims; i++){
                    position[i] = other.position[i];
                }
            }
            if (other.velocity != NULL){
                for (i = 0; i < velDims; i++){
                    velocity[i] = other.velocity[i];
                }
            }
            if (other.continuousParams != NULL){
                for (i = 0; i < numContinuous; i++){
                    continuousParams[i] = other.continuousParams[i];
                }
            }
            if (other.discreteParams != NULL){
                for (i = 0; i < numDiscrete; i++){
                    discreteParams[i] = other.discreteParams[i];
                }
            }
        }
        return *this;
    }
    
    Snapshot(const Snapshot &other){
        if (this != &other){
            this->posDims = other.posDims;
            this->velDims = other.velDims;
            this->numContinuous = other.numContinuous;
            this->numDiscrete = other.numDiscrete;
            
            int i;
            this->indices = other.indices;
            position = new Array<float>[posDims];
            for (i = 0; i < posDims; i++){
                position[i] = other.position[i];
            }
            velocity = new Array<float>[velDims];
            for (i = 0; i < velDims; i++){
                velocity[i] = other.velocity[i];
            }
            continuousParams = new Array<float>[numContinuous];
            for (i = 0; i < numContinuous; i++){
                continuousParams[i] = other.continuousParams[i];
            }
            discreteParams = new Array<float>[numDiscrete];
            for (i = 0; i < numDiscrete; i++){
                discreteParams[i] = other.discreteParams[i];
            }
        }
    }
};


#endif /* SNAPSHOT_H */

