#include "snapshot.h"

float LSFf_temp(float x, float alpha){
    float result;
    if (x > 0){
        result = logf((x+alpha)/alpha);
    }
    else if (x < 0){
        result = logf(alpha/(alpha-x));
    }
    else{
        result = 0;
    }
    return result;
}

float LSFf_inv_temp(float x, float alpha){
    float result;
    if (x > 0){
        result = alpha*(expf(x)-1);
    }
    else if (x < 0){
        result = alpha*(1-expf(-x));
    }
    else{
        result = 0;
    }
    return result;
}

Array<float> LSF_Array(const Array<float>& a, float alpha){
    Array<float> temp(a.numElements, a.maxElements);
    for (int i = 0; i < a.numElements; i++){
        temp[i] = LSFf_temp(a[i], alpha);
    }
    return temp;
}
Array<float> LSF_Inv_Array(const Array<float>& a, float alpha){
    Array<float> temp(a.numElements, a.maxElements);
    for (int i = 0; i < a.numElements; i++){
        temp[i] = LSFf_inv_temp(a[i], alpha);
    }
    return temp;
}

