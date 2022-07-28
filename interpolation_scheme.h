/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   interpolation_scheme.h
 * Author: Marshall
 *
 * Created on June 7, 2022, 4:54 PM
 */

#ifndef INTERPOLATION_SCHEME_H
#define INTERPOLATION_SCHEME_H

#include "auxiliary_functions.h"
extern "C" {
#include "pwl_interp_2d_scattered.h"
}

class InterpolationScheme {
public:
    InterpolationScheme(BackgroundEfield bgField = BackgroundEfield::none) : bgField{bgField}
    {
        setReferencePosition();
    };
    ~InterpolationScheme();

    float r0 = 0.0;
    float E0 = 0.0;

    double a = 0.0;
    double nu0 = 0.0;
    double V = 0.0;

    vec3<float> referencePosition = 0.0f;

    void changeBackgroundField(BackgroundEfield bgField);

    void initializeInterpolationScheme(const char* folderName, const char* fieldsName,
            const char* trianglesName);

    vec3<float> Efield(float x, float y, float z);

    BackgroundEfield bgField;
private:
    int m_nodeNum = 0;
    int m_elementNum = 0; //number of elements (triangles)

    int* m_elementNeighbor = nullptr; //connectivity matrix for the triangles of the interpolation scheme
    int* m_triangle = nullptr; //connectivity matrix for the nodes of the interpolation scheme
    double* m_nodeXY = nullptr; //interpolation scheme coordinate data
    double* m_zEr = nullptr; //interpolation scheme Er data
    double* m_zEz = nullptr; //interpolation scheme Ez data

    vec3<float> computeEfieldSingleEmitter(float x, float y, float z);
    void initializeSingleEmitterInterpolationScheme(const char* folderName,
            const char* fieldsName, const char* trianglesName);
    void fillInterpolationScheme(char* fieldsName);
    void fillTriangleScheme(char* triangleName);
    void freeInterpolationScheme();

    void setReferencePosition();
};

#endif /* INTERPOLATION_SCHEME_H */

