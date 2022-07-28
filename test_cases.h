/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   test_cases.h
 * Author: Marshall
 *
 * Created on July 6, 2022, 11:52 AM
 */

#ifndef TEST_CASES_H
#define TEST_CASES_H

#include "particles.h"
#include "force_calculation_methods.h"
#include "integration_schemes.h"

void compareForceCalculationMethodsAcceleration(float dt, float simTime,
        int numItersStored, int numThreads, const char* folderName,
        const char* infoFileName, const char* fieldsName, 
        const char* trianglesName);

void compareForceCalculationMethodsEvolution(float dt, float simTime,
        int numItersStored, int numThreads, const char* folderName,
        const char* infoFileName, const char* fieldsName, 
        const char* trianglesName);

void compareForceCalculationMethodsEndSim(float dt, float simTime,
        int numItersStored, int numThreads, const char* folderName,
        const char* infoFileName, const char* fieldsName, 
        const char* trianglesName);

#endif /* TEST_CASES_H */

