/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   simulation_parameters.h
 * Author: Marshall
 *
 * Created on June 7, 2022, 8:30 AM
 */

/**
 * Initializes and stores various simulation parameters
 * 
 * 
 */

#ifndef SIMULATION_PARAMETERS_H
#define SIMULATION_PARAMETERS_H

#include <cstddef>
#include "constants.h"
#include "structures.h"
#include "auxiliary_functions.h"

//class SimulationParameters {
//public:
//    SimulationParameters();
//    ~SimulationParameters();
//    
//    void initializeInterpolationScheme(const char* fieldsName, const char* trianglesName);
//    void readInputInfoFile(const char* infoFileName);
//    void setTimeStepInfo();
//    void createOutputInfoFile();
//     
//    int m_numParticles = -1;
//    int m_numItersStored = 100;
//    int m_numThreads = 1;
//    int m_numIters;
//    int m_interval;
//    int m_numSymmetryLevels = 0;
//    
//    float m_timeStep;
//    float m_simTime;
//    float m_I0 = 0.0;
//    float m_v0;
//    float m_E0;
//    float m_r0;
//    float m_worldDimR;
//    float m_worldDimZ;
//    float m_geom[6];
//    float m_dN;
//    
//    const char* m_folder;
//    const char* m_timerName;
//    
//    vec3<int> symmetryAxes;
//    vec3<float> domainSize;
//    
//    SaveDataOptions m_saveDataToggle;
//    InitializationOptions m_initInfo;
//    BackgroundEField m_bgField;
//    InterpolationScheme m_interpScheme;
//    FILE** m_saveDataFiles = nullptr;
//    TimerOptions m_timerInfo;
//    SymmetryType m_symType;
//    
//    Compression<float> compression;
//    
//    
//};

class SimulationParameters {
public:
    SimulationParameters(){
        timeStep = 0.0;
        simTime = 0.0;
//        printf("\nWarning: Simulation Parameters class uninitialized.\n");
//        printf("Please manually initialize member variables.\n");
    };
    SimulationParameters(float timeStep, float simTime, int numItersStored, 
            int numThreads, const char* folder) : timeStep{timeStep}, 
            simTime{simTime}, numItersStored{numItersStored}, numThreads{numThreads},
            folder{folder} {};
    
    ~SimulationParameters();
    
    void setDefaultInjectionTimes();
    void setInjectionTimes(float* startTimes, float* stopTimes, int size);
    void initialize(float timeStep, float simTime, int numItersStored,
            int numThreads, const char* folder){
        this->timeStep = timeStep;
        this->simTime = simTime;
        this->numItersStored = numItersStored;
        this->numThreads = numThreads;
        this->folder = folder;
    };
    void readInputInfoFile(const char* infoFileName);
    void setTimeStepInfo();
    void createOutputInfoFile();
    void initializeDataStorageFiles(const char* additionalName = "");

    int numStartingParticles = 0;
    int numParticles = -1;
    int numItersStored = 100;
    int numThreads = 1;
    int numIters;
    int interval;
    int numSymmetryLevels = 0;
    int sizeInjectionTimes = 0;

    float timeStep;
    float simTime;
    float I0 = 0.0;
    float v0;
    float E0;
    float r0;
    float worldDimR;
    float worldDimZ;
    float geom[6];
    float dN;
    float a;
    float nu0;
    float V;
    
    float totalInjectingTime;
    
    float* startInjectionTimes = nullptr;
    float* stopInjectionTimes = nullptr;

    const char* folder;
    const char* timerName = nullptr;

    vec3<int> symmetryAxes;
    vec3<float> domainSize;

    SymmetryType symType = SymmetryType::none;
    
//    BoundaryConditionOptions boundaryConditionType = BoundaryConditionOptions::none;
//    vec4<float> spatialBoundaryConditions = 0.0f;
//    float temporalBoundaryCondition = 0.0f;
    BoundaryConditionParameters boundaryConditions;
//    CrashConditionOptions crashConditions = CrashConditionOptions::none;
    CrashParameters crashParameters;
    
    //    Compression<float> compression;
    SaveDataOptions saveDataToggle = SaveDataOptions::newOriginal;
    TimerOptions timerInfo = TimerOptions::off;
    FILE** saveDataFiles = nullptr;
    int numSaveDataFiles;
    
    FILE* eventsFile = nullptr;
    FILE* neutralsFile = nullptr;

private:
    void freeInjectionTimes(){
        if (startInjectionTimes != nullptr) delete[] startInjectionTimes;
        if (stopInjectionTimes != nullptr) delete[] stopInjectionTimes;
        startInjectionTimes = nullptr;
        stopInjectionTimes = nullptr;
    }
    FILE* initializeWritingStateFiles(const char* fileName);
    FILE* initializeWritingEventsFiles(const char* fileName);
    FILE* initializeWritingNeutralsFiles(const char* fileName);
};

#endif /* SIMULATION_PARAMETERS_H */

