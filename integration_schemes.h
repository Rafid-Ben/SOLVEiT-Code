/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   integration_schemes.h
 * Author: Marshall
 *
 * Created on June 13, 2022, 3:03 PM
 */

#ifndef INTEGRATION_SCHEMES_H
#define INTEGRATION_SCHEMES_H

#include "particles.h"
#include "force_calculation_methods.h"
#include "snapshot.h"

class Integrator {
public:
    virtual void runSimulation() = 0;
    Integrator(ParticleSystem& particleSystem, ForceCalculationMethod& forceMethod);
    ~Integrator();
    void setComparisons(Array<ForceCalculationMethod*>& forceMethodComparisons);
    void toggleEnergy(int n, int initSteps, int normalSteps, const char* folderName);
    void initializeEnergyFile(const char* folderName);
    void setSaveDataRate(SaveDataRateOptions saveDataRateOptions){
        saveDataRate = saveDataRateOptions;
    }
    void saveFullCurrentConditions(FILE* filename);
protected:
    void incrementFirstChargedParticle(int dN);
    void incrementLastChargedParticle(int dN);
    void setFirstChargedParticle(int dN);
    void setLastChargedParticle(int dN);
    void resetFirstChargedParticle();
    void saveData();
    void updateAcceleration();
    void compareAcceleration(int iter);
    void finishIter(int iter);
    void finishInterval();
    void finishFullSim();
    void initializeSimulation();
    void initializeEnergyData();
    void updateEnergyInjectionData(int templastParticlePrev, int particleDelta);
    void updateEnergyIntervalData();
    
    ParticleSystem& particleSystem;
    ForceCalculationMethod& forceMethod;
    
    Timer fullSimTimer;
    Timer intervalTimer;
    Timer iterTimer;
    
    Array<ForceCalculationMethod*> forceMethodComparisons;
    Array<Timer> fullSimTimerComparisons;
    Array<Timer> intervalTimerComparisons;
    Array<Timer> iterTimerComparisons;
    
    float prevTime = 0.0f;
    float currentTime = 0.0f;
    int firstParticle = 0;
    int lastParticle = 0;
    vec4<float>* comparisonAcceleration = nullptr;
    
    float tempPosTime;
    float tempVelTime;
    int energyToggle = 0;
    int intervalCount = 0;
    SaveDataRateOptions saveDataRate = SaveDataRateOptions::constantRate;
    
    vec3<float>* initialPositions = nullptr;
    
    int numEnergySteps_initialization;
    int numEnergySteps;
    const char* energySaveName;
    FILE* energyFile = nullptr;
    
    float KE_inj = 0.0f;
    float PE_inj = 0.0f;
    float FE_inj = 0.0f;
    float KE = 0.0f;
    float PE = 0.0f;
    float FE = 0.0f;
};

class LeapfrogDKD : public Integrator {
public:
    LeapfrogDKD(ParticleSystem& particleSystem, ForceCalculationMethod& forceMethod) : 
    Integrator{particleSystem, forceMethod}
    {
        energySaveName = "_LeapfrogDKD.bin";
    };

    void runSimulation();
private:
    void drift(int& tempLastParticlePrev, int& tempLastParticle, float dt);
    void driftNew(int& tempLastParticlePrev, int& tempLastParticle, float dt);
    void kick(int iter, float dt);
};

class LeapfrogDKDSubsteps : public Integrator {
public:
    LeapfrogDKDSubsteps(ParticleSystem& particleSystem, ForceCalculationMethod& forceMethod,
            ForceCalculationMethod& substepMethod) : 
    Integrator{particleSystem, forceMethod}, substepMethod{substepMethod}
    {
        energySaveName = "_LeapfrogDKDSubsteps.bin";
    };
    
    LeapfrogDKDSubsteps(ParticleSystem& particleSystem, ForceCalculationMethod& forceMethod) :
    Integrator{particleSystem, forceMethod}, substepMethod{forceMethod}
    {
        energySaveName = "_LeapfrogDKDSubsteps.bin";
    };
    
    ~LeapfrogDKDSubsteps(){
        freeBoundaryConditions();
    }
    
    void runSimulation();
    
//    float substep_dt;
//    float substep_time;
    
//    float* substep_dt_array = nullptr;
//    BoundaryConditionParameters* boundaryConditions = nullptr;
    SubstepBoundaryConditions* substepBoundaryConditions = nullptr;
    
    void setSubstepTimes(float* substepTimes, BoundaryConditionParameters* bc, int size);
    void setSubstepTimes(SubstepBoundaryConditions* substepBC, int size);
    void setDefaultSubstepTimes();
    
private:
    ForceCalculationMethod& substepMethod;
    
    int maxLevel = 0;
    
    void freeBoundaryConditions();
    
    void drift(int iter, int& tempLastParticlePrev, int& tempLastParticle, 
            int level, int toggle);
    void driftNew(int iter, int& tempLastParticlePrev, int& tempLastParticle, 
            int level, int toggle);
    void kick(int iter, int& tempLastParticlePrev, int& tempLastParticle, 
            int level, int topLevel);
    
    void variableTimestepKick();
    
    void initializeParticleLevels();
    void checkConditions(int index, int level);
    void checkSpatialConditions(int index, int level);
    void checkTemporalConditions(int index, int level);
    void crossoverParticles(int level);
    void setParticleStatus(int level, int topLevel);
    void updateAccelerationSubsteps();
    
    int* particleLevels = nullptr;
};

#endif /* INTEGRATION_SCHEMES_H */

