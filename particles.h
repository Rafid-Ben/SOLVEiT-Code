/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   particles.h
 * Author: Marshall
 *
 * Created on June 7, 2022, 3:55 PM
 */

#ifndef PARTICLES_H
#define PARTICLES_H

#include <cstddef>
#include "interpolation_scheme.h"
#include "simulation_parameters.h"
#include "randomizer.h"

namespace ParticleSpeciesFunction {
    inline int getCharge(ParticleSpecies particleSpecies){
        switch (particleSpecies){
            case (ParticleSpecies::monomerPositive):
                return PhysicalConstant::CHARGE;
            case (ParticleSpecies::dimerPositive):
                return PhysicalConstant::CHARGE;
            case (ParticleSpecies::trimerPositive):
                return PhysicalConstant::CHARGE;
            case (ParticleSpecies::droplet):
                // No implementation for now
                return 0;
            case (ParticleSpecies::neutral):
                return 0;
            case (ParticleSpecies::monomerNegative):
                return -PhysicalConstant::CHARGE;
            case (ParticleSpecies::dimerNegative):
                return -PhysicalConstant::CHARGE;
            case (ParticleSpecies::trimerNegative):
                return -PhysicalConstant::CHARGE;
        }
    }

    inline int getTemperature(ParticleSpecies particleSpecies){
        switch (particleSpecies){
            case (ParticleSpecies::monomerPositive):
                return PhysicalConstant::TEMPERATURE;
            case (ParticleSpecies::dimerPositive):
                return PhysicalConstant::TEMPERATURE;
            case (ParticleSpecies::trimerPositive):
                return PhysicalConstant::TEMPERATURE;
            case (ParticleSpecies::droplet):
                // No implementation for now
                return -1;
            case (ParticleSpecies::neutral):
                return PhysicalConstant::TEMPERATURE;
            case (ParticleSpecies::monomerNegative):
                return PhysicalConstant::TEMPERATURE;
            case (ParticleSpecies::dimerNegative):
                return PhysicalConstant::TEMPERATURE;
            case (ParticleSpecies::trimerNegative):
                return PhysicalConstant::TEMPERATURE;
        }
    }
    
    inline long double getMass(ParticleSpecies particleSpecies) {
        switch (particleSpecies) {
            case (ParticleSpecies::monomerPositive):
                return PhysicalConstant::EMI_MASS;
            case (ParticleSpecies::dimerPositive):
                return (2 * PhysicalConstant::EMI_MASS + PhysicalConstant::BF4_MASS);
            case (ParticleSpecies::trimerPositive):
                return (3 * PhysicalConstant::EMI_MASS + 2 * PhysicalConstant::BF4_MASS);
            case (ParticleSpecies::droplet):
                // No implementation for now
                return -1;
            case (ParticleSpecies::neutral):
                return (PhysicalConstant::EMI_MASS + PhysicalConstant::BF4_MASS);
            case (ParticleSpecies::monomerNegative):
                return (PhysicalConstant::BF4_MASS);
            case (ParticleSpecies::dimerNegative):
                return (PhysicalConstant::EMI_MASS + 2 * PhysicalConstant::BF4_MASS);
            case (ParticleSpecies::trimerNegative):
                return (2 * PhysicalConstant::EMI_MASS + 3 * PhysicalConstant::BF4_MASS);
        }
    }
}

class ChargedParticles {
public:
//    ChargedParticles();
    ~ChargedParticles();
    int firstParticle = 0;
    int lastParticle = 0;
    int maxNeutrals = 0;
    int maxNumParticles = 0;
    
    vec4<float>* position = nullptr;
    vec4<float>* velocity = nullptr;
    vec4<float>* acceleration = nullptr;
    float* temperature = nullptr;
    float* timeFlag = nullptr;
    float* timeInjected = nullptr;
    int* timesFragmented = nullptr;
    int* crash = nullptr;
    
    vec3<float>* Efield1 = nullptr;
    vec3<float>* Efield2 = nullptr;
    
    float* injectionKE = nullptr;
    float* injectionPE = nullptr;
    float* fieldEnergy = nullptr;
    
    void initialize(Initialization initialization, ParticleFractions particleFractions,
            const char* folderName = "");
    
    void propagatePosition(float dt);
    void propagateVelocity(float dt);
    
    void propagatePositionNew(float tstart, float dt, 
            CrashParameters crashParameters, FILE* eventsFile);
    void propagateVelocityNew(float tstart, float dt);
    
    void setLastParticle(int n);
    void setFirstParticle(int n);
    void incrementFirstParticle(int dN);
    void incrementLastParticle(int dN);
    
    void zeroAcceleration();
    
    void allocateMemory(int numParticles);
    void allocateMemory_Energy();
    void freeMemory();
    
    float getCharge(int index);
    float getTemperature(int index);
    long double getMass(int index);
    
    float getSpeciesCharge(int index);
    float getSpeciesTemperature(int index);
    long double getSpeciesMass(int index);
    
    void printNumParticles();
    void printPosInfo(int n);
    void printVelInfo(int n);
    void printAccelInfo(int n);
    void printTimeInjected(int n);
    
    int findMostRecentInjected(float time);
    
    void writeToEventsFile(int index, float time, FILE* eventsFile);
        
private:
//    float* Efield(int i);
//    float* Efield(float x, float y, float z);
    
    void initSingleEmitter(ParticleFractions particleFractions, const char* folderName);
    void initParaxialRay(ParticleFractions particleFractions);
    
    void bodyErrorIndexOutOfBounds(int lineNumber, const char* fileName);
};

class NeutralParticles {
public:
//    NeutralParticles();
    ~NeutralParticles();
    int firstParticle = 0;
    int lastParticle = 0;
    int maxNumParticles = 0;
    
    void initialize(Initialization initialization, int numParticles);
    
    vec3<float>* position = nullptr;
    vec3<float>* velocity = nullptr;
    float* timeCreated = nullptr;
    int* particleIndex = nullptr;
    
    void propagatePosition(float dt);
    void propagateVelocity(float dt);
    
    void propagatePositionNew(float tstart, float dt,
            CrashParameters crashParameters, FILE* eventsFile);
    void propagateVelocityNew(float tstart, float dt);
    
    void freeMemory();
    void allocateMemory(int numParticles);
    
//    int firstParticle = 0;
//    int lastParticle = 0;
    
};

class ParticleSystem {
public:
    ChargedParticles chargedParticles;
    NeutralParticles neutralParticles;
    
    InterpolationScheme interpScheme;
    
    SimulationParameters params;
    
//    void propagatePosition(float dt);
//    void propagateVelocity(float dt);
    
    void propagatePositionNew(float tstart, float dt);
    void propagateVelocityNew(float tstart, float dt);
    
    void finishIteration();
    
    void enforceSymmetryConditions();
    void enforceBoundaryConditions(float currentTime);
    void resetFirstChargedParticle();
    
    void setInjectionTiming(InjectionTimingOptions injectOptions,
            ParticleFractions particleFractions);
    void initialize(Initialization initialization, InjectionTimingOptions injectOptions,
            const char* folderName = "");
    
    void initializeBackgroundEField(BackgroundEfield bgField, const char* folderName,
            const char* fieldsName, const char* trianglesName);
    
    void fragment();
    void computeBackgroundAcceleration();
//    SymmetryType symType;
//    vec3<int> symmetryAxes;
//    vec3<float> domainSize;
    void initializeFragmentationModel(FragmentationModel fragmentationModel);
    
    void writeInitialConditionsFile();
    
    void computeChargedInjectionKineticEnergy(int n);
    float computeChargedKineticEnergy(int n);
    void computeChargedInjectionElectricPotentialEnergy(int n);
    float computeChargedElectricPotentialEnergy(int n);
    void computeChargedInjectionFieldEnergy(int index, int numSteps);
    void updateChargedFieldEnergy(int index, int numSteps, vec3<float> initPos);
    
private:
    void enforceRotationalSymmetry();
    void rotateToFirstQuadrant(float& posX, float& posY, float& velX, float& velY);
    void enforceTranslationalSymmetry();
    void translateToSymmetryDomain(float& posX, float& posY, float& posZ);
    
    void enforceSpatialBoundaryConditions(float currentTime);
    int checkOutOfSpatialBounds(float& x, float& y, float& z);
    void enforceTemporalBoundaryConditions(float currentTime);
    int checkOutOfTemporalBounds(float& timeInjected, float currentTime);
    
    void writeToEventsFile(int index, float time, ParticleType particleType);
        
    
    FragmentationModel fragmentationModel = FragmentationModel::none;
};

#endif /* PARTICLES_H */

