/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   main.cpp
 * Author: Marshall
 *
 * Created on June 6, 2022, 8:19 PM
 */

#include <cstdlib>

#include "particles.h"
#include "force_calculation_methods.h"
#include "integration_schemes.h"
#include "test_cases.h"

/*
 * 
 */
int main(int argc, char** argv) {
    
    // Call example test functions
    
    // Run basic simulation
    
    // Run more advanced simulation
    
    // Compare methods (L2-norm)
    
    // Compare methods (propagation)
    
    // Compare integrators
    
    const char* folderName = "EFIELD_DATABASE/1.9/"; // Folder name
    const char* infoFileName = "info_323.txt";
    const char* fieldsName = "fields_323.csv";
    const char* trianglesName = "triangles_323.csv";
//    
//    
//    
//    InterpolationScheme interpScheme(BackgroundEfield::singleEmitter);
//    interpScheme.initializeInterpolationScheme(folderName, fieldsName, trianglesName);
//    
//    vec3<float> ptemp = interpScheme.Efield(5, 10, 10);
//    
//    interpScheme.changeBackgroundField(BackgroundEfield::none);
//    
//    vec3<float> ptemp2 = interpScheme.Efield(5, 10, 10);
//    
//    interpScheme.changeBackgroundField(BackgroundEfield::singleEmitter);
//    interpScheme.initializeInterpolationScheme(folderName, fieldsName, trianglesName);
//    interpScheme.initializeInterpolationScheme(folderName, fieldsName, trianglesName);
//    
//    vec3<float> ptemp3 = interpScheme.Efield(5, 10, 10);
    
    float dt = 1e-12;
    float simTime = 1e-10;
    int numItersStored = 100;
    int numThreads = 1;
    
    compareForceCalculationMethodsAcceleration(dt, simTime, numItersStored,
            numThreads, folderName, infoFileName, fieldsName, trianglesName);
    
    
//    compareForceCalculationMethodsEndSim(dt, simTime, numItersStored, numThreads,
//            folderName, infoFileName, fieldsName, trianglesName);
    
    
    
//    ParticleSystem particleSystem;
//    particleSystem.params.initialize(dt, simTime, numItersStored, numThreads,
//            folderName);
//    particleSystem.params.readInputInfoFile(infoFileName);
//    particleSystem.params.createOutputInfoFile();
//    particleSystem.params.initializeDataStorageFiles();
//    particleSystem.params.timerInfo = TimerOptions::minimized;
//    
////    particleSystem.params.boundaryConditions.boundaryConditionType = BoundaryConditionOptions::space;
////    particleSystem.params.boundaryConditions.temporalBoundaryCondition = 5e-10;
////    particleSystem.params.boundaryConditions.spatialBoundaryConditions.x = 2e-6;
////    particleSystem.params.boundaryConditions.spatialBoundaryConditions.y = 2e-6;
////    particleSystem.params.boundaryConditions.spatialBoundaryConditions.z = 15e-6;
////    particleSystem.params.boundaryConditions.spatialBoundaryConditions.w = 2.5e-6;
//    
////    particleSystem.params.crashParameters.crashOptions = CrashConditionOptions::extractorGrid;
////    particleSystem.params.crashParameters.r_distance = 1e-6;
////    particleSystem.params.crashParameters.z_distance = 5e-6;
////    particleSystem.params.crashParameters.z_width = 10e-6;
//    
////    float startTimes[] = {0, 5e-10};
////    float stopTimes[] = {2.5e-10, 7.5e-10};
////    
////    particleSystem.params.setInjectionTimes(startTimes, stopTimes, 2);
//    
//    particleSystem.initialize(Initialization::singleEmitter, 
//            InjectionTimingOptions::constantRate, folderName);
//    
//    particleSystem.initializeBackgroundEField(BackgroundEfield::singleEmitter, 
//            folderName, fieldsName, trianglesName);
//    
//    particleSystem.initializeFragmentationModel(FragmentationModel::none);
//    
////    for (int i = 0; i < 10; i++){
////        particleSystem.chargedParticles.printPosInfo(i);
////    }
////    
////    for (int i = 0; i < 10; i++){
////        particleSystem.chargedParticles.printVelInfo(i);
////    }
////    
////    for (int i = 0; i < 10; i++){
////        particleSystem.chargedParticles.printTimeInjected(i);
////    }
////    
////    particleSystem.chargedParticles.setFirstParticle(0);
////    particleSystem.chargedParticles.setLastParticle(2);
////    particleSystem.chargedParticles.zeroAcceleration();
//    
//    DirectForceCalculation df(&particleSystem);
////    BarnesHutForceCalculation bh(particleSystem);
//
////    df.setFirstParticle(0);
////    df.setLastParticle(2);
////    
////    df.updateAcceleration();
////    
////    for (int i = 0; i < 10; i++){
////        particleSystem.chargedParticles.printAccelInfo(i);
////    }
//    
////    LeapfrogDKD leap(particleSystem, df);
////    
////    leap.setSaveDataRate(SaveDataRateOptions::constantRate_and_lastSnapshots);
////    leap.toggleEnergy(true, 100000, 10000, folderName);
////    
////    leap.runSimulation();
//    
//    int numSubsteps = 1;
//    float substepTimes[] = {1.5e-13};
//    BoundaryConditionParameters* bc = new BoundaryConditionParameters[numSubsteps];
//    for (int i = 0; i < numSubsteps; i++){
//        bc[i].boundaryConditionType = BoundaryConditionOptions::space_and_time;
//        if (i == 0){
//            bc[i].temporalBoundaryCondition = 50.0e-12;
//            bc[i].spatialBoundaryConditions.z = 2.0e-6;
//        }
//        if (i == 1){
//            bc[i].temporalBoundaryCondition = 5e-13;
//            bc[i].spatialBoundaryConditions.w = 5e-6;
//        }
//        if (i == 2){
//            bc[i].temporalBoundaryCondition = 7e-13;
//            bc[i].spatialBoundaryConditions.z = 10e-6;
//        }
//        if (i == 3){
//            bc[i].temporalBoundaryCondition = 6e-13;
//            bc[i].spatialBoundaryConditions.y = 1e-6;
//        }
//        if (i == 4){
//            bc[i].temporalBoundaryCondition = 1e-12;
//            bc[i].spatialBoundaryConditions.w = 4e-6;
//        }
//    }
////    
////    SubstepBoundaryConditions* substepBC = new SubstepBoundaryConditions[numSubsteps];
////    
////    for (int i = 0; i < numSubsteps; i++){
////        substepBC[i].dt = substepTimes[i];
////        substepBC[i].bc = bc[i];
////    }
////    
////    for (int i = 0; i < numSubsteps; i++){
////        printf("%d: %.9g\n",i,substepBC[i].dt);
////    }
////    printf("\n");
////    std::sort(substepBC, substepBC + numSubsteps);
////    for (int i = 0; i < numSubsteps; i++){
////        printf("%d: %.9g\n",i,substepBC[i].dt);
////    }
//    
////    BarnesHutForceCalculation bh(&particleSystem);
//    MultipoleMethodConstants c;
//    MultipoleMethodVariables v;
//    FmmType t = FmmType::fmm;
//    MultipoleMethodForceCalculation fmm(&particleSystem, &c, &v, t);
//    Array<ForceCalculationMethod*> fca(1);
//    fca[0] = &fmm;
//    
////    LeapfrogDKDSubsteps leap(particleSystem, df);
////    leap.setSubstepTimes(substepTimes, bc, numSubsteps);
//    
//    LeapfrogDKD leap(particleSystem, df);
//    leap.setComparisons(fca);
//    
//    for (int j = 0; j < 10; j++){
//        particleSystem.chargedParticles.printPosInfo(j);
//    }
//    printf("\n");
//    
//    leap.runSimulation();
//    
//    for (int j = 0; j < 10; j++){
//        particleSystem.chargedParticles.printPosInfo(j);
//    }
//    printf("\n");
//    
//    delete[] bc;

    return 0;
}

