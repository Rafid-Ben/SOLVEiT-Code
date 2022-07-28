/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   force_calculation_methods.h
 * Author: Marshall
 *
 * Created on June 11, 2022, 4:50 PM
 */

#ifndef FORCE_CALCULATION_METHODS_H
#define FORCE_CALCULATION_METHODS_H

#include "constants.h"
#include "structures.h"
#include "particles.h"
#include "lattice_point.h"
#include "timer.h"
#include "cell.h"
#include "multipole_method_functions.h"

#include "direct_force_method.h"
#include "domain_bounds_functions.h"

class ForceCalculationMethod {
public:
    virtual void updateAcceleration() = 0;
    void setNumThreads(int n = 1);
    
    ForceCalculationMethod(){};
    ForceCalculationMethod(ParticleSystem* particleSystem) : 
    particleSystem{particleSystem} {};
    
    ~ForceCalculationMethod();
    void printNumThreads();
    void setFirstParticle(int n);
    void setLastParticle(int n);
    void incrementLastParticle(int dN);
    void incrementFirstParticle(int dN);
    void setGrid(Grid* grid);
    
    void initializeTimingFile(const char* folderName);
    void initializeL2NormFile(const char* folderName);
    void initializeFullL2NormFile(const char* folderName);
    
    void writeL2NormInfo(int numParticles, float L2norm, float maxL2norm,
            int maxL2normIndex);
    
    virtual void increaseLevel() = 0;
    
    const char* timerName = "";
    
    FILE* timeInfoFile = nullptr;
    FILE* L2NormFile = nullptr;
    FILE* fullL2NormFile = nullptr;
    
    int compareRate = 1;
    float previousTotalTimes[FMMConstant::windowSize];
    int crossoverPoints[FMMConstant::crossoverSize];
    ForceCalculationMethod* olderSibling = nullptr;
    
    virtual void setSaveName() = 0;
    virtual void setTimerName() = 0;
protected:
    void allocateThreadMemory_Base();
    void deallocateThreadMemory_Base();
    virtual void allocateThreadMemory() = 0;
    virtual void deallocateThreadMemory() = 0;
    
    ParticleSystem* particleSystem = nullptr;
    
    int numThreads = 1;
    int firstParticle = 0;
    int lastParticle = 0;
    
    const char* saveName = "";
    
    Grid* grid = nullptr;
    
    vec4<int>* rangevalues = nullptr;
#ifdef _PTHREAD_H
    pthread_t* threads = nullptr;
#endif
};


class NoParticleParticleForceCalculation : public ForceCalculationMethod{
public:
    void updateAcceleration() { particleSystem->chargedParticles.zeroAcceleration(); };
    void increaseLevel() {};
    NoParticleParticleForceCalculation(ParticleSystem* particleSystem) : 
    ForceCalculationMethod{particleSystem}
    {
        setTimerName();
        setSaveName();
    };
    void setSaveName(){
        if (saveName[0] != '\0'){
            saveName = "_nop2p.bin";
        }
    }
    void setTimerName(){
        if (timerName[0] != '\0'){
            timerName = "No Particle-Particle Force Method";
        }
    }
};


class DirectForceCalculation;

typedef struct {
    DirectForceCalculation* d;
    int threadNum;
} directForceStruct;

class DirectForceCalculation : public ForceCalculationMethod {
public:
    void updateAcceleration();
    void increaseLevel(){};
    DirectForceCalculation(){};
    DirectForceCalculation(ParticleSystem* particleSystem);
    ~DirectForceCalculation();
    
    void setVectorizedToggle(int toggle);
    
    void setSaveName(){
        if (saveName[0] == '\0'){
            if (vectorizedToggle){
                switch (particleSystem->params.symType){
                    case (SymmetryType::none):
                        saveName = "_directVectorized.bin";
                        break;
                    case (SymmetryType::rotational):
                        saveName = "_directRotationalVectorized.bin";
                        break;
                    case (SymmetryType::translational):
                        saveName = "_directTranslationalVectorized.bin";
                        break;
                }
            }
            else {
                switch (particleSystem->params.symType){
                    case (SymmetryType::none):
                        saveName = "_direct.bin";
                        break;
                    case (SymmetryType::rotational):
                        saveName = "_directRotational.bin";
                        break;
                    case (SymmetryType::translational):
                        saveName = "_directTranslational.bin";
                        break;
                }
            }
        }
    }
    void setTimerName(){
        if (timerName[0] == '\0'){
            if (vectorizedToggle){
                switch (particleSystem->params.symType){
                    case (SymmetryType::none):
                        timerName = "Direct Force Method Vectorized";
                        break;
                    case (SymmetryType::rotational):
                        timerName = "Rotational Direct Force Method Vectorized";
                        break;
                    case (SymmetryType::translational):
                        timerName = "Translational Direct Force Method Vectorized";
                        break;
                }
            }
            else{
                switch (particleSystem->params.symType){
                    case (SymmetryType::none):
                        timerName = "Direct Force Method";
                        break;
                    case (SymmetryType::rotational):
                        timerName = "Rotational Direct Force Method";
                        break;
                    case (SymmetryType::translational):
                        timerName = "Translational Direct Force Method";
                        break;
                }
            }
        }
    }

protected:
    static void* directKernel(void* dF);
    void directKernelCalculation(int threadNum = 0);
    
    static void* directKernelRotation(void* dF);
    void directKernelRotationCalculation(int threadNum = 0);
    
    static void* directKernelTranslation(void* dF);
    void directKernelTranslationCalculation(int threadNum = 0);
    
    
    void allocateThreadMemory();
    void deallocateThreadMemory();
    
    int vectorizedToggle = false;
    
    Jpdata *jptc = nullptr;

#ifdef _PTHREAD_H
    directForceStruct* dF = nullptr;
#endif
};


class BarnesHutForceCalculation;

typedef struct {
    BarnesHutForceCalculation* d;
    int threadNum;
} barnesHutForceStruct;

class BarnesHutForceCalculation : public ForceCalculationMethod {
public:
//    BarnesHutForceCalculation& operator=(const BarnesHutForceCalculation& rhs){
//        printf("Assignment operator.\n");
//        return *this;
//    }
    void updateAcceleration();
    void increaseLevel(){};
    BarnesHutForceCalculation(){};
    BarnesHutForceCalculation(ParticleSystem* particleSystem);
    ~BarnesHutForceCalculation();
    
    void setBoxDimVals(vec3<int> boxDimVals);
    void setDomainToggle(DomainToggle toggle);
    void setSubcellBlockSize(int subcellBlockSize);
    void setTheta(float tempTheta);
    void setOutlierToggle(int toggle);
    
    void setSaveName(){
        if (saveName[0] == '\0'){
            char bufferTemp[512];
            snprintf(bufferTemp, sizeof(char) * 512, "_%.0ftheta", 1000*theta);
            switch (particleSystem->params.symType){
                case (SymmetryType::none):
                    saveName = concat(bufferTemp, "_barnesHut.bin");
                    break;
                case (SymmetryType::rotational):
                    saveName = concat(bufferTemp, "_barnesHutRotational.bin");
                    break;
                case (SymmetryType::translational):
                    saveName = concat(bufferTemp, "_barnesHutTranslational.bin");
                    break;
            }
        }
    };
    void setTimerName(){
        if (timerName[0] == '\0'){
            char bufferTemp[512];
            snprintf(bufferTemp, sizeof(char) * 512, " (Theta = %f)", theta);
            switch (particleSystem->params.symType){
                case (SymmetryType::none):
                    timerName = concat("Barnes Hut Method", bufferTemp);
                    break;
                case (SymmetryType::rotational):
                    timerName = concat("Rotational Barnes Hut Method", bufferTemp);
                    break;
                case (SymmetryType::translational):
                    timerName = concat("Translational Barnes Hut Method", bufferTemp);
                    break;
            }
        }
    };

protected:
    float computeEffectiveDistance(int index, Cell* cell);
    float computeEffectiveDistance(vec3<float> pos, Cell* cell);
    int locateSubcell(Cell* cell, int index);
    void addToCell(Cell* cell, int index);
    void generateTree();
    void computeCellProperties(Cell* cell);
    void computeForceFromCell(Cell* cell, int index, vec3<float>& accel, float deff);
    void computeForceFromTree(Cell* cell, int index, vec3<float>& accel);
    void computeWorldLimits();
    
    static void* barnesHutKernel(void* dF);
    void barnesHutKernelCalculation(int threadNum = 0);
    
    
    float computeEffectiveDistanceRotations(vec3<float> pos, Cell* cell);
    void computeForceFromCellRotations(Cell* cell, vec3<float> pos,
            vec3<float>& accel, float deff, int cb, int sb);
    void computeForceFromTreeRotations(Cell* cell, vec3<float> pos, int index, 
            vec3<float>& accel, int cb, int sb, int quadrant);
    static void* rotationalBarnesHutKernel(void* dF);
    void rotationalBarnesHutKernelCalculation(int threadNum = 0);
    
    
    void computeForceFromCellTranslations(Cell* cell, vec3<float> pos,
            vec3<float>& accel, float deff);
    void computeForceFromTreeTranslations(Cell* cell, vec3<float> pos, int index,
            vec3<float>& accel, int centerCheck);
    static void* translationalBarnesHutKernel(void* dF);
    void translationalBarnesHutKernelCalculation(int threadNum = 0);
    void translationalBarnesHut(LatticePoint lattice, int threadNum = 0);
    
    void allocateThreadMemory();
    void deallocateThreadMemory();
    
#ifdef _PTHREAD_H
    barnesHutForceStruct* dF = nullptr;
#endif
    
    int* outlierCheck = nullptr;
    int numOutliers;
    int outlierToggle = false;
    
    int forceBoxDims = false;
    vec3<int> forcedBoxDimVals = {1, 1, 1};
    DomainToggle domainToggle = DomainToggle::cubicCells;
    
    int subcellBlock = 2;
    
    vec3<float> dims;
    vec3<float> roots;
    vec3<int> domainDimVals;
    Cell* rootCell = nullptr;
    float theta = BHConstant::THETA;
};


class MultipoleMethodConstants {
public:
    MultipoleMethodConstants();
    ~MultipoleMethodConstants();
private:
    void precalc();

    float* factorial;
    std::complex<double>* Ynm;
    std::complex<double>*** Dnm;
    int* rotationLink;
    float* Anm;
    float* anm;

    friend class MultipoleMethodForceCalculation;
};

class MultipoleMethodVariables {
public:
    MultipoleMethodVariables(){};
    ~MultipoleMethodVariables();
private:
    void allocateMMvariables(int maxNumParticles, int numBoxIndexLeaf, int numBoxIndexTotal,
            int numBoxIndexFull, int numExpansions, int maxLevel, int numLevelSC,
            int numDomains, int numOutlierDomains);
    
    void allocateDomainVariables(int numLevelSC, int numDomains, int avgPPB);
    
    void allocateOutlierVariables(int numOutlierDomains, int numOutliers, int outlierNumBoxIndexLeaf);
    void reallocateTempParticles(int tempIndex, int additionalParticles);
    
    void allocateDomainLinkage();
    
    void deallocateAll();
    
    int maxNumBoxIndexLeaf = 0;
    int maxNumBoxIndexTotal = 0;
    int maxNumBoxIndexFull = 0;
    int maxNumExpansions = 0;
    int maxMaxLevel = 0;
    int maxNumLevelSC = 0;
    int maxNumDomains = 0;
    int maxNumOutlierDomains = 0;
    int maxMaxNumParticles = 0;
    int maxNumOutliers = 0;
    int maxOutlierNumBoxIndexLeaf = 0;
    
    // Can be allocated together - requires info from reallocate() inputs
    int* permutation = nullptr;
    int* sortValue = nullptr;
    int* sortIndex = nullptr;
    int* sortValueBuffer = nullptr;
    int* outlierCheck = nullptr;
    int* domainMortonIndex = nullptr;
    
    int** particleOffset = nullptr;
    int* levelOffset = nullptr;
    int* numInteraction = nullptr;
    int*** interactionList = nullptr;
    int* boxOffsetStart = nullptr;
    int* boxOffsetEnd = nullptr;
    int* sortIndexBuffer = nullptr;
//    std::complex<double>(*Lnm)[FMMConstant::numCoefficients]; // local expansion coefficients
//    std::complex<double>(*LnmOld)[FMMConstant::numCoefficients]; // Lnm from previous level
//    std::complex<double>(*Mnm)[FMMConstant::numCoefficients]; // multipole expansion coefficients
    std::complex<double> **Lnm = nullptr;
    std::complex<double> **LnmOld = nullptr;
    std::complex<double> **Mnm = nullptr;
    int* domainBoxIndexMask = nullptr;
    int* domainBoxIndexFull = nullptr;
    int*** domainOffset = nullptr;
    int** outlierParticleOffset = nullptr;
    int*** outlierDomainOffset = nullptr;
    int* neo = nullptr;
    
    // Can be allocated from numDomains, avgPPB, and 
    vec3<int>* dimValsSC = nullptr;
    vec3<float>* boxMinSC = nullptr;
    int** tempParticles = nullptr;
    int** tempParticlesPerDomain = nullptr;
    vec3<int>* offsetSC = nullptr;
    int* numDomainsLevel = nullptr;
    vec3<int>* boxMinIndexSC = nullptr;
    
    // Might make the most sense to allocate this separately
    domainInfo** domainLinkage = nullptr;
    
    // Can be allocated with outlier info
    outlierDomainInfo* outlierDomainLinkage = nullptr;
    outlierInfo* outlierList = nullptr;
    int* outlierBoxIndexFull = nullptr;
    
    
    
    friend class MultipoleMethodForceCalculation;
};

class MultipoleMethodForceCalculation;

typedef struct {
    MultipoleMethodForceCalculation* d;
    int threadNum;
} multipoleMethodForceStruct;

class MultipoleMethodForceCalculation : public ForceCalculationMethod {
public:
    void updateAcceleration();
    void increaseLevel(){if (forcedMaxLevel < 6) forcedMaxLevel++;};
    MultipoleMethodForceCalculation(){};
    MultipoleMethodForceCalculation(ParticleSystem* particleSystem,
            MultipoleMethodConstants* c, MultipoleMethodVariables* v, FmmType fType);
    ~MultipoleMethodForceCalculation();
    
    void setNumExpansions(int p);
    void setBoxDimVals(vec3<int> boxDimVals);
    void setVectorizedToggle(int toggle);
    void setMaxLevel(int maxLevel);
    
    void setSaveName(){
        if (saveName[0] == '\0'){
            char bufferTemp[512];
            snprintf(bufferTemp, sizeof(char) * 512, "_%dexpansions_%dm", numExpansions,
                    forcedMaxLevel);
            switch (particleSystem->params.symType){
                case (SymmetryType::none):
                    switch (fmmType){
                        case (FmmType::treecode):
                            saveName = concat(bufferTemp, "_treecodeMultipoleMethod.bin");
                            break;
                        case (FmmType::treecodeOutliers):
                            saveName = concat(bufferTemp, "_treecodeMultipoleMethodOutliers.bin");
                            break;
                        case (FmmType::fmm):
                            saveName = concat(bufferTemp, "_fmm.bin");
                            break;
                        case (FmmType::fmmOutliers):
                            saveName = concat(bufferTemp, "_fmmOutliers.bin");
                            break;
                    }
                    break;
                case (SymmetryType::rotational):
                    switch (fmmType){
                        case (FmmType::treecode):
                            saveName = concat(bufferTemp, "_treecodeMultipoleMethod.bin");
                            break;
                        case (FmmType::treecodeOutliers):
                            saveName = concat(bufferTemp, "_treecodeMultipoleMethodOutliers.bin");
                            break;
                        case (FmmType::fmm):
                            saveName = concat(bufferTemp, "_fmm.bin");
                            break;
                        case (FmmType::fmmOutliers):
                            saveName = concat(bufferTemp, "_fmmOutliers.bin");
                            break;
                    }
                    break;
                case (SymmetryType::translational):
                    switch (fmmType){
                        case (FmmType::treecode):
                            saveName = concat(bufferTemp, "_treecodeMultipoleMethod.bin");
                            break;
                        case (FmmType::treecodeOutliers):
                            saveName = concat(bufferTemp, "_treecodeMultipoleMethodOutliers.bin");
                            break;
                        case (FmmType::fmm):
                            saveName = concat(bufferTemp, "_fmm.bin");
                            break;
                        case (FmmType::fmmOutliers):
                            saveName = concat(bufferTemp, "_fmmOutliers.bin");
                            break;
                    }
                    break;
            }
        }
    };
    void setTimerName(){
        if (timerName[0] == '\0'){
            char bufferTemp[512];
            snprintf(bufferTemp, sizeof(char) * 512, " (p = %d, forcedMaxLevel = %d)", 
                    numExpansions, forcedMaxLevel);
            switch (particleSystem->params.symType){
                case (SymmetryType::none):
                    switch (fmmType){
                        case (FmmType::treecode):
                            timerName = concat("Treecode Multipole Method", bufferTemp);
                            break;
                        case (FmmType::treecodeOutliers):
                            timerName = concat("Treecode Multipole Method (Outliers)", bufferTemp);
                            break;
                        case (FmmType::fmm):
                            timerName = concat("FMM", bufferTemp);
                            break;
                        case (FmmType::fmmOutliers):
                            timerName = concat("FMM (Outliers)", bufferTemp);
                            break;
                    }
                    break;
                case (SymmetryType::rotational):
                    switch (fmmType){
                        case (FmmType::treecode):
                            timerName = concat("Rotational Treecode Multipole Method", bufferTemp);
                            break;
                        case (FmmType::treecodeOutliers):
                            timerName = concat("Rotational Treecode Multipole Method (Outliers)", bufferTemp);
                            break;
                        case (FmmType::fmm):
                            timerName = concat("Rotational FMM", bufferTemp);
                            break;
                        case (FmmType::fmmOutliers):
                            timerName = concat("Rotational FMM (Outliers)", bufferTemp);
                            break;
                    }
                    break;
                case (SymmetryType::translational):
                    switch (fmmType){
                        case (FmmType::treecode):
                            timerName = concat("Translational Treecode Multipole Method", bufferTemp);
                            break;
                        case (FmmType::treecodeOutliers):
                            timerName = concat("Translational Treecode Multipole Method (Outliers)", bufferTemp);
                            break;
                        case (FmmType::fmm):
                            timerName = concat("Translational FMM", bufferTemp);
                            break;
                        case (FmmType::fmmOutliers):
                            timerName = concat("Translational FMM (Outliers)", bufferTemp);
                            break;
                    }
                    break;
            }
        }
    };

protected:
    int numExpansions = FMMConstant::maxNumExpansions;
    int numCoefficients;
    int forcedMaxLevel = -2;
    int forceBoxDims = -1;
    vec3<int> forcedBoxDimVals = {1, 1, 1};
    
    void updateAccelerationNoSym();
    void updateAccelerationRotational();
    void updateAccelerationTranslational();
    
    void setDomainSizeRotations();
    void setOptimumLevelRotations();
    void getInteractionListRotations(int numBoxIndex, int numLevel,
            FmmInteractionType interactionType);
    
    void setDomainSizeTranslations();
    void setOptimumLevelTranslations();
    void presortMortonTranslations();
    void mortonTranslations();
    void getBoxDataTranslations(int& numBoxIndex, int toggle);
    void getBoxDataOfParentTranslations(int& numBoxIndex, int numLevel, int toggle);
    void getBoxIndexMaskTranslations(int numBoxIndex, int numLevel, int toggle);
    void getInteractionListTranslations(int numBoxIndex, int numLevel, FmmInteractionType interactionType);
    void sortParticlesTranslations();
    void countNonEmptyBoxesTranslations();
    
    void createDomainMortonLink();
    void setDomainSize();
    void setOptimumLevel();
    void presortMorton();
    void morton();
    void sort(int domainNum, int toggle);
    void sortParticles();
    void unsortParticles();
    void countNonEmptyBoxes();
    void getBoxData(int& numBoxIndex, int toggle);
    void getBoxDataOfParent(int& numBoxIndex, int numLevel, int toggle);
    void getBoxIndexMask(int numBoxIndex, int numLevel, int toggle);
    void getInteractionList(int numBoxIndex, int numLevel, FmmInteractionType interactionType);
    int findDomain(int numBoxIndex, int numLevel, int minIndex = 0);
    int findMortonIndex(int tempBoxIndex, int numLevel, int minIndex = 0);
    void rotationMorton(vec3<int> boxIndex3D, int& forwardIndex, int& backwardIndex, 
            int numLevel = 1, int boxSize = 11);
    void sortOutliers(int domainNum, int toggle);
    void rotation(std::complex<double>* Cnm, std::complex<double>* CnmOut, std::complex<double>** Dnm);
    
    void direct();
    static void* direct_threaded(void* dF);
    void direct_kernel(int threadNum = 0);
    
    static void* direct_threaded_rotations(void* dF);
    void direct_kernel_rotations(int threadNum = 0);
    
    static void* direct_threaded_translations(void* dF);
    void direct_kernel_translations(int threadNum = 0);

    void p2p(int numBoxIndex, int toggle);
    static void* p2p_threaded(void* dF);
    void p2p_kernel(int threadNum = 0);
    
    static void* p2p_threaded_rotations(void* dF);
    void p2p_kernel_rotations(int threadNum = 0);
    
    static void* p2p_threaded_translations(void* dF);
    void p2p_kernel_translations(int threadNum = 0);

    void p2m(int numBoxIndex, int toggle);
    static void* p2m_threaded(void* dF);
    void p2m_kernel(int threadNum = 0);
    
    static void* p2m_threaded_translations(void* dF);
    void p2m_kernel_translations(int threadNum = 0);

    void m2m(int numBoxIndex, int numBoxIndexOld, int numLevel, int toggle);
    static void* m2m_threaded(void* dF);
    void m2m_kernel(int threadNum = 0);
    
    static void* m2m_threaded_rotations(void* dF);
    void m2m_kernel_rotations(int threadNum = 0);
    
    static void* m2m_threaded_translations(void* dF);
    void m2m_kernel_translations(int threadNum = 0);

    void m2l(int numBoxIndex, int numLevel, int toggle, int topLevel = 0);
    static void* m2l_threaded(void* dF);
    void m2l_kernel(int threadNum = 0);
    
    static void* m2l_threaded_rotations(void* dF);
    void m2l_kernel_rotations(int threadNum = 0);
    
    static void* m2l_threaded_translations(void* dF);
    void m2l_kernel_translations(int threadNum = 0);

    void l2l(int numBoxIndex, int numLevel, int toggle);
    static void* l2l_threaded(void* dF);
    void l2l_kernel(int threadNum = 0);
    
    static void* l2l_threaded_translations(void* dF);
    void l2l_kernel_translations(int threadNum = 0);
    
    void l2p(int numBoxIndex, int toggle);
    static void* l2p_threaded(void* dF);
    void l2p_kernel(int threadNum = 0);
    
    static void* l2p_threaded_translations(void* dF);
    void l2p_kernel_translations(int threadNum = 0);

    void m2p(int numBoxIndex, int numLevel, int toggle);
    static void* m2p_threaded(void* dF);
    void m2p_kernel(int threadNum = 0);
    
    static void* m2p_threaded_rotations(void* dF);
    void m2p_kernel_rotations(int threadNum = 0);
    
    static void* m2p_threaded_translations(void* dF);
    void m2p_kernel_translations(int threadNum = 0);

    static void* sm2m_threaded(void* dF);
    void sm2m_kernel(int threadNum = 0);
    
    static void* sm2m_threaded_translations(void* dF);
    void sm2m_kernel_translations(int threadNum = 0);

    static void* sm2l_threaded(void* dF);
    void sm2l_kernel(int threadNum = 0);
    
    static void* sm2l_threaded_translations(void* dF);
    void sm2l_kernel_translations(int threadNum = 0);

    static void* sl2l_threaded(void* dF);
    void sl2l_kernel(int threadNum = 0);

    static void* sm2p_threaded(void* dF);
    void sm2p_kernel(int threadNum = 0);
    
    static void* sm2p_threaded_translations(void* dF);
    void sm2p_kernel_translations(int threadNum = 0);

    static void* sp2p_threaded(void* dF);
    void sp2p_kernel(int threadNum = 0);
    
    static void* sp2p_threaded_rotations(void* dF);
    void sp2p_kernel_rotations(int threadNum = 0);
    
    static void* sp2p_threaded_translations(void* dF);
    void sp2p_kernel_translations(int threadNum = 0);

    static void* sp2m_threaded(void* dF);
    void sp2m_kernel(int threadNum = 0);

    static void* sl2p_threaded(void* dF);
    void sl2p_kernel(int threadNum = 0);

    void outlier_p2p();

    void p2outlier_p(int numBoxIndex, int toggle);
    static void* p2outlier_p_threaded(void* dF);
    void p2outlier_p_kernel(int threadNum = 0);
    
    static void* p2outlier_p_threaded_rotations(void* dF);
    void p2outlier_p_kernel_rotations(int threadNum = 0);
    
    static void* p2outlier_p_threaded_translations(void* dF);
    void p2outlier_p_kernel_translations(int threadNum = 0);

    static void* sp2outlier_p_threaded(void* dF);
    void sp2outlier_p_kernel(int threadNum = 0);
    
    static void* sp2outlier_p_threaded_translations(void* dF);
    void sp2outlier_p_kernel_translations(int threadNum = 0);

    void m2outlier_p(int numBoxIndex, int numLevel, int toggle);
    static void* m2outlier_p_threaded(void* dF);
    void m2outlier_p_kernel(int threadNum = 0);
    
    static void* m2outlier_p_threaded_rotations(void* dF);
    void m2outlier_p_kernel_rotations(int threadNum = 0);
    
    static void* m2outlier_p_threaded_translations(void* dF);
    void m2outlier_p_kernel_translations(int threadNum = 0);
    
    static void* sm2outlier_p_threaded(void* dF);
    void sm2outlier_p_kernel(int threadNum = 0);
    
    static void* sm2outlier_p_threaded_translations(void* dF);
    void sm2outlier_p_kernel_translations(int threadNum = 0);

    MultipoleMethodConstants* constants = nullptr;
    MultipoleMethodVariables* variables = nullptr;
    FmmType fmmType;

    int maxLevel;
    int numBoxIndexFull;
    int numBoxIndexLeaf;
    int numBoxIndexTotal;
    float rootBoxSize;
    vec3<float> boxMin;
    int boxDimVals[3];
    int domainNumBoxIndexFull;
    int numDomains;
    int numLevelSC;
    int numOutliers;
    int numOutlierDomains;
    int outlierNumBoxIndex;
    int outlierNumBoxIndexLeaf;
    int numBoxIndexClass;

    FILE* boxInfo = nullptr;
        
    int treeOrFMM;
    int outlierToggle;
        
    multipoleMethodTiming timing;
    
    void allocateThreadMemory();
    void deallocateThreadMemory();

    int vectorizedToggle = false;
    Jpdata *jptc = nullptr;
#ifdef _PTHREAD_H
    multipoleMethodForceStruct* dF = nullptr;
    pthread_mutex_t mutex_fmm = PTHREAD_MUTEX_INITIALIZER;
#endif
};

#endif /* FORCE_CALCULATION_METHODS_H */

