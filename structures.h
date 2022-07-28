/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   structures.h
 * Author: Marshall
 *
 * Created on June 8, 2022, 2:50 PM
 */

#ifndef STRUCTURES_H
#define STRUCTURES_H
#include <cassert>
#include <cstdio>
#include <cmath>
#include <algorithm>
#include "timer.h"

template<typename T>
class vec3 {
public:
    T x;
    T y;
    T z;
    
    vec3 (): x{0}, y{0}, z{0} {}
    vec3(const T x, const T y, const T z) : x{x}, y{y}, z{z} {}
    vec3(const T a[3]) : x{a[0]}, y{a[1]}, z{a[2]} {}
    vec3(const T c) : x{c}, y{c}, z{c} {}
    vec3<T>& operator=(double s) {x = s; y = s; z = s; return (*this);}
    vec3<T>& operator+=(vec3<T> o) {x += o.x; y += o.y; z += o.z; return (*this);}
    vec3<T>& operator+=(double s) {x += s; y += s; z += s; return (*this);}
    vec3<T>& operator*=(vec3<T> o) {x *= o.x; y *= o.y; z *= o.z; return (*this);}
    vec3<T>& operator*=(double s) {x *= s; y *= s; z *= s; return (*this);}
    vec3<T>& operator/=(vec3<T> o) {x /= o.x; y /= o.y; z /= o.z; return (*this);}
    vec3<T>& operator/=(double s) {x /= s; y /= s; z /= s; return (*this);}
    vec3<T>& operator-=(vec3<T> o) {x -= o.x; y -= o.y; z -= o.z; return (*this);}
    vec3<T>& operator-=(double s) {x -= s; y -= s; z -= s; return (*this);}
    
    vec3<T>& operator=(const vec3<T>& other){
        if (this != &other){
            this->x = other.x;
            this->y = other.y;
            this->z = other.z;
        }
        return *this;
    }
    
    T& operator[] (int i) {
        assert(i >= 0 && i < 3);
        if (i == 0){
            return this->x;
        }
        else if (i == 1){
            return this->y;
        }
        else {
            return this->z;
        }
    }
    T& operator[] (int i) const{
        assert(i >= 0 && i < 3);
        if (i == 0){
            return this->x;
        }
        else if (i == 1){
            return this->y;
        }
        else {
            return this->z;
        }
    }
};

template<typename T>
vec3<T> operator+(const vec3<T>& a, const vec3<T>& b){
    return vec3<T> (a.x + b.x, a.y + b.y, a.z + b.z);
}
template<typename T>
vec3<T> operator*(const vec3<T>& a, const vec3<T>& b){
    return vec3<T> (a.x * b.x, a.y * b.y, a.z * b.z);
}
template<typename T>
vec3<T> operator/(const vec3<T>& a, const vec3<T>& b){
    return vec3<T> (a.x / b.x, a.y / b.y, a.z / b.z);
}
template<typename T>
vec3<T> operator-(const vec3<T>& a, const vec3<T>& b){
    return vec3<T> (a.x - b.x, a.y - b.y, a.z - b.z);
}

template<typename T, typename S>
vec3<T> operator+(const vec3<T>& a, const vec3<S>& b){
    return vec3<T> (a.x + b.x, a.y + b.y, a.z + b.z);
}
template<typename T, typename S>
vec3<T> operator*(const vec3<T>& a, const vec3<S>& b){
    return vec3<T> (a.x * b.x, a.y * b.y, a.z * b.z);
}
template<typename T, typename S>
vec3<T> operator/(const vec3<T>& a, const vec3<S>& b){
    return vec3<T> (a.x / b.x, a.y / b.y, a.z / b.z);
}
template<typename T, typename S>
vec3<T> operator-(const vec3<T>& a, const vec3<S>& b){
    return vec3<T> (a.x - b.x, a.y - b.y, a.z - b.z);
}

template<typename T>
vec3<T> operator*(const vec3<T> &a, T s){
    return vec3<T>(a.x * s, a.y * s, a.z * s);
}
template<typename T>
vec3<T> operator*(T s, const vec3<T> & a){
    return vec3<T>(a.x * s, a.y * s, a.z * s);
}
template<typename T>
vec3<T> operator+(const vec3<T> &a, T s){
    return vec3<T>(a.x + s, a.y + s, a.z + s);
}
template<typename T>
vec3<T> operator+(T s, const vec3<T> & a){
    return vec3<T>(a.x + s, a.y + s, a.z + s);
}
template<typename T>
vec3<T> operator-(const vec3<T> &a, T s){
    return vec3<T>(a.x - s, a.y - s, a.z - s);
}
template<typename T>
vec3<T> operator-(T s, const vec3<T> & a){
    return vec3<T>(s - a.x, s - a.y, s - a.z);
}
template<typename T>
vec3<T> operator/(const vec3<T> &a, T s){
    return vec3<T>(a.x / s, a.y / s, a.z / s);
}
template<typename T>
vec3<T> operator/(T s, const vec3<T> & a){
    return vec3<T>(s / a.x, s / a.y, s / a.z);
}

template<typename T, typename S>
vec3<T> operator*(const vec3<T> &a, S s){
    return vec3<T>(a.x * s, a.y * s, a.z * s);
}
template<typename T, typename S>
vec3<T> operator*(S s, const vec3<T> & a){
    return vec3<T>(a.x * s, a.y * s, a.z * s);
}
template<typename T, typename S>
vec3<T> operator+(const vec3<T> &a, S s){
    return vec3<T>(a.x + s, a.y + s, a.z + s);
}
template<typename T, typename S>
vec3<T> operator+(S s, const vec3<T> & a){
    return vec3<T>(a.x + s, a.y + s, a.z + s);
}
template<typename T, typename S>
vec3<T> operator-(const vec3<T> &a, S s){
    return vec3<T>(a.x - s, a.y - s, a.z - s);
}
template<typename T, typename S>
vec3<T> operator-(S s, const vec3<T> & a){
    return vec3<T>(s - a.x, s - a.y, s - a.z);
}
template<typename T, typename S>
vec3<T> operator/(const vec3<T> &a, S s){
    return vec3<T>(a.x / s, a.y / s, a.z / s);
}
template<typename T, typename S>
vec3<T> operator/(S s, const vec3<T> & a){
    return vec3<T>(s / a.x, s / a.y, s / a.z);
}

template<typename T>
class vec4 {
public:
    T x;
    T y;
    T z;
    T w;
    
    vec4 (): x{0}, y{0}, z{0}, w{0} {}
    vec4(const T x, const T y, const T z, const T w) : x{x}, y{y}, z{z}, w{w} {}
    vec4(const T a[4]) : x{a[0]}, y{a[1]}, z{a[2]}, w{a[3]} {}
    vec4(const T c) : x{c}, y{c}, z{c}, w{c} {}
};




enum class BackgroundEfield {
    singleEmitter,
    none
};

enum class SymmetryType {
    rotational,
    translational,
    none
};

enum class ParticleStatus {
    inbound = 1,
    crashed = -1,
    outOfBounds = 0,
    affectsAcceleration = 2
};

enum class Initialization {
    paraxialRay,
    singleEmitter,
    none
};

enum class InjectionTimingOptions {
    constantRate,
    exponentialDistribution
};

enum class ParticleSpecies {
    monomerPositive = 0,
    dimerPositive = 1,
    trimerPositive = 2,
    droplet = 3,
    neutral = 4,
    monomerNegative = 5,
    dimerNegative = 6,
    trimerNegative = 7
};

enum class ParticleProportionOptions {
    byNumber,
    byCharge
};

enum class FragmentationModel {
    none
};

enum class SaveDataOptions {
    original,
    newOriginal
};

enum class TimerOptions {
    verbose,
    minimized,
    endOnly,
    off
};

enum class BoundaryConditionOptions {
    none,
    space,
    time,
    space_and_time
};

enum class DomainToggle {
    cubicCells,
    rectangularPrismCells
};

enum class FmmInteractionType {
    P2P,
    P2P_Supercell,
    P2P_Outliers,
    P2P_SupercellOutliers,
    M2L_SubcellTopLevel,
    M2L_SubcellTopLevelOutliers,
    M2L_Subcell,
    M2L_SubcellOutliers,
    M2L_Supercell,
    M2L_SupercellOutliers,
    M2L_SupercellTopLevel,
    M2L_SupercellTopLevelOutliers
};

enum class FmmType {
    treecode,
    treecodeOutliers,
    fmm,
    fmmOutliers
};

class BoundaryConditionParameters {
public:
    BoundaryConditionParameters() {};
    BoundaryConditionOptions boundaryConditionType = BoundaryConditionOptions::none;
    float temporalBoundaryCondition = 0.0f;
    vec4<float> spatialBoundaryConditions = 0.0f;
    
    BoundaryConditionParameters(const BoundaryConditionParameters& bc){
        boundaryConditionType = bc.boundaryConditionType;
        temporalBoundaryCondition = bc.temporalBoundaryCondition;
        spatialBoundaryConditions = bc.spatialBoundaryConditions;
    }
    
    BoundaryConditionParameters& operator=(const BoundaryConditionParameters& bc){
        boundaryConditionType = bc.boundaryConditionType;
        temporalBoundaryCondition = bc.temporalBoundaryCondition;
        spatialBoundaryConditions = bc.spatialBoundaryConditions;
        return *this;
    }
    
    void setMostStringentBoundaryConditionType(){
        switch (boundaryConditionType){
            case (BoundaryConditionOptions::none):
                temporalBoundaryCondition = 0.0f;
                spatialBoundaryConditions = 0.0f;
                break;
            case (BoundaryConditionOptions::space):
                temporalBoundaryCondition = 0.0f;
                break;
            case (BoundaryConditionOptions::time):
                spatialBoundaryConditions = 0.0f;
                break;
            case (BoundaryConditionOptions::space_and_time):
                break;
        }
        
        boundaryConditionType = BoundaryConditionOptions::space_and_time;
    }
    
    void chooseMinimumBoundaryConditionParameters(const BoundaryConditionParameters& bc){
        
        if (temporalBoundaryCondition == 0 && bc.temporalBoundaryCondition == 0){
            temporalBoundaryCondition = 0.0f;
        }
        else {
            temporalBoundaryCondition = std::min(temporalBoundaryCondition, bc.temporalBoundaryCondition);
        }
        
        if (spatialBoundaryConditions.x == 0 && bc.spatialBoundaryConditions.x == 0){
            spatialBoundaryConditions.x = 0.0f;
        }
        else {
            spatialBoundaryConditions.x = std::min(spatialBoundaryConditions.x, bc.spatialBoundaryConditions.x);
        }
        
        if (spatialBoundaryConditions.y == 0 && bc.spatialBoundaryConditions.y == 0){
            spatialBoundaryConditions.y = 0.0f;
        }
        else {
            spatialBoundaryConditions.y = std::min(spatialBoundaryConditions.y, bc.spatialBoundaryConditions.y);
        }
        
        if (spatialBoundaryConditions.z == 0 && bc.spatialBoundaryConditions.z == 0){
            spatialBoundaryConditions.z = 0.0f;
        }
        else {
            spatialBoundaryConditions.z = std::min(spatialBoundaryConditions.z, bc.spatialBoundaryConditions.z);
        }
        
        if (spatialBoundaryConditions.w == 0 && bc.spatialBoundaryConditions.w == 0){
            spatialBoundaryConditions.w = 0.0f;
        }
        else {
            spatialBoundaryConditions.w = std::min(spatialBoundaryConditions.w, bc.spatialBoundaryConditions.w);
        }
    }
    
    void chooseMaximumBoundaryConditionParameters(const BoundaryConditionParameters& bc){
        if (temporalBoundaryCondition == 0 || bc.temporalBoundaryCondition == 0){
            temporalBoundaryCondition = 0.0f;
        }
        else {
            temporalBoundaryCondition = std::max(temporalBoundaryCondition, bc.temporalBoundaryCondition);
        }
        
        if (spatialBoundaryConditions.x == 0 || bc.spatialBoundaryConditions.x == 0){
            spatialBoundaryConditions.x = 0.0f;
        }
        else {
            spatialBoundaryConditions.x = std::max(spatialBoundaryConditions.x, bc.spatialBoundaryConditions.x);
        }
        
        if (spatialBoundaryConditions.y == 0 || bc.spatialBoundaryConditions.y == 0){
            spatialBoundaryConditions.y = 0.0f;
        }
        else {
            spatialBoundaryConditions.y = std::max(spatialBoundaryConditions.y, bc.spatialBoundaryConditions.y);
        }
        
        if (spatialBoundaryConditions.z == 0 || bc.spatialBoundaryConditions.z == 0){
            spatialBoundaryConditions.z = 0.0f;
        }
        else {
            spatialBoundaryConditions.z = std::max(spatialBoundaryConditions.z, bc.spatialBoundaryConditions.z);
        }
        
        if (spatialBoundaryConditions.w == 0 || bc.spatialBoundaryConditions.w == 0){
            spatialBoundaryConditions.w = 0.0f;
        }
        else {
            spatialBoundaryConditions.w = std::max(spatialBoundaryConditions.w, bc.spatialBoundaryConditions.w);
        }
    }
};

class SubstepBoundaryConditions {
public:
    SubstepBoundaryConditions() {};
    float dt = 0.0f;
    BoundaryConditionParameters bc;
    
    SubstepBoundaryConditions(const SubstepBoundaryConditions& substepBC){
        dt = substepBC.dt;
        bc = substepBC.bc;
    }
    
    SubstepBoundaryConditions& operator=(const SubstepBoundaryConditions& substepBC){
        dt = substepBC.dt;
        bc = substepBC.bc;
        return *this;
    }
    
    bool operator < (SubstepBoundaryConditions const & rhs) const {
        printf("dt = %.9g, rhs.dt = %.9g\n",dt, rhs.dt);
        printf("%d\n", (int)(dt < rhs.dt));
        return (dt < rhs.dt);
    }
    
    bool operator > (SubstepBoundaryConditions const & rhs) const {
        return (dt > rhs.dt);
    }
};


enum class CrashConditionOptions {
    none,
    extractorGrid
};

enum class SaveDataRateOptions {
    constantRate,
    lastSnapshots,
    constantRate_and_lastSnapshots
};

enum class ParticleType{
    charged,
    neutral
};

class CrashParameters {
public:
    CrashConditionOptions crashOptions = CrashConditionOptions::none;
    float z_distance = 0.0; // Distance from zero to bottom of grid
    float z_width = 0.0; // Grid width
    float r_distance = 0.0; // Radial distance from centerline
    int checkCrash(vec3<float> initialPos, vec3<float>& finalPos){
        switch (crashOptions){
            case (CrashConditionOptions::none):
                return false;
                break;
            case (CrashConditionOptions::extractorGrid):
                return checkExtractorGridCrash(initialPos, finalPos);
                break;
        }
        return false;
    }
private:
    int checkExtractorGridCrash(vec3<float> initialPos, vec3<float>& finalPos){
        float initial_r = sqrt(initialPos.x * initialPos.x + initialPos.y * initialPos.y);
        float final_r = sqrt(finalPos.x * finalPos.x + finalPos.y * finalPos.y);
        float m = (finalPos.z - initialPos.z) / (final_r - initial_r);
        float b = initialPos.z - m * initial_r;
        
        float min_r = std::min(initial_r, final_r);
        float max_r = std::max(initial_r, final_r);
        float min_z = std::min(initialPos.z, finalPos.z);
        float max_z = std::max(initialPos.z, finalPos.z);
        
        float temp_r, temp_z;
        
        
        // Check intersection with bottom of grid
        if (z_distance >= min_z && z_distance <= max_z){
            temp_r = (z_distance - b) / m;
//            if (temp >= min_r && temp <= max_r && temp >= r_distance){
            if (temp_r >= std::max(min_r, r_distance) && temp_r <= max_r){
                temp_z = z_distance;
                finalPos.z = temp_z;
                return true; // Crashed
            }
        }
        
        // Check intersection with top of grid
        if (z_distance + z_width >= min_z && z_distance + z_width <= max_z){
            temp_r = (z_distance + z_width - b) / m;
//            if (temp >= min_r && temp <= max_r && temp >= r_distance){
            if (temp_r >= std::max(min_r, r_distance) && temp_r <= max_r){
                temp_z = z_distance + z_width;
                finalPos.z = temp_z;
                return true; // Crashed
            }
        }
        
        // Check intersection with side of grid
        if (r_distance >= min_r && r_distance <= max_r){
            temp_z = m * r_distance + b;
//            if (temp >= z_distance && temp <= z_distance + z_width && 
//                    temp >= min_z && temp <= max_z){
            if (temp_z >= std::max(min_z, z_distance) && 
                    temp_z <= std::min(max_z, z_distance + z_width)){
                temp_r = r_distance;
                finalPos.z = temp_z;
                return true; // Crashed
            }
        }
        
        return false;
    }
};

class ParticleFractions {
public:
    ParticleFractions(int maxNumSpecies) : maxNumSpecies{maxNumSpecies} {
        allocateMemory();
    };
    ParticleFractions(Initialization initialization){
        switch (initialization){
            case (Initialization::none):
                maxNumSpecies = 0;
                allocateMemory();
                break;
            case (Initialization::paraxialRay):
                maxNumSpecies = 1;
                allocateMemory();
                species[0] = ParticleSpecies::monomerPositive;
                proportion[0] = 1.0;
                break;
            case (Initialization::singleEmitter):
                maxNumSpecies = 1;
                allocateMemory();
                species[0] = ParticleSpecies::monomerPositive;
                proportion[0] = 1.0;
                break;
        }
        generateCDF();
    }
    ~ParticleFractions(){
        freeMemory();
    }
    
    ParticleFractions(const ParticleFractions& p){
        int i;
        this->maxNumSpecies = p.maxNumSpecies;
        this->species = new ParticleSpecies[p.maxNumSpecies];
        this->proportion = new float[p.maxNumSpecies];
        this->cdf = new float[p.maxNumSpecies];
        this->proportionType = p.proportionType;
        for (i = 0; i < p.maxNumSpecies; i++){
            this->species[i] = p.species[i];
            this->proportion[i] = p.proportion[i];
            this->cdf[i] = p.cdf[i];
        }
    }
    
    void generateCDF(){
        if (cdf != nullptr) delete[] cdf;
        cdf = new float[maxNumSpecies + 1];
        cdf[0] = 0.0;
        int i;
        for (i = 0; i < maxNumSpecies; i++){
            cdf[i + 1] = cdf[i] + proportion[i];
        }
    }
    
    ParticleSpecies selectSpecies(float randVal){
        int maxIndex = maxNumSpecies - 1;
        int minIndex = 0;
        int guessIndex;
        int indexFound = 0;
        
        while (indexFound == 0){
            guessIndex = (minIndex + maxIndex) / 2;
            if (randVal >= cdf[guessIndex + 1]){
                minIndex = guessIndex + 1;
            }
            else if (randVal < cdf[guessIndex]){
                maxIndex = guessIndex;
            }
            else {
                indexFound = 1;
            }
            
            if (maxIndex < minIndex){
                indexFound = 1;
            }
        }
        
        return species[guessIndex];
        
    }
    
    void allocateMemory(){
        freeMemory();
        species = new ParticleSpecies[maxNumSpecies];
        proportion = new float[maxNumSpecies];
    }
    
    void freeMemory(){
        if (species != nullptr) delete[] species;
        if (proportion != nullptr) delete[] proportion;
        if (cdf != nullptr) delete[] cdf;
        species = nullptr;
        proportion = nullptr;
        cdf = nullptr;
    }
    
    
    
    int maxNumSpecies = 0;
    ParticleSpecies* species = nullptr;
    float* proportion = nullptr;
    float* cdf = nullptr;
    ParticleProportionOptions proportionType = ParticleProportionOptions::byNumber;
};


typedef struct {
    float setup_time = 0;
    float interactionList_time = 0;
    float direct_time = 0;
    float p2p_time = 0;
    float p2m_time = 0;
    float m2m_time = 0;
    float m2l_time = 0;
    float l2l_time = 0;
    float l2p_time = 0;
    float m2p_time = 0;
    float sm2m_time = 0;
    float sm2l_time = 0;
    float sl2l_time = 0;
    float sm2p_time = 0;
    float sp2p_time = 0;
    float sp2m_time = 0;
    float sl2p_time = 0;
    float outlier_p2p_time = 0;
    float p2outlier_p_time = 0;
    float sp2outlier_p_time = 0;
    float m2outlier_p_time = 0;
    float sm2outlier_p_time = 0;
    float total_time = 0;
    
    Timer timer;
    Timer total_timer;
    void writeTiming(FILE* filename){
        if (filename != nullptr){
            fwrite(&setup_time, sizeof(*&setup_time), 1, filename);
            fwrite(&interactionList_time, sizeof(*&interactionList_time), 1, filename);
            fwrite(&direct_time, sizeof(*&direct_time), 1, filename);
            fwrite(&p2p_time, sizeof(*&p2p_time), 1, filename);
            fwrite(&p2m_time, sizeof(*&p2m_time), 1, filename);
            fwrite(&m2m_time, sizeof(*&m2m_time), 1, filename);
            fwrite(&m2l_time, sizeof(*&m2l_time), 1, filename);
            fwrite(&l2l_time, sizeof(*&l2l_time), 1, filename);
            fwrite(&l2p_time, sizeof(*&l2p_time), 1, filename);
            fwrite(&m2p_time, sizeof(*&m2p_time), 1, filename);
            fwrite(&sm2m_time, sizeof(*&sm2m_time), 1, filename);
            fwrite(&sm2l_time, sizeof(*&sm2l_time), 1, filename);
            fwrite(&sl2l_time, sizeof(*&sl2l_time), 1, filename);
            fwrite(&sm2p_time, sizeof(*&sm2p_time), 1, filename);
            fwrite(&sp2p_time, sizeof(*&sp2p_time), 1, filename);
            fwrite(&sp2m_time, sizeof(*&sp2m_time), 1, filename);
            fwrite(&sl2p_time, sizeof(*&sl2p_time), 1, filename);
            fwrite(&outlier_p2p_time, sizeof(*&outlier_p2p_time), 1, filename);
            fwrite(&p2outlier_p_time, sizeof(*&p2outlier_p_time), 1, filename);
            fwrite(&sp2outlier_p_time, sizeof(*&sp2outlier_p_time), 1, filename);
            fwrite(&m2outlier_p_time, sizeof(*&m2outlier_p_time), 1, filename);
            fwrite(&sm2outlier_p_time, sizeof(*&sm2outlier_p_time), 1, filename);
            fwrite(&total_time, sizeof(*&total_time), 1, filename);
        }
    }
    void resetTiming(){
        timer.resetTimer();
        total_timer.resetTimer();
        setup_time = 0;
        interactionList_time = 0;
        direct_time = 0;
        p2p_time = 0;
        p2m_time = 0;
        m2m_time = 0;
        m2l_time = 0;
        l2l_time = 0;
        l2p_time = 0;
        m2p_time = 0;
        sm2m_time = 0;
        sm2l_time = 0;
        sl2l_time = 0;
        sm2p_time = 0;
        sp2p_time = 0;
        sp2m_time = 0;
        sl2p_time = 0;
        outlier_p2p_time = 0;
        p2outlier_p_time = 0;
        sp2outlier_p_time = 0;
        m2outlier_p_time = 0;
        sm2outlier_p_time = 0;
        total_time = 0;
    }
} multipoleMethodTiming;

typedef struct {int mortonIndex, firstParticle, lastParticle, parent;} domainInfo;
typedef struct {int firstParticle, lastParticle; unsigned int mortonIndex;} outlierDomainInfo;
typedef struct {int index; unsigned int signedMorton;} outlierInfo;

#endif /* STRUCTURES_H */

