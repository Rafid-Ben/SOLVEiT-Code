/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   randomizer.h
 * Author: Marshall
 *
 * Created on January 22, 2022, 4:13 PM
 */

#ifndef RANDOMIZER_H
#define RANDOMIZER_H

#include <random>

class Randomizer {
public:
    Randomizer(): mt_gen{std::random_device()()}, rnd_dist{0, 1.0} {}
    Randomizer(unsigned int seed): mt_gen{seed}, rnd_dist{0, 1.0} {}
    double operator() () {return rnd_dist(mt_gen);}
protected:
    std::mt19937 mt_gen;
    std::uniform_real_distribution<double> rnd_dist;
};

class RandomizerNew {
public:
    RandomizerNew() : mt_gen{std::random_device()()}, rnd_dist{0, 1.0} {}
    RandomizerNew(unsigned int seed): mt_gen{seed}, rnd_dist{0, 1.0} {}
    double uniformRealDistribution() {return rnd_dist(mt_gen);}
    double exponentialDistribution() {return exp_dist(mt_gen);}
    void setExponentialDistributionLambda(double lambda){
        std::exponential_distribution<double>::param_type new_lambda(lambda);
        exp_dist.param(new_lambda);
    };
protected:
    std::mt19937 mt_gen;
    std::uniform_real_distribution<double> rnd_dist;
    std::exponential_distribution<double> exp_dist;
};

extern Randomizer randomizer;
extern RandomizerNew randomizerNew;

#endif /* RANDOMIZER_H */

