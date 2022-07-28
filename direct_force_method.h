/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   direct_force_method.h
 * Author: Marshall
 *
 * Created on June 12, 2022, 5:11 PM
 */

#ifndef DIRECT_FORCE_METHOD_H
#define DIRECT_FORCE_METHOD_H

#include <cassert>
#include "vectorized_functions.h"

#define AX      "%xmm0"
#define AY      "%xmm1"
#define AZ      "%xmm2"
#define QI      "%xmm3"

#define XJ      "%xmm4"
#define YJ      "%xmm5"
#define ZJ      "%xmm6"
#define QJ      "%xmm7"

#define RINV    "%xmm8"
#define X2      "%xmm9"
#define Y2      "%xmm10"
#define Z2      "%xmm11"

#define XI      "%xmm12"
#define YI      "%xmm13"
#define ZI      "%xmm14"
#define R2      "%xmm15"

namespace DirectForceFunction {
//    void direct(int i_min, int i_max, int j_min, int j_max, 
//            ChargedParticles& particles);
    
    void setJpdata(int& tempNumParticles, ChargedParticles& particles, Jpdata* jptc);
    
    void direct(vec4<int> rangevalues, ChargedParticles& particles);
        
    void direct_vectorized(vec4<int> rangevalues, ChargedParticles& particles, Jpdata *jptc);
    
    inline void direct_vectorized_kernel(Ipdata &ipdata, Fodata &fodata, Jpdata *jpdata,
            int nj, int offset = 0){
        int j;

        assert(((unsigned long)jpdata & 15) == 0);
        assert(((unsigned long)&ipdata & 15) == 0);
        assert(((unsigned long)&fodata & 15) == 0);

        XORPS(AX, AX);                      // AX = {0, 0, 0, 0};
        XORPS(AY, AY);                      // AY = {0, 0, 0, 0};
        XORPS(AZ, AZ);                      // AZ = {0, 0, 0, 0};

        LOADPS(ipdata.x, XI);               // XI = {ipdata.x[0], ipdata.x[1], ipdata.x[2], ipdata.x[3]};
        LOADPS(ipdata.y, YI);               // YI = {ipdata.y[0], ipdata.y[1], ipdata.y[2], ipdata.y[3]};
        LOADPS(ipdata.z, ZI);               // ZI = {ipdata.z[0], ipdata.z[1], ipdata.z[2], ipdata.z[3]};

        LOADPS(fodata.eps2, R2);            // R2 = {eps2, eps2, eps2, eps2};
        MOVAPS(R2, QI);

    //    LOADPS(jpdata[0], QJ);              // QJ = {jpdata[j].x, jpdata[j].y, jpdata[j].z, jpdata[j].q};
        LOADPS(jpdata[offset], QJ);
        MOVAPS(QJ, X2);                     // X2 = QJ;
        MOVAPS(QJ, Y2);                     // Y2 = QJ;
        MOVAPS(QJ, Z2);                     // Z2 = QJ;

        BCAST0(X2);                         // X2 = {jpdata[j].x, jpdata[j].x, jpdata[j].x, jpdata[j].x};
        BCAST1(Y2);                         // Y2 = {jpdata[j].y, jpdata[j].y, jpdata[j].y, jpdata[j].y};
        BCAST2(Z2);                         // Z2 = {jpdata[j].z, jpdata[j].z, jpdata[j].z, jpdata[j].z};
        BCAST3(QJ);                         // QJ = {jpdata[j].q, jpdata[j].q, jpdata[j].q, jpdata[j].q};

        SUBPS(XI, X2);                      // X2 = X2 - XI;
        SUBPS(YI, Y2);                      // Y2 = Y2 - YI;
        SUBPS(ZI, Z2);                      // Z2 = Z2 - ZI;

        MOVAPS(X2, XJ);                     // XJ = X2;
        MOVAPS(Y2, YJ);                     // YJ = Y2;
        MOVAPS(Z2, ZJ);                     // ZJ = Z2;

        MULPS(X2, X2);                      // X2 = X2 * X2;
        MULPS(Y2, Y2);                      // Y2 = Y2 * Y2;
        MULPS(Z2, Z2);                      // Z2 = Z2 * Z2;

        ADDPS(X2, R2);                      // R2 = R2 + X2;
        ADDPS(Y2, R2);                      // R2 = R2 + Y2;
        ADDPS(Z2, R2);                      // R2 = R2 + Z2;

    //    LOADPS(jpdata[1], X2);              // X2 = {jpdata[j].x, jpdata[j].y, jpdata[j].z, jpdata[j].q};
        LOADPS(jpdata[offset + 1], X2);
        MOVAPS(X2, Y2);                     // Y2 = X2;
        MOVAPS(X2, Z2);                     // Z2 = X2;

    //    for (j = 1; j <= nj; j++){
        for (j = offset + 1; j <= nj + offset; j++){

            CMPNEQPS(R2, QI);               // QI = (R2 != QI);

            RSQRTPS(R2, RINV);              // RINV = rsqrt(R2);
            ANDPS(QI, RINV);                // RINV = (RINV & QI);

            LOADPS(fodata.eps2, R2);
            MOVAPS(R2, QI);
            BCAST0(X2);
            BCAST1(Y2);
            BCAST2(Z2);
            SUBPS(XI, X2);                  // X2 = X2 - XI;
            SUBPS(YI, Y2);                  // Y2 = Y2 - YI;
            SUBPS(ZI, Z2);                  // Z2 = Z2 - ZI;

            MULPS(RINV, QJ);                // QJ = QJ * RINV;
            MULPS(RINV, RINV);              // RINV = RINV * RINV;
            MULPS(QJ, RINV);                // RINV = QJ * RINV;
            LOADPS(jpdata[j], QJ);
            BCAST3(QJ);

            MULPS(RINV, XJ);                // XJ = XJ * RINV;

            SUBPS(XJ, AX);                  // AX = AX - XJ;
            MOVAPS(X2, XJ);
            MULPS(X2, X2);
            ADDPS(X2, R2);
            LOADPS(jpdata[j + 1], X2);

            MULPS(RINV, YJ);                // YJ = YJ * RINV;
            SUBPS(YJ, AY);                  // AY = AY - YJ;
            MOVAPS(Y2, YJ);
            MULPS(Y2, Y2);
            ADDPS(Y2, R2);
            MOVAPS(X2, Y2);

            MULPS(RINV, ZJ);                // ZJ = ZJ * RINV;
            SUBPS(ZJ, AZ);                  // AZ = AZ - ZJ;
            MOVAPS(Z2, ZJ);
            MULPS(Z2, Z2);
            ADDPS(Z2, R2);
            MOVAPS(X2, Z2);
        }

        LOADPS(ipdata.q, QI);               // QI = {ipdata.q[0], ipdata.q[1], ipdata.q[2], ipdata.q[3]};
        MULPS(QI, AX);                      // AX = AX * QI;
        MULPS(QI, AY);                      // AY = AY * QI;
        MULPS(QI, AZ);                      // AZ = AZ * QI;

        STORPS(AX, fodata.ax);              // fodata.ax = AX;
        STORPS(AY, fodata.ay);              // fodata.ay = AY;
        STORPS(AZ, fodata.az);              // fodata.az = AZ;
    }
    
    void direct_rotational(int quadrant, vec4<int> rangevalues, ChargedParticles& particles);
    
    void direct_vectorized_rotational(int quadrant, vec4<int> rangevalues, ChargedParticles& particles, Jpdata *jptc);
    
//    void direct_translational(volatile const vec3<float> domainOffset, vec4<int> rangevalues, ChargedParticles& particles);
    
//    void direct_vectorized_translational(volatile const vec3<float> domainOffset, vec4<int> rangevalues, ChargedParticles& particles, Jpdata *jptc);
    
    void direct_translational(LatticePoint lattice, vec3<float> domainSize,
            vec4<int> rangevalues, ChargedParticles& particles);
    
    void direct_vectorized_translational(LatticePoint lattice, vec3<float> domainSize,
            vec4<int> rangevalues, ChargedParticles& particles, Jpdata* jptc);
}

#endif /* DIRECT_FORCE_METHOD_H */

