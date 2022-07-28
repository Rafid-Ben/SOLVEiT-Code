/*
 * To change this license header, choose License Headers in Project Properties.
 * To change this template file, choose Tools | Templates
 * and open the template in the editor.
 */

/* 
 * File:   vectorized_functions.h
 * Author: Marshall
 *
 * Created on February 26, 2022, 7:40 PM
 */

#ifndef VECTORIZED_FUNCTIONS_H
#define VECTORIZED_FUNCTIONS_H

#define XORPS(a, b) asm("xorps " a "," b);
#define LOADPS(mem, reg) asm("movaps %0, %" reg:: "m" (mem));
#define STORPS(reg, mem) asm("movaps %" reg " , %0" :: "m" (mem));
#define MOVAPS(src, dst) asm("movaps " src "," dst);
#define MOVQ(src, dst) asm("movq " src "," dst);
#define BCAST0(reg) asm("shufps $0x00, " reg "," reg);
#define BCAST1(reg) asm("shufps $0x55, " reg "," reg);
#define BCAST2(reg) asm("shufps $0xaa, " reg "," reg);
#define BCAST3(reg) asm("shufps $0xff, " reg "," reg);
#define MULPS(src, dst) asm("mulps " src "," dst);
#define ADDPS(src, dst) asm("addps " src "," dst);
#define SUBPS(src, dst) asm("subps " src "," dst);
#define RSQRTPS(src, dst) asm("rsqrtps " src "," dst);
#define MOVHLPS(src, dst) asm("movhlps " src "," dst);
#define DEBUGPS(reg)
#define CMPUNORDPS(src, dst) asm("cmpunordps " src "," dst);
#define CMPORDPS(src, dst) asm("cmpordps " src "," dst);
#define CMPEQPS(src, dst) asm("cmpeqps " src "," dst);
#define CMPNEQPS(src, dst) asm("cmpneqps " src "," dst);
#define ANDNPS(src, dst) asm("andnps " src "," dst);
#define ANDPS(src, dst) asm("andps " src "," dst);

#define ALIGN16 __attribute__ ((aligned(16)))

typedef double v2df __attribute__ ((vector_size(16)));
typedef float  v4sf __attribute__ ((vector_size(16)));
typedef int    v4si __attribute__ ((vector_size(16)));
typedef short  v8hi __attribute__ ((vector_size(16)));

typedef struct iptdata{
    float x[4];
    float y[4];
    float z[4];
    float q[4];
} Ipdata ALIGN16;

typedef struct fodata{
    float ax[4];
    float ay[4];
    float az[4];
    float eps2[4];
} Fodata ALIGN16;

typedef struct jpdata{
    float x, y, z, q;
} Jpdata ALIGN16;

#endif /* VECTORIZED_FUNCTIONS_H */

