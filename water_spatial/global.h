
/*************************************************************************/
/*                                                                       */
/*  Copyright (c) 1994 Stanford University                               */
/*                                                                       */
/*  All rights reserved.                                                 */
/*                                                                       */
/*  Permission is given to use, copy, and modify this software for any   */
/*  non-commercial purpose as long as this copyright notice is not       */
/*  removed.  All other uses, including redistribution in whole or in    */
/*  part, are forbidden without prior written permission.                */
/*                                                                       */
/*  This software is provided with absolutely no warranty and no         */
/*  support.                                                             */
/*                                                                       */
/*************************************************************************/

/*  This file contains the declaration of the GlobalMemory 
structure and the maximum number of molecules allowed 
by the program. */

#ifndef __GLOBAL_H__
#define __GLOBAL_H__
#define MAX_THREADS 16
#define MAX_BOX_PER_SIDE 100
#define MAX_NS_CUBED 9*9*9
#define MAXMOLS	1728

#include "include/src/tm_basics.h"
typedef tm_type<double> tm_double; // AICI !!!!!

struct GlobalMemory {
    int (IOLock);
    int (IndexLock);
    int (IntrafVirLock);
    int (InterfVirLock);
    int (KinetiSumLock);
    int (PotengSumLock);
    int (start);
    int (InterfBar);
    int (PotengBar);
    int MolLock[MAXMOLS];
    int Index;
    tm_double VIR;
    tm_double SUM[3];
    tm_double POTA, POTR, POTRF;
    //double VIR;
    //double SUM[3];
    //double POTA, POTR, POTRF;
    unsigned long createstart,createend,computestart,computeend;
    unsigned long trackstart, trackend, tracktime;
    unsigned long intrastart, intraend, intratime;
    unsigned long interstart, interend, intertime;
};

extern struct GlobalMemory *gl;

/* controlling parallism */


/* externally defined function primitives */
extern void INTERF(int DEST, unsigned threadID);
//extern void UPDATE_FORCES(struct link *lint_ptr, int DEST, double XL[], double YL[], double ZL[], double FF[]);
extern void POTENG(tm_double * POTA, tm_double * POTR, tm_double * PTRF, unsigned threadID);
extern void CSHIFT(double XA[], double XB[], double XMA, double XMB, double XL[], double BOXH, double BOXL);
extern void INTRAF(unsigned threadID);
extern void INITIA();
extern void CNSTNT(int N, double *C);
extern void BNDRY(unsigned threadID);
extern void KINETI(tm_double SUM[3], double HMAS, double OMAS, unsigned threadID);
extern void SYSCNS();
extern double MDMAIN(int NFSV, int NFRST, int NSTEP, int NRST, int NPRINT, int NSAVE, int LKT, int NORD1, unsigned threadID);
extern void PREDIC(double C[], int NOR1, unsigned threadID);
extern void CORREC( int NOR1, unsigned threadID);


#endif //__GLOBAL_H__
