#ifndef __WATER_MDDATA_HEADER__
#define __WATER_MDDATA_HEADER__

#include "data.h"
#include "global.h"

//typedef double vm_type[3];

typedef struct mol_dummy {
	tm_double VM[3];
	tm_double F[MXOD2][NDIR][NATOM];
} molecule_type;

typedef struct non_tm_mol_dummy {
	double VM[3];
	double F[MXOD2][NDIR][NATOM];
} non_tm_molecule_type;


extern molecule_type *VAR;

// external data

extern double  TLC[100], ELPST,TKIN[MAX_PROCS],TVIR[MAX_PROCS],TTMV[MAX_PROCS],FPOT,FKIN;
extern int IX[3*MXOD2+1],IRST,NVAR,NXYZ,NXV,IXF,IYF,IZF,IMY,IMZ;
extern int StartMol[MAX_PROCS+1];
extern int StartPos[MAX_PROCS+1];
extern int MolsPerProc;
extern double tempXTT[MAX_PROCS];
void CSHIFT(tm_double XA[], tm_double XB[],double XMA, double XMB, double XL[], double BOXH, double BOXL);
void PREDIC(double C[], int NOR1, int threadID);
void INTRAF(tm_double * VIR, int threadID);
void CORREC(double PCC[], int NOR1, int threadID);
void BNDRY(int threadID);
void KINETI(int NMOL, tm_double SUM[], double HMAS, double OMAS, int threadID);
void POTENG(tm_double * POTA, tm_double * POTR, tm_double * PTRF, int threadID);
void INTERF(int, tm_double *, int);
void SYSCNS(void);
void INITIA(FILE * nfmc);
double  MDMAIN(int NFSV, int NFRST, int NSTEP, int NRST, int NPRINT, int NSAVE, int LKT, int NORD1, int threadID);
void CNSTNT(int N, double *C);

#endif //__WATER_MDDATA_HEADER__

