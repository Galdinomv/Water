#ifndef __WATER_GLOBAL_HEADER__
#define __WATER_GLOBAL_HEADER__

#define MAXMOLS	1728
#include "include/src/tm_basics.h"

typedef tm_type<int> tm_int;
typedef tm_type<double> tm_double;

struct GlobalMemory {
	//int IntrafVirLock;
    	//int InterfVirLock;
    	//int KinetiSumLock;
	//int PotengSumLock;
	//int MolLock[MAXMOLS];
    	int start;	// start     has just 1 value assigned : "0"
	int InterfBar;	// InterfBar has just 1 value assigned : "1"
	int PotengBar;	// PotengBar has just 1 value assigned : "2"
	tm_double VIR;
	tm_double SUM[3];
	tm_double POTA, POTR, POTRF;
};

extern struct GlobalMemory *gl;


#endif //__WATER_GLOBAL_HEADER__
