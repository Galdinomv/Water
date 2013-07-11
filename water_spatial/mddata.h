
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

  /* this file contains the declarations of the main data
  structure types used by the program */
  
#define BOTH 2
#include "global.h"
#include "parameters.h"
struct link;
struct list_of_boxes;

typedef struct mol_dummy {
  double VM[3];
  double F[MXOD2][NDIR][NATOM];
} molecule_type;

typedef struct link {
  molecule_type mol;
  struct link * next_mol;
} link_type;

struct tm_link;
typedef struct tm_link* p_tm_link;
typedef tm_type<p_tm_link> tm_p_tm_link; // AICI !!!!!

typedef struct tm_link {
  molecule_type mol;
  tm_p_tm_link next_mol;
} tm_link_type;

typedef struct box_dummy {
  struct link * list;
  //int (boxlock);
} box_type;


typedef struct tm_box_dummy {
  tm_p_tm_link list;
  //int (boxlock);
} tm_box_type;
//extern box_type ***BOX;

typedef struct array_dummy {
  int box[NDIR][BOTH];
} first_last_array;

//extern first_last_array **start_end;

typedef struct list_of_boxes {
  int coord[3];
  struct list_of_boxes * next_box;    
} box_list;

//extern box_list **my_boxes;

extern double  TLC[100], FPOT, FKIN;
extern int IX[3*MXOD2+1], IRST,NVAR,NXYZ,NXV,IXF,IYF,IZF,IMY,IMZ;

extern int NumProcs;
extern int NumBoxes;

struct GlobalSharedMemory {
  first_last_array start_end[MAX_THREADS][MAX_THREADS];
  box_list *my_boxes[MAX_THREADS];
  tm_box_type BOX[MAX_BOX_PER_SIDE][MAX_BOX_PER_SIDE][MAX_BOX_PER_SIDE];
  tm_link_type GLOBAL_NS[MAX_NS_CUBED];
  int global_ns_count;
  box_list GLOBAL_BOXES[MAX_BOX_PER_SIDE*MAX_THREADS*3];
  int global_boxes_count;
};

extern struct GlobalSharedMemory *gsm;
extern double tempXTT;

#define start_end gsm->start_end
#define my_boxes gsm->my_boxes
#define BOX gsm->BOX
