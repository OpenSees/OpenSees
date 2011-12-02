/*
 * File:  FeStructs.h
 * ==================
 * Data structures for finite element program.
 * This is changed from the original program by only including
 * the define of OFFDBLK.
 * ----------------------
 * Written by:
 * Jun Peng (junpeng@stanford.edu)
 * Prof. Kincho H. Law
 * Stanford University
 */

#ifndef FeStructs_h
#define FeStructs_h

struct offdblk
{
	int             row    ;
	int             beg    ;
	struct offdblk  *bnext ;
	struct offdblk  *next  ;
	double          *nz    ;
};

typedef  struct offdblk  OFFDBLK;

#endif

/************************************************************/

/*
 * data structure
 *  struct - pointer to pointer to nonzeros in L
 *  struct - pointer to data structure containing nonzeros in L
 *  struct .row - row number of i th row segment
 *  struct .beg - begining column of i th row segment
 *  struct .bnext - pointer to next row segment in this column block
 *  struct .nz   - pointer to coefficients of L in rowsegment i
 *  struct .next - pointer to next row segment with the same row number 
 *		   this row segment if it exists, other wise the row segment
 *		   with the next highest row number
 *
 *  begblk       - pointer to the first structure in each block 
 *  xblk         - an integer array containint the begining column/row
 *		   of each block
 *  penv         - pointers to pointers of the nonzeros in the 
 *		   profile next to the diagonal.
 *  penv[i]      - pointer to the nonzeros of L on the profile
 *		   structure next to teh diagonal
 *  diag[i]      - the i the diagonal of L
 */



