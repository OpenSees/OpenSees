/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */
                                                                        
// $Revision: 1.2 $
// $Date: 2003-10-15 00:38:07 $
// $Source: /usr/local/cvs/OpenSees/EXAMPLES/ShadowTruss/ShadowActorTruss.h,v $

// Written: fmk
// Revision: A
//
// Description: This file contains the integer codes used in ShadowTruss
// and the ActorTruss classes.
//
// What: "@(#) ShadowActorTruss.h, revA"
static const int ShadowActorTruss_setMaterial        = 1;
static const int ShadowActorTruss_setDomain          = 2;
static const int ShadowActorTruss_commitState        = 3;
static const int ShadowActorTruss_revertToLastCommit = 4;
static const int ShadowActorTruss_revertToStart      = 5;
static const int ShadowActorTruss_update             = 6;
static const int ShadowActorTruss_getTangentStiff    = 7;
static const int ShadowActorTruss_getInitialStiff    = 8;
static const int ShadowActorTruss_getResistingForce  = 9;
static const int ShadowActorTruss_DIE                = 10;
