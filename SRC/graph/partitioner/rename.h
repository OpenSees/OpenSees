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
                                                                        
// $Revision: 1.1.1.1 $
// $Date: 2000-09-15 08:23:21 $
// $Source: /usr/local/cvs/OpenSees/SRC/graph/partitioner/rename.h,v $
                                                                        
                                                                        
/*
 * rename.h
 *
 * This file contains renaming functions for the metis library
 *
 * Started 8/22/95
 * George
 *
 * $Id: rename.h,v 1.1.1.1 2000-09-15 08:23:21 fmk Exp $
 *
 */

#ifndef _RENAME_H_
#define _RENAME_H_

/* entrypoint.h */
#define ConvertGraph		__ConvertGraph
#define CleanUpRootGraph	__CleanUpRootGraph


/* balpart.c */
#define FastInitBalance		__FastInitBalance
#define FastBalance		__FastBalance
#define FastBalance2		__FastBalance2

/* bucketlist.c */
#define initbucket		__initbucket
#define resetbucket		__resetbucket
#define freebucket		__freebucket
#define Add2Part		__Add2Part
#define UpdatePart		__UpdatePart
#define GetMaxGainVtx		__GetMaxGainVtx
#define SeeMaxGainVtx		__SeeMaxGainVtx
#define PrintPlusPart		__PrintPlusPart
#define PrintPartGains		__PrintPartGains

/* coarsen.c */
#define Coarsen			__Coarsen
#define KwayCoarsen		__KwayCoarsen
#define SelectMatching		__SelectMatching
#define KwaySelectMatching	__KwaySelectMatching

/* fm.c */
#define FMR_Refine		__FMR_Refine
#define BFMR_Refine		__BFMR_Refine
#define BFMR_Refine_Weighted	__BFMR_Refine_Weighted
#define BFMR_Refine_EqWgt	__BFMR_Refine_EqWgt
#define Greedy_Refine		__Greedy_Refine
#define printwhere		__printwhere
#define CheckBndSize		__CheckBndSize


/* htable.c */
#define CreateHTable		__CreateHTable
#define SelHTSize		__SelHTSize
#define AddHTable		__AddHTable
#define DelHTable		__DelHTable
#define IncreaseHTable		__IncreaseHTable


/* initpart.c */
#define InitPartition		__InitPartition
#define GGPPartition		__GGPPartition
#define GGGPPartition		__GGGPPartition
#define EigPartition		__EigPartition
#define inccompeinz		__inccompeinz
#define ComputeCut		__ComputeCut


/* kwaygreedy.c */
#define KWay_RefineGreedy	__KWay_RefineGreedy
#define KWayUpdateDegrees	__KWayUpdateDegrees
#define KWayUpdateVtxDegrees	__KWayUpdateVtxDegrees


/* kwaypart.c */
#define KWayPart		__KWayPart
#define KWayRefine		__KWayRefine
#define KWayComputePartitionParams	__KWayComputePartitionParams
#define KWayProjectPartition		__KWayProjectPartition
#define KWayCheckDegrees		__KWayCheckDegrees


/* lanczos.c */
#define lanczos			__lanczos
#define givens			__givens
#define pivot			__pivot
#define matvec			__matvec


/* list.c */
#define addthisnode		__addthisnode
#define delthisnode		__delthisnode


/* match.h */
#define RM_Match		__RM_Match
#define RM_Match_W		__RM_Match_W
#define HEM_Match		__HEM_Match
#define HEM_Match_W		__HEM_Match_W
#define LEM_Match		__LEM_Match
#define LEM_Match_W		__LEM_Match_W
#define HCM_Match		__HCM_Match
#define HCM_Match_W		__HCM_Match_W
#define MHEM_Match		__MHEM_Match
#define MHEM_Match_W		__MHEM_Match_W
#define SRM_Match		__SRM_Match
#define SHEM_Match		__SHEM_Match
#define SMHEM_Match		__SMHEM_Match
#define CreateCoarseGraph	__CreateCoarseGraph
#define mergevertices		__mergevertices


/* memory.c */
#define AllocatePools		__AllocatePools
#define FreePools		__FreePools
#define CreateGraph		__CreateGraph
#define FreeRootGraph		__FreeRootGraph
#define FreeGraph		__FreeGraph
#define InitGraph		__InitGraph
#define ResetPools		__ResetPools
#define GetEdgePool		__GetEdgePool
#define SetEdgePool		__SetEdgePool
#define FreeEdgePool		__FreeEdgePool
#define EdgePoolSizeLeft	__EdgePoolSizeLeft
#define GetnExtDegrees		__GetnExtDegrees
#define ResetExtDegrees		__ResetExtDegrees
#define icoremalloc		__icoremalloc
#define icorefree		__icorefree

/* mincover.c */
#define MinCover		__MinCover
#define MinCover_Augment	__MinCover_Augment
#define MinCover_Decompose	__MinCover_Decompose
#define MinCover_ColDFS		__MinCover_ColDFS
#define MinCover_RowDFS		__MinCover_RowDFS


/* mlevelorder.c */
#define MultiLevelOrder		__MultiLevelOrder
#define MLND			__MLND
#define SplitGraphOrder		__SplitGraphOrder
#define SimpleOrder		__SimpleOrder
#define MDOrder			__MDOrder


/* mlevelpart.c */
#define MultiLevelPart		__MultiLevelPart
#define RMLB			__RMLB
#define SplitGraphPart		__SplitGraphPart
#define SplitGraphPart1_2	__SplitGraphPart1_2


/* mmd.c */
#define genmmd			__genmmd
#define mmdelm			__mmdelm
#define mmdint			__mmdint
#define mmdnum			__mmdnum
#define mmdupd			__mmdupd


/* refine.c */
#define Refine			__Refine
#define ComputePartitionParams	__ComputePartitionParams
#define ProjectPartition	__ProjectPartition


/* separator.c */
#define FindMinCovNodeSeparator	__FindMinCovNodeSeparator


/* smbfactor.c */
#define ComputeFillIn		__ComputeFillIn
#define smbfactor		__smbfactor
#define ComputeElTree		__ComputeElTree


/* smbfct.c */
#define smbfct			__smbfct


/* stat.c */
#define PrintGraphMMM		__PrintGraphMMM
#define PrintPartResults	__PrintPartResults
#define PrintOrderResults	__PrintOrderResults
#define CalcNodeOpc		__CalcNodeOpc
#define CalcParOpc		__CalcParOpc
#define ComputePartBalance	__ComputePartBalance


/* util.c */
#define InitRandom		__InitRandom
#define PermuteGraphRandom	__PermuteGraphRandom
#define RandomPermute		__RandomPermute
#define RelDiff			__RelDiff
#define CheckDegrees		__CheckDegrees
#define SortKeyValueNodesDec	__SortKeyValueNodesDec
#define DecKeyValueCmp		__DecKeyValueCmp


/* Renaming for GKlib.c */
#define iasum		__iasum
#define icopy		__icopy
#define iset		__iset
#define iamax		__iamax
#define iamin		__iamin
#define sasum		__sasum
#define dasum 		__dasum
#define daxpy		__daxpy
#define dcopy		__dcopy
#define ddot		__ddot
#define dnrm2		__dnrm2
#define dscal		__dscal
#define dswap		__dswap
#define GKfopen		__GKfopen
#define GKfclose	__GKfclose
#define imalloc		__imalloc
#define ismalloc	__ismalloc
#define GKmalloc	__GKmalloc
#define GKfree		__GKfree
#define iincsort	__iincsort
#define idecsort	__idecsort
#define errexit		__errexit

#ifndef METISLIB
#define cleartimer	__cleartimer
#define starttimer	__startimer
#define stoptimer	__stoptimer
#define printtimer	__printitmer
#define gettimer	__gettimer
#define seconds		__seconds
#else
#define cleartimer(a)		;
#define starttimer(a)		;
#define stoptimer(a)		;
#define printtimer(a, b)	;
#define gettimer(a)		;
#define seconds(a)		;
#endif


#endif
