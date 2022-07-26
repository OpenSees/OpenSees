/*****************************************************************************
/
/ SPACE (SPArse Cholesky Elimination) Library: const.h
/
/ author        J"urgen Schulze, University of Paderborn
/ created       99sep14
/
/ This file contains constant definitions
/
******************************************************************************/

/* matrix types */
#define GRID		0
#define MESH		1
#define TORUS		2
#define HB              3

/* graph types */
#define UNWEIGHTED      0
#define WEIGHTED        1

/* type of ordering */
#define MINIMUM_PRIORITY       0
#define INCOMPLETE_ND          1
#define MULTISECTION           2
#define TRISTAGE_MULTISECTION  3

/* fill-reducing node selection strategies */
#define AMD        0
#define AMF        1
#define AMMF       2
#define AMIND      3

/* node selection strategies for generating the domain decompositions */
#define QMD    0
#define QMRDV  1
#define QRAND  2

/* default options for SPACE */
#define SPACE_ORDTYPE          MULTISECTION
#define SPACE_NODE_SELECTION1  AMMF
#define SPACE_NODE_SELECTION2  AMMF
#define SPACE_NODE_SELECTION3  QMRDV
#define SPACE_DOMAIN_SIZE      200
#define SPACE_MSGLVL           2
#define SPACE_ETREE_NONZ       256
#define SPACE_ETREE_BAL        5
#define SPACE_MASK_OFFSET      2

/* misc. constants */
#define TRUE		1
#define FALSE		0
#define ERR             -1
#define NOERR           0
#define MAX_LINE_LEN    255
#define MAX_INT         ((1<<30)-1)
#define MAX_FLOAT       1e31
#define EPS             0.001

/* constants used in color array */
/* these constants are also used as an index (do not change) */
#define GRAY       0
#define BLACK      1
#define WHITE      2

/* constants for the Dulmage-Mendelsohn decomposition (dmflags) */
/* these constants are also used as an index (do not change) */
#define SI         0     /* node e X is reachable via exposed node e X */
#define SX         1     /* node e X is reachable via exposed node e Y */
#define SR         2     /* SR = X - (SI u SX) */
#define BI         3     /* node e Y is reachable via exposed node e Y */
#define BX         4     /* node e Y is reachable via exposed node e X */
#define BR         5     /* BR = Y - (BI u BX) */

/* size/indices of option array (do not change) */
#define ORD_OPTION_SLOTS        7

#define OPTION_ORDTYPE          0
#define OPTION_NODE_SELECTION1  1
#define OPTION_NODE_SELECTION2  2
#define OPTION_NODE_SELECTION3  3
#define OPTION_DOMAIN_SIZE      4
#define OPTION_MSGLVL           5
#define OPTION_ETREE_NONZ       6

/* size/indices for timing array in ordering computation */
#define ORD_TIME_SLOTS    12

#define TIME_COMPRESS     0    /*  0. TIME_COMPRESS                */
#define TIME_MS           1    /*  1. TIME_MS                      */
#define TIME_MULTILEVEL   2    /*     1.1 TIME_MULTILEVEL          */
#define TIME_INITDOMDEC   3    /*         1.1.1 TIME_INITDOMDEC    */
#define TIME_COARSEDOMDEC 4    /*         1.1.2 TIME_COARSEDOMDEC  */
#define TIME_INITSEP      5    /*         1.1.3 TIME_INITSEP       */
#define TIME_REFINESEP    6    /*         1.1.4 TIME_REFINESEP     */
#define TIME_SMOOTH       7    /*     1.2 TIME_SMOOTH              */
#define TIME_BOTTOMUP     8    /*  2. TIME_BOTTOMUP                */
#define TIME_UPDADJNCY    9    /*     2.1 TIME_UPDADJNCY           */
#define TIME_FINDINODES  10    /*     2.2 TIME_FINDINODES          */
#define TIME_UPDSCORE    11    /*     2.3 TIME_UPDSCORE            */

/* size/indices for timing array in sequential numerical factorization */
#define NUMFAC_TIME_SLOTS 4
 
#define TIME_INITFRONT    0
#define TIME_EXADD        1
#define TIME_KERNEL       2
#define TIME_INITUPD      3

/* size/indices for timing array in parallel numerical factorization */
#define NUMFACPAR_TIME_SLOTS 9    
 
#define TIME_INITFRONT       0    
#define TIME_EXADD           1    
#define TIME_KERNEL          2
#define TIME_INITUPD         3
#define TIME_EXCHANGE        4
#define TIME_INITFRONTPAR    5
#define TIME_EXADDPAR        6
#define TIME_KERNELPAR       7
#define TIME_INITUPDPAR      8

/* size/indices for timing array in parallel kernel */
#define KERNELPAR_TIME_SLOTS 4

#define TIME_PIVOT           0
#define TIME_PIVOT_WAIT      1
#define TIME_CMOD            2
#define TIME_CMOD_WAIT       3
