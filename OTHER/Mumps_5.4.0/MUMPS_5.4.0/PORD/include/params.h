/*****************************************************************************
/
/ SPACE (SPArse Cholesky Elimination) Library: params.h
/
/ author        J"urgen Schulze, University of Paderborn
/ created       99sep14
/
/ This file contains parameter definitions
/
******************************************************************************/

/* default parameters */
#define MAX_BAD_FLIPS         100   /* interrupt/stop FM */
#define COMPRESS_FRACTION     0.75  /* node reduction in compressed graph */
#define MIN_NODES             100   /* stop recursive separator construction */
#define DEFAULT_SEPS          31    /* default number of separators */
#define MAX_SEPS              255   /* max. number of separators */
#define MIN_DOMAINS           100   /* min. number of domains in a decomp. */
#define MAX_COARSENING_STEPS  10    /* max. number of generated dom. decomp. */

