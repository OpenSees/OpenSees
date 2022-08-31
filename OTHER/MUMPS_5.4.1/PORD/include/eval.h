/*****************************************************************************
/
/ PORD Ordering Library: eval.h
/
/ author        J"urgen Schulze, University of Paderborn
/ created       99mar30
/
/ This file contains the definition of various separator evaluation functions
/
******************************************************************************/

#define F          eval1  /* default separator evaluation function */


  /* --------------------------------------------------------------------- */
  /* SEPARATOR EVALUATION FUNCTION 1                                       */
  /* Size of domains W and B is allowed to differ TOLERANCE * 100 percent. */
  /* Within this tolerance the difference is not penalized and only the    */
  /* size of the separator is returned. Additionally, the mantissa of the  */
  /* returned value is set to (max-min)/max.                               */
  /* --------------------------------------------------------------------- */

#define TOL1       0.50   /* tolerated imbalance induced by bisector */
#define PEN1       100    /* penalty in case of higher imbalance */

#define eval1(S, B, W) \
         S + PEN1 * max(0, max(W,B) * (1-TOL1) - min(W,B)) \
           + (FLOAT)(max(W,B)-min(W,B)) / (FLOAT)max(W,B)

  /* --------------------------------------------------------------------- */
  /* SEPARATOR EVALUATION FUNCTION 2                                       */
  /* Ashcraft and Liu (Using domain decomposition to find graph bisectors) */
  /* --------------------------------------------------------------------- */

#define alpha      0.1
#define TOL2       0.70
#define PEN2       100

#define eval2(S, B, W) \
         S * (1 + alpha * ((FLOAT)max(W,B)/(FLOAT)max(1,min(W,B)))) \
         + PEN2 * max(0, max(W,B) * (1-TOL2) - min(W,B))

  /* --------------------------------------------------------------------- */
  /* SEPARATOR EVALUATION FUNCTION 3                                       */
  /* Ashcraft and Liu (Generalized nested dissection:some recent progress) */
  /* --------------------------------------------------------------------- */

#define alpha2     0.33

#define eval3(S, B, W) \
         S * S + alpha2 * (max(W,B)-min(W,B)) * (max(W,B)-min(W,B))

