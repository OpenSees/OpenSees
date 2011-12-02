
#include <GKlib.h>

/*****************************************************************          
**********     SMBFCT ..... SYMBOLIC FACTORIZATION       ********* 
******************************************************************          
*   PURPOSE - THIS ROUTINE PERFORMS SYMBOLIC FACTORIZATION               
*   ON A PERMUTED LINEAR SYSTEM AND IT ALSO SETS UP THE               
*   COMPRESSED DATA STRUCTURE FOR THE SYSTEM.                         
*
*   INPUT PARAMETERS -                                               
*      NEQNS - NUMBER OF EQUATIONS.                                 
*      (XADJ, ADJNCY) - THE ADJACENCY STRUCTURE.                   
*      (PERM, INVP) - THE PERMUTATION VECTOR AND ITS INVERSE.     
*
*   UPDATED PARAMETERS -                                         
*      MAXSUB - SIZE OF THE SUBSCRIPT ARRAY NZSUB.  ON RETURN,  
*             IT CONTAINS THE NUMBER OF SUBSCRIPTS USED        
*
*   OUTPUT PARAMETERS -                                       
*      XLNZ - INDEX INTO THE NONZERO STORAGE VECTOR LNZ.   
*      (XNZSUB, NZSUB) - THE COMPRESSED SUBSCRIPT VECTORS. 
*      MAXLNZ - THE NUMBER OF NONZEROS FOUND.             
*
*******************************************************************/
int smbfct(int neqns, int *xadj, int *adjncy, int *perm, int *invp, int *xlnz, int *maxlnz, 
	int *xnzsub, int *nzsub, int *maxsub)
{
  /* Local variables */
  int node, rchm, mrgk, lmax, i, j, k, m, nabor, nzbeg, nzend;
  int kxsub, jstop, jstrt, mrkflg, inz, knz, flag;
  int *mrglnk, *marker, *rchlnk;

  rchlnk = imalloc(neqns+1, "smbfct: rchlnk");
  marker = ismalloc(neqns+1, 0, "smbfct: marker");
  mrglnk = ismalloc(neqns+1, 0, "smbfct: mgrlnk");

  /* Parameter adjustments */
  --marker;
  --mrglnk;
  --rchlnk;
  --nzsub;
  --xnzsub;
  --xlnz;
  --invp;
  --perm;
  --adjncy;
  --xadj;

  /* Function Body */
  flag = 0;
  nzbeg = 1;
  nzend = 0;
  xlnz[1] = 1;

  /* FOR EACH COLUMN KNZ COUNTS THE NUMBER OF NONZEROS IN COLUMN K ACCUMULATED IN RCHLNK. */
  for (k = 1; k <= neqns; ++k) {
    knz = 0;
    mrgk = mrglnk[k];
    mrkflg = 0;
    marker[k] = k;
    if (mrgk != 0) 
      marker[k] = marker[mrgk];
    xnzsub[k] = nzend;
    node = perm[k];

    if (xadj[node] >= xadj[node+1]) {
      xlnz[k+1] = xlnz[k];
      continue;
    }

    /* USE RCHLNK TO LINK THROUGH THE STRUCTURE OF A(*,K) BELOW DIAGONAL */
    rchlnk[k] = neqns+1;
    for (j=xadj[node]; j<xadj[node+1]; j++) {
      nabor = invp[adjncy[j]];
      if (nabor <= k) 
        continue;
      rchm = k;

      do {
        m = rchm;
        rchm = rchlnk[m];
      } while (rchm <= nabor); 

      knz++;
      rchlnk[m] = nabor;
      rchlnk[nabor] = rchm;
      if (marker[nabor] != marker[k]) 
        mrkflg = 1;
    }

    /* TEST FOR MASS SYMBOLIC ELIMINATION */
    lmax = 0;
    if (mrkflg != 0 || mrgk == 0 || mrglnk[mrgk] != 0) 
      goto L350;
    xnzsub[k] = xnzsub[mrgk] + 1;
    knz = xlnz[mrgk + 1] - (xlnz[mrgk] + 1);
    goto L1400;


    /* LINK THROUGH EACH COLUMN I THAT AFFECTS L(*,K) */
L350:
    i = k;
    while ((i = mrglnk[i]) != 0) {
      inz = xlnz[i+1] - (xlnz[i]+1);
      jstrt = xnzsub[i] + 1;
      jstop = xnzsub[i] + inz;

      if (inz > lmax) { 
        lmax = inz;
        xnzsub[k] = jstrt;
      }

      /* MERGE STRUCTURE OF L(*,I) IN NZSUB INTO RCHLNK. */ 
      rchm = k;
      for (j = jstrt; j <= jstop; ++j) {
        nabor = nzsub[j];
        do {
          m = rchm;
          rchm = rchlnk[m];
        } while (rchm < nabor);

        if (rchm != nabor) {
          knz++;
          rchlnk[m] = nabor;
          rchlnk[nabor] = rchm;
          rchm = nabor;
        }
      }
    }

    /* CHECK IF SUBSCRIPTS DUPLICATE THOSE OF ANOTHER COLUMN */
    if (knz == lmax) 
      goto L1400;

    /* OR IF TAIL OF K-1ST COLUMN MATCHES HEAD OF KTH */
    if (nzbeg > nzend) 
      goto L1200;

    i = rchlnk[k];
    for (jstrt = nzbeg; jstrt <= nzend; ++jstrt) {
      if (nzsub[jstrt] < i) 
        continue;

      if (nzsub[jstrt] == i) 
        goto L1000;
      else 
        goto L1200;
    }
    goto L1200;

L1000:
    xnzsub[k] = jstrt;
    for (j = jstrt; j <= nzend; ++j) {
      if (nzsub[j] != i) 
        goto L1200;
      
      i = rchlnk[i];
      if (i > neqns) 
        goto L1400;
    }
    nzend = jstrt - 1;

    /* COPY THE STRUCTURE OF L(*,K) FROM RCHLNK TO THE DATA STRUCTURE (XNZSUB, NZSUB) */
L1200:
    nzbeg = nzend + 1;
    nzend += knz;

    if (nzend > *maxsub) {
      flag = 1; /* Out of memory */
      break;
    }

    i = k;
    for (j=nzbeg; j<=nzend; ++j) {
      i = rchlnk[i];
      nzsub[j] = i;
      marker[i] = k;
    }
    xnzsub[k] = nzbeg;
    marker[k] = k;

    /*
     * UPDATE THE VECTOR MRGLNK.  NOTE COLUMN L(*,K) JUST FOUND   
     * IS REQUIRED TO DETERMINE COLUMN L(*,J), WHERE              
     * L(J,K) IS THE FIRST NONZERO IN L(*,K) BELOW DIAGONAL.      
     */
L1400:
    if (knz > 1) { 
      kxsub = xnzsub[k];
      i = nzsub[kxsub];
      mrglnk[k] = mrglnk[i];
      mrglnk[i] = k;
    }

    xlnz[k + 1] = xlnz[k] + knz;
  }

  if (flag == 0) {
    *maxlnz = xlnz[neqns] - 1;
    *maxsub = xnzsub[neqns];
    xnzsub[neqns + 1] = xnzsub[neqns];
  }

  marker++;
  mrglnk++;
  rchlnk++;
  nzsub++;
  xnzsub++;
  xlnz++;
  invp++;
  perm++;
  adjncy++;
  xadj++;
  GKfree(rchlnk, mrglnk, marker, -1);

  return flag;
  
} 

