/* 
 * This is a program to try the symbolic factorization
 */

#include <stdio.h>
#include "FeStructs.h"
#include "globalVars.h"


main()
{

  int neq, more;
  int LSPARSE;
  int *xadj, *adjncy;
  int i; 

  printf("Please input LSPARSE\n");
  scanf("%d", &LSPARSE);
 
  neq = 6;
  more = 30;
  xadj = (int *)calloc(neq+1, sizeof(int));
  adjncy = (int *)calloc(more, sizeof(int));

  xadj[0] = 0;
  xadj[1] = 5; 
  xadj[2] = 10;
  xadj[3] = 15;
  xadj[4] = 20;
  xadj[5] = 25;
  xadj[6] = 30;
 
  adjncy[0] = 1;
  adjncy[1] = 2;
  adjncy[2] = 3;
  adjncy[3] = 4;
  adjncy[4] = 5;
  adjncy[5] = 0;
  adjncy[6] = 2;
  adjncy[7] = 3;
  adjncy[8] = 4;
  adjncy[9] = 5;
  adjncy[10] = 0;
  adjncy[11] = 1;
  adjncy[12] = 3;
  adjncy[13] = 4;
  adjncy[14] = 5;
  adjncy[15] = 0;
  adjncy[16] = 1;
  adjncy[17] = 2;
  adjncy[18] = 4;
  adjncy[19] = 5;
  adjncy[20] = 0;
  adjncy[21] = 1;
  adjncy[22] = 2;
  adjncy[23] = 3;
  adjncy[24] = 5;
  adjncy[25] = 0;
  adjncy[26] = 1;
  adjncy[27] = 2;
  adjncy[28] = 3;
  adjncy[29] = 4;

 
 printf("Begin symbolic factorization.\n");
 nblks = symFactorization(xadj, adjncy, neq, LSPARSE); 

printf("\nnblks is: %d ", nblks); 
printf("\n\nThe matrix invp is:\n");
  for (i=0; i<neq+1; i++)
  {
    printf("%d  ", invp[i]);
  }
}

