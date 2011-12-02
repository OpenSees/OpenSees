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
// $Date: 2000-09-15 08:23:25 $
// $Source: /usr/local/cvs/OpenSees/SRC/remote/remote.c,v $
                                                                        
                                                                        
/* File: ~/bin/remote.c
**
** Written: fmk 10/96
** Rev: 
**
** Purpose: To start a process running on a remote machine
**
** usage: remote machine program host_Inet_Address
**
** machine and program being character strings.
*/

#include <stdio.h>
#include <string.h>

char  remotecmd[400];	

main(int argv, char **argc)
{

    int chilpid,i,j;

    /* now get remote program running */
    
    strcpy(remotecmd,"ssh  ");

    for (i =1; i <argv; i++) {
	strcat(remotecmd,argc[i]);
	strcat(remotecmd," ");	
    }
    fprintf(stderr,"%s\n",remotecmd);
    /* This is an expensive hack to make this work */

    if ( (chilpid = fork()) < 0) {  /* FORK !! - hopefully done */
	fprintf(stderr,"REMOTE: could not fork");
	exit(-1);
    }

    if (chilpid == 0) {
       system(remotecmd);
       exit(0);
   }

    exit(0);
}
