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
                                                                       
// Created: Pedro Arduino, UW, 11.2011
//
// Description: This file contains the class definition for ManzariDafalias.

#ifndef ManzariDafalias_h
#define ManzariDafalias_h

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <NDMaterial.h>
#include <Matrix.h>
#include <Vector.h>

class ManzariDafalias : public NDMaterial
{
  public:

    // full constructor
    ManzariDafalias(int tag, int classTag, double massDen);
    // null constructor
    ManzariDafalias();
    
    // destructor
    ~ManzariDafalias();
 
    NDMaterial *getCopy(const char *type);

    int commitState(void);
    int revertToLastCommit(void);
    int revertToStart(void);

    NDMaterial *getCopy(void);
    const char *getType(void) const;
    int getOrder(void) const;

    Response *setResponse (const char **argv, int argc, OPS_Stream &output);
    int getResponse (int responseID, Information &matInformation);

    int sendSelf(int commitTag, Channel &theChannel);  
    int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker); 

    void Print(OPS_Stream &s, int flag =0);

	int setParameter(const char **argv, int argc, Parameter &param);
  	int updateParameter(int responseID, Information &eleInformation);

	// send mass density to element in dynamic analysis
	double getRho(void) {return massDen;};

  protected:

  // memeber variables and functions

};
#endif
