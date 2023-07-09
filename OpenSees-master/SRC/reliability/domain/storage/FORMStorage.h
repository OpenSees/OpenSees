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
                                                        
                                                                        
// Written: Kevin Mackie 
//
// Description: This file contains the FORMStorage class interface

#ifndef FORMStorage_h
#define FORMStorage_h

#include <ReliabilityStorage.h>
#include <Information.h>
#include <Vector.h>


class FORMStorage : public ReliabilityStorage
{
 public:
    FORMStorage();  
    ~FORMStorage();
  
    const char *getClassType(void) const;
    
    int setVariable(const char *variable, Information &);
    int getVariable(const char *variable, Information &);

private:
    Vector *alpha;
    Vector *uStar;
    Vector *xStar;
    Vector *gradientU;
    Vector *gradientX;
    double beta;
    double pf;
    double firstCurvature;
    
};

#endif
