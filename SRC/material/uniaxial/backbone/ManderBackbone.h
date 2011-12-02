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

// $Date: 2009-01-08 22:00:17 $                                                                  

// $Source: /usr/local/cvs/OpenSees/SRC/material/uniaxial/backbone/ManderBackbone.h,v $                                                                



// Written: MHS

// Created: Mar 2001

//

// Description: This file contains the implementation of 

// ManderBackbone, which the concrete backbone function given

// by Mander, Priestly, and Park (1988)



#ifndef ManderBackbone_h

#define ManderBackbone_h



#include <HystereticBackbone.h>

#include <Vector.h>



class ManderBackbone : public HystereticBackbone

{

 public:

  ManderBackbone(int tag, double fc, double epsc, double Ec);

  ManderBackbone();

  ~ManderBackbone();

  

  double getStress(double strain);

  double getTangent(double strain);

  double getEnergy(double strain);

  

  double getYieldStrain(void);

  

  HystereticBackbone *getCopy(void);

  

  void Print(OPS_Stream &s, int flag = 0);

  

  int setVariable(char *argv);

  int getVariable(int varID, double &theValue);

  

  int sendSelf(int commitTag, Channel &theChannel);  

  int recvSelf(int commitTag, Channel &theChannel, 

	       FEM_ObjectBroker &theBroker);    

  

 protected:

  

 private:

  double fpc;

  double epsc;

  double Ec;

};



#endif

