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
                                                                        
// $Revision: 1.4 $
// $Date: 2007-10-17 22:11:35 $
// $Source: /usr/local/cvs/OpenSees/SRC/domain/load/Beam2dTempLoad.h,v $
                                                                        
#ifndef Beam2dTempLoad_h
#define Beam2dTempLoad_h

// Written: Scott R. Hamilton	15 July 2002 

// Purpose: This file contains the class definition for Beam2dTempLoad, 
// and is modeleed after Beam2dPointLoad by fmk.

#include <ElementalLoad.h>

class Beam2dTempLoad : public ElementalLoad
{
  public:
  // Constructors based on 4, 2, 1 or 0 temperature changes given
  Beam2dTempLoad(int tag, 
		 double Ttop1, double Tbot1, double Ttop2, 
		 double Tbot2, 
		 int theElementTag);

  Beam2dTempLoad(int tag, 
		 double Tuniform, 
		 int theElementTag);

  Beam2dTempLoad(int tag, 
		 double Ttop, double Tbot, 
		 int theElementTag);

  Beam2dTempLoad(int tag, int theElementTag);

  Beam2dTempLoad();    

  ~Beam2dTempLoad();
  
  const Vector &getData(int &type, double loadFactor);
  
  int sendSelf(int commitTag, Channel &theChannel);  
  int recvSelf(int commitTag, Channel &theChannel,  
	       FEM_ObjectBroker &theBroker);
  void Print(OPS_Stream &s, int flag =0);       
  
 protected:
  
 private:
  double Ttop1;	      // Temp change at top node 1 end of member	
  double Tbot1;       // Temp change at bottom node 1 end of member
  double Ttop2;       // Temp change at top node 2 end of member
  double Tbot2;	      // Temp change at bottom node 2 end of member	
  static Vector data; // data for temp loads
};

#endif

