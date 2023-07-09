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
  
                             
// $Revision: 1.3 $
// $Date: 2008-04-14 22:38:26 $
// $Source: /usr/local/cvs/OpenSees/SRC/damage/DamageModel.h,v $
                                                                        
#ifndef DamageModel_h
#define DamageModel_h         
                                                               
// Written: Arash Altoontash, Gregory Deierlein
// Created: 08/02
// Revision: AA
//
// Description: This file contains the class definition for 
// damage model Damage modelis a base class and 
// thus no objects of it's type can be instantiated. It has pure virtual 
// functions which must be implemented in it's derived classes. 

#include <DomainComponent.h>
#include <MovableObject.h>
#include <TaggedObject.h>
#include <Vector.h>

class Response;
class DamageResponse;

enum DamageType {
	NotSpecified,
	Force,
	Deformation,
	PlasticDefo,
	TotalEnergy,
	PlasticEnergy,
};

class DamageModel :  public TaggedObject, public MovableObject
{
 public:
    DamageModel(int tag, int classTag);    
    virtual ~DamageModel();

    virtual int setTrial(const Vector &trialVector) = 0;
    virtual double getDamage (void) = 0;
    virtual double getPosDamage (void) = 0;
    virtual double getNegDamage (void) = 0;
    
    virtual int commitState (void) = 0;
    virtual int revertToLastCommit (void) = 0;    
    virtual int revertToStart (void) = 0;        
    
    virtual DamageModel *getCopy (void) = 0;
    
    virtual Response *setResponse(const char **argv, int argc, OPS_Stream &theOutputStream);
    virtual int getResponse(int responseID, Information &info);
    
    virtual int sendSelf(int commitTag, Channel &theChannel) = 0;  
    virtual int recvSelf(int commitTag, Channel &theChannel, 
			 FEM_ObjectBroker &theBroker) = 0;
    virtual void Print(OPS_Stream &s, int flag =0) =0;
    
  protected:
    
  private:

};

extern bool OPS_addDamageModel(DamageModel *newComponent);
extern DamageModel *OPS_getDamageModel(int tag);
extern void OPS_clearAllDamageModel(void);

#endif
