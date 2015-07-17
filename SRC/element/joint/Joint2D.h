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

// $Revision: 1.11 $
// $Date: 2010-04-23 22:53:56 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/joint/Joint2D.h,v $

// Written: Arash & GGD
// Created: 03/02
// Revision: Arash

// Joint2D.h: interface for the Joint2d class.
//
//////////////////////////////////////////////////////////////////////

#ifndef Joint2D_h
#define Joint2D_h

#include <bool.h>
#include <Matrix.h>
#include <Vector.h>
#include <Element.h>
#include <ID.h>
#include <Domain.h>

class Node;
class UniaxialMaterial;
class Response;
class DamageModel;

class Joint2D : public Element  
{
public:
  Joint2D();

  Joint2D(int tag, int nd1, int nd2, int nd3, int nd4, int IntNodeTag,
	UniaxialMaterial &spring1,
	UniaxialMaterial &spring2,
	UniaxialMaterial &spring3,
	UniaxialMaterial &spring4,
	UniaxialMaterial &springC,
	Domain *theDomain,
	int LrgDisp);
  
    Joint2D(int tag, int nd1, int nd2, int nd3, int nd4, int IntNodeTag,
	UniaxialMaterial &spring1,
	UniaxialMaterial &spring2,
	UniaxialMaterial &spring3,
	UniaxialMaterial &spring4,
	UniaxialMaterial &springC,
	Domain *theDomain,
	int LrgDisp,
	DamageModel &dmg1,
	DamageModel &dmg2,
	DamageModel &dmg3,
	DamageModel &dmg4,
	DamageModel &dmgC);

  ~Joint2D();

  const char *getClassType(void) const {return "Joint2D";};
  
  // methods dealing with domain
  int	getNumExternalNodes(void) const;
  const	ID &getExternalNodes(void);
  Node **getNodePtrs(void);
  int	getNumDOF(void);
  
  void	setDomain(Domain *theDomain);  
  bool	isSubdomain(void) { return false; } ;
	
  // methods dealing with committed state and update
  int update(void);
  int commitState(void);
  int revertToLastCommit(void);
  int revertToStart(void);
  
  // methods to return the current linearized stiffness,
  // damping and mass matrices
  const	Matrix &getTangentStiff(void);
  const Matrix &getInitialStiff(void);   
  const	Matrix &getDamp(void);
  const	Matrix &getMass(void);
	
  // methods for returning and applying loads
  //virtual Vector &getUVLoadVector(double q1, double q2);
  void	zeroLoad(void); 
  int addLoad(ElementalLoad *theLoad, double loadFactor);
  int addInertiaLoadToUnbalance(const Vector &accel);

  const	Vector &getResistingForce(void);
  const	Vector &getResistingForceIncInertia(void);     

  // method for graphics
  int displaySelf(Renderer &, int mode, float fact, const char **displayModes=0, int numModes=0);
	
  // method for obtaining information specific to an element
  Response* setResponse(const char **argv, int argc, OPS_Stream &s);
  int getResponse(int responseID, Information &eleInformation);
  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
  void Print(OPS_Stream &s, int flag =0);

    // AddingSensitivity:BEGIN //////////////////////////////////////////
    int	  addInertiaLoadSensitivityToUnbalance(const Vector &accel, bool tag);
    int   setParameter(const char **argv, int argc, Parameter &param);
    const Vector & getResistingForceSensitivity(int gradNumber);
    const Matrix & getKiSensitivity(int gradNumber);
    const Matrix & getMassSensitivity(int gradNumber);
    int   commitSensitivity(int gradNumber, int numGrads);
    // AddingSensitivity:END ///////////////////////////////////////////

 protected:
  int 	addMP_Joint(Domain *theDomain, int RnodeID, int CnodeID, int MainDOF, int FixedEnd, int LrgDispFlag );   

 private:
  UniaxialMaterial *theSprings[5]; 
  DamageModel *theDamages[5];
  ID		ExternalNodes;
  ID		InternalConstraints;
  int       fixedEnd[5];
  Node		*theNodes[5];
  Domain	*TheDomain;
  int		numDof, nodeDbTag, dofDbTag;
  static	Matrix K;
  static	Vector V;

  // AddingSensitivity:BEGIN //////////////////////////////////////////
  int parameterID;
  Vector *theLoadSens;
  // AddingSensitivity:END ///////////////////////////////////////////
};

#endif
