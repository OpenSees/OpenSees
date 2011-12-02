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
// $Date: 2006-08-04 22:22:37 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/joint/Joint3D.h,v $

// Written: Arash Altoontash, Gregory Deierlein
// Created: 04/03
// Revision: Arash

// Joint3D.h: interface for the Joint3D class.
//
//////////////////////////////////////////////////////////////////////

#ifndef Joint3D_h
#define Joint3D_h

#include <bool.h>
#include <Matrix.h>
#include <Vector.h>
#include <Element.h>
#include <ID.h>
#include <Domain.h>

class Node;
class UniaxialMaterial;
class Response;

class Joint3D : public Element  
{
public:
  Joint3D();

  Joint3D(int tag, int nd1, int nd2, int nd3, int nd4, int nd5, int nd6, int IntNodeTag,
	  UniaxialMaterial &springx,
	  UniaxialMaterial &springy,
	  UniaxialMaterial &springz,
	  Domain *theDomain,
	  int LrgDisp);
  
  
  ~Joint3D();

  const char *getClassType(void) const {return "Joint3D";};
  
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
  int	displaySelf(Renderer &theViewer, int displayMode, float fact);  
	
  // method for obtaining information specific to an element
  Response* setResponse(const char **argv, int argc, Information &eleInformation, OPS_Stream &s);
  int getResponse(int responseID, Information &eleInformation);
  int sendSelf(int commitTag, Channel &theChannel);
  int recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker);
  void Print(OPS_Stream &s, int flag =0);


 protected:
  int 	addMP_Joint(Domain *theDomain, int mpNum, int RetNodeID, int ConNodeID,
		    int RotNodeID, int Rdof, int DspNodeID, int Ddof, 
		    int LrgDispFlag );   
  
 private:
  UniaxialMaterial *theSprings[3]; 
  ID		ExternalNodes;
  ID		InternalConstraints;
  Node		*theNodes[7];
  Domain	*TheDomain;
  int		numDof, nodeDbTag, dofDbTag;
  static	Matrix K;
  static	Vector V;
};

#endif
