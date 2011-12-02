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
// $Date: 2000-09-15 08:23:17 $
// $Source: /usr/local/cvs/OpenSees/SRC/coordTransformation/LinearCrdTransf3d.cpp,v $
                                                                        
                                                                        
//
// Written: Remo Magalhaes de Souza (rmsouza@ce.berkeley.edu)
// Created: 04/2000
// Revision: A
// 
// Purpose: This file contains the implementation for the 
// LinearCrdTransf3d class. LinearCrdTransf3d is a linear
// transformation for a spatial frame between the global 
// and basic coordinate systems


#include <Vector.h>
#include <Matrix.h>
#include <Node.h>
#include <Channel.h>

#include <LinearCrdTransf3d.h>

// initialize static variables
Vector LinearCrdTransf3d::ul(12); 
Matrix LinearCrdTransf3d::Tlg(12,12); 
Matrix LinearCrdTransf3d::Tbl(6,12); 

 
// constructor:
LinearCrdTransf3d::LinearCrdTransf3d(int tag, const Vector &vecInLocXZPlane,
                                     const Vector &rigJntOffsetI, const Vector &rigJntOffsetJ,
				     int PDeltaFlag):
  CrdTransf3d(tag, CRDTR_TAG_LinearCrdTransf3d),
  nodeIPtr(0), nodeJPtr(0),
  vAxis(3), nodeIOffset(3), nodeJOffset(3), pDeltaFlag(PDeltaFlag), xAxis(3),
  Rlj(3,3), L(0)
{
   // check vector that defines local xz plane
   if (&vecInLocXZPlane == 0 || vecInLocXZPlane.Size() != 3 )
   {
      cerr << "LinearCrdTransf3d::LinearCrdTransf3d:  Vector that defines local xz plane is invalid\n";
      cerr << "Size must be 3\n. Using (0,0,1)";      
      vAxis(0) = 0;       vAxis(1) = 0;      vAxis(2) = 1;
   }
   else
     vAxis = vecInLocXZPlane;
       
   // check rigid joint offset for node I
   if (&rigJntOffsetI == 0 || rigJntOffsetI.Size() != 3 )
   {
      cerr << "LinearCrdTransf3d::LinearCrdTransf3d:  Invalid rigid joint offset vector for node I\n";
      cerr << "Size must be 3\n";      
      nodeIOffset.Zero();      
   }
   else
     nodeIOffset = rigJntOffsetI;
   
   // check rigid joint offset for node J
   if (&rigJntOffsetJ == 0 || rigJntOffsetJ.Size() != 3 )
   {
      cerr << "LinearCrdTransf3d::LinearCrdTransf3d:  Invalid rigid joint offset vector for node J\n";
      cerr << "Size must be 3\n";      
      nodeJOffset.Zero(); 
   }
   else
     nodeJOffset = rigJntOffsetJ;
}



// constructor:
// invoked by a FEM_ObjectBroker, recvSelf() needs to be invoked on this object.
LinearCrdTransf3d::LinearCrdTransf3d():
  CrdTransf3d(0, CRDTR_TAG_LinearCrdTransf3d),
  nodeIPtr(0), nodeJPtr(0),
  vAxis(3), nodeIOffset(3), nodeJOffset(3), pDeltaFlag(0), xAxis(3),
  Rlj(3,3), L(0)
{
}



// destructor:
LinearCrdTransf3d::~LinearCrdTransf3d() 
{
}

int
LinearCrdTransf3d::commitState(void)
{
   // linear transformation - nothing to commit
   return 0;
}


int
LinearCrdTransf3d::revertToLastCommit(void)
{
   this->update();
   return 0;
}


int
LinearCrdTransf3d::revertToStart(void)
{
   this->update();
   return 0;
}


int LinearCrdTransf3d::initialize(Node *nodeIPointer, Node *nodeJPointer)
{       
   int error;

   nodeIPtr = nodeIPointer;
   nodeJPtr = nodeJPointer;

   if ((!nodeIPtr) || (!nodeJPtr))
   {
      cerr << "\nLinearCrdTransf3d::initialize";
      cerr << "\ninvalid pointers to the element nodes\n";
      return -1;
   }
       
   // get element length and orientation
   if ((error = this->computeElemtLengthAndOrient()))
      return error;
      
   // get 3by3 rotation matrix (the columns of which are the element local axes)
   if ((error = this->getLocalAxes()))
      return error;
   
   // compute transformation matrices
   // this->update();     // ************* check

   return 0;
}


int
LinearCrdTransf3d::update(void)
{       
   // get matrix that transforms displacements from global coordinates
   // to the joint offsets position (eliminating rigid joint effect)    
   static Matrix Tjg(12,12);
   this->getTransfMatrixJointGlobal(Tjg);

   // get matrix that transforms displacements at the joint offsets to 
   // local coordinates (at the ends of the flexible part of the element)
   static Matrix Tlj(12,12);
   this->getTransfMatrixLocalJoint(Tlj);     // OPTIMIZE LATER

   //cerr << "Rlj: " << Rlj;
   //cerr << "Tjg: " << Tjg;
   
   // compute composite matrix that transforms displacements from the original
   // joint position to local coordinates    
   Tlg.addMatrixProduct(0.0, Tlj, Tjg, 1.0);     // Tlg = Tlj * Tjg; // OPTIMIZE LATER

   // get matrix that transforms displacements from local 
   // coordinates to basic system without rigid body modes
   this->getTransfMatrixBasicLocal(Tbl);

   //cerr << "Tbl: " << Tbl;
   //cerr << "Tlg: " << Tlg;

   return 0;
}


int 
LinearCrdTransf3d::computeElemtLengthAndOrient()
{
   // element projection

   static Vector dx(3);

   dx = (nodeJPtr->getCrds() + nodeJOffset) - (nodeIPtr->getCrds() + nodeIOffset);  
	    
   // calculate the element length

   L = dx.Norm();

   if (L == 0.0) 
   {
      cerr << "\nLinearCrdTransf3d::computeElemtLengthAndOrien: 0 length\n";
      return -2;  
   }

   // calculate the element local x axis components (direction cossines)
   // wrt to the global coordinates 
   xAxis = dx/L;

   return 0;
}



int 
LinearCrdTransf3d::getLocalAxes (void)
{
   // calculate the cross-product y = v * x   

   static Vector yAxis(3), zAxis(3);
   
   yAxis(0) = vAxis(1)*xAxis(2) - vAxis(2)*xAxis(1);
   yAxis(1) = vAxis(2)*xAxis(0) - vAxis(0)*xAxis(2);
   yAxis(2) = vAxis(0)*xAxis(1) - vAxis(1)*xAxis(0);
   
   double ynorm = yAxis.Norm();
 
   if (ynorm == 0)
   {
      cerr << "\nLinearCrdTransf3d::getElementLengthAndOrientation";
      cerr << "\nvector v that defines plane xz is parallel to x axis\n";
      return -3;
   }

   yAxis /= ynorm;
     
   // calculate the cross-product z = x * y 

   zAxis(0) = xAxis(1)*yAxis(2) - xAxis(2)*yAxis(1);
   zAxis(1) = xAxis(2)*yAxis(0) - xAxis(0)*yAxis(2);
   zAxis(2) = xAxis(0)*yAxis(1) - xAxis(1)*yAxis(0);

   for (int i=0; i < 3; i++)
   {
      Rlj(0,i) = xAxis(i);
      Rlj(1,i) = yAxis(i);
      Rlj(2,i) = zAxis(i);
   }      
       
   return 0;
}


double 
LinearCrdTransf3d::getInitialLength(void)
{
   return L;
}


double 
LinearCrdTransf3d::getDeformedLength(void)
{
   return L;
}


const Vector &
LinearCrdTransf3d::getBasicTrialDisp (void)
{
   // determine global displacements
   const Vector &disp1 = nodeIPtr->getTrialDisp();
   const Vector &disp2 = nodeJPtr->getTrialDisp();

   static Vector ug(12);
   for (int i = 0; i < 6; i++)
   {
      ug(i)   = disp1(i);
      ug(i+6) = disp2(i);
   }

   // transform global end displacements to local coordinates
   // static Vector ul(12);      // total displacements

   ul.addMatrixVector(0.0, Tlg,  ug, 1.0);       //  ul = Tlg *  ug;

   // eliminate rigid body modes from the displacements of the 
   // local system
   static Vector ub(6);
   ub.addMatrixVector(0.0, Tbl,  ul, 1.0);       //  ub = Tbl *  ul

   return ub;
}


const Vector &
LinearCrdTransf3d::getBasicIncrDisp (void)
{
   // determine global displacement increments wrt last converged state
   const Vector &dispIncr1 = nodeIPtr->getIncrDisp();
   const Vector &dispIncr2 = nodeJPtr->getIncrDisp();

   static Vector dug(12);
   for (int i = 0; i < 6; i++)
   {
      dug(i)   = dispIncr1(i);
      dug(i+6) = dispIncr2(i);
   }

   // transform global end displacement increments to local coordinates

   static Vector dul(12);     // displ. increments from last converged state
   
   dul.addMatrixVector(0.0, Tlg, dug, 1.0);       // dul = Tlg * dug;

   // eliminate rigid body modes from the displacement increments of the 
   // local system
   static Vector dub(6);
   dub.addMatrixVector(0.0, Tbl, dul, 1.0);       // dub = Tbl * dul;

   return dub;
}


const Vector &
LinearCrdTransf3d::getBasicIncrDeltaDisp(void)
{
   // determine global displacement increments wrt last iteration
   const Vector &dispIncr1 = nodeIPtr->getIncrDeltaDisp();
   const Vector &dispIncr2 = nodeJPtr->getIncrDeltaDisp();

   static Vector Dug(12);
   for (int i = 0; i < 6; i++)
   {
      Dug(i)   = dispIncr1(i);
      Dug(i+6) = dispIncr2(i);
   }

   // transform global end displacement increments to local coordinates
   static Vector Dul(12);     // displ. increments from last iteration
   
   Dul.addMatrixVector(0.0, Tlg, Dug, 1.0);       // Dul = Tlg * Dug;

   // eliminate rigid body modes from the displacement increments of the 
   // local system
   static Vector Dub(6);
   Dub.addMatrixVector(0.0, Tbl, Dul, 1.0);       // Dub = Tbl * Dul;

   return Dub;
}


void 
LinearCrdTransf3d::getTransfMatrixJointGlobal(Matrix &Tjg)
{
   // setup transformation matrix that considers rigid joint offsets

   Tjg.Zero();

   for (int i = 0; i < 12; i++)
      Tjg(i,i) = 1;

   // cerr << "getTrans rig= " << rigJntOffset << endl;

   Tjg(1, 3) = -nodeIOffset(2);
   Tjg(2, 3) =  nodeIOffset(1);
   Tjg(0, 4) =  nodeIOffset(2);
   Tjg(2, 4) = -nodeIOffset(0);
   Tjg(0, 5) = -nodeIOffset(1);
   Tjg(1, 5) =  nodeIOffset(0);

   Tjg(7, 9) = -nodeJOffset(2);
   Tjg(8, 9) =  nodeJOffset(1);
   Tjg(6,10) =  nodeJOffset(2);
   Tjg(8,10) = -nodeJOffset(0);
   Tjg(6,11) = -nodeJOffset(1);
   Tjg(7,11) =  nodeJOffset(0);
}


void
LinearCrdTransf3d::getTransfMatrixLocalJoint (Matrix &Tlj) 
{
   // setup transformation matrix
   Tlj.Zero();

   // put rotation matrix R on the 4 diagonal submatrices of Tlj
   int i, j, k;
   for (k = 0; k < 12; k += 3)
      for (j = 0; j < 3; j++)
         for (i = 0; i < 3; i++)
            Tlj(k+i,k+j) = Rlj(i,j);
}


void 
LinearCrdTransf3d::getTransfMatrixBasicLocal(Matrix &Tbl)
{
   Tbl.Zero();
   
   double oneOverL = 1/L;
    
   Tbl(0, 0) = -1;            
   Tbl(0, 6) =  1;
   Tbl(1, 1) =  oneOverL;
   Tbl(1, 5) =  1;
   Tbl(1, 7) = -oneOverL;
   Tbl(2, 1) =  oneOverL;
   Tbl(2, 7) = -oneOverL;
   Tbl(2,11) =  1;
   Tbl(3, 2) = -oneOverL;
   Tbl(3, 4) =  1;
   Tbl(3, 8) =  oneOverL;
   Tbl(4, 2) = -oneOverL;
   Tbl(4, 8) =  oneOverL;
   Tbl(4,10) =  1;
   Tbl(5, 3) = -1;  
   Tbl(5, 9) =  1;
}


const Vector &
LinearCrdTransf3d::getGlobalResistingForce(const Vector &pb, const Vector &unifLoad)
{
   // transform resisting forces from the basic system to local coordinates
   static Vector pl(12);
   pl.addMatrixTransposeVector(0.0, Tbl, pb, 1.0);    // pl = Tbl ^ pb;

   // add end forces due to uniforme distributed loads to the system with rigid body modes
   pl(0) -= unifLoad(0)*L;
   pl(1) -= unifLoad(1)*L/2;
   pl(7) -= unifLoad(1)*L/2;
   pl(2) -= unifLoad(2)*L/2;
   pl(8) -= unifLoad(2)*L/2;
         
   // include P-Delta effects
   if (pDeltaFlag)
   {
       double NoverL = pb(0)/L;
       // cerr << "force NoverL: " << NoverL;
       
       pl(1) += (ul(1)-ul(7))*NoverL;
       pl(7) -= (ul(1)-ul(7))*NoverL;
       pl(2) += (ul(2)-ul(8))*NoverL;
       pl(8) -= (ul(2)-ul(8))*NoverL; 
   }
   
   // transform resisting forces  from local to global coordinates
   static Vector pg(12);
   pg.addMatrixTransposeVector(0.0, Tlg, pl, 1.0);    // pg = Tlg ^ pl;

   return pg;
}
  


const Matrix &
LinearCrdTransf3d::getGlobalStiffMatrix (const Matrix &kb, const Vector &pb)
{
   // transform tangent stiffness matrix from the basic system to local coordinates
   static Matrix kl(12,12);

   //cerr << "Tbl: " << Tbl;
   
   //cerr << "Tlg: " << Tlg;
   
   kl.addMatrixTripleProduct(0.0, Tbl, kb, 1.0);      // kl = Tbl ^ kb * Tbl;

   //cerr << "kl: " << kl;
   // include P-Delta effects
   if (pDeltaFlag)
   {
       double NoverL = pb(0)/L;
       // cerr << "stiffness  NoverL: " << NoverL;
           
       kl(1,1) += NoverL;
       kl(2,2) += NoverL;
       kl(7,7) += NoverL; 
       kl(8,8) += NoverL;
       kl(1,7) -= NoverL; 
       kl(7,1) -= NoverL;
       kl(2,8) -= NoverL;
       kl(8,2) -= NoverL; 
   }

   // transform tangent  stiffness matrix from local to global coordinates
   static Matrix kg(12,12);
   kg.addMatrixTripleProduct(0.0, Tlg, kl,  1.0);     // kg =  Tlj ^ kl * Tlj;

   return kg;
}



CrdTransf3d *
LinearCrdTransf3d::getCopy(void)
{
  // create a new instance of LinearCrdTransf3d 

  LinearCrdTransf3d *theCopy = new LinearCrdTransf3d (this->getTag(), vAxis, nodeIOffset, nodeJOffset, pDeltaFlag);
  
  if (!theCopy)
  {
     g3ErrorHandler->fatal("LinearCrdTransf3d::getCopy() - out of memory creating copy");
     return 0;
  }    
    
  theCopy->nodeIPtr = nodeIPtr;
  theCopy->nodeJPtr = nodeJPtr;
  theCopy->xAxis = xAxis;
  theCopy->Rlj = Rlj;

  theCopy->L = L;
  
  return theCopy;
}



int 
LinearCrdTransf3d::sendSelf(int cTag, Channel &theChannel)
{
	int res = 0;

	static Vector data(21);

	data(0)  = this->getTag();
	data(1)  = pDeltaFlag;
	data(2)  = nodeIOffset(0);
	data(3)  = nodeIOffset(1);
	data(4)  = nodeIOffset(2);
	data(5)  = nodeJOffset(0);
	data(6)  = nodeJOffset(1);
	data(7) = nodeJOffset(2);
	data(8) = vAxis(0);
	data(9) = vAxis(1);
	data(10) = vAxis(2);
	data(11) = L;

	// Can't just send vAxis because node pointers are needed
	int loc = 12;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			data(loc++) = Rlj(i,j);

	res += theChannel.sendVector(this->getDbTag(), cTag, data);
	if (res < 0) {
		g3ErrorHandler->warning("%s - failed to send Vector",
			"LinearCrdTransf3d::sendSelf");
		return res;
	}

    return res;
}

    

int 
LinearCrdTransf3d::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	int res = 0;

	static Vector data(21);

	res += theChannel.recvVector(this->getDbTag(), cTag, data);
	if (res < 0) {
		g3ErrorHandler->warning("%s - failed to receive Vector",
			"LinearCrdTransf3d::recvSelf");
		return res;
	}

	this->setTag((int)data(0));
	pDeltaFlag = (int)data(1);
	nodeIOffset(0) = data(2);
	nodeIOffset(1) = data(3);
	nodeIOffset(2) = data(4);
	nodeJOffset(0) = data(5);
	nodeJOffset(1) = data(6);
	nodeJOffset(2) = data(7);
	vAxis(0) = data(8);
	vAxis(1) = data(9);
	vAxis(2) = data(10);
	L = data(11);

	int loc = 12;
	for (int i = 0; i < 3; i++)
		for (int j = 0; j < 3; j++)
			Rlj(i,j) = data(loc++);

	return res;
}
 	
const Vector &
LinearCrdTransf3d::getPointGlobalCoordFromLocal(const Vector &xl)
{
   static Vector xg(3);

   xg = nodeIPtr->getCrds() + nodeIOffset;
   // xg = xg + RljT*xl
   xg.addMatrixTransposeVector(1.0, Rlj, xl, 1.0);
     
   return xg;  
}

    
const Vector &
LinearCrdTransf3d::getPointGlobalDisplFromBasic (double xi, const Vector &uxb)
{
   // determine global displacements
   const Vector &disp1 = nodeIPtr->getTrialDisp();
   const Vector &disp2 = nodeJPtr->getTrialDisp();

   static Vector ug(12);
   for (int i = 0; i < 6; i++)
   {
      ug(i)   = disp1(i);
      ug(i+6) = disp2(i);
   }

   // transform global end displacements to local coordinates
   static Vector ul(12);      // total displacements

   ul.addMatrixVector(0.0, Tlg,  ug, 1.0);       //  ul = Tlg *  ug;
   
   // compute displacements at point xi, in local coordinates
   static Vector uxl(3),  uxg(3);

   uxl(0) = uxb(0) +        ul(0);
   uxl(1) = uxb(1) + (1-xi)*ul(1) + xi*ul(7);
   uxl(2) = uxb(2) + (1-xi)*ul(2) + xi*ul(8);
  
   // rotate displacements to global coordinates
   // uxg = RljT*uxl
   uxg.addMatrixTransposeVector(0.0, Rlj, uxl, 1.0);
     
   return uxg;  
}





void
LinearCrdTransf3d::Print(ostream &s, int flag)
{
   s << "\nCrdTransf: " << this->getTag() << " Type: LinearCrdTransf3d";
   s << "\tvAxis: " << vAxis;
   s << "\tnodeI Offset: " << nodeIOffset;
   s << "\tnodeJ Offset: " << nodeJOffset;
   s << "\tpDelta flag: "  << pDeltaFlag;
}


