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
// $Source: /usr/local/cvs/OpenSees/SRC/coordTransformation/LinearCrdTransf2d.cpp,v $
                                                                        
                                                                        
// File: ~/crdTransf/LinearCrdTransf2d.C
//
// Written: Remo Magalhaes de Souza (rmsouza@ce.berkeley.edu)
// Created: 04/2000
// Revision: A
// 
// Purpose: This file contains the implementation for the 
// LinearCrdTransf2d class. LinearCrdTransf2d is a linear
// transformation for a planar frame between the global 
// and basic coordinate systems


#include <Vector.h>
#include <Matrix.h>
#include <Node.h>
#include <Channel.h>

#include <iomanip.h>

#include <LinearCrdTransf2d.h>

// initialize static variables
Vector LinearCrdTransf2d::ul(6); 
Matrix LinearCrdTransf2d::Tlg(6,6); 
Matrix LinearCrdTransf2d::Tbl(3,6);

 
// constructor:
LinearCrdTransf2d::LinearCrdTransf2d(int tag, const Vector &rigJntOffset1, const Vector &rigJntOffset2, 
				     int PDeltaFlag):
  CrdTransf2d(tag, CRDTR_TAG_LinearCrdTransf2d),
  nodeIPtr(0), nodeJPtr(0),
  nodeIOffset(2), nodeJOffset(2), pDeltaFlag(PDeltaFlag), 
  cosTheta(0), sinTheta(0), L(0)
{
   // check rigid joint offset for node I
   if (&rigJntOffset1 == 0 || rigJntOffset1.Size() != 2 )
   {
      cerr << "LinearCrdTransf2d::LinearCrdTransf2d:  Invalid rigid joint offset vector for node I\n";
      cerr << "Size must be 2\n";      
      nodeIOffset.Zero();      
   }
   else
     nodeIOffset = rigJntOffset1;
   
   // check rigid joint offset for node J
   if (&rigJntOffset2 == 0 || rigJntOffset2.Size() != 2 )
   {
      cerr << "LinearCrdTransf2d::LinearCrdTransf2d:  Invalid rigid joint offset vector for node J\n";
      cerr << "Size must be 2\n";      
      nodeJOffset.Zero(); 
   }
   else
     nodeJOffset = rigJntOffset2;
}



 
// constructor:
// invoked by a FEM_ObjectBroker, recvSelf() needs to be invoked on this object.
LinearCrdTransf2d::LinearCrdTransf2d():
  CrdTransf2d(0, CRDTR_TAG_LinearCrdTransf2d),
  nodeIPtr(0), nodeJPtr(0),
  nodeIOffset(2), nodeJOffset(2), pDeltaFlag(0), 
  cosTheta(0), sinTheta(0), L(0)
{

}



// destructor:
LinearCrdTransf2d::~LinearCrdTransf2d() 
{
}


int
LinearCrdTransf2d::commitState(void)
{
   // linear transformation - nothing to commit
   return 0;
}


int
LinearCrdTransf2d::revertToLastCommit(void)
{
   this->update();
   return 0;
}


int
LinearCrdTransf2d::revertToStart(void)
{
   this->update();
   return 0;
}


int 
LinearCrdTransf2d::initialize(Node *nodeIPointer, Node *nodeJPointer)
{       
   //cerr << setiosflags(ios::scientific);
   //cerr << setiosflags(ios::showpos);
   //cerr << setprecision(16);
  
   int error;

   nodeIPtr = nodeIPointer;
   nodeJPtr = nodeJPointer;

   if ((!nodeIPtr) || (!nodeJPtr))
   {
      cerr << "\nLinearCrdTransf2d::initialize";
      cerr << "\ninvalid pointers to the element nodes\n";
      return -1;
   }
       
   // get element length and orientation
   if ((error = this->computeElemtLengthAndOrient()))
      return error;
      
   return 0;
}


int
LinearCrdTransf2d::update(void)
{       
   // get matrix that transforms displacements from global coordinates
   // to the joint offsets position (eliminating rigid joint effect)    
   static Matrix Tjg(6,6);
   this->getTransfMatrixJointGlobal(Tjg);

   // get matrix that transforms displacements at the joint offsets to 
   // local coordinates (at the ends of the flexible part of the element)
   static Matrix Tlj(6,6);
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
LinearCrdTransf2d::computeElemtLengthAndOrient()
{
   // element projection

   static Vector dx(2);

   dx = (nodeJPtr->getCrds() + nodeJOffset) - (nodeIPtr->getCrds() + nodeIOffset);  
	    
   // calculate the element length

   L = dx.Norm();

   if (L == 0.0) 
   {
      cerr << "\nLinearCrdTransf2d::computeElemtLengthAndOrien: 0 length\n";
      return -2;  
   }

   // calculate the element local x axis components (direction cossines)
   // wrt to the global coordinates 
   cosTheta = dx(0)/L;
   sinTheta = dx(1)/L;
   
   return 0;
}





double 
LinearCrdTransf2d::getInitialLength(void)
{
   return L;
}


double 
LinearCrdTransf2d::getDeformedLength(void)
{
   return L;
}


const Vector &
LinearCrdTransf2d::getBasicTrialDisp (void)
{
   // determine global displacements
   const Vector &disp1 = nodeIPtr->getTrialDisp();
   const Vector &disp2 = nodeJPtr->getTrialDisp();

   static Vector ug(6);
   for (int i = 0; i < 3; i++)
   {
      ug(i)   = disp1(i);
      ug(i+3) = disp2(i);
   }

   // transform global end displacements to local coordinates
   // static Vector ul(6);      // total displacements

   ul.addMatrixVector(0.0, Tlg,  ug, 1.0);       //  ul = Tlg *  ug;

   // eliminate rigid body modes from the displacements of the 
   // local system
   static Vector ub(3);
   ub.addMatrixVector(0.0, Tbl,  ul, 1.0);       //  ub = Tbl *  ul

   return ub;
}


const Vector &
LinearCrdTransf2d::getBasicIncrDisp (void)
{
   // determine global displacement increments wrt last converged state
   const Vector &dispIncr1 = nodeIPtr->getIncrDisp();
   const Vector &dispIncr2 = nodeJPtr->getIncrDisp();

   static Vector dug(6);
   for (int i = 0; i < 3; i++)
   {
      dug(i)   = dispIncr1(i);
      dug(i+3) = dispIncr2(i);
   }

   // transform global end displacement increments to local coordinates

   static Vector dul(6);     // displ. increments from last converged state
   
   dul.addMatrixVector(0.0, Tlg, dug, 1.0);       // dul = Tlg * dug;

   // eliminate rigid body modes from the displacement increments of the 
   // local system
   static Vector dub(3);
   dub.addMatrixVector(0.0, Tbl, dul, 1.0);       // dub = Tbl * dul;

   return dub;
}


const Vector &
LinearCrdTransf2d::getBasicIncrDeltaDisp(void)
{
   // determine global displacement increments wrt last iteration
   const Vector &dispIncr1 = nodeIPtr->getIncrDeltaDisp();
   const Vector &dispIncr2 = nodeJPtr->getIncrDeltaDisp();

   static Vector Dug(6);
   for (int i = 0; i < 3; i++)
   {
      Dug(i)   = dispIncr1(i);
      Dug(i+3) = dispIncr2(i);
   }

   // transform global end displacement increments to local coordinates
   static Vector Dul(6);     // displ. increments from last iteration
   
   Dul.addMatrixVector(0.0, Tlg, Dug, 1.0);       // Dul = Tlg * Dug;

   // eliminate rigid body modes from the displacement increments of the 
   // local system
   static Vector Dub(3);
   Dub.addMatrixVector(0.0, Tbl, Dul, 1.0);       // Dub = Tbl * Dul;

   return Dub;
}


void 
LinearCrdTransf2d::getTransfMatrixJointGlobal(Matrix &Tjg)
{
   // setup transformation matrix that considers rigid joint offsets

   Tjg.Zero();

   for (int i = 0; i < 6; i++)
      Tjg(i,i) = 1;

   // cerr << "getTrans rig= " << rigJntOffset << endl;

   Tjg(0,2) = -nodeIOffset(1);
   Tjg(1,2) =  nodeIOffset(0);

   Tjg(3,5) = -nodeJOffset(1);
   Tjg(4,5) =  nodeJOffset(0); 
}


void
LinearCrdTransf2d::getTransfMatrixLocalJoint (Matrix &Tlj) 
{
   // setup transformation matrix
   Tlj.Zero();

   Tlj(0,0) = Tlj(3,3) =  cosTheta;           
   Tlj(0,1) = Tlj(3,4) =  sinTheta;
   Tlj(1,0) = Tlj(4,3) = -sinTheta;
   Tlj(1,1) = Tlj(4,4) =  cosTheta;
   Tlj(2,2) = Tlj(5,5) =  1.0;

}


void 
LinearCrdTransf2d::getTransfMatrixBasicLocal(Matrix &Tbl)
{
   Tbl.Zero();
   
   double oneOverL = 1/L;

   Tbl(0,0) = -1;            
   Tbl(0,3) =  1;
   Tbl(1,1) =  oneOverL;
   Tbl(1,4) = -oneOverL;
   Tbl(1,2) =  1;
   Tbl(2,1) =  oneOverL;
   Tbl(2,5) =  1;
   Tbl(2,4) = -oneOverL;
}


const Vector &
LinearCrdTransf2d::getGlobalResistingForce(const Vector &pb, const Vector &unifLoad)
{
   // transform resisting forces from the basic system to local coordinates
   static Vector pl(6);
   pl.addMatrixTransposeVector(0.0, Tbl, pb, 1.0);    // pl = Tbl ^ pb;

   // add end forces due to uniforme distributed loads to the system with rigid body modes
   pl(0) -= unifLoad(0)*L;
   pl(1) -= unifLoad(1)*L/2;
   pl(4) -= unifLoad(1)*L/2;
         
   // include P-Delta effects
   if (pDeltaFlag)
   {
       double NoverL = pb(0)/L;
             
       pl(1) += (ul(1) - ul(4)) * NoverL;
       pl(4) -= (ul(1) - ul(4)) * NoverL;

       /*
       double coeff = pb(0)/(30*L);
       double L2 = L*L;
       static Matrix kg(6,6);
       cerr << " N = " << pb(0);
       
       kg.Zero();
       
       kg(1,1) =  36*coeff;
       kg(1,2) =  3*L*coeff;
       kg(1,4) = -36*coeff;
       kg(1,5) =  3*L*coeff;

       kg(2,1) =  3*L*coeff;
       kg(2,2) =  4*L2*coeff;
       kg(2,4) = -3*L*coeff;
       kg(2,5) = -L2*coeff;

       kg(4,1) = -36*coeff;
       kg(4,2) = -3*L*coeff;
       kg(4,4) =  36*coeff;
       kg(4,5) = -3*L*coeff;

       kg(5,1) =  3*L*coeff;
       kg(5,2) = -L2*coeff; 
       kg(5,4) = -3*L*coeff;
       kg(5,5) =  4*L2*coeff;

       pl.addMatrixVector(1.0, kg, ul, 1.0); 
       */
   }
   
   // transform resisting forces  from local to global coordinates
   static Vector pg(6);
   pg.addMatrixTransposeVector(0.0, Tlg, pl, 1.0);    // pg = Tlg ^ pl;

   return pg;
}
  


const Matrix &
LinearCrdTransf2d::getGlobalStiffMatrix (const Matrix &kb, const Vector &pb)
{
   // transform tangent stiffness matrix from the basic system to local coordinates
   static Matrix kl(6,6);

   //cerr << "Tbl: " << Tbl;
   
   //cerr << "Tlg: " << Tlg;
   
   kl.addMatrixTripleProduct(0.0, Tbl, kb, 1.0);      // kl = Tbl ^ kb * Tbl;

   //cerr << "kl: " << kl;
   // include P-Delta effects
   if (pDeltaFlag)
   {
       double NoverL = pb(0)/L;
       //cerr << "stiffness  NoverL: " << NoverL;
           
       kl(1,1) += NoverL;
       kl(4,4) += NoverL;
       kl(1,4) -= NoverL;
       kl(4,1) -= NoverL;
   
       /*
       double coeff = pb(0)/(30*L);
       double L2 = L*L;
       //cerr << "stiffness  NoverL: " << NoverL;

       kl(1,1) += 36*coeff;
       kl(1,2) += 3*L*coeff;
       kl(1,4) -= 36*coeff;
       kl(1,5) += 3*L*coeff;

       kl(2,1) += 3*L*coeff;
       kl(2,2) += 4*L2*coeff;
       kl(2,4) -= 3*L*coeff;
       kl(2,5) -= L2*coeff;

       kl(4,1) -= 36*coeff;
       kl(4,2) -= 3*L*coeff;
       kl(4,4) += 36*coeff;    
       kl(4,5) -= 3*L*coeff;

       kl(5,1) += 3*L*coeff;
       kl(5,2) -= L2*coeff; 
       kl(5,4) -= 3*L*coeff;
       kl(5,5) += 4*L2*coeff;    
       */
   }   
   
   // transform tangent  stiffness matrix from local to global coordinates
   static Matrix kg(6,6);
   kg.addMatrixTripleProduct(0.0, Tlg, kl,  1.0);     // kg =  Tlj ^ kl * Tlj;

   return kg;
}
  



CrdTransf2d *
LinearCrdTransf2d::getCopy(void)
{
  // create a new instance of LinearCrdTransf2d 

  LinearCrdTransf2d *theCopy = new LinearCrdTransf2d (this->getTag(), nodeIOffset, nodeJOffset, pDeltaFlag);
  
  if (!theCopy)
  {
     g3ErrorHandler->fatal("LinearCrdTransf2d::getCopy() - out of memory creating copy");
     return 0;
  }    
    
  theCopy->nodeIPtr = nodeIPtr;
  theCopy->nodeJPtr = nodeJPtr;
  theCopy->cosTheta = cosTheta;
  theCopy->sinTheta = sinTheta;
  theCopy->L = L;
  
  return theCopy;
}


int 
LinearCrdTransf2d::sendSelf(int cTag, Channel &theChannel)
{
	int res = 0;

	static Vector data(9);

	data(0) = this->getTag();
	data(1) = pDeltaFlag;
	data(2) = nodeIOffset(0);
	data(3) = nodeIOffset(1);
	data(4) = nodeJOffset(0);
	data(5) = nodeJOffset(1);
	data(6) = L;
	data(7) = cosTheta;
	data(8) = sinTheta;

	res += theChannel.sendVector(this->getDbTag(), cTag, data);
	if (res < 0) {
		g3ErrorHandler->warning("%s - failed to send Vector",
			"LinearCrdTransf2d::sendSelf");
		return res;
	}

    return res;
}

    

int 
LinearCrdTransf2d::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	int res = 0;

	static Vector data(9);

	res += theChannel.recvVector(this->getDbTag(), cTag, data);
	if (res < 0) {
		g3ErrorHandler->warning("%s - failed to receive Vector",
			"LinearCrdTransf2d::recvSelf");
		return res;
	}

	this->setTag((int)data(0));
	pDeltaFlag = (int)data(1);
	nodeIOffset(0) = data(2);
	nodeIOffset(1) = data(3);
	nodeJOffset(0) = data(4);
	nodeJOffset(1) = data(5);
	L = data(6);
	cosTheta = data(7);
	sinTheta = data(8);

    return res;
}
 	

const Vector &
LinearCrdTransf2d::getPointGlobalCoordFromLocal(const Vector &xl)
{
   static Vector xg(2);
   static Matrix Rlj(2,2);

   xg = nodeIPtr->getCrds() + nodeIOffset;

   Rlj(0,0) =  cosTheta;
   Rlj(0,1) =  sinTheta;
   Rlj(1,0) = -sinTheta;
   Rlj(1,1) =  cosTheta;
   
   // xg = xg + RljT*xl
   xg.addMatrixTransposeVector(1.0, Rlj, xl, 1.0);
     
   return xg;  
}

    
const Vector &
LinearCrdTransf2d::getPointGlobalDisplFromBasic (double xi, const Vector &uxb)
{
   // determine global displacements
   const Vector &disp1 = nodeIPtr->getTrialDisp();
   const Vector &disp2 = nodeJPtr->getTrialDisp();

   static Vector ug(6);
   for (int i = 0; i < 3; i++)
   {
      ug(i)   = disp1(i);
      ug(i+3) = disp2(i);
   }

   // transform global end displacements to local coordinates
   static Vector ul(6);      // total displacements

   ul.addMatrixVector(0.0, Tlg,  ug, 1.0);       //  ul = Tlg *  ug;
   
   // compute displacements at point xi, in local coordinates
   static Vector uxl(2),  uxg(2);

   uxl(0) = uxb(0) +        ul(0);
   uxl(1) = uxb(1) + (1-xi)*ul(1) + xi*ul(4);

   // rotate displacements to global coordinates
   // uxg = RljT*uxl
   static Matrix Rlj(2,2);
  
   Rlj(0,0) =  cosTheta;
   Rlj(0,1) =  sinTheta;
   Rlj(1,0) = -sinTheta;
   Rlj(1,1) =  cosTheta;
   
   uxg.addMatrixTransposeVector(0.0, Rlj, uxl, 1.0);
     
   return uxg;  
}




void
LinearCrdTransf2d::Print(ostream &s, int flag)
{
   s << "\nCrdTransf: " << this->getTag() << " Type: LinearCrdTransf2d";
   s << "\tnodeI Offset: " << nodeIOffset;
   s << "\tnodeJ Offset: " << nodeJOffset;
   s << "\tPDelta flag: "  << pDeltaFlag;
}








