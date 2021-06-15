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

// $Revision: 1.19 $
// $Date: 2009-08-25 22:32:08 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/Element.h,v $


// Written: jaabell
// Created: 06/2017
// Revision: A
//
// What: "@(#) Element.h, revA"
/*

Element is based on the following references by Carlos Felippa

Membrane and drilling parts:
(i)   Membrane Triangles with Corner Drilling Freedoms Part I: The EFF Element.
        Ken Alvin, Horacio M. de la Fuente, Bjorn Haugen, & Carlos A. Felippa.
        August 1991,
        University of Colorado, Boulder,
        Report No. CU-CDDC-91-24
(ii)  Membrane Triangles with Corner Drilling Freedoms Part II: The ANDES Element.
        Carlos A. Felippa & Carmello Militello,
        August 1991,
        University of Colorado, Boulder,
        Report No. CU-CDDC-91-24b
(iii) Membrane Triangles with Corner Drilling Freedoms Part III: Implementation and performance evaluation.
        Carlos A. Felippa & Scott Alexander,
        August 1991,
        University of Colorado, Boulder,
        Report No. CU-CDDC-91-24c
(iv)  A Study of Optimal Triangles with Drilling Freedoms.
        Carlos A. Felippa
        February 2003
        Report No. CU-CAS-03-02

For the bending component
(v)   The First ANDES Elements: 9-DOF Plate Bending Trianges
        Carmello Militello & Carlos A. Felippa
        December 1989
        Report No. CU-CSSC-89-22
(vi)  Chapter 32 of Felipp's Lecture Notes
      Finite element templates for bending

Obtainable as of July 22 2012 at http://www.colorado.edu/engineering/CAS/Felippa.d/FelippaHome.d/Home.html

These documents all mirror published works in indexed journals.
*/


#include <iomanip>
#include <iostream>

using namespace std;

#include <Vector.h>
#include <cmath>
#include <limits>
#include <Renderer.h>

#include <elementAPI.h>
#include <ElementResponse.h>
#include "ShellANDeS.h"


// Initialize static variables
double ShellANDeS::alpha_membrane = 1.5;          //Cfr. Reference (iii)
Vector ShellANDeS::beta_membrane(10);
unsigned int ShellANDeS::number_of_three_node_andes_membrane = 0;
Matrix ShellANDeS::Mq(9, 9);

// Using the optimal parameter set found in reference (iv) page 14

unsigned int ShellANDeS::number_of_three_node_andes_shells = 0;

Vector calculate_cross_product(const Vector& a, const Vector& b)
{
     Vector a_cross_b(3); // Store the result here

     if ( (a.Size() != 3) || (b.Size() != 3) )
     {
          opserr << "Error: calculate_cross_product only defined for 3x1 vectors.\n";
          exit(-1);
     }

     a_cross_b(0) =   a(1) * b(2) - b(1) * a(2);
     a_cross_b(1) = - a(0) * b(2) + b(0) * a(2);
     a_cross_b(2) =   a(0) * b(1) - b(0) * a(1);

     return a_cross_b;
}


void *
OPS_ShellANDeS(void)
{

     Element *theElement = 0;

     int numArgs = OPS_GetNumRemainingInputArgs();

     if (numArgs < 6)
     {
          opserr << "Want: element ShellANDeS $tag $iNode $jNode $kNode $thick $E $nu $rho";
          return 0;
     }

     int iData[4];
     int numData = 4;
     if (OPS_GetIntInput(&numData, iData) != 0)
     {
          opserr << "WARNING invalid integer tag: element ShellANDeS \n";
          return 0;
     }

     double dData[11];
     numArgs = OPS_GetNumRemainingInputArgs();
     if (OPS_GetDoubleInput(&numArgs, dData) != 0)
     {
          opserr << "WARNING invalid double thickness: element ShellANDeS \n";
          return 0;
     }
     if (numArgs == 4)
     {
          theElement = new ShellANDeS(iData[0], iData[1], iData[2], iData[3], dData[0], dData[1], dData[2], dData[3]);
     }
     else if (numArgs == 11)
     {
          theElement = new ShellANDeS(iData[0], iData[1], iData[2], iData[3],
                                      dData[0], dData[1], dData[2], dData[3], dData[4], dData[5],
                                      dData[6], dData[7], dData[8], dData[9], dData[10]);
     }


     return theElement;
}


//===================================================================================
// Element Constructors / Destructors
//===================================================================================

ShellANDeS::ShellANDeS()
     :
     Element(0, ELE_TAG_ShellANDeS ),
     connectedExternalNodes(3),
     K(18, 18), M(18, 18),
     P(18), Q(18), bf(3),
     is_stiffness_calculated(false),
     is_mass_calculated(false),
     thickness(0),
     xl1(3), xl2(3), xl3(3), x0(3), T_lg(3, 3),
     rho(0),
     E_planestress(3, 3), initialized_disps(false)
{
     // zero node pointers
     for (int i = 0; i < 3; i++)
     {
          theNodes[i] = 0;
     }

}


ShellANDeS::ShellANDeS(int element_number,
                       int node_numb_1, int node_numb_2, int node_numb_3, double t,
                       double E_, double nu_, double rho_)
     :
     Element(element_number, ELE_TAG_ShellANDeS ),
     connectedExternalNodes(3),
     K(18, 18), M(18, 18),
     P(18), Q(18), bf(3),
     is_stiffness_calculated(false),
     is_mass_calculated(false),
     thickness(t),
     xl1(3), xl2(3), xl3(3), x0(3), T_lg(3, 3),
     rho(rho_),
     E_planestress(3, 3), initialized_disps(false)
{

     double M  = E_ / ( 1.0 - nu_ * nu_ ) ; //membrane modulus
     double G  =  0.5 * E_ / ( 1.0 + nu_ ) ; //shear modulus

     // G *= thickness ;  //multiply by thickness
     // M *= thickness ;

     mE11 = M;
     mE22 = M;
     mE33 = G;
     mE12 = nu_ * M;
     mE13 = 0;
     mE23 = 0;

     n1 = n2 = n3 = 0.0;

     // Set connected external node IDs
     connectedExternalNodes(0) = node_numb_1;
     connectedExternalNodes(1) = node_numb_2;
     connectedExternalNodes(2) = node_numb_3;

     // zero node pointers
     for (int i = 0; i < 3; i++)
     {
          theNodes[i] = 0;
     }

     initializeBetaArrays();
}



ShellANDeS::ShellANDeS(int element_number, int node_numb_1, int node_numb_2, int node_numb_3, double t, double E11, double E22,
                       double E33, double E12, double E13, double E23, double n1_, double n2_, double n3_, double rho_)
     :
     Element(element_number, ELE_TAG_ShellANDeS ),
     connectedExternalNodes(3),
     K(18, 18), M(18, 18),
     P(18), Q(18), bf(3),
     is_stiffness_calculated(false),
     is_mass_calculated(false),
     thickness(t),
     xl1(3), xl2(3), xl3(3), x0(3), T_lg(3, 3),
     rho(rho_),
     E_planestress(3, 3),initialized_disps(false)
{

     // G *= thickness ;  //multiply by thickness
     // M *= thickness ;

     mE11 = E11;
     mE22 = E22;
     mE33 = E33;
     mE12 = E12;
     mE13 = E13;
     mE23 = E23;

     n1 = n1_;
     n2 = n2_;
     n3 = n3_;

     // Set connected external node IDs
     connectedExternalNodes(0) = node_numb_1;
     connectedExternalNodes(1) = node_numb_2;
     connectedExternalNodes(2) = node_numb_3;

     // zero node pointers
     for (int i = 0; i < 3; i++)
     {
          theNodes[i] = 0;
     }

     initializeBetaArrays();
}



ShellANDeS::~ShellANDeS()
{


}

//===================================================================================
// Element interface
//===================================================================================

int ShellANDeS::getNumExternalNodes() const
{
     return 3;
}

const ID &ShellANDeS::getExternalNodes()
{
     return connectedExternalNodes;
}

Node **ShellANDeS::getNodePtrs(void)
{
     return theNodes;
}

int ShellANDeS::getNumDOF()
{
     return 18;
}

void ShellANDeS::setDomain(Domain *theDomain)
{
     if (theDomain == 0)
     {
          theNodes[0] = 0;
          theNodes[1] = 0;
          theNodes[2] = 0;
     }
     else
     {
          int Nd1 = connectedExternalNodes(0);
          int Nd2 = connectedExternalNodes(1);
          int Nd3 = connectedExternalNodes(2);

          theNodes[0] = theDomain->getNode(Nd1);
          theNodes[1] = theDomain->getNode(Nd2);
          theNodes[2] = theDomain->getNode(Nd3);

          if (theNodes[0] == 0 || theNodes[1] == 0 || theNodes[2] == 0 )
          {
               opserr << "FATAL ERROR ShellANDeS (tag: " << this->getTag() << "), node not found in domain\n";
               exit(-1);
          }

          int dofNd1 = theNodes[0]->getNumberDOF();
          int dofNd2 = theNodes[1]->getNumberDOF();
          int dofNd3 = theNodes[2]->getNumberDOF();

          if (dofNd1 != 6 || dofNd2 != 6 || dofNd3 != 6)
          {
               opserr << "FATAL ERROR ShellANDeS (tag: " << this->getTag() << "), has differing number of DOFs at its nodes\n";
               exit(-1);
          }

          this->DomainComponent::setDomain(theDomain);
          initializeGeometry(n1, n2, n3);

          if(!initialized_disps)
          {
               for (int node = 0; node < 3; node++)
               {
                    const Vector &disps = theNodes[node]->getTrialDisp();
                    disp_init[node][0] = disps(0);
                    disp_init[node][1] = disps(1);
                    disp_init[node][2] = disps(2);
                    disp_init[node][3] = disps(3);
                    disp_init[node][4] = disps(4);
                    disp_init[node][5] = disps(5);
                    // opserr << disps << endln;
               }
               initialized_disps=true;
          }
     }
}

int ShellANDeS::commitState ()
{
     return 0;
}

int ShellANDeS::revertToLastCommit ()
{
     return 0;
}

int ShellANDeS::revertToStart ()
{
     return 0;
}

int ShellANDeS::update(void)
{
     return 0;
}

const Matrix &ShellANDeS::getTangentStiff ()
{
     // Since this is a linear element, the stiffness is calculated only once.
     if (!is_stiffness_calculated)
     {
          // Form completestiffness
          Matrix Kb(18, 18);       // Bending stiffness
          Matrix Km(18, 18);       // Membrane stiffness

          //Get matrices
          // Kb = bending_element -> getTangentStiff();
          Kb = getBendingTangentStiffness();
          Km = getMembraneTangentStiffness();
          // Km = membrane_element -> getTangentStiff();

          K = Kb + Km;

          is_stiffness_calculated = true;
     }

     return K;
}

const Matrix &ShellANDeS::getInitialStiff()
{
     return getTangentStiff();
}

const Matrix &ShellANDeS::getMass ()
{
     if (!is_mass_calculated)
     {
          // Form completestiffness
          Matrix Mb(18, 18);       // Bending mass
          Matrix Mm(18, 18);       // Membrane mass

          //Get matrices
          Mb = getBendingMass();
          Mm = getMembraneMass();

          M = Mb + Mm;

          is_mass_calculated = true;
     }


     return M;
}

void ShellANDeS::zeroLoad ()
{
     Q.Zero();
}

int ShellANDeS::addLoad(ElementalLoad *theLoad, double loadFactor)
{


     int type;
     const Vector &data = theLoad->getData(type, loadFactor);

     if (type == LOAD_TAG_SelfWeight)
     {

          Vector Fbody = this->getBodyForce(loadFactor, data);

          // opserr << "   loadfactor = " << loadFactor << endln;
          // opserr << "   data       = " << data << endln;
          // opserr << "   Fbody      = " << Fbody << endln;

          Q.addVector(1.0, Fbody, 1.0);

     }
     //    else if (type == LOAD_TAG_BrickSurfaceLoad8Node)
     //    {
     //        Vector Fsurface = this->getSurfaceForce(loadFactor, data);
     //
     //        Q.addVector(1.0, Fsurface, 1.0);
     //
     //    }
     else
     {
          opserr << "ShellANDeS::addLoad() - addLoad " << this->getTag() << ",load type " << type << "unknown\n";
          return -1;
     }



     return 0;
}

int ShellANDeS::addInertiaLoadToUnbalance(const Vector &accel)
{
     // Get R * accel from the nodes
     const Vector &Raccel1 = theNodes[0]->getRV(accel);
     const Vector &Raccel2 = theNodes[1]->getRV(accel);
     const Vector &Raccel3 = theNodes[2]->getRV(accel);

     if (6 != Raccel1.Size() || 6 != Raccel2.Size() || 6 != Raccel3.Size()  )
     {
          // Xiaoyan changed 2 to 3 and added Racce15-18  09/27/00
          opserr << "ShellANDeS::addInertiaLoadToUnbalance matrix and vector sizes are incompatable\n";
          return -1;
     }

     static Vector ra(18);  // Changed form 8 to 24(3*8)  Xiaoyan 09/27/00

     ra.Zero();

     ra(0)  =  Raccel1(0);
     ra(1)  =  Raccel1(1);
     ra(2)  =  Raccel1(2);
     ra(3)  =  Raccel1(3);
     ra(4)  =  Raccel1(4);
     ra(5)  =  Raccel1(5);

     ra(6)  =  Raccel2(0);
     ra(7)  =  Raccel2(1);
     ra(8)  =  Raccel2(2);
     ra(9)  =  Raccel2(3);
     ra(10) =  Raccel2(4);
     ra(11) =  Raccel2(5);

     ra(12) =  Raccel3(0);
     ra(13) =  Raccel3(1);
     ra(14) =  Raccel3(2);
     ra(15) =  Raccel3(3);
     ra(16) =  Raccel3(4);
     ra(17) =  Raccel3(5);

     Q.addMatrixVector(1.0, M, ra, -1.0);
     //Q.addMatrixVector(1.0, M, ra, 1.0);

     return 0;
}

const Vector &ShellANDeS::getResistingForce ()
{
     P.Zero();

     static Vector NodalDisplacements(18);
     static Vector disp_i(6);

     NodalDisplacements.Zero();
     disp_i.Zero();

     for (int node = 0; node < 3; node++)
     {
          Node *node_i = theNodes[node];
          disp_i = node_i->getDisp();
          //disp_i += node_i->getIncrDeltaDisp();
          disp_i += node_i->getIncrDisp();

          NodalDisplacements(6 * node + 0) = disp_i(0) - disp_init[node][0];;
          NodalDisplacements(6 * node + 1) = disp_i(1) - disp_init[node][1];;
          NodalDisplacements(6 * node + 2) = disp_i(2) - disp_init[node][2];;
          NodalDisplacements(6 * node + 3) = disp_i(3) - disp_init[node][3];;
          NodalDisplacements(6 * node + 4) = disp_i(4) - disp_init[node][4];;
          NodalDisplacements(6 * node + 5) = disp_i(5) - disp_init[node][5];;
     }

     P = K * NodalDisplacements;
     // Add nodal equivalent loads from surface or body loads
     P += Q;

     return P;
}

const Vector &ShellANDeS::getResistingForceIncInertia ()
{
     static Vector NodalDisplacements(18);
     static Vector disp_i(6);
     static Vector NodalAccelerations(18);
     static Vector accel_i(6);

     NodalDisplacements.Zero();
     disp_i.Zero();
     NodalAccelerations.Zero();
     accel_i.Zero();

     for (int node = 0; node < 3; node++)
     {
          Node *node_i = theNodes[node];
          disp_i = node_i->getDisp();
          //disp_i += node_i->getIncrDeltaDisp();
          disp_i += node_i->getIncrDisp();

          accel_i = node_i->getTrialAccel();

          NodalDisplacements(6 * node + 0) = disp_i(0) - disp_init[node][0];;
          NodalDisplacements(6 * node + 1) = disp_i(1) - disp_init[node][1];;
          NodalDisplacements(6 * node + 2) = disp_i(2) - disp_init[node][2];;
          NodalDisplacements(6 * node + 3) = disp_i(3) - disp_init[node][3];;
          NodalDisplacements(6 * node + 4) = disp_i(4) - disp_init[node][4];;
          NodalDisplacements(6 * node + 5) = disp_i(5) - disp_init[node][5];;

          NodalAccelerations(6 * node + 0) = accel_i(0);
          NodalAccelerations(6 * node + 1) = accel_i(1);
          NodalAccelerations(6 * node + 2) = accel_i(2);
          NodalAccelerations(6 * node + 3) = accel_i(3);
          NodalAccelerations(6 * node + 4) = accel_i(4);
          NodalAccelerations(6 * node + 5) = accel_i(5);
     }

     P = K * NodalDisplacements + this->getMass() * NodalAccelerations;

     // Add nodal equivalent loads from surface or body loads
     P += Q;

     return P;
}


int ShellANDeS::sendAndCheckVector(int dataTag , int commitTag , Vector& v , Channel & channel , const std::string &name)
{
     int res = channel.sendVector(dataTag, commitTag, v);
     if (res < 0) {
         opserr << "WARNING LysmerTriangle::sendAndCheckVector() - " << this->getTag() << " failed to send " << name.c_str() << endln;
         return -1;
     }  
     return res;
}

int ShellANDeS::recvAndCheckVector(int dataTag , int commitTag , Vector& v , Channel & channel , const std::string &name)
{
     int res = channel.recvVector(dataTag, commitTag, v);
     if (res < 0) {
         opserr << "WARNING LysmerTriangle::recvAndCheckVector() - " << this->getTag() << " failed to receive " << name.c_str() << endln;
         return -1;
     }  
     return res;
}

int ShellANDeS::sendAndCheckMatrix(int dataTag , int commitTag , Matrix& v , Channel & channel , const std::string &name)
{
     int res = channel.sendMatrix(dataTag, commitTag, v);
     if (res < 0) {
         opserr << "WARNING LysmerTriangle::sendAndCheckMatrix() - " << this->getTag() << " failed to send " << name.c_str() << endln;
         return -1;
     }  
     return res;
}

int ShellANDeS::recvAndCheckMatrix(int dataTag , int commitTag , Matrix& v , Channel & channel , const std::string &name)
{
     int res = channel.recvMatrix(dataTag, commitTag, v);
     if (res < 0) {
         opserr << "WARNING LysmerTriangle::recvAndCheckMatrix() - " << this->getTag() << " failed to receive " << name.c_str() << endln;
         return -1;
     }  
     return res;
}


int ShellANDeS::sendAndCheckID(int dataTag     , int commitTag , ID& v     , Channel & channel , const std::string &name)
{
     int res = channel.sendID(dataTag, commitTag, v);
     if (res < 0) {
         opserr << "WARNING LysmerTriangle::sendAndCheckID() - " << this->getTag() << " failed to send " << name.c_str() << endln;
         return -1;
     }  
     return res;
}

int ShellANDeS::recvAndCheckID(int dataTag     , int commitTag , ID& v     , Channel & channel , const std::string &name)
{
     int res = channel.recvID(dataTag, commitTag, v);
     if (res < 0) {
         opserr <<" WARNING LysmerTriangle::recvAndCheckID() - " << this->getTag() << " failed to receive " << name.c_str() << endln;
         return -1;
     }  
     return res;
}


int ShellANDeS::sendSelf (int commitTag, Channel &theChannel)
{
     int pos = 0;
     int dataTag = this->getDbTag();

     static ID idata(4);
     idata(pos++) = this->getTag();
     idata(pos++) = connectedExternalNodes(0);
     idata(pos++) = connectedExternalNodes(1);
     idata(pos++) = connectedExternalNodes(2);
     
     sendAndCheckID(commitTag, dataTag, idata, theChannel, "idata");


     static Vector ddata(21+18);
     pos = 0;
     ddata(pos++) = thickness;
     ddata(pos++) = Area;
     ddata(pos++) = x12;
     ddata(pos++) = x23;
     ddata(pos++) = x31;
     ddata(pos++) = y12;
     ddata(pos++) = y23;
     ddata(pos++) = y31;
     ddata(pos++) = rho;
     ddata(pos++) = mE11;
     ddata(pos++) = mE22;
     ddata(pos++) = mE33;
     ddata(pos++) = mE12;
     ddata(pos++) = mE13;
     ddata(pos++) = mE23;
     ddata(pos++) = n1;
     ddata(pos++) = n2;
     ddata(pos++) = n3;
     ddata(pos++) = alpha_membrane;
     ddata(pos++) = beta0;

     for (int i = 0; i < 18; ++i)
          ddata(pos++) = disp_init[i/6][i%6];
     
     ddata(pos++) = initialized_disps;

     sendAndCheckVector(commitTag, dataTag, ddata, theChannel, "ddata");
     sendAndCheckVector(commitTag, dataTag, P, theChannel, "P");
     sendAndCheckVector(commitTag, dataTag, Q, theChannel, "Q");
     sendAndCheckVector(commitTag, dataTag, bf, theChannel, "bf");
     sendAndCheckVector(commitTag, dataTag, xl1, theChannel, "xl1");
     sendAndCheckVector(commitTag, dataTag, xl2, theChannel, "xl2");
     sendAndCheckVector(commitTag, dataTag, xl3, theChannel, "xl3");
     sendAndCheckVector(commitTag, dataTag, x0, theChannel, "x0");
     sendAndCheckMatrix(commitTag, dataTag, T_lg, theChannel, "T_lg");
     sendAndCheckMatrix(commitTag, dataTag, E_planestress, theChannel, "E_planestress");

     return 0;
}

int ShellANDeS::recvSelf (int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
     int pos = 0;
     int dataTag = this->getDbTag();

     static ID idata(4);
     recvAndCheckID(commitTag, dataTag, idata, theChannel, "idata");
     
     this->setTag(idata(pos++));
     connectedExternalNodes(0) = idata(pos++);
     connectedExternalNodes(1) = idata(pos++);
     connectedExternalNodes(2) = idata(pos++);
     


     static Vector ddata(21+18);
     recvAndCheckVector(commitTag, dataTag, ddata, theChannel, "ddata");
     pos = 0;
     thickness = ddata(pos++);
     Area = ddata(pos++);
     x12 = ddata(pos++);
     x23 = ddata(pos++);
     x31 = ddata(pos++);
     y12 = ddata(pos++);
     y23 = ddata(pos++);
     y31 = ddata(pos++);
     rho = ddata(pos++);
     mE11 = ddata(pos++);
     mE22 = ddata(pos++);
     mE33 = ddata(pos++);
     mE12 = ddata(pos++);
     mE13 = ddata(pos++);
     mE23 = ddata(pos++);
     n1 = ddata(pos++);
     n2 = ddata(pos++);
     n3 = ddata(pos++);
     alpha_membrane = ddata(pos++);
     beta0 = ddata(pos++);
     for (int i = 0; i < 18; ++i)
          disp_init[i/6][i%6] = ddata(pos++);
     initialized_disps = ddata(pos++);

     recvAndCheckVector(commitTag, dataTag, P, theChannel, "P");
     recvAndCheckVector(commitTag, dataTag, Q, theChannel, "Q");
     recvAndCheckVector(commitTag, dataTag, bf, theChannel, "bf");
     recvAndCheckVector(commitTag, dataTag, xl1, theChannel, "xl1");
     recvAndCheckVector(commitTag, dataTag, xl2, theChannel, "xl2");
     recvAndCheckVector(commitTag, dataTag, xl3, theChannel, "xl3");
     recvAndCheckVector(commitTag, dataTag, x0, theChannel, "x0");
     recvAndCheckMatrix(commitTag, dataTag, T_lg, theChannel, "T_lg");
     recvAndCheckMatrix(commitTag, dataTag, E_planestress, theChannel, "E_planestress");


     return 0;
}

void ShellANDeS::Print(OPS_Stream &s, int flag)
{
     if (flag == OPS_PRINT_CURRENTSTATE) {
          s << "\nShell ANDeS ----- tag = " << this->getTag() << endln;
          s << "       connectedExternalNodes = " << connectedExternalNodes;
          s << "       thickness = " << thickness << endln;
          s << "       xl1 = " << xl1;
          s << "       xl2 = " << xl2;
          s << "       xl3 = " << xl3;
          s << "       x0 = " << x0;
          s << "       Area = " << Area << endln;
          s << "       x12 = " << x12 << endln;
          s << "       x23 = " << x23 << endln;
          s << "       x31 = " << x31 << endln;
          s << "       y12 = " << y12 << endln;
          s << "       y23 = " << y23 << endln;
          s << "       y31 = " << y31 << endln;
          s << "       mE11 = " << mE11 << endln;
          s << "       mE22 = " << mE22 << endln;
          s << "       mE33 = " << mE33 << endln;
          s << "       mE12 = " << mE12 << endln;
          s << "       mE13 = " << mE13 << endln;
          s << "       mE23 = " << mE23 << endln;
          s << "       rho = " << rho << endln;
     }

     if (flag == OPS_PRINT_PRINTMODEL_JSON) {
          s << "\t\t\t{";
          s << "\"name\": " << this->getTag() << ", ";
          s << "\"type\": \"ShellANDeS\", ";
          s << "\"nodes\": [" << connectedExternalNodes(0) << ", ";
          s << connectedExternalNodes(1) << ", ";
          s << connectedExternalNodes(2) << "], ";
          s << "\"thickness\": " << thickness << ", ";
          s << "       mE11 = " << mE11 << endln;
          s << "       mE22 = " << mE22 << endln;
          s << "       mE33 = " << mE33 << endln;
          s << "       mE12 = " << mE12 << endln;
          s << "       mE13 = " << mE13 << endln;
          s << "       mE23 = " << mE23 << endln;
          s << "\"masspervolume\": " << rho << "\"}";
     }
}

Response* ShellANDeS::setResponse (const char** argv, int argc, OPS_Stream& theHandler)
{
     if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "forces") == 0)
     {
          return new ElementResponse(this, 1, P);
     }
     else if (strcmp(argv[0], "stiff") == 0 || strcmp(argv[0], "stiffness") == 0)
     {
          return new ElementResponse(this, 5, K);
     }
     else if (strcmp(argv[0], "moments") == 0 || strcmp(argv[0], "stresses") == 0)
     {
          return new ElementResponse(this, 1313, Vector(3));
     }
     else
     {
          return 0;  // Return a null pointer
     }
}

int ShellANDeS::getResponse (int responseID, Information& eleInformation)
{
     if (responseID == 1) //forces
     {
          return eleInformation.setVector(P);
     }
     else if (responseID == 5) //stiffness
     {
          return eleInformation.setMatrix(K);
     }
     else if (responseID == 1313) //bending-moment
     {
          return eleInformation.setVector(get_bending_moment_field());
     }
     else
     {
          return -1;
     }
}

Matrix ShellANDeS::returnMass(void)
{
     if (!is_mass_calculated)
     {
          return getMass();
     }
     else
     {
          return M;
     }
}

void ShellANDeS::useThisCoordinateSystem(Vector e1, Vector e2, Vector e3)
{

     //Local-to-global transformation matrix
     for (int i = 0; i < 3; i++)
     {
          T_lg(i, 0) = e1(i);
          T_lg(i, 1) = e2(i);
          T_lg(i, 2) = e3(i);
     }

     //Matrix T_gl(3,3);
     //T_lg.Invert(T_gl);  //Compute global -> local transformation matrix

     // Project global nodal coordinates into local coordinates
     T_lg.Solve(theNodes[0]->getCrds() - x0, xl1);  // Solve T_lg * xl1 = x1
     T_lg.Solve(theNodes[1]->getCrds() - x0, xl2);  // Solve T_lg * xl2 = x2
     T_lg.Solve(theNodes[2]->getCrds() - x0, xl3);  // Solve T_lg * xl3 = x3

     // Node differences
     x12 = xl1(0) - xl2(0);
     x23 = xl2(0) - xl3(0);
     x31 = xl3(0) - xl1(0);
     y12 = xl1(1) - xl2(1);
     y23 = xl2(1) - xl3(1);
     y31 = xl3(1) - xl1(1);

}

double ShellANDeS::getArea() const
{
     return Area;
}

bool ShellANDeS::gotMass() const
{

     return is_mass_calculated;
}

const Vector &ShellANDeS::getBodyForce(double loadFactor, const Vector &data)
{
     static Vector bforce(18);
     static Vector ba(18);
     static Vector bfx(3);

     bforce.Zero();
     ba.Zero();
     bfx.Zero();

     //    Vector bf(3);
     bf(0) = data(0);
     bf(1) = data(1);
     bf(2) = data(2);


     bfx = bf * loadFactor;

     ba(0) =  bfx(0);
     ba(1) =  bfx(1);
     ba(2) =  bfx(2);
     ba(3) =  0.0;
     ba(4) =  0.0;
     ba(5) =  0.0;
     ba(6) =  bfx(0);
     ba(7) =  bfx(1);
     ba(8) =  bfx(2);
     ba(9) =  0.0;
     ba(10) = 0.0;
     ba(11) = 0.0;
     ba(12) = bfx(0);
     ba(13) = bfx(1);
     ba(14) = bfx(2);
     ba(15) = 0.0;
     ba(16) = 0.0;
     ba(17) = 0.0;

     //Form equivalent body force
     this->getMass();

     bforce.addMatrixVector(0.0, M, ba, 1.0);

     return bforce;
}



const Vector & ShellANDeS::get_bending_moment_field()
{
     static Vector m(3);
     static Vector disps_global(18);
     static Vector disps_local(18);
     // static Vector disps_local_bending(9);
     static Vector d1(6); 
     static Vector d2(6); 
     static Vector d3(6); 
     static Matrix TLG(18, 18);             // Local-to-global transformation matrix
     static Matrix L(3, 9);
     static Matrix EL(3, 9);

     m.Zero();
     disps_global.Zero();
     disps_local.Zero();
     // disps_local_bending.Zero();

     TLG.Zero();
     L.Zero();
     EL.Zero();

     d1 = theNodes[0]->getTrialDisp();
     d2 = theNodes[1]->getTrialDisp();
     d3 = theNodes[2]->getTrialDisp();

     // Form TLG
     int I, J;
     static Matrix T_gl(3, 3);
     T_gl.Zero();
     T_lg.Invert(T_gl);
     int offset = 0;

     for (int n = 0; n < 6; n++)
     {
          I = 0;

          for (int i = 0; i < 3; i++)
          {
               J = 0;

               for (int j = 0; j < 3; j++)
               {
                    TLG(I + offset, J + offset) = T_gl(i, j);
                    J++;
               }

               I++;
          }
          offset += 3;
     }

     disps_global(0)  = d1(0) - disp_init[0][0];
     disps_global(1)  = d1(1) - disp_init[0][1];
     disps_global(2)  = d1(2) - disp_init[0][2];
     disps_global(3)  = d1(3) - disp_init[0][3];
     disps_global(4)  = d1(4) - disp_init[0][4];
     disps_global(5)  = d1(5) - disp_init[0][5];
     disps_global(6)  = d2(0) - disp_init[1][0];
     disps_global(7)  = d2(1) - disp_init[1][1];
     disps_global(8)  = d2(2) - disp_init[1][2];
     disps_global(9)  = d2(3) - disp_init[1][3];
     disps_global(10) = d2(4) - disp_init[1][4];
     disps_global(11) = d2(5) - disp_init[1][5];
     disps_global(12) = d3(0) - disp_init[2][0];
     disps_global(13) = d3(1) - disp_init[2][1];
     disps_global(14) = d3(2) - disp_init[2][2];
     disps_global(15) = d3(3) - disp_init[2][3];
     disps_global(16) = d3(4) - disp_init[2][4];
     disps_global(17) = d3(5) - disp_init[2][5];


     disps_local.addMatrixVector(1.0, TLG, disps_global, 1.0);

     // disps_local_bending(0) = disps_local(2);
     // disps_local_bending(1) = disps_local(3);
     // disps_local_bending(2) = disps_local(4);
     // disps_local_bending(3) = disps_local(8);
     // disps_local_bending(4) = disps_local(9);
     // disps_local_bending(5) = disps_local(10);
     // disps_local_bending(6) = disps_local(14);
     // disps_local_bending(7) = disps_local(15);
     // disps_local_bending(8) = disps_local(16);


     // Geometric data
     double x21 = -x12;
     double y21 = -y12;
     double x32 = -x23;
     double y32 = -y23;
     double x13 = -x31;
     double y13 = -y31;

     double l12 = sqrt(x12 * x12 + y12 * y12);
     double l23 = sqrt(x23 * x23 + y23 * y23);
     double l31 = sqrt(x31 * x31 + y31 * y31);

     double c12 = x21 / l12;
     double c23 = x32 / l23;
     double c31 = x13 / l31;
     double s12 = y21 / l12;
     double s23 = y32 / l23;
     double s31 = y13 / l31;



     L(0, 0) = -c12 * s12 + c31 * s31;
     L(1, 0) = -c31 * s31 + c12 * s12;
     L(2, 0) = (s31 * s31 - c31 * c31) - (s12 * s12 - c12 * c12);
     L(0, 1) = (s12 * s12 * x12 + s31 * s31 * x31) / 2.0;
     L(1, 1) = (c12 * c12 * x12 + c31 * c31 * x31) / 2.0;
     L(2, 1) = c12 * c12 * y21 + c31 * c31 * y13;
     L(0, 2) = (s12 * s12 * y12 + s31 * s31 * y31) / 2.0;
     L(1, 2) = (c12 * c12 * y12 + c31 * c31 * y31) / 2.0;
     L(2, 2) = s12 * s12 * x21 + s31 * s31 * x13;
     L(0, 3) = -c23 * s23 + c12 * s12;
     L(1, 3) = -c12 * s12 + c23 * s23;
     L(2, 3) = (s12 * s12 - c12 * c12) - (s23 * s23 - c23 * c23);
     L(0, 4) = (s12 * s12 * x12 + s23 * s23 * x23) / 2.0;
     L(1, 4) = (c12 * c12 * x12 + c23 * c23 * x23) / 2.0;
     L(2, 4) = c12 * c12 * y21 + c23 * c23 * y32;
     L(0, 5) = (s12 * s12 * y12 + s23 * s23 * y23) / 2.0;
     L(1, 5) = (c12 * c12 * y12 + c23 * c23 * y23) / 2.0;
     L(2, 5) = s12 * s12 * x21 + s23 * s23 * x32;
     L(0, 6) = -c31 * s31 + c23 * s23;
     L(1, 6) = -c23 * s23 + c31 * s31;
     L(2, 6) = (s23 * s23 - c23 * c23) - (s31 * s31 - c31 * c31);
     L(0, 7) = (s23 * s23 * x23 + s31 * s31 * x31) / 2.0;
     L(1, 7) = (c23 * c23 * x23 + c31 * c31 * x31) / 2.0;
     L(2, 7) = c23 * c23 * y32 + c31 * c31 * y13;
     L(0, 8) = (s23 * s23 * y23 + s31 * s31 * y31) / 2.0;
     L(1, 8) = (c23 * c23 * y23 + c31 * c31 * y31) / 2.0;
     L(2, 8) = s23 * s23 * x32 + s31 * s31 * x13;

     EL.Zero();
     EL.addMatrixProduct(thickness * thickness * thickness / 12 / Area, E_planestress, L, 1.0);

     m.Zero();
     m.addMatrixVector(1.0, EL , disps_local, 1.0);


     // High-order part
     // double lam_12 =



     return m;
}




// Vector ShellANDeS::get_bending_moment_field()
// {

//     Vector m(3);
//     Vector disps(9);

//     Vector d1 = theNodes[0]->getTrialDisp();
//     Vector d2 = theNodes[1]->getTrialDisp();
//     Vector d3 = theNodes[2]->getTrialDisp();

//     disps(0) = d1(0);
//     disps(1) = d1(1);
//     disps(2) = d1(2);
//     disps(3) = d2(0);
//     disps(4) = d2(1);
//     disps(5) = d2(2);
//     disps(6) = d3(0);
//     disps(7) = d3(1);
//     disps(8) = d3(2);

//     // Geometric data
//     double x21 = -x12;
//     double y21 = -y12;
//     double x32 = -x23;
//     double y32 = -y23;
//     double x13 = -x31;
//     double y13 = -y31;

//     double l12 = sqrt(x12 * x12 + y12 * y12);
//     double l23 = sqrt(x23 * x23 + y23 * y23);
//     double l31 = sqrt(x31 * x31 + y31 * y31);

//     double c12 = x21 / l12;
//     double c23 = x32 / l23;
//     double c31 = x13 / l31;
//     double s12 = y21 / l12;
//     double s23 = y32 / l23;
//     double s31 = y13 / l31;


//     Matrix L(3, 9);
//     Matrix EL(3, 9);

//     L(0, 0) = -c12 * s12 + c31 * s31;
//     L(1, 0) = -c31 * s31 + c12 * s12;
//     L(2, 0) = (s31 * s31 - c31 * c31) - (s12 * s12 - c12 * c12);
//     L(0, 1) = (s12 * s12 * x12 + s31 * s31 * x31) / 2.0;
//     L(1, 1) = (c12 * c12 * x12 + c31 * c31 * x31) / 2.0;
//     L(2, 1) = c12 * c12 * y21 + c31 * c31 * y13;
//     L(0, 2) = (s12 * s12 * y12 + s31 * s31 * y31) / 2.0;
//     L(1, 2) = (c12 * c12 * y12 + c31 * c31 * y31) / 2.0;
//     L(2, 2) = s12 * s12 * x21 + s31 * s31 * x13;
//     L(0, 3) = -c23 * s23 + c12 * s12;
//     L(1, 3) = -c12 * s12 + c23 * s23;
//     L(2, 3) = (s12 * s12 - c12 * c12) - (s23 * s23 - c23 * c23);
//     L(0, 4) = (s12 * s12 * x12 + s23 * s23 * x23) / 2.0;
//     L(1, 4) = (c12 * c12 * x12 + c23 * c23 * x23) / 2.0;
//     L(2, 4) = c12 * c12 * y21 + c23 * c23 * y32;
//     L(0, 5) = (s12 * s12 * y12 + s23 * s23 * y23) / 2.0;
//     L(1, 5) = (c12 * c12 * y12 + c23 * c23 * y23) / 2.0;
//     L(2, 5) = s12 * s12 * x21 + s23 * s23 * x32;
//     L(0, 6) = -c31 * s31 + c23 * s23;
//     L(1, 6) = -c23 * s23 + c31 * s31;
//     L(2, 6) = (s23 * s23 - c23 * c23) - (s31 * s31 - c31 * c31);
//     L(0, 7) = (s23 * s23 * x23 + s31 * s31 * x31) / 2.0;
//     L(1, 7) = (c23 * c23 * x23 + c31 * c31 * x31) / 2.0;
//     L(2, 7) = c23 * c23 * y32 + c31 * c31 * y13;
//     L(0, 8) = (s23 * s23 * y23 + s31 * s31 * y31) / 2.0;
//     L(1, 8) = (c23 * c23 * y23 + c31 * c31 * y31) / 2.0;
//     L(2, 8) = s23 * s23 * x32 + s31 * s31 * x13;

//     EL.Zero();
//     EL.addMatrixProduct(1.0, E_planestress * (thickness * thickness * thickness / 12 / Area), L, 1.0);

//     m.Zero();
//     m.addMatrixVector(1.0, EL , disps, 1.0);


//     return m;
// }

const Matrix &ShellANDeS::getMembraneMass ()
{

     double x13 = -x31;
     double y13 = -y31;
     double alpha = 1.5;

     static Matrix TLG(18, 18);             // Local-to-global transformation matrix
     static Matrix Ml(9, 9);
     static Matrix Mlocal(18, 18);

     TLG.Zero();
     Ml.Zero();
     Mlocal.Zero();

     Ml(0, 0) = 30;
     Ml(0, 1) = 0;
     Ml(0, 2) = 3 * alpha * y12 + 3 * alpha * y13;
     Ml(0, 3) = 15;
     Ml(0, 4) = 0;
     Ml(0, 5) = -3 * alpha * y12 + 3 * alpha * y23 / 2;
     Ml(0, 6) = 15;
     Ml(0, 7) = 0;
     Ml(0, 8) = -3 * alpha * y13 - 3 * alpha * y23 / 2;
     Ml(1, 0) = 0;
     Ml(1, 1) = 30;
     Ml(1, 2) = -3 * alpha * x12 - 3 * alpha * x13;
     Ml(1, 3) = 0;
     Ml(1, 4) = 15;
     Ml(1, 5) = 3 * alpha * x12 - 3 * alpha * x23 / 2;
     Ml(1, 6) = 0;
     Ml(1, 7) = 15;
     Ml(1, 8) = 3 * alpha * x13 + 3 * alpha * x23 / 2;
     Ml(2, 0) = 3 * alpha * y12 + 3 * alpha * y13;
     Ml(2, 1) = -3 * alpha * x12 - 3 * alpha * x13;
     Ml(2, 2) = -45 * alpha * x12 * (-alpha * x12 / 45 - alpha * x13 / 90) / 2 - 45 * alpha * x13 * (-alpha * x12 / 90 - alpha * x13 / 45) / 2 + 45 * alpha * y12 * (alpha * y12 / 45 + alpha * y13 / 90) / 2 + 45 * alpha * y13 * (alpha * y12 / 90 + alpha * y13 / 45) / 2;
     Ml(2, 3) = 3 * alpha * y12 + 3 * alpha * y13 / 2;
     Ml(2, 4) = -3 * alpha * x12 - 3 * alpha * x13 / 2;
     Ml(2, 5) = 45 * alpha * x12 * (-alpha * x12 / 45 - alpha * x13 / 90) / 2 - 45 * alpha * x23 * (-alpha * x12 / 90 - alpha * x13 / 90) / 2 - 45 * alpha * y12 * (alpha * y12 / 45 + alpha * y13 / 90) / 2 + 45 * alpha * y23 * (alpha * y12 / 90 + alpha * y13 / 90) / 2;
     Ml(2, 6) = 3 * alpha * y12 / 2 + 3 * alpha * y13;
     Ml(2, 7) = -3 * alpha * x12 / 2 - 3 * alpha * x13;
     Ml(2, 8) = 45 * alpha * x13 * (-alpha * x12 / 90 - alpha * x13 / 45) / 2 + 45 * alpha * x23 * (-alpha * x12 / 90 - alpha * x13 / 90) / 2 - 45 * alpha * y13 * (alpha * y12 / 90 + alpha * y13 / 45) / 2 - 45 * alpha * y23 * (alpha * y12 / 90 + alpha * y13 / 90) / 2;
     Ml(3, 0) = 15;
     Ml(3, 1) = 0;
     Ml(3, 2) = 3 * alpha * y12 + 3 * alpha * y13 / 2;
     Ml(3, 3) = 30;
     Ml(3, 4) = 0;
     Ml(3, 5) = -3 * alpha * y12 + 3 * alpha * y23;
     Ml(3, 6) = 15;
     Ml(3, 7) = 0;
     Ml(3, 8) = -3 * alpha * y13 / 2 - 3 * alpha * y23;
     Ml(4, 0) = 0;
     Ml(4, 1) = 15;
     Ml(4, 2) = -3 * alpha * x12 - 3 * alpha * x13 / 2;
     Ml(4, 3) = 0;
     Ml(4, 4) = 30;
     Ml(4, 5) = 3 * alpha * x12 - 3 * alpha * x23;
     Ml(4, 6) = 0;
     Ml(4, 7) = 15;
     Ml(4, 8) = 3 * alpha * x13 / 2 + 3 * alpha * x23;
     Ml(5, 0) = -3 * alpha * y12 + 3 * alpha * y23 / 2;
     Ml(5, 1) = 3 * alpha * x12 - 3 * alpha * x23 / 2;
     Ml(5, 2) = -45 * alpha * x12 * (alpha * x12 / 45 - alpha * x23 / 90) / 2 - 45 * alpha * x13 * (alpha * x12 / 90 - alpha * x23 / 90) / 2 + 45 * alpha * y12 * (-alpha * y12 / 45 + alpha * y23 / 90) / 2 + 45 * alpha * y13 * (-alpha * y12 / 90 + alpha * y23 / 90) / 2;
     Ml(5, 3) = -3 * alpha * y12 + 3 * alpha * y23;
     Ml(5, 4) = 3 * alpha * x12 - 3 * alpha * x23;
     Ml(5, 5) = 45 * alpha * x12 * (alpha * x12 / 45 - alpha * x23 / 90) / 2 - 45 * alpha * x23 * (alpha * x12 / 90 - alpha * x23 / 45) / 2 - 45 * alpha * y12 * (-alpha * y12 / 45 + alpha * y23 / 90) / 2 + 45 * alpha * y23 * (-alpha * y12 / 90 + alpha * y23 / 45) / 2;
     Ml(5, 6) = -3 * alpha * y12 / 2 + 3 * alpha * y23;
     Ml(5, 7) = 3 * alpha * x12 / 2 - 3 * alpha * x23;
     Ml(5, 8) = 45 * alpha * x13 * (alpha * x12 / 90 - alpha * x23 / 90) / 2 + 45 * alpha * x23 * (alpha * x12 / 90 - alpha * x23 / 45) / 2 - 45 * alpha * y13 * (-alpha * y12 / 90 + alpha * y23 / 90) / 2 - 45 * alpha * y23 * (-alpha * y12 / 90 + alpha * y23 / 45) / 2;
     Ml(6, 0) = 15;
     Ml(6, 1) = 0;
     Ml(6, 2) = 3 * alpha * y12 / 2 + 3 * alpha * y13;
     Ml(6, 3) = 15;
     Ml(6, 4) = 0;
     Ml(6, 5) = -3 * alpha * y12 / 2 + 3 * alpha * y23;
     Ml(6, 6) = 30;
     Ml(6, 7) = 0;
     Ml(6, 8) = -3 * alpha * y13 - 3 * alpha * y23;
     Ml(7, 0) = 0;
     Ml(7, 1) = 15;
     Ml(7, 2) = -3 * alpha * x12 / 2 - 3 * alpha * x13;
     Ml(7, 3) = 0;
     Ml(7, 4) = 15;
     Ml(7, 5) = 3 * alpha * x12 / 2 - 3 * alpha * x23;
     Ml(7, 6) = 0;
     Ml(7, 7) = 30;
     Ml(7, 8) = 3 * alpha * x13 + 3 * alpha * x23;
     Ml(8, 0) = -3 * alpha * y13 - 3 * alpha * y23 / 2;
     Ml(8, 1) = 3 * alpha * x13 + 3 * alpha * x23 / 2;
     Ml(8, 2) = -45 * alpha * x12 * (alpha * x13 / 90 + alpha * x23 / 90) / 2 - 45 * alpha * x13 * (alpha * x13 / 45 + alpha * x23 / 90) / 2 + 45 * alpha * y12 * (-alpha * y13 / 90 - alpha * y23 / 90) / 2 + 45 * alpha * y13 * (-alpha * y13 / 45 - alpha * y23 / 90) / 2;
     Ml(8, 3) = -3 * alpha * y13 / 2 - 3 * alpha * y23;
     Ml(8, 4) = 3 * alpha * x13 / 2 + 3 * alpha * x23;
     Ml(8, 5) = 45 * alpha * x12 * (alpha * x13 / 90 + alpha * x23 / 90) / 2 - 45 * alpha * x23 * (alpha * x13 / 90 + alpha * x23 / 45) / 2 - 45 * alpha * y12 * (-alpha * y13 / 90 - alpha * y23 / 90) / 2 + 45 * alpha * y23 * (-alpha * y13 / 90 - alpha * y23 / 45) / 2;
     Ml(8, 6) = -3 * alpha * y13 - 3 * alpha * y23;
     Ml(8, 7) = 3 * alpha * x13 + 3 * alpha * x23;
     Ml(8, 8) = 45 * alpha * x13 * (alpha * x13 / 45 + alpha * x23 / 90) / 2 + 45 * alpha * x23 * (alpha * x13 / 90 + alpha * x23 / 45) / 2 - 45 * alpha * y13 * (-alpha * y13 / 45 - alpha * y23 / 90) / 2 - 45 * alpha * y23 * (-alpha * y13 / 90 - alpha * y23 / 45) / 2;


     // Expand matrices to local complete stiffness
     //   DOFS   --->        ux1,  uy1, uz1, rx1, ry1, rz1  |  ux2, uy2, uz2, rx2, ry2, rz2  |  ux3, uy3, uz3, rx3, ry3, rz3
     int membrane_dofs[9] = {  0,   1,                  5,      6,   7,                 11,     12,  13,                 17 };
     //int  bending_dofs[9] = {            2,   3,   4,                     8,   9,  10,                    14,  15,  16      };

     double weight_over_180 = thickness * Area * rho / 180.0;

     int I = 0;
     int J = 0;

     for (int i = 0; i < 9; i++)
     {
          for (int j = 0; j < 9; j++)
          {
               // Expand membrane mass to local membrane (& drilling) DOFS.
               I = membrane_dofs[i];
               J = membrane_dofs[j];
               Mlocal(I, J) = weight_over_180 * Ml(i, j);
          }
     }

     // Transform K from local to global, K_glob = [TLG]^T * K * [TLG]
     // Where [T_lg] is blockwise diagonal arrangement of 6 T_lg matrices.
     //
     //           / T_lg     0     0     0 \                                 .
     // [TLG] =  |     0  T_lg     0     0  |                                .
     //          |     0     0  ...      0  |                                .
     //           \    0     0     0  T_lg / 18 x 18                         .
     //

     // Form TLG
     TLG.Zero();
     static Matrix T_gl(3, 3);
     T_gl.Zero();
     T_lg.Invert(T_gl);
     int offset = 0;

     for (int n = 0; n < 6; n++)
     {
          I = 0;

          for (int i = 0; i < 3; i++)
          {
               J = 0;

               for (int j = 0; j < 3; j++)
               {
                    TLG(I + offset, J + offset) = T_gl(i, j);
                    J++;
               }

               I++;
          }

          offset += 3;
     }

     // Compute total mass
     static Matrix Mass(18, 18);
     Mass.Zero();
     Mass.addMatrixTripleProduct(1.0, TLG, Mlocal, 1.0);

     return Mass;
}


const Matrix &ShellANDeS::getBendingMass ()
{

     double x_1 = this->xl1[0];
     double y_1 = this->xl1[1];

     double x_2 = this->xl2[0];
     double y_2 = this->xl2[1];

     double x_3 = this->xl3[0];
     double y_3 = this->xl3[1];


     static Matrix Gp(9, 9);
     static Matrix Hp(9, 9);

     Gp.Zero();
     Hp.Zero();

     Gp(0, 0) = 2;
     Gp(0, 1) = 1;
     Gp(0, 2) = 0;
     Gp(0, 3) = 1;
     Gp(0, 4) = 0;
     Gp(0, 5) = 1;
     Gp(0, 6) = 1;
     Gp(0, 7) = 0;
     Gp(0, 8) = -1;

     Gp(1, 0) = (-x_1 - x_2 + 2 * x_3) / (x_1 * y_2 - x_1 * y_3 - x_2 * y_1 + x_2 * y_3 + x_3 * y_1 - x_3 * y_2);
     Gp(1, 1) = (2 * x_1 - x_2 - x_3) / (x_1 * y_2 - x_1 * y_3 - x_2 * y_1 + x_2 * y_3 + x_3 * y_1 - x_3 * y_2);
     Gp(1, 2) = (-x_1 + 2 * x_2 - x_3) / (x_1 * y_2 - x_1 * y_3 - x_2 * y_1 + x_2 * y_3 + x_3 * y_1 - x_3 * y_2);
     Gp(1, 3) = 2 * (-x_1 - x_2 + 2 * x_3) / (x_1 * y_2 - x_1 * y_3 - x_2 * y_1 + x_2 * y_3 + x_3 * y_1 - x_3 * y_2);
     Gp(1, 4) = 0;
     Gp(1, 5) = 2 * (x_1 - 2 * x_2 + x_3) / (x_1 * y_2 - x_1 * y_3 - x_2 * y_1 + x_2 * y_3 + x_3 * y_1 - x_3 * y_2);
     Gp(1, 6) = 3 * (-x_1 - x_2 + 2 * x_3) / (x_1 * y_2 - x_1 * y_3 - x_2 * y_1 + x_2 * y_3 + x_3 * y_1 - x_3 * y_2);
     Gp(1, 7) = 0;
     Gp(1, 8) = 3 * (-x_1 + 2 * x_2 - x_3) / (x_1 * y_2 - x_1 * y_3 - x_2 * y_1 + x_2 * y_3 + x_3 * y_1 - x_3 * y_2);

     Gp(2, 0) = (-y_1 - y_2 + 2 * y_3) / (x_1 * y_2 - x_1 * y_3 - x_2 * y_1 + x_2 * y_3 + x_3 * y_1 - x_3 * y_2);
     Gp(2, 1) = (2 * y_1 - y_2 - y_3) / (x_1 * y_2 - x_1 * y_3 - x_2 * y_1 + x_2 * y_3 + x_3 * y_1 - x_3 * y_2);
     Gp(2, 2) = (-y_1 + 2 * y_2 - y_3) / (x_1 * y_2 - x_1 * y_3 - x_2 * y_1 + x_2 * y_3 + x_3 * y_1 - x_3 * y_2);
     Gp(2, 3) = 2 * (-y_1 - y_2 + 2 * y_3) / (x_1 * y_2 - x_1 * y_3 - x_2 * y_1 + x_2 * y_3 + x_3 * y_1 - x_3 * y_2);
     Gp(2, 4) = 0;
     Gp(2, 5) = 2 * (y_1 - 2 * y_2 + y_3) / (x_1 * y_2 - x_1 * y_3 - x_2 * y_1 + x_2 * y_3 + x_3 * y_1 - x_3 * y_2);
     Gp(2, 6) = 3 * (-y_1 - y_2 + 2 * y_3) / (x_1 * y_2 - x_1 * y_3 - x_2 * y_1 + x_2 * y_3 + x_3 * y_1 - x_3 * y_2);
     Gp(2, 7) = 0;
     Gp(2, 8) = 3 * (-y_1 + 2 * y_2 - y_3) / (x_1 * y_2 - x_1 * y_3 - x_2 * y_1 + x_2 * y_3 + x_3 * y_1 - x_3 * y_2);

     Gp(3, 0) = 0;
     Gp(3, 1) = 2;
     Gp(3, 2) = 1;
     Gp(3, 3) = 1;
     Gp(3, 4) = 1;
     Gp(3, 5) = 0;
     Gp(3, 6) = -1;
     Gp(3, 7) = 1;
     Gp(3, 8) = 0;

     Gp(4, 0) = (-x_1 - x_2 + 2 * x_3) / (x_1 * y_2 - x_1 * y_3 - x_2 * y_1 + x_2 * y_3 + x_3 * y_1 - x_3 * y_2);
     Gp(4, 1) = (2 * x_1 - x_2 - x_3) / (x_1 * y_2 - x_1 * y_3 - x_2 * y_1 + x_2 * y_3 + x_3 * y_1 - x_3 * y_2);
     Gp(4, 2) = (-x_1 + 2 * x_2 - x_3) / (x_1 * y_2 - x_1 * y_3 - x_2 * y_1 + x_2 * y_3 + x_3 * y_1 - x_3 * y_2);
     Gp(4, 3) = 2 * (x_1 + x_2 - 2 * x_3) / (x_1 * y_2 - x_1 * y_3 - x_2 * y_1 + x_2 * y_3 + x_3 * y_1 - x_3 * y_2);
     Gp(4, 4) = 2 * (2 * x_1 - x_2 - x_3) / (x_1 * y_2 - x_1 * y_3 - x_2 * y_1 + x_2 * y_3 + x_3 * y_1 - x_3 * y_2);
     Gp(4, 5) = 0;
     Gp(4, 6) = 3 * (-x_1 - x_2 + 2 * x_3) / (x_1 * y_2 - x_1 * y_3 - x_2 * y_1 + x_2 * y_3 + x_3 * y_1 - x_3 * y_2);
     Gp(4, 7) = 3 * (2 * x_1 - x_2 - x_3) / (x_1 * y_2 - x_1 * y_3 - x_2 * y_1 + x_2 * y_3 + x_3 * y_1 - x_3 * y_2);
     Gp(4, 8) = 0;

     Gp(5, 0) = (-y_1 - y_2 + 2 * y_3) / (x_1 * y_2 - x_1 * y_3 - x_2 * y_1 + x_2 * y_3 + x_3 * y_1 - x_3 * y_2);
     Gp(5, 1) = (2 * y_1 - y_2 - y_3) / (x_1 * y_2 - x_1 * y_3 - x_2 * y_1 + x_2 * y_3 + x_3 * y_1 - x_3 * y_2);
     Gp(5, 2) = (-y_1 + 2 * y_2 - y_3) / (x_1 * y_2 - x_1 * y_3 - x_2 * y_1 + x_2 * y_3 + x_3 * y_1 - x_3 * y_2);
     Gp(5, 3) = 2 * (y_1 + y_2 - 2 * y_3) / (x_1 * y_2 - x_1 * y_3 - x_2 * y_1 + x_2 * y_3 + x_3 * y_1 - x_3 * y_2);
     Gp(5, 4) = 2 * (2 * y_1 - y_2 - y_3) / (x_1 * y_2 - x_1 * y_3 - x_2 * y_1 + x_2 * y_3 + x_3 * y_1 - x_3 * y_2);
     Gp(5, 5) = 0;
     Gp(5, 6) = 3 * (-y_1 - y_2 + 2 * y_3) / (x_1 * y_2 - x_1 * y_3 - x_2 * y_1 + x_2 * y_3 + x_3 * y_1 - x_3 * y_2);
     Gp(5, 7) = 3 * (2 * y_1 - y_2 - y_3) / (x_1 * y_2 - x_1 * y_3 - x_2 * y_1 + x_2 * y_3 + x_3 * y_1 - x_3 * y_2);
     Gp(5, 8) = 0;

     Gp(6, 0) = 1;
     Gp(6, 1) = 0;
     Gp(6, 2) = 2;
     Gp(6, 3) = 0;
     Gp(6, 4) = 1;
     Gp(6, 5) = 1;
     Gp(6, 6) = 0;
     Gp(6, 7) = -1;
     Gp(6, 8) = 1;

     Gp(7, 0) = (-x_1 - x_2 + 2 * x_3) / (x_1 * y_2 - x_1 * y_3 - x_2 * y_1 + x_2 * y_3 + x_3 * y_1 - x_3 * y_2);
     Gp(7, 1) = (2 * x_1 - x_2 - x_3) / (x_1 * y_2 - x_1 * y_3 - x_2 * y_1 + x_2 * y_3 + x_3 * y_1 - x_3 * y_2);
     Gp(7, 2) = (-x_1 + 2 * x_2 - x_3) / (x_1 * y_2 - x_1 * y_3 - x_2 * y_1 + x_2 * y_3 + x_3 * y_1 - x_3 * y_2);
     Gp(7, 3) = 0;
     Gp(7, 4) = 2 * (-2 * x_1 + x_2 + x_3) / (x_1 * y_2 - x_1 * y_3 - x_2 * y_1 + x_2 * y_3 + x_3 * y_1 - x_3 * y_2);
     Gp(7, 5) = 2 * (-x_1 + 2 * x_2 - x_3) / (x_1 * y_2 - x_1 * y_3 - x_2 * y_1 + x_2 * y_3 + x_3 * y_1 - x_3 * y_2);
     Gp(7, 6) = 0;
     Gp(7, 7) = 3 * (2 * x_1 - x_2 - x_3) / (x_1 * y_2 - x_1 * y_3 - x_2 * y_1 + x_2 * y_3 + x_3 * y_1 - x_3 * y_2);
     Gp(7, 8) = 3 * (-x_1 + 2 * x_2 - x_3) / (x_1 * y_2 - x_1 * y_3 - x_2 * y_1 + x_2 * y_3 + x_3 * y_1 - x_3 * y_2);

     Gp(8, 0) = (-y_1 - y_2 + 2 * y_3) / (x_1 * y_2 - x_1 * y_3 - x_2 * y_1 + x_2 * y_3 + x_3 * y_1 - x_3 * y_2);
     Gp(8, 1) = (2 * y_1 - y_2 - y_3) / (x_1 * y_2 - x_1 * y_3 - x_2 * y_1 + x_2 * y_3 + x_3 * y_1 - x_3 * y_2);
     Gp(8, 2) = (-y_1 + 2 * y_2 - y_3) / (x_1 * y_2 - x_1 * y_3 - x_2 * y_1 + x_2 * y_3 + x_3 * y_1 - x_3 * y_2);
     Gp(8, 3) = 0;
     Gp(8, 4) = 2 * (-2 * y_1 + y_2 + y_3) / (x_1 * y_2 - x_1 * y_3 - x_2 * y_1 + x_2 * y_3 + x_3 * y_1 - x_3 * y_2);
     Gp(8, 5) = 2 * (-y_1 + 2 * y_2 - y_3) / (x_1 * y_2 - x_1 * y_3 - x_2 * y_1 + x_2 * y_3 + x_3 * y_1 - x_3 * y_2);
     Gp(8, 6) = 0;
     Gp(8, 7) = 3 * (2 * y_1 - y_2 - y_3) / (x_1 * y_2 - x_1 * y_3 - x_2 * y_1 + x_2 * y_3 + x_3 * y_1 - x_3 * y_2);
     Gp(8, 8) = 3 * (-y_1 + 2 * y_2 - y_3) / (x_1 * y_2 - x_1 * y_3 - x_2 * y_1 + x_2 * y_3 + x_3 * y_1 - x_3 * y_2);

     Gp.Invert(Hp);

     // Get local matrix
     static Matrix Ml(9, 9);
     Ml.Zero();
     Ml.addMatrixTripleProduct(1.0, Hp, Mq, rho * Area * thickness);


     //Expand local matrix to local complete dofs.
     static Matrix Mlocal(18, 18);
     Mlocal.Zero();

     // Expand matrices to local complete mass
     //   DOFS   --->          ux1, uy1, uz1, rx1, ry1, rz1  |  ux2, uy2, uz2, rx2, ry2, rz2  |  ux3, uy3, uz3, rx3, ry3, rz3
     //int membrane_dofs[9] = {  0,   1,                  5,      6,   7,                 11,     12,  13,                 17 };
     int  bending_dofs[9] = {            2,   3,   4,                     8,   9,  10,                    14,  15,  16      };
     int I, J;

     for (int i = 0; i < 9; i++)
     {
          for (int j = 0; j < 9; j++)
          {
               // Expand bending mass to local bending DOFS.
               I = bending_dofs[i];
               J = bending_dofs[j];
               Mlocal(I, J) = Ml(i, j);
          }
     }

     // Form TLG Transform local -> global
     static Matrix TLG(18, 18);
     TLG.Zero();
     static Matrix T_gl(3, 3);
     T_gl.Zero();
     T_lg.Invert(T_gl);
     int offset = 0;

     for (int n = 0; n < 6; n++)
     {
          I = 0;

          for (int i = 0; i < 3; i++)
          {
               J = 0;

               for (int j = 0; j < 3; j++)
               {
                    TLG(I + offset, J + offset) = T_gl(i, j);
                    J++;
               }

               I++;
          }

          offset += 3;
     }

     // Compute total stiffness as M = TLG^T * Mlocal *TLG
     static Matrix Mass(18, 18);
     Mass.Zero();
     Mass.addMatrixTripleProduct(1.0, TLG, Mlocal, 1.0);

     return Mass;
}



const Matrix &ShellANDeS::getMembraneTangentStiffness ()
{

     int I, J;

     // Form completestiffness
     static Matrix Kb_membrane(9, 9);       // Basic membrane stiffness
     static Matrix Kh_membrane(9, 9);       // High-order membrane stiffness
     static Matrix Klocal(18, 18);          // Total (membrane + bending) stiffness in local coordinates
     static Matrix TLG(18, 18);             // Local-to-global transformation matrix

     Kb_membrane.Zero();
     Kh_membrane.Zero();
     Klocal.Zero();
     TLG.Zero();

     //Get matrices
     Kb_membrane = getMembraneBasicStiffness();
     Kh_membrane = getMembraneHighOrderStiffness();
     //Kb_bending  = getBendingBasicStiffness();
     //Kh_bending  = getBendingHighOrderStiffness();

     // Expand matrices to local complete stiffness
     //   DOFS   --->        ux1, uy1, uz1, rx1, ry1, rz1  |  ux2, uy2, uz2, rx2, ry2, rz2  |  ux3, uy3, uz3, rx3, ry3, rz3
     int membrane_dofs[9] = {  0,   1,                  5,      6,   7,                 11,     12,  13,                 17 };
     //int  bending_dofs[9] = {            2,   3,   4,                     8,   9,  10,                    14,  15,  16      };

     for (int i = 0; i < 9; i++)
     {
          for (int j = 0; j < 9; j++)
          {
               // Expand membrane stiffness to local membrane (& drilling) DOFS.
               I = membrane_dofs[i];
               J = membrane_dofs[j];
               Klocal(I, J) = Kb_membrane(i, j) + Kh_membrane(i, j);

               // Expand bending stiffness to local bending DOFS.
               //I = bending_dofs[i];
               //J = bending_dofs[j];
               //Klocal(I,J) = Kb_bending(i,j)/Area + Area*Kh_bending(i,j);
          }
     }

     // Transform K from local to global, K_glob = [TLG]^T * K * [TLG]
     // Where [T_lg] is blockwise diagonal arrangement of 6 T_lg matrices.
     //
     //           / T_lg     0     0     0 \                               .
     // [TLG] =  |     0  T_lg     0     0  |                              .
     //          |     0     0  ...      0  |                              .
     //           \    0     0     0  T_lg / 18x18                         .
     //

     // Form TLG
     TLG.Zero();
     static Matrix T_gl(3, 3);
     T_gl.Zero();
     T_lg.Invert(T_gl);
     int offset = 0;

     for (int n = 0; n < 6; n++)
     {
          I = 0;

          for (int i = 0; i < 3; i++)
          {
               J = 0;

               for (int j = 0; j < 3; j++)
               {
                    TLG(I + offset, J + offset) = T_gl(i, j);
                    J++;
               }

               I++;
          }

          offset += 3;
     }

     // Compute total stiffness as K = TLG^T * Klocal *TLG
     K.Zero();

     K.addMatrixTripleProduct(1.0, TLG, Klocal, 1.0);

     return K;
}


const Matrix &ShellANDeS::getBendingTangentStiffness ()
{
     int I, J;

     // Form completestiffness
     static Matrix Kb_bending(9, 9);        // Basic bending stiffness
     static Matrix Kh_bending(9, 9);        // High-order bending stiffness
     static Matrix Klocal(18, 18);          // Total (membrane + bending) stiffness in local coordinates
     static Matrix TLG(18, 18);             // Local-to-global transformation matrix

     Kb_bending.Zero();
     Kh_bending.Zero();
     Klocal.Zero();
     TLG.Zero();

     //Get matrices
     Kb_bending  = getBendingBasicStiffness();
     Kh_bending  = getBendingHighOrderStiffness();
     // Expand matrices to local complete stiffness
     //   DOFS   --->        ux1, uy1, uz1, rx1, ry1, rz1  |  ux2, uy2, uz2, rx2, ry2, rz2  |  ux3, uy3, uz3, rx3, ry3, rz3
     //int membrane_dofs[9] = {  0,   1,                  5,      6,   7,                 11,     12,  13,                 17 };
     int  bending_dofs[9] = {            2,   3,   4,                     8,   9,  10,                    14,  15,  16      };

     for (int i = 0; i < 9; i++)
     {
          for (int j = 0; j < 9; j++)
          {
               // Expand bending stiffness to local bending DOFS.
               I = bending_dofs[i];
               J = bending_dofs[j];
               Klocal(I, J) = Kb_bending(i, j) + Kh_bending(i, j);
          }
     }

     // Transform K from local to global, K_glob = [TLG]^T * K * [TLG]
     // Where [T_lg] is blockwise diagonal arrangement of 6 T_lg matrices.
     //
     //           / T_lg     0     0     0 \                               .
     // [TLG] =  |     0  T_lg     0     0  |                              .
     //          |     0     0  ...      0  |                              .
     //           \    0     0     0  T_lg / 18x18                         .
     //

     // Form TLG
     TLG.Zero();
     static Matrix T_gl(3, 3);
     T_gl.Zero();
     T_lg.Invert(T_gl);
     int offset = 0;

     for (int n = 0; n < 6; n++)
     {
          I = 0;

          for (int i = 0; i < 3; i++)
          {
               J = 0;

               for (int j = 0; j < 3; j++)
               {
                    TLG(I + offset, J + offset) = T_gl(i, j);
                    J++;
               }

               I++;
          }

          offset += 3;
     }

     // Compute total stiffness as K = TLG^T * Klocal *TLG
     K.Zero();
     K.addMatrixTripleProduct(1.0, TLG, Klocal, 1.0);

     return K;
}


//===================================================================================
// Bending
//===================================================================================

Matrix ShellANDeS::getBendingBasicStiffness()
{
     //    // Original implementation
     //    Matrix L(3,9), D(3,3), Kb(9,9);
     //
     //    L = getBendingForceLumpingMatrix();
     //    D = E_planestress*(thickness*thickness*thickness/12);
     //
     //    Kb.Zero();
     //    Kb.addMatrixTripleProduct(1.0, L, D, 1.0);
     //
     //    return Kb / Area;
     //


     // =====================================================================================
     // Other implementation (ported from code by Sasa Stosic)
     // =====================================================================================

     // Jose: Original comments in Serbian were translated with Google Translate for pseudoclarity.
     //
     //
     //
     // Racuna i nalazi lokalnu matricu 9x9 u kojoj su pomeranja svih 9
     //
     // matrica C je matrica krutosti preseka sa osu X kao osu koordinatnog sistema
     //
     // C je takva da je treca vrsta pomnozena sa 2 ( ovo je znacajno jer je kod mene ova vrsta
     // deljena sa 2 jer vezujem sigma-epsilon a ovde ocekuju sigma-gamma
     //    -------------------------------------------------------------------
     //    Izrazi se mogu naci u matrci (56) na strani 233 rada na osnovu koga su izvedeni
     //    ---------------------------------------------------------------------
     //
     // [[ Google translate:  Account and finds the local 9x9 matrix in which the displacement of 9
     //
     // Matrix stiffness matrix C is the intersection with the X axis as the axis of the coordinate system
     //
     // C is such that a third species multiplied by 2 (this is important to me because this kind of
     // Shared by 2 because bind sigma-epsilon and here expect sigma-gamma
     //     ------------------------------------------------ -------------------
     //     The terms can be found in matrci (56) on page 233 of whom are based on the derived
     //     -------------------------------------------------- -------------------  ]]
     //
     //
     // Ovo bi trebalo da bude bolje ali se kod primera PLATE_10x10_B javlja
     // glupost ( nisam uspeo da nadjem zasto ) pa se vracam na linearno koje
     // do sad nije davalo gluposti,
     //
     //
     // [[ Google translate:   This is supposed to be better but in the example occurs PLATE_10x10_B
     // Stupid (I could not find out why) and returns to the linear
     // So far gave no nonsense, ]]



     // Geometric data
     double x21 = -x12;
     double y21 = -y12;
     double x32 = -x23;
     double y32 = -y23;
     double x13 = -x31;
     double y13 = -y31;

     double l12 = sqrt(x12 * x12 + y12 * y12);
     double l23 = sqrt(x23 * x23 + y23 * y23);
     double l31 = sqrt(x31 * x31 + y31 * y31);

     double c12 = x21 / l12;
     double c23 = x32 / l23;
     double c31 = x13 / l31;
     double s12 = y21 / l12;
     double s23 = y32 / l23;
     double s31 = y13 / l31;


     static Matrix L(3, 9);
     static Matrix Kb(9, 9);

     L.Zero();
     Kb.Zero();

     L(0, 0) = -c12 * s12 + c31 * s31;
     L(1, 0) = -c31 * s31 + c12 * s12;
     L(2, 0) = (s31 * s31 - c31 * c31) - (s12 * s12 - c12 * c12);

     L(0, 1) = (s12 * s12 * x12 + s31 * s31 * x31) / 2.0;
     L(1, 1) = (c12 * c12 * x12 + c31 * c31 * x31) / 2.0;
     L(2, 1) = c12 * c12 * y21 + c31 * c31 * y13;

     L(0, 2) = (s12 * s12 * y12 + s31 * s31 * y31) / 2.0;
     L(1, 2) = (c12 * c12 * y12 + c31 * c31 * y31) / 2.0;
     L(2, 2) = s12 * s12 * x21 + s31 * s31 * x13;
     //--------------------------------------------------
     L(0, 3) = -c23 * s23 + c12 * s12;
     L(1, 3) = -c12 * s12 + c23 * s23;
     L(2, 3) = (s12 * s12 - c12 * c12) - (s23 * s23 - c23 * c23);

     L(0, 4) = (s12 * s12 * x12 + s23 * s23 * x23) / 2.0;
     L(1, 4) = (c12 * c12 * x12 + c23 * c23 * x23) / 2.0;
     L(2, 4) = c12 * c12 * y21 + c23 * c23 * y32;

     L(0, 5) = (s12 * s12 * y12 + s23 * s23 * y23) / 2.0;
     //  L(1,5) = (c12*c12*y12 + c23*c23*y32)/2.0;
     L(1, 5) = (c12 * c12 * y12 + c23 * c23 * y23) / 2.0; //- Izmena 11-05-2002 [[ Google translate: changed 11-05-2002]]
     L(2, 5) = s12 * s12 * x21 + s23 * s23 * x32;
     //--------------------------------------------------
     L(0, 6) = -c31 * s31 + c23 * s23;
     L(1, 6) = -c23 * s23 + c31 * s31;
     L(2, 6) = (s23 * s23 - c23 * c23) - (s31 * s31 - c31 * c31);

     L(0, 7) = (s23 * s23 * x23 + s31 * s31 * x31) / 2.0;
     L(1, 7) = (c23 * c23 * x23 + c31 * c31 * x31) / 2.0;
     // Izmena 16-jan-2003 ovo je bilo kod mene :
     // L(7,2,7) = c23*c23*y32 + c31*c31*y31; // L(8,3) = c23*c23*y23 + c31*c31*y13;
     //         // Bilo je pa nije valjalo i treba y32
     // A u Radimpex-u je bilo:
     L(2, 7) = c23 * c23 * y32 + c31 * c31 * y13; //# odnosno permutovan je index y31 u y13
     // Treba proveriti navedeni primer PLATE_10x10_B jer je kod mene radilo
     // A sad je postavljeno kao u RadImpexu.

     L(0, 8) = (s23 * s23 * y23 + s31 * s31 * y31) / 2.0;
     L(1, 8) = (c23 * c23 * y23 + c31 * c31 * y31) / 2.0;
     L(2, 8) = s23 * s23 * x32 + s31 * s31 * x13;

     // U odnosu na izraze izvedene u Mathematici u file-u LMatrica.nb, svi clanovi
     // u vrstama 2,3,5,6,8 i 9 imaju negativne vrednosti odnosno treba proveriti
     // konvenciju o znaku momenata savijanja.
     //
     //  Lt = L.transpose();
     //
     // Modifikuje matricu C:
     // C(1,3) = C(1,3)/2 ; C(2,3) = C(2,3)/2 ; C(3,3)  =  C(3,3)/2;
     // To sad radi rutina koja poziva
     //
     //
     //    C1 = C *(1./Ar) ;
     //    Kb = C1*Lt;
     //    Kb = L*Kb;
     // cout<<"Kb = " <<Kb.OutString();
     // cout<<"det Kb = ";
     // Kb.eigenvalues().print();


     Kb.Zero();
     calculate_E_planestress_and_beta0();

     Kb.addMatrixTripleProduct(1.0, L, E_planestress * thickness * thickness * thickness / 12 / Area, 1.0);

     return Kb;
}

Matrix ShellANDeS::getBendingHighOrderStiffness()
{
     static Matrix Kh(9, 9);
     Kh.Zero();

     // Geometric data
     double x21 = -x12;
     double y21 = -y12;
     double x32 = -x23;
     double y32 = -y23;
     double x13 = -x31;
     double y13 = -y31;

     double l12 = sqrt(x12 * x12 + y12 * y12);
     double l23 = sqrt(x23 * x23 + y23 * y23);
     double l31 = sqrt(x31 * x31 + y31 * y31);
     double k_A = 4.0 * Area * Area;

     static Matrix T(3, 3);
     T.Zero();

     T(0, 0) = y23 * y13 / k_A;
     T(0, 1) = y31 * y21 / k_A;
     T(0, 2) = y12 * y32 / k_A;
     T(1, 0) = x23 * x13 / k_A;
     T(1, 1) = x31 * x21 / k_A;
     T(1, 2) = x12 * x32 / k_A;
     T(2, 0) = (y23 * x31 + x32 * y13) / k_A;
     T(2, 1) = (y31 * x12 + x13 * y21) / k_A;
     T(2, 2) = (y12 * x23 + x21 * y32) / k_A;

     static Matrix D(3, 3);
     D.Zero();

     D.addMatrixTripleProduct(1.0, T, E_planestress * thickness * thickness * thickness / 12, 1.0);

     // D = T.transpose()*(C*T);
     // D = T.inverse()*(C*T);
     // Ovo posle kad proradi prebaci preko promenjivih da brze radi
     // [[ Google translate: This starts to work after you get over the variables that work fast ]]

     double La12 = (x12 * x13 + y12 * y13) / pow(l12, 2);
     double La23 = (x23 * x21 + y23 * y21) / pow(l23, 2);
     double La31 = (x31 * x32 + y31 * y32) / pow(l31, 2);

     double B11 = 2.0 * (La12 * La12 - La12 + 1);
     double B22 = 2.0 * (La23 * La23 - La23 + 1);
     double B33 = 2.0 * (La31 * La31 - La31 + 1);

     double B12 = (2.0 - La12) * La23 - La12 - 1.0;
     double B23 = (2.0 - La23) * La31 - La23 - 1.0;
     double B13 = (2.0 - La31) * La12 - La31 - 1.0;

     double r11 = B11 * D(0, 0);
     double r12 = B12 * D(0, 1);
     double r13 = B13 * D(0, 2);
     double r22 = B22 * D(1, 1);
     double r23 = B23 * D(1, 2);
     double r33 = B33 * D(2, 2);

     //matrix_m Kh(9,9,0.0);

     Kh(0, 0) = 4.*(r33 - 2 * r13 + r11);
     Kh(0, 1) = 2.*((r11 - r13) * y21 + (r13 - r33) * y13);
     Kh(0, 2) = 2.*((r11 - r13) * x12 + (r13 - r33) * x31);
     Kh(0, 3) = 4.*( - r23 + r13 + r12 - r11);
     Kh(0, 4) = 2.*((r12 - r23) * y32 + (r11 - r13) * y21);
     Kh(0, 5) = 2.*((r12 - r23) * x23 + (r11 - r13) * x12);
     Kh(0, 6) = 4.*( - r33 + r23 + r13 - r12);
     Kh(0, 7) = 2.*((r12 - r23) * y32 + (r13 - r33) * y13);
     Kh(0, 8) = 2.*((r12 - r23) * x23 + (r13 - r33) * x31);

     Kh(1, 0) = Kh(0, 1);
     Kh(1, 1) = r11 * pow(y21, 2) + 2 * r13 * y13 * y21 + r33 * pow(y13, 2);
     Kh(1, 2) = (r11 * x21 + r13 * x13) * y12 + (r13 * x21 + r33 * x13) * y31;
     Kh(1, 3) = 2.*((r12 - r11) * y21 + (r23 - r13) * y13);
     Kh(1, 4) = (r12 * y21 + r23 * y13) * y32 + r11 * pow(y21, 2) + r13 * y13 * y21;
     Kh(1, 5) = (r12 * x32 + r11 * x21) * y12 + (r23 * x32 + r13 * x21) * y31;
     Kh(1, 6) = 2.*((r13 - r12) * y21 + (r33 - r23) * y13);
     Kh(1, 7) = (r12 * y21 + r23 * y13) * y32 + r13 * y13 * y21 + r33 * y13 * y13;
     Kh(1, 8) = (r12 * x32 + r13 * x13) * y12 + (r23 * x32 + r33 * x13) * y31; // Ovde je bila greska !!! 6 - 04 - 2002 [[ Google translate: There was a mistake! 6 - 04-2002 ]]

     Kh(2, 0) = Kh(0, 2);
     Kh(2, 1) = Kh(1, 2);
     Kh(2, 2) = r11 * pow(x21, 2) + 2.*r13 * x13 * x21 + r33 * pow(x13, 2);
     Kh(2, 3) = 2.*((r11 - r12) * x21 + (r13 - r23) * x13);
     Kh(2, 4) = (r12 * x21 + r23 * x13) * y23 + (r11 * x21 + r13 * x13) * y12;
     Kh(2, 5) = (r12 * x21 + r23 * x13) * x32 + r11 * pow(x21, 2) + r13 * x13 * x21;
     Kh(2, 6) = 2.*((r12 - r13) * x21 + (r23 - r33) * x13);
     Kh(2, 7) = (r12 * x21 + r23 * x13) * y23 + (r13 * x21 + r33 * x13) * y31;
     Kh(2, 8) = (r12 * x21 + r23 * x13) * x32 + (r13 * x21 + r33 * x13) * x13;

     Kh(3, 0) = Kh(0, 3);
     Kh(3, 1) = Kh(1, 3) ;
     Kh(3, 2) = Kh(2, 3) ;
     Kh(3, 3) = 4.0 * (r22 - 2 * r12 + r11);
     Kh(3, 4) = 2.*((r22 - r12) * y32 + (r12 - r11) * y21);
     Kh(3, 5) = 2.*((r22 - r12) * x23 + (r12 - r11) * x12);
     Kh(3, 6) = 4.0 * (r23 - r22 - r13 + r12);
     Kh(3, 7) = 2.*((r22 - r12) * y32 + (r23 - r13) * y13);
     Kh(3, 8) = 2.*((r22 - r12) * x23 + (r23 - r13) * x31);

     Kh(4, 0) = Kh(0, 4);
     Kh(4, 1) = Kh(1, 4);
     Kh(4, 2) = Kh(2, 4);
     Kh(4, 3) = Kh(3, 4);
     Kh(4, 4) = r22 * pow(y32, 2) + 2.*r12 * y21 * y32 + r11 * pow(y21, 2);
     Kh(4, 5) = (r22 * x32 + r12 * x21) * y23 + (r12 * x32 + r11 * x21) * y12;
     Kh(4, 6) = 2.*((r23 - r22) * y32 + (r13 - r12) * y21);
     Kh(4, 7) = r22 * pow(y32, 2) + (r12 * y21 + r23 * y13) * y32 + r13 * y13 * y21;
     Kh(4, 8) = (r22 * x32 + r23 * x13) * y23 + (r12 * x32 + r13 * x13) * y12;

     Kh(5, 0) = Kh(0, 5);
     Kh(5, 1) = Kh(1, 5);
     Kh(5, 2) = Kh(2, 5);
     Kh(5, 3) = Kh(3, 5);
     Kh(5, 4) = Kh(4, 5);
     Kh(5, 5) = r22 * pow(x32, 2) + 2.*r12 * x21 * x32 + r11 * pow(x21, 2);
     Kh(5, 6) = 2.*((r22 - r23) * x32 + (r12 - r13) * x21);
     Kh(5, 7) = (r22 * x32 + r12 * x21) * y23 + (r23 * x32 + r13 * x21) * y31;
     Kh(5, 8) = r22 * pow(x32, 2) + (r12 * x21 + r23 * x13) * x32 + r13 * x13 * x21;

     Kh(6, 0) = Kh(0, 6);
     Kh(6, 1) = Kh(1, 6);
     Kh(6, 2) = Kh(2, 6);
     Kh(6, 3) = Kh(3, 6);
     Kh(6, 4) = Kh(4, 6);
     Kh(6, 5) = Kh(5, 6);
     Kh(6, 6) = 4.0 * (r33 - 2 * r23 + r22);
     Kh(6, 7) = 2.*((r23 - r22) * y32 + (r33 - r23) * y13);
     Kh(6, 8) = 2.*((r23 - r22) * x23 + (r33 - r23) * x31);

     Kh(7, 0) = Kh(0, 7);
     Kh(7, 1) = Kh(1, 7);
     Kh(7, 2) = Kh(2, 7);
     Kh(7, 3) = Kh(3, 7);
     Kh(7, 4) = Kh(4, 7);
     Kh(7, 5) = Kh(5, 7);
     Kh(7, 6) = Kh(6, 7);
     Kh(7, 7) = r22 * pow(y32, 2) + 2.*r23 * y13 * y32 + r33 * pow(y13, 2);
     Kh(7, 8) = (r22 * x32 + r23 * x13) * y23 + (r23 * x32 + r33 * x13) * y31;

     Kh(8, 0) = Kh(0, 8);
     Kh(8, 1) = Kh(1, 8);
     Kh(8, 2) = Kh(2, 8);
     Kh(8, 3) = Kh(3, 8);
     Kh(8, 4) = Kh(4, 8);
     Kh(8, 5) = Kh(5, 8);
     Kh(8, 6) = Kh(6, 8);
     Kh(8, 7) = Kh(7, 8);
     Kh(8, 8) = r22 * pow(x32, 2) + 2.*r23 * x13 * x32 + r33 * pow(x13, 2);

     return Kh * Area;
}


void ShellANDeS::initializeMq()
{
     // This matrix is the same for all instances of the element.
     // Appendix A of Ref. (vii) extended to the mass matrix.
     // Generated from mass_matrix.py
     Mq(0, 0) =  1960;
     Mq(0, 1) =  1540;
     Mq(0, 2) =  1540;
     Mq(0, 3) =  280;
     Mq(0, 4) =  224;
     Mq(0, 5) =  336;
     Mq(0, 6) =  112;
     Mq(0, 7) =  -56;
     Mq(0, 8) =  -56;
     Mq(1, 0) =  1540;
     Mq(1, 1) =  1960;
     Mq(1, 2) =  1540;
     Mq(1, 3) =  336;
     Mq(1, 4) =  280;
     Mq(1, 5) =  224;
     Mq(1, 6) =  -56;
     Mq(1, 7) =  112;
     Mq(1, 8) =  -56;
     Mq(2, 0) =  1540;
     Mq(2, 1) =  1540;
     Mq(2, 2) =  1960;
     Mq(2, 3) =  224;
     Mq(2, 4) =  336;
     Mq(2, 5) =  280;
     Mq(2, 6) =  -56;
     Mq(2, 7) =  -56;
     Mq(2, 8) =  112;
     Mq(3, 0) =  280;
     Mq(3, 1) =  336;
     Mq(3, 2) =  224;
     Mq(3, 3) =  112;
     Mq(3, 4) =  56;
     Mq(3, 5) =  56;
     Mq(3, 6) =  0;
     Mq(3, 7) =  32;
     Mq(3, 8) =  -32;
     Mq(4, 0) =  224;
     Mq(4, 1) =  280;
     Mq(4, 2) =  336;
     Mq(4, 3) =  56;
     Mq(4, 4) =  112;
     Mq(4, 5) =  56;
     Mq(4, 6) =  -32;
     Mq(4, 7) =  0;
     Mq(4, 8) =  32;
     Mq(5, 0) =  336;
     Mq(5, 1) =  224;
     Mq(5, 2) =  280;
     Mq(5, 3) =  56;
     Mq(5, 4) =  56;
     Mq(5, 5) =  112;
     Mq(5, 6) =  32;
     Mq(5, 7) =  -32;
     Mq(5, 8) =  0;
     Mq(6, 0) =  112;
     Mq(6, 1) =  -56;
     Mq(6, 2) =  -56;
     Mq(6, 3) =  0;
     Mq(6, 4) =  -32;
     Mq(6, 5) =  32;
     Mq(6, 6) =  60;
     Mq(6, 7) =  -27;
     Mq(6, 8) =  -27;
     Mq(7, 0) =  -56;
     Mq(7, 1) =  112;
     Mq(7, 2) =  -56;
     Mq(7, 3) =  32;
     Mq(7, 4) =  0;
     Mq(7, 5) =  -32;
     Mq(7, 6) =  -27;
     Mq(7, 7) =  60;
     Mq(7, 8) =  -27;
     Mq(8, 0) =  -56;
     Mq(8, 1) =  -56;
     Mq(8, 2) =  112;
     Mq(8, 3) =  -32;
     Mq(8, 4) =  32;
     Mq(8, 5) =  0;
     Mq(8, 6) =  -27;
     Mq(8, 7) =  -27;
     Mq(8, 8) =  60;

     Mq = Mq / 1680.;
}


Matrix ShellANDeS::getMembraneForceLumpingMatrix()
{
     static Matrix L(3, 9); // This is the transpose of Felippa's article
     L.Zero();

     // Filling each row
     L(0, 0) = y23;
     L(2, 0) = -x23;
     L(1, 1) = -x23;
     L(2, 1) = y23;
     L(0, 2) = alpha_membrane * y23 * (-y31 + y12) / 6;
     L(1, 2) = alpha_membrane * (-x23) * (x31 - x12) / 6;
     L(2, 2) = alpha_membrane * (-x31 * y31 + x12 * y12) / 3;

     L(0, 3) = y31;
     L(2, 3) = -x31;
     L(1, 4) = -x31;
     L(2, 4) = y31;
     L(0, 5) = alpha_membrane * y31 * (-y12 + y23) / 6;
     L(1, 5) = alpha_membrane * (-x31) * (x12 - x23) / 6;
     L(2, 5) = alpha_membrane * (-x12 * y12 + x23 * y23) / 3;

     L(0, 6) = y12;
     L(2, 6) = -x12;
     L(1, 7) = -x12;
     L(2, 7) = y12;
     L(0, 8) = alpha_membrane * y12 * (-y23 + y31) / 6;
     L(1, 8) = alpha_membrane * (-x12) * (x23 - x31) / 6;
     L(2, 8) = alpha_membrane * (-x23 * y23 + x31 * y31) / 3;

     return 0.5 * thickness * L;
}

Matrix ShellANDeS::getMembraneHierarchicalProjectionMatrix()
{
     static Matrix T_thu(3, 9);
     T_thu.Zero();

     double x21 = -x12;
     double y21 = -y12;
     double x32 = -x23;
     double y32 = -y23;
     double x13 = -x31;
     double y13 = -y31;


     T_thu(0, 0) = x32 / (4 * Area);
     T_thu(1, 0) = x32 / (4 * Area);
     T_thu(2, 0) = x32 / (4 * Area);

     T_thu(0, 1) = y32 / (4 * Area);
     T_thu(1, 1) = y32 / (4 * Area);
     T_thu(2, 1) = y32 / (4 * Area);

     T_thu(0, 2) = 1;


     T_thu(0, 3) = x13 / (4 * Area);
     T_thu(1, 3) = x13 / (4 * Area);
     T_thu(2, 3) = x13 / (4 * Area);

     T_thu(0, 4) = y13 / (4 * Area);
     T_thu(1, 4) = y13 / (4 * Area);
     T_thu(2, 4) = y13 / (4 * Area);

     T_thu(1, 5) = 1;


     T_thu(0, 6) = x21 / (4 * Area);
     T_thu(1, 6) = x21 / (4 * Area);
     T_thu(2, 6) = x21 / (4 * Area);

     T_thu(0, 7) = y21 / (4 * Area);
     T_thu(1, 7) = y21 / (4 * Area);
     T_thu(2, 7) = y21 / (4 * Area);

     T_thu(2, 8) = 1;

     return T_thu;
}

Matrix ShellANDeS::getMembraneNaturalStrainProjectionMatrix()
{
     static Matrix T_e(3, 3);
     T_e.Zero();

     double x21 = -x12;
     double y21 = -y12;
     double x32 = -x23;
     double y32 = -y23;
     double x13 = -x31;
     double y13 = -y31;


     // From Felippa's code
     double ll21 =       x21 * x21 + y21 * y21;
     double ll32 =       x32 * x32 + y32 * y32;
     double ll13 =       x13 * x13 + y13 * y13;
     double tfac =       1.0 / (4.0 * Area * Area);

     T_e(0, 0) =    tfac * y23 * y13 * ll21;
     T_e(0, 1) =    tfac * y31 * y21 * ll32;
     T_e(0, 2) =    tfac * y12 * y32 * ll13;
     T_e(1, 0) =    tfac * x23 * x13 * ll21;
     T_e(1, 1) =    tfac * x31 * x21 * ll32;
     T_e(1, 2) =    tfac * x12 * x32 * ll13;
     T_e(2, 0) =    tfac * (y23 * x31 + x32 * y13) * ll21;
     T_e(2, 1) =    tfac * (y31 * x12 + x13 * y21) * ll32;
     T_e(2, 2) =    tfac * (y12 * x23 + x21 * y32) * ll13;

     return T_e;
}

Matrix ShellANDeS::getMembraneBasicStiffness()
{
     //Initialize basic stiffness
     static Matrix Kb(9, 9);
     Kb.Zero();

     //Condense E into plane stress matrix
     calculate_E_planestress_and_beta0();

     // Force lumping matrix
     static Matrix L(3, 9);
     L.Zero();
     L = getMembraneForceLumpingMatrix();

     double Vol = thickness * Area;

     Kb.Zero();
     Kb.addMatrixTripleProduct(1.0, L, E_planestress, 1.0);

     return Kb / Vol;
}

Matrix ShellANDeS::getMembraneHighOrderStiffness()
{
     static Matrix Kh(9, 9);
     static Matrix Ktheta(3, 3);
     static Matrix Q1(3, 3), Q2(3, 3), Q3(3, 3), Q4(3, 3), Q5(3, 3), Q6(3, 3);

     Kh.Zero();
     Ktheta.Zero();
     Q1.Zero();
     Q2.Zero();
     Q3.Zero();
     Q4.Zero();
     Q5.Zero();
     Q6.Zero();


     double c12 = (2 * Area) / (3 * ( x12 * x12 + y12 * y12 ));
     double c23 = (2 * Area) / (3 * ( x23 * x23 + y23 * y23 ));
     double c31 = (2 * Area) / (3 * ( x31 * x31 + y31 * y31 ));

     // Q1
     Q1(0, 0) = beta_membrane(1) * c12;
     Q1(0, 1) = beta_membrane(2) * c12;
     Q1(0, 2) = beta_membrane(3) * c12;

     Q1(1, 0) = beta_membrane(4) * c23;
     Q1(1, 1) = beta_membrane(5) * c23;
     Q1(1, 2) = beta_membrane(6) * c23;

     Q1(2, 0) = beta_membrane(7) * c31;
     Q1(2, 1) = beta_membrane(8) * c31;
     Q1(2, 2) = beta_membrane(9) * c31;

     // Q2
     Q2(0, 0) = beta_membrane(9) * c12;
     Q2(0, 1) = beta_membrane(7) * c12;
     Q2(0, 2) = beta_membrane(8) * c12;

     Q2(1, 0) = beta_membrane(3) * c23;
     Q2(1, 1) = beta_membrane(1) * c23;
     Q2(1, 2) = beta_membrane(2) * c23;

     Q2(2, 0) = beta_membrane(6) * c31;
     Q2(2, 1) = beta_membrane(4) * c31;
     Q2(2, 2) = beta_membrane(5) * c31;

     // Q3
     Q3(0, 0) = beta_membrane(5) * c12;
     Q3(0, 1) = beta_membrane(6) * c12;
     Q3(0, 2) = beta_membrane(4) * c12;

     Q3(1, 0) = beta_membrane(8) * c23;
     Q3(1, 1) = beta_membrane(9) * c23;
     Q3(1, 2) = beta_membrane(7) * c23;

     Q3(2, 0) = beta_membrane(2) * c31;
     Q3(2, 1) = beta_membrane(3) * c31;
     Q3(2, 2) = beta_membrane(1) * c31;

     // Q4 Q5 and Q6
     Q4 = 0.5 * (Q1 + Q2);
     Q5 = 0.5 * (Q2 + Q3);
     Q6 = 0.5 * (Q3 + Q1);

     //Transpose Qs


     Matrix T_e = getMembraneNaturalStrainProjectionMatrix();
     static Matrix Enat(3, 3);
     Enat.Zero();
     Enat.addMatrixTripleProduct(1.0, T_e, E_planestress, 1.0);  // (T_e)^T*E_planestress*T_e


     // Form high order matrix in Natural coordinates

     //Felippa's unrolled loop implementation
     // this is very slow in C++ because the array accesses are done through function
     // calling instead of pointers. It can be changed..... though....
     double x21 = -x12;
     double y21 = -y12;
     double x32 = -x23;
     double y32 = -y23;
     double x13 = -x31;
     double y13 = -y31;


     // This code calculates the Kh
     double  sfac =    0.75 * beta0 * Area;

     double  a4 =     Enat(0, 0) * Q4(0, 0) + Enat(0, 1) * Q4(1, 0) + Enat(0, 2) * Q4(2, 0);
     double  b4 =     Enat(1, 0) * Q4(0, 0) + Enat(1, 1) * Q4(1, 0) + Enat(1, 2) * Q4(2, 0);
     double  c4 =     Enat(2, 0) * Q4(0, 0) + Enat(2, 1) * Q4(1, 0) + Enat(2, 2) * Q4(2, 0);
     double  a5 =     Enat(0, 0) * Q5(0, 0) + Enat(0, 1) * Q5(1, 0) + Enat(0, 2) * Q5(2, 0);
     double  b5 =     Enat(1, 0) * Q5(0, 0) + Enat(1, 1) * Q5(1, 0) + Enat(1, 2) * Q5(2, 0);
     double  c5 =     Enat(2, 0) * Q5(0, 0) + Enat(2, 1) * Q5(1, 0) + Enat(2, 2) * Q5(2, 0);
     double  a6 =     Enat(0, 0) * Q6(0, 0) + Enat(0, 1) * Q6(1, 0) + Enat(0, 2) * Q6(2, 0);
     double  b6 =     Enat(1, 0) * Q6(0, 0) + Enat(1, 1) * Q6(1, 0) + Enat(1, 2) * Q6(2, 0);
     double  c6 =     Enat(2, 0) * Q6(0, 0) + Enat(2, 1) * Q6(1, 0) + Enat(2, 2) * Q6(2, 0);

     double  S11 =   (Q4(0, 0) * a4 + Q4(1, 0) * b4 + Q4(2, 0) * c4 + Q5(0, 0) * a5 + Q5(1, 0) * b5
                      + Q5(2, 0) * c5 + Q6(0, 0) * a6 + Q6(1, 0) * b6 + Q6(2, 0) * c6) * sfac;
     double  S12 =   (Q4(0, 1) * a4 + Q4(1, 1) * b4 + Q4(2, 1) * c4 + Q5(0, 1) * a5 + Q5(1, 1) * b5
                      + Q5(2, 1) * c5 + Q6(0, 1) * a6 + Q6(1, 1) * b6 + Q6(2, 1) * c6) * sfac;
     double  S13 =   (Q4(0, 2) * a4 + Q4(1, 2) * b4 + Q4(2, 2) * c4 + Q5(0, 2) * a5 + Q5(1, 2) * b5
                      + Q5(2, 2) * c5 + Q6(0, 2) * a6 + Q6(1, 2) * b6 + Q6(2, 2) * c6) * sfac;

     a4 =     Enat(0, 0) * Q4(0, 1) + Enat(0, 1) * Q4(1, 1) + Enat(0, 2) * Q4(2, 1);
     b4 =     Enat(1, 0) * Q4(0, 1) + Enat(1, 1) * Q4(1, 1) + Enat(1, 2) * Q4(2, 1);
     c4 =     Enat(2, 0) * Q4(0, 1) + Enat(2, 1) * Q4(1, 1) + Enat(2, 2) * Q4(2, 1);
     a5 =     Enat(0, 0) * Q5(0, 1) + Enat(0, 1) * Q5(1, 1) + Enat(0, 2) * Q5(2, 1);
     b5 =     Enat(1, 0) * Q5(0, 1) + Enat(1, 1) * Q5(1, 1) + Enat(1, 2) * Q5(2, 1);
     c5 =     Enat(2, 0) * Q5(0, 1) + Enat(2, 1) * Q5(1, 1) + Enat(2, 2) * Q5(2, 1);
     a6 =     Enat(0, 0) * Q6(0, 1) + Enat(0, 1) * Q6(1, 1) + Enat(0, 2) * Q6(2, 1);
     b6 =     Enat(1, 0) * Q6(0, 1) + Enat(1, 1) * Q6(1, 1) + Enat(1, 2) * Q6(2, 1);
     c6 =     Enat(2, 0) * Q6(0, 1) + Enat(2, 1) * Q6(1, 1) + Enat(2, 2) * Q6(2, 1);

     double  S22 =   (Q4(0, 1) * a4 + Q4(1, 1) * b4 + Q4(2, 1) * c4 + Q5(0, 1) * a5 + Q5(1, 1) * b5
                      + Q5(2, 1) * c5 + Q6(0, 1) * a6 + Q6(1, 1) * b6 + Q6(2, 1) * c6) * sfac;
     double   S23 =   (Q4(0, 2) * a4 + Q4(1, 2) * b4 + Q4(2, 2) * c4 + Q5(0, 2) * a5 + Q5(1, 2) * b5
                       + Q5(2, 2) * c5 + Q6(0, 2) * a6 + Q6(1, 2) * b6 + Q6(2, 2) * c6) * sfac;

     a4 =     Enat(0, 0) * Q4(0, 2) + Enat(0, 1) * Q4(1, 2) + Enat(0, 2) * Q4(2, 2);
     b4 =     Enat(1, 0) * Q4(0, 2) + Enat(1, 1) * Q4(1, 2) + Enat(1, 2) * Q4(2, 2);
     c4 =     Enat(2, 0) * Q4(0, 2) + Enat(2, 1) * Q4(1, 2) + Enat(2, 2) * Q4(2, 2);
     a5 =     Enat(0, 0) * Q5(0, 2) + Enat(0, 1) * Q5(1, 2) + Enat(0, 2) * Q5(2, 2);
     b5 =     Enat(1, 0) * Q5(0, 2) + Enat(1, 1) * Q5(1, 2) + Enat(1, 2) * Q5(2, 2);
     c5 =     Enat(2, 0) * Q5(0, 2) + Enat(2, 1) * Q5(1, 2) + Enat(2, 2) * Q5(2, 2);
     a6 =     Enat(0, 0) * Q6(0, 2) + Enat(0, 1) * Q6(1, 2) + Enat(0, 2) * Q6(2, 2);
     b6 =     Enat(1, 0) * Q6(0, 2) + Enat(1, 1) * Q6(1, 2) + Enat(1, 2) * Q6(2, 2);
     c6 =     Enat(2, 0) * Q6(0, 2) + Enat(2, 1) * Q6(1, 2) + Enat(2, 2) * Q6(2, 2);

     double  S33 =   (Q4(0, 2) * a4 + Q4(1, 2) * b4 + Q4(2, 2) * c4 + Q5(0, 2) * a5 + Q5(1, 2) * b5
                      + Q5(2, 2) * c5 + Q6(0, 2) * a6 + Q6(1, 2) * b6 + Q6(2, 2) * c6) * sfac;

     double  Ssum1 =    (S11 + S12 + S13) / (4 * Area);
     double  Ssum2 =    (S12 + S22 + S23) / (4 * Area);
     double  Ssum3 =    (S13 + S23 + S33) / (4 * Area);
     double  Ssum123 =  (Ssum1 + Ssum2 + Ssum3) / (4 * Area);

     Kh(0, 0) =  Ssum123 * x32 * x32;
     Kh(0, 1) =  Ssum123 * x32 * y32;
     Kh(0, 2) =  Ssum1 * x32;
     Kh(0, 3) =  Ssum123 * x13 * x32;
     Kh(0, 4) =  Ssum123 * x32 * y13;
     Kh(0, 5) =  Ssum2 * x32;
     Kh(0, 6) =  Ssum123 * x21 * x32;
     Kh(0, 7) =  Ssum123 * x32 * y21;
     Kh(0, 8) =  Ssum3 * x32;
     Kh(1, 1) =  Ssum123 * y32 * y32;
     Kh(1, 2) =  Ssum1 * y32;
     Kh(1, 3) =  Ssum123 * x13 * y32;
     Kh(1, 4) =  Ssum123 * y13 * y32;
     Kh(1, 5) =  Ssum2 * y32;
     Kh(1, 6) =  Ssum123 * x21 * y32;
     Kh(1, 7) =  Ssum123 * y21 * y32;
     Kh(1, 8) =  Ssum3 * y32;
     Kh(2, 2) =  S11;
     Kh(2, 3) =  Ssum1 * x13;
     Kh(2, 4) =  Ssum1 * y13;
     Kh(2, 5) =  S12;
     Kh(2, 6) =  Ssum1 * x21;
     Kh(2, 7) =  Ssum1 * y21;
     Kh(2, 8) =  S13;
     Kh(3, 3) =  Ssum123 * x13 * x13;
     Kh(3, 4) =  Ssum123 * x13 * y13;
     Kh(3, 5) =  Ssum2 * x13;
     Kh(3, 6) =  Ssum123 * x13 * x21;
     Kh(3, 7) =  Ssum123 * x13 * y21;
     Kh(3, 8) =  Ssum3 * x13;
     Kh(4, 4) =  Ssum123 * y13 * y13;
     Kh(4, 5) =  Ssum2 * y13;
     Kh(4, 6) =  Ssum123 * x21 * y13;
     Kh(4, 7) =  Ssum123 * y13 * y21;
     Kh(4, 8) =  Ssum3 * y13;
     Kh(5, 5) =  S22;
     Kh(5, 6) =  Ssum2 * x21;
     Kh(5, 7) =  Ssum2 * y21;
     Kh(5, 8) =  S23;
     Kh(6, 6) =  Ssum123 * x21 * x21;
     Kh(6, 7) =  Ssum123 * x21 * y21;
     Kh(6, 8) =  Ssum3 * x21;
     Kh(7, 7) =  Ssum123 * y21 * y21         ;
     Kh(7, 8) =  Ssum3 * y21;
     Kh(8, 8) =  S33;

     // Fill in the rest of the matrix by symmetry
     for (int i = 1; i < 9; i++ )
     {
          for (int j = 0; j < i; j++ )
          {
               Kh(i, j) = Kh(j, i);
          }
     }

     return Kh * thickness;
}


void ShellANDeS::calculate_E_planestress_and_beta0()
{
     // This condenses the E matrix which holds the elastic coefficients
     // into the plane-stress E matrix.
     // Matrix Etang(6, 6);
     // // E = (*theMaterial)->getTangent();

     // //Static condensation
     // Matrix A(3, 3);
     // Matrix B(3, 3);
     // Matrix C(3, 3);

     // int keep_these_dofs[3] = {0, 1, 3};
     // int condensed_these_dofs[3] = {2, 4, 5};

     // for (int i = 0; i < 3; i++)
     // {
     //     for (int j = 0; j < 3; j++)
     //     {
     //         A(i, j) = Etang(keep_these_dofs[i], keep_these_dofs[j]);
     //         B(i, j) = Etang(keep_these_dofs[i], condensed_these_dofs[j]);
     //         C(i, j) = Etang(condensed_these_dofs[i], condensed_these_dofs[j]);
     //     }
     // }

     // Matrix C_inv(3, 3);
     // C.Invert(C_inv); // Invert C and place result in C_inv

     // E_planestress.Zero();

     // // Compute condensed matrix
     // // E_planestress = A - B*(C^-1)*B^T
     // E_planestress = A;

     // for (int i = 0; i < 3; i++)
     // {
     //     for (int j = 0; j < 3; j++)
     //     {
     //         for (int m = 0; m < 3; m ++)
     //         {
     //             for (int n = 0; n < 3; n ++)
     //             {
     //                 E_planestress(i, j) -= B(i, m) * C_inv(m, n) * B(j, n);
     //             }
     //         }
     //     }
     // }

     // Calculate beta0
     // From Carlos Felippa's fortran code:

     // double M  = E / ( 1.0 - nu*nu ) ; //membrane modulus
     // double G  =  0.5 * E / ( 1.0 + nu ) ; //shear modulus

     // // G *= thickness ;  //multiply by thickness
     // // M *= thickness ;

     // E_planestress.Zero() ;

     E_planestress(0, 0) = mE11;
     E_planestress(1, 1) = mE22;
     E_planestress(2, 2) = mE33;

     E_planestress(1, 0) = mE12;
     E_planestress(0, 1) = mE12;

     E_planestress(2, 0) = mE13;
     E_planestress(0, 2) = mE13;

     E_planestress(2, 1) = mE23;
     E_planestress(1, 2) = mE23;


     double E11 =  mE11;  //E_planestress(0, 0);
     double E22 =  mE22;  //E_planestress(1, 1);
     double E33 =  mE33;  //E_planestress(2, 2);
     double E12 =  mE12;  //E_planestress(0, 1);
     double E13 =  mE13;  //E_planestress(0, 2);
     double E23 =  mE23;  //E_planestress(1, 2);

     double Edet = E11 * E22 * E33 + 2 * E12 * E13 * E23 - E11 * pow(E23, 2) - E22 * pow(E13, 2) - E33 * pow(E12, 2);

     double E11C11 = (-5 * E11 * pow(E12, 2) - 6 * pow(E12, 3) - 3 * E11 * pow(E13, 2) + 14 * E12 * pow(E13, 2) +
                      5 * pow(E11, 2) * E22 + 6 * E11 * E12 * E22 - 5 * pow(E12, 2) * E22 - 75 * pow(E13, 2) * E22 +
                      5 * E11 * pow(E22, 2) - 14 * E11 * E13 * E23 + 92 * E12 * E13 * E23 - 14 * E13 * E22 * E23 -
                      75 * E11 * pow(E23, 2) + 14 * E12 * pow(E23 , 2) - 3 * E22 * pow(E23, 2) + (3 * pow(E11, 2) + 82 * E11 * E22 +
                                3 * pow(E22, 2) - 4 * (6 * pow(E12, 2) + 5 * pow(E13, 2) - 6 * E13 * E23 + 5 * pow(E23, 2))) * E33 +
                      4 * (5 * E11 - 6 * E12 + 5 * E22) * pow(E33, 2)) / (128 * Edet);

     beta0 = fmax(2.0 / E11C11 - 1.5, 0.01);
}

void ShellANDeS::initializeGeometry(double n1, double n2, double n3)
{
     //Local unit vectors
     static Vector e1(3), e2(3), e3(3);
     e1.Zero();
     e2.Zero();
     e3.Zero();
     //node coordinates
     static Vector x1(3);
     static Vector x2(3);
     static Vector x3(3);
     x1 = theNodes[0]->getCrds();
     x2 = theNodes[1]->getCrds();
     x3 = theNodes[2]->getCrds();

     //Centroid
     x0 = (x1 + x2 + x3) / 3;
     
     // Local x axis points paralell to side 1-2
     e1 = x2 - x1;
     e1.Normalize();

     // Local z axis is normal to the triangle
     e3 = calculate_cross_product(x2 - x1, x3 - x1);

     Area = 0.5 * e3.Norm();
     e3.Normalize();

     if (Area < 0)
     {
          cout << "ThreeNodeAndesMembrane::initializeGeometry() -> Element # " << getTag() << " has A < 0!! " << endl;
     }

     //Local y axis is orthogonal to e1 and e3.
     e2 = calculate_cross_product(e3, e1);
     e2.Normalize();

     //If local e1 axis is given, then use it instead but keep the normal!!
     if (n1 != 0 || n2 != 0 || n3 != 0)
     {
          static Vector nn(3);
          static Vector temp1(3);
          static Vector temp2(3);
          nn(0) = n1;
          nn(1) = n2;
          nn(2) = n3;

          // Make sure that e1 is in the plane of the element
          double a = (e1 ^ nn);
          double b = (e2 ^ nn);
          
          temp1 = a*e1;
          temp2 = b*e2;
          
          e1 = temp1 + temp2;
          e1.Normalize();

          // Recompute e2
          e2 = calculate_cross_product(e3, e1);
          e2.Normalize();
     }


     useThisCoordinateSystem(e1, e2, e3);

     initializeMq();

}


void ShellANDeS::initializeBetaArrays()
{
     // We only want to initialize these numbers once and for all elements of this class
     if (number_of_three_node_andes_membrane == 0)
     {

          // Membrane
          beta_membrane.Zero();
          beta_membrane(0) = 0;
          beta_membrane(1) = 1;
          beta_membrane(2) = 2;
          beta_membrane(3) = 1;
          beta_membrane(4) = 0;
          beta_membrane(5) = 1;
          beta_membrane(6) = -1;
          beta_membrane(7) = -1;
          beta_membrane(8) = -1;
          beta_membrane(9) = -2;
     }

}

int
ShellANDeS::displaySelf(Renderer& theViewer, int displayMode, float fact, const char** displayModes, int numModes)
{
    // get the end point display coords
    static Vector v1(3);
    static Vector v2(3);
    static Vector v3(3);
    theNodes[0]->getDisplayCrds(v1, fact, displayMode);
    theNodes[1]->getDisplayCrds(v2, fact, displayMode);
    theNodes[2]->getDisplayCrds(v3, fact, displayMode);

    // place values in coords matrix
    static Matrix coords(3, 3);
    for (int i = 0; i < 3; i++) {
        coords(0, i) = v1(i);
        coords(1, i) = v2(i);
        coords(2, i) = v3(i);
    }

    // basic colors
    static Vector values(3);
    for (int i = 0; i < 3; i++)
        values(i) = 0.0;

    // draw the polygon
    return theViewer.drawPolygon(coords, values, this->getTag());
}