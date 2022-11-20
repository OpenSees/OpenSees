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

// $Revision: 1.10 $
// $Date: 2011/03/10 22:51:21 $
// $Source: /usr/local/cvs/OpenSees/SRC/element/shell/IGAKLShell_BendingStrip.cpp,v $

// Original implementation: Ed "C++" Love
// Reimplementation: Leopoldo Tesser, Diego A. Talledo, Véronique Le Corvec
//
// Bathe MITC 4 four node shell element with membrane and drill
// Ref: Dvorkin,Bathe, A continuum mechanics based four node shell
//      element for general nonlinear analysis,
//      Eng.Comput.,1,77-88,1984

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

//STD library
#include <float.h>
#include <vector>

#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <SectionForceDeformation.h>
#include <Domain.h>
#include <ErrorHandler.h>
#include <IGAKLShell_BendingStrip.h>
#include <R3vectors.h>
#include <Renderer.h>
#include <ElementResponse.h>
#include <ElementalLoad.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>
#include <map>

#define min(a,b) ( (a)<(b) ? (a):(b) )
#define max(a,b) ( (a)>(b) ? (a):(b) )

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#include "gaussQuadrature.h"
#include "R3vectors.h"

#include <set>
using namespace std;
template<typename T>
std::vector<T> arange(T start, T stop, T step = 1) {
  std::vector<T> values;
  for (T value = start; value < stop; value += step)
    values.push_back(value);
  return values;
}

static int numIGAKLShell_BendingStrip = 0;


//static data
Matrix*  IGAKLShell_BendingStrip::stiff = 0;
Matrix*  IGAKLShell_BendingStrip::mass = 0;
Vector*  IGAKLShell_BendingStrip::resid = 0;
// Vector*  IGAKLShell_BendingStrip::load = 0;




//null constructor
IGAKLShell_BendingStrip::IGAKLShell_BendingStrip( ) :
  Element( 0, ELE_TAG_IGAKLShell_BendingStrip ),
  connectedExternalNodes(4)
{

}


//*********************************************************************
//full constructor
IGAKLShell_BendingStrip::IGAKLShell_BendingStrip( int tag,
                        IGASurfacePatch *myPatch_,
                        const ID& nodes,
                        int ngauss_,
                        const Vector& xiE_,
                        const Vector& etaE_,
                        const ID& matTags):
  Element( tag, ELE_TAG_IGAKLShell_BendingStrip ),
  ngauss(ngauss_),
  myPatch(myPatch_),
  xiE(xiE_),
  etaE(etaE_),
  connectedExternalNodes(nodes)
{
  if (numIGAKLShell_BendingStrip == 0) {
    // opserr << "Using IGAKLShell_BendingStrip - Developed by: Felipe Elgueta and Jose A. Abell (www.joseabell.com)\n";
    numIGAKLShell_BendingStrip++;
  }
  // ngauss  = quadorder * quadorder;
  int nLayers = myPatch->getNLayers();

  // ngauss=36;

  quadPoint = new Matrix(ngauss, 2);
  quadWeight = new Vector(ngauss);

  ID PQ = myPatch -> getOrders();
  int P = PQ(0);
  int Q = PQ(1);

  // P=5;
  // Q=5;

  gaussQuad2dNurbs(P + 1, Q + 1, quadPoint, quadWeight);
  // gaussQuad2dNurbs(P + 1, P + 1, quadPoint, quadWeight);
  // gaussQuad2dNurbs(4, 4, quadPoint, quadWeight);

  // opserr << "*quadWeight = " << *quadWeight << endln;
  // opserr << "*quadPoint = " << *quadPoint << endln;

  P = PQ(0);
  Q = PQ(1);



  materialPointers = new NDMaterial** [ngauss]; //Double array
  for (int i = 0; i < ngauss; ++i)
  {
    materialPointers[i] = new NDMaterial* [nLayers];
  }


  for (int gp = 0 ;  gp < ngauss; gp++ )
  {
    for (int capa = 0; capa < nLayers; capa++)
    {
      NDMaterial* theReferenceMaterial = OPS_getNDMaterial(myPatch->getMatTag(capa)); // Pointer to NDMaterial
      NDMaterial* newmat = theReferenceMaterial->getCopy( );  // Copy of pointer to NDMaterial
      materialPointers[gp][capa] =  newmat;   // llama new

      if (materialPointers[gp][capa] == 0) {
        opserr << "ShellMITC4::constructor - failed to get a material of type: ShellSection\n";
      } //end if
    }
  } //end for gp

  applyLoad = 0;

  appliedB[0] = 0.0;
  appliedB[1] = 0.0;
  appliedB[2] = 0.0;

  pressure = 0;

  load = 0; // null pointer

}

//******************************************************************

//destructor
IGAKLShell_BendingStrip::~IGAKLShell_BendingStrip( )
{
  // borrar todos los punteros. ARREGLAR PARA ARREGLO 2D
  int nLayers = myPatch->getNLayers();
  // opserr << "quadorder = " << quadorder << endln;
  // int i ;
  // for ( i = 0 ;  i < 4; i++ ) {

  //   delete materialPointers[i] ;
  //   materialPointers[i] = 0 ;

  // } //end for i
  for (int gp = 0; gp < ngauss; ++gp)
  {
    for (int capa = 0; capa < nLayers; ++capa)
    {
      delete materialPointers[gp][capa];
      materialPointers[gp][capa] = 0;
    }
  }
  // if (Ki != 0)
  //   delete Ki;


  if (load != 0)
    delete load;

  // if (Ki != 0)
  //   delete Ki;

}
//**************************************************************************


//set domain
void  IGAKLShell_BendingStrip::setDomain( Domain *theDomain )
{
  int i;

  int Nnodes = connectedExternalNodes.Size();
  int NDOF = 3 * Nnodes;

  // Initialiaze class-wide static data
  if (numIGAKLShell_BendingStrip == 1)
  {
    stiff = new Matrix(NDOF, NDOF);
    mass = new Matrix(NDOF, NDOF);
    resid = new Vector(NDOF);
    // load = new Vector(NDOF);
  }

  //node pointers
  nodePointers = new Node* [connectedExternalNodes.Size()];

  for ( i = 0; i < connectedExternalNodes.Size(); i++ ) {
    nodePointers[i] = theDomain->getNode( connectedExternalNodes(i) ) ;
    if (nodePointers[i] == 0) {
      opserr << "IGAKLShell_BendingStrip::setDomain - no node " << connectedExternalNodes(i);
      opserr << " exists in the model\n";
    }
  }

  this->DomainComponent::setDomain(theDomain);
}


//get the number of external nodes
int  IGAKLShell_BendingStrip::getNumExternalNodes( ) const
{
  return connectedExternalNodes.Size() ;
}


//return connected external nodes
const ID&  IGAKLShell_BendingStrip::getExternalNodes( )
{
  return connectedExternalNodes ;
}


Node **
IGAKLShell_BendingStrip::getNodePtrs(void)
{
  return nodePointers;
}

//return number of dofs
int  IGAKLShell_BendingStrip::getNumDOF( )
{
  return 3 * connectedExternalNodes.Size() ;
}


//commit state
int  IGAKLShell_BendingStrip::commitState( )
{
  int success = 0 ;

  // call element commitState to do any base class stuff
  if ((success = this->Element::commitState()) != 0) {
    opserr << "IGAKLShell_BendingStrip::commitState () - failed in base class";
  }



  for (int gp = 0; gp < ngauss; gp++ )
    for (int capa = 0; capa < myPatch->getNLayers(); ++capa)
    {
      success += materialPointers[gp][capa]->commitState( ) ;
    }

  return success ;
}



//revert to last commit
int  IGAKLShell_BendingStrip::revertToLastCommit( )
{
  int success = 0 ;

  for (int gp = 0; gp < ngauss; gp++ )
    for (int capa = 0; capa < myPatch->getNLayers(); ++capa)
    {
      success += materialPointers[gp][capa]->revertToLastCommit( ) ;
    }

  return success ;
}


//revert to start
int  IGAKLShell_BendingStrip::revertToStart( )
{
  int success = 0 ;

  for (int gp = 0; gp < ngauss; gp++ )
    for (int capa = 0; capa < myPatch->getNLayers(); ++capa)
    {
      success += materialPointers[gp][capa]->revertToStart( ) ;
    }

  return success ;
}

//print out element data
void  IGAKLShell_BendingStrip::Print( OPS_Stream &s, int flag )
{

}

Response*
IGAKLShell_BendingStrip::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  Response *theResponse = 0;

  // output.tag("ElementOutput");
  // output.attr("eleType", "IGAKLShell_BendingStrip");
  // output.attr("eleTag", this->getTag());
  // int numNodes = this->getNumExternalNodes();
  // const ID &nodes = this->getExternalNodes();
  // static char nodeData[32];

  // for (int i = 0; i < numNodes; i++) {
  //   sprintf(nodeData, "node%d", i + 1);
  //   output.attr(nodeData, nodes(i));
  // }

  // if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "forces") == 0 ||
  //     strcmp(argv[0], "globalForce") == 0 || strcmp(argv[0], "globalForces") == 0) {
  //   const Vector &force = this->getResistingForce();
  //   int size = force.Size();
  //   for (int i = 0; i < size; i++) {
  //     sprintf(nodeData, "P%d", i + 1);
  //     output.tag("ResponseType", nodeData);
  //   }
  //   theResponse = new ElementResponse(this, 1, this->getResistingForce());
  // }

  // else if (strcmp(argv[0], "material") == 0 || strcmp(argv[0], "Material") == 0) {
  //   if (argc < 2) {
  //     opserr << "IGAKLShell_BendingStrip::setResponse() - need to specify more data\n";
  //     return 0;
  //   }
  //   int pointNum = atoi(argv[1]);
  //   if (pointNum > 0 && pointNum <= ngauss) {

  //     output.tag("GaussPoint");
  //     output.attr("number", pointNum);
  //     output.attr("eta", sg[pointNum - 1]);
  //     output.attr("neta", tg[pointNum - 1]);

  //     theResponse =  materialPointers[pointNum - 1]->setResponse(&argv[2], argc - 2, output);

  //     output.endTag();
  //   }

  // } else if (strcmp(argv[0], "stresses") == 0) {

  //   for (int i = 0; i < ngauss; i++) {
  //     output.tag("GaussPoint");
  //     output.attr("number", i + 1);
  //     output.attr("eta", sg[i]);
  //     output.attr("neta", tg[i]);

  //     output.tag("SectionForceDeformation");
  //     output.attr("classType", materialPointers[i]->getClassTag());
  //     output.attr("tag", materialPointers[i]->getTag());

  //     output.tag("ResponseType", "p11");
  //     output.tag("ResponseType", "p22");
  //     output.tag("ResponseType", "p1212");
  //     output.tag("ResponseType", "m11");
  //     output.tag("ResponseType", "m22");
  //     output.tag("ResponseType", "m12");
  //     output.tag("ResponseType", "q1");
  //     output.tag("ResponseType", "q2");

  //     output.endTag(); // GaussPoint
  //     output.endTag(); // NdMaterialOutput
  //   }

  //   theResponse =  new ElementResponse(this, 2, Vector(32));
  // }

  // else if (strcmp(argv[0], "strains") == 0) {

  //   for (int i = 0; i < ngauss; i++) {
  //     output.tag("GaussPoint");
  //     output.attr("number", i + 1);
  //     output.attr("eta", sg[i]);
  //     output.attr("neta", tg[i]);

  //     output.tag("SectionForceDeformation");
  //     output.attr("classType", materialPointers[i]->getClassTag());
  //     output.attr("tag", materialPointers[i]->getTag());

  //     output.tag("ResponseType", "eps11");
  //     output.tag("ResponseType", "eps22");
  //     output.tag("ResponseType", "gamma12");
  //     output.tag("ResponseType", "theta11");
  //     output.tag("ResponseType", "theta22");
  //     output.tag("ResponseType", "theta33");
  //     output.tag("ResponseType", "gamma13");
  //     output.tag("ResponseType", "gamma23");

  //     output.endTag(); // GaussPoint
  //     output.endTag(); // NdMaterialOutput
  //   }

  //   theResponse =  new ElementResponse(this, 3, Vector(32));
  // }

  // output.endTag();
  return theResponse;
}

int
IGAKLShell_BendingStrip::getResponse(int responseID, Information &eleInfo)
{
  // int cnt = 0;
  // static Vector stresses(32);
  // static Vector strains(32);

  // switch (responseID) {
  // case 1: // global forces
  //   return eleInfo.setVector(this->getResistingForce());
  //   break;

  // case 2: // stresses
  //   for (int i = 0; i < ngauss; i++) {

  //     // Get material stress response
  //     const Vector &sigma = materialPointers[i]->getStressResultant();
  //     stresses(cnt) = sigma(0);
  //     stresses(cnt + 1) = sigma(1);
  //     stresses(cnt + 2) = sigma(2);
  //     stresses(cnt + 3) = sigma(3);
  //     stresses(cnt + 4) = sigma(4);
  //     stresses(cnt + 5) = sigma(5);
  //     stresses(cnt + 6) = sigma(6);
  //     stresses(cnt + 7) = sigma(7);
  //     cnt += 8;
  //   }
  //   return eleInfo.setVector(stresses);
  //   break;
  // case 3: //strain
  //   for (int i = 0; i < ngauss; i++) {

  //     // Get section deformation
  //     const Vector &deformation = materialPointers[i]->getSectionDeformation();
  //     strains(cnt) = deformation(0);
  //     strains(cnt + 1) = deformation(1);
  //     strains(cnt + 2) = deformation(2);
  //     strains(cnt + 3) = deformation(3);
  //     strains(cnt + 4) = deformation(4);
  //     strains(cnt + 5) = deformation(5);
  //     strains(cnt + 6) = deformation(6);
  //     strains(cnt + 7) = deformation(7);
  //     cnt += 8;
  //   }
  //   return eleInfo.setVector(strains);
  //   break;
  // default:
  //   return -1;
  // }
  // cnt = 0;
  return 0;
}


//return stiffness matrix
const Matrix&  IGAKLShell_BendingStrip::getTangentStiff( )
{
  // opserr << "IGAKLShell_BendingStrip::getTangentStiff - eleTag" << this->getTag() << " called! " << endln;

  // opserr << "connectedExternalNodes = " << connectedExternalNodes << endln;

  // opserr << "Element number:  = " << this->getTag() << endln;

  // opserr << "xiE = " << xiE << endln;
  // opserr << "etaE = " << etaE << endln;

  int tang_flag = 1 ; //get the tangent

  //do tangent and residual here
  formResidAndTangent( tang_flag ) ;
  // opserr << "(*stiff) = " << (*stiff) << endln;

  return *stiff ;



  // stiff->Zero();

  // bool nonLinearGeometry = myPatch->getAnalysisType();

  // // opserr << "xiE = " << xiE << endln;
  // // opserr << "etaE = " << etaE << endln;

  // float wt;
  // float ptU;
  // float ptV;
  // double xi;
  // double eta;

  // int noFuncs = myPatch->getNoFuncs();

  // Vector R(noFuncs);
  // Vector dRdxi(noFuncs);
  // Vector dRdeta(noFuncs);
  // Vector dR2dxi(noFuncs);
  // Vector dR2deta(noFuncs);
  // Vector dR2dxideta(noFuncs);

  // double J2; //Jacobian for the parametric mapping

  // Vector point(3);
  // Vector point_disp(3);
  // Matrix pts(connectedExternalNodes.Size(), 3);
  // Matrix pts_d(connectedExternalNodes.Size(), 3);
  // // Get pts span matrix for jacobian
  // for (int i = 0; i < connectedExternalNodes.Size(); ++i)
  // {
  //   point = nodePointers[i]->getCrds();
  //   point_disp = nodePointers[i]->getTrialDisp(); //getTrialDisp
  //   // points_disp = nodePointers[i]->getDisp(); //getTrialDisp
  //   for (int j = 0; j < 3; ++j)
  //   {
  //     pts(i, j) = point(j);
  //     if (isnan(point_disp(j)))
  //     {
  //       opserr << "Nan found on getTangentStiff = " << endln;
  //     }
  //     if (nonLinearGeometry)
  //     {
  //       pts_d(i, j) = pts(i, j) + point_disp(j);
  //     }
  //     else
  //     {
  //       pts_d(i, j) = pts(i, j);
  //     }
  //   }
  // }
  // // opserr << "connectedExternalNodes = " << connectedExternalNodes << endln;
  // // opserr << "pts = " << pts << endln;

  // // opserr << "pts = " << pts << endln;

  // // Loop over integrations points
  // for (int gp = 0; gp < ngauss; ++gp)
  // {
  //   wt = (*quadWeight)(gp);
  //   ptU = (*quadPoint)(gp, 0);
  //   ptV = (*quadPoint)(gp, 1);
  //   xi = myPatch->parent2ParametricSpace(xiE, ptU);
  //   eta = myPatch->parent2ParametricSpace(etaE, ptV);
  //   J2 = 0.5 * (xiE(1) - xiE(0)) * 0.5 * (etaE(1) - etaE(0));
  //   // opserr << "xi = " << xi << endln;
  //   // opserr << "eta = " << eta << endln;
  //   R.Zero();
  //   dRdxi.Zero();
  //   dRdeta.Zero();
  //   dR2dxi.Zero();
  //   dR2deta.Zero();
  //   dR2dxideta.Zero();
  //   myPatch->Nurbs2DBasis2ndDers(xi, eta, R, dRdxi, dRdeta, dR2dxi, dR2deta, dR2dxideta);

  //   // opserr << "R = " << R << endln;
  //   // opserr << "dRdxi = " << dRdxi << endln;
  //   // opserr << "dRdeta = " << dRdeta << endln;
  //   // opserr << "dR2dxi = " << dR2dxi << endln;
  //   // opserr << "dR2deta = " << dR2deta << endln;
  //   // opserr << "dR2dxideta = " << dR2dxideta << endln;




  //   // Get first and second order jacobians
  //   Matrix dr(2, noFuncs);
  //   Matrix dr2(3, noFuncs);


  //   for (int i = 0; i < noFuncs; ++i)
  //   {
  //     dr(0, i) = dRdxi(i);
  //     dr(1, i) = dRdeta(i);

  //     dr2(0, i) = dR2dxi(i);
  //     dr2(1, i) = dR2deta(i);
  //     dr2(2, i) = dR2dxideta(i);
  //   }

  //   //compute the jacobian of physical and parameter domain mapping
  //   // then the derivative w.r.t spatial physical coordinates, current configuration


  //   Matrix jacob = dr * pts;  //G: Covariant base vectors, needed in the reference configuration
  //   Matrix jacob2 = dr2 * pts; //H: Hessian matrix, needed in the reference configuration


  //   // !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  //   // !    DESCRIPTION OF THE VARIABLES !
  //   // !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

  //   // ! g: matrix with the 3 components of the in plane vectors g1, g2
  //   // ! dR: first derivatives of the shape functions
  //   // ! ddR: second derivatives of the shape functions
  //   // ! dn: first variation of the normal vector g3
  //   // ! lg3: length of the normal vector g3

  //   // ! E_cu: strain vector [E11,E22,E12] referred to the curvilinear coord. sys
  //   // ! E_ca: strain vector [E11,E22,2*E12] referred to cartesian coord. sys
  //   // ! K_cu: curvature vector [K11,K22,K12] referred to the curvilinear coord. sys
  //   // ! K_ca: curvature vector [K11,K22,2*K12] referred to the cartesian coord. sys

  //   // ! dE_cu: first variation of the membrane strain in local coord. system
  //   // ! dE_ca: first variation of the membrane strain in global coord. system
  //   // ! dK_cu: first variation of the curvature in local coord. system
  //   // ! dK_ca: first variation of the curvature in global coord. system

  //   // ! ddE_cu: second variation of the membrane strain in local coord. system
  //   // ! ddE_ca: second variation of the membrane strain in global coord. system
  //   // ! ddK_cu: second variation of the curvature in local coord. system
  //   // ! ddK_ca: second variation of the curvature in global coord. system

  //   Vector E_cu(3);
  //   Vector K_cu(3);

  //   Vector dE_cu(3);
  //   Matrix dE_ca(3, 3 * noFuncs);
  //   Vector dK_cu(3);
  //   Matrix dK_ca(3, 3 * noFuncs);

  //   Matrix ddR(3, noFuncs);
  //   for (int i = 0; i < noFuncs; ++i)
  //   {
  //     ddR(0, i) = dR2dxi(i);
  //     ddR(1, i) = dR2deta(i);
  //     ddR(2, i) = dR2dxideta(i);
  //   }


  //   Matrix G = dr  * pts;
  //   Matrix H = ddR * pts;

  //   Matrix g = dr  * pts_d;
  //   Matrix h = ddR * pts_d;

  //   G = transpose(2, 3, G); // Transposing because it is needed in cols drdxi, drdeta
  //   g = transpose(2, 3, g);
  //   H = transpose(3, 3, H);
  //   h = transpose(3, 3, h);

  //   // Shell geo reference parameters and temporal
  //   Vector g3_tmp(3);
  //   Vector n_tmp(3);
  //   double dA; // = lG3
  //   Matrix Gab_r(2, 2); // = Gab
  //   Vector Bv_r(3); // = Bv
  //   Matrix Tb(3, 3); // = T_Gcon_E
  //   Matrix T_E_G_tmp(3, 3); // = T_Gcon_E
  //   Matrix T_G_E_tmp(3, 3); // = T_Gcon_E

  //   // Shell geo deformed parameters and temporal
  //   Vector g3(3); // = g3
  //   Vector n(3); // = n
  //   double lg3; // = dA
  //   Matrix gab(2, 2); // = Gab
  //   Vector bv(3); // = Bv
  //   Matrix Tb_tmp(3, 3);

  //   // Call shell geometry functions
  //   shellGeo(G , H , g3_tmp , dA  , n_tmp , Gab_r , Bv_r , Tb, T_E_G_tmp, T_G_E_tmp);
  //   shellGeo(g , h , g3     , lg3 , n     , gab   , bv   , Tb_tmp, T_E_G_tmp, T_G_E_tmp);

  //   //
  //   // Compute the membrane strains and the curvatures according to the eqns 1.161 - 1.162 of Computational FSI (Yuri's book)
  //   //

  //   //strain vector [E11,E22,E12] referred to curvilinear coor sys
  //   E_cu(0) = 0.5 * (gab(0, 0) - Gab_r(0, 0));
  //   E_cu(1) = 0.5 * (gab(1, 1) - Gab_r(1, 1));
  //   E_cu(2) = 0.5 * (gab(0, 1) - Gab_r(0, 1));

  //   //curvature vector [K11,K22,K12] referred to curvilinear coor sys
  //   K_cu(0) = -(bv(0) - Bv_r(0));
  //   K_cu(1) = -(bv(1) - Bv_r(1));
  //   K_cu(2) = -(bv(2) - Bv_r(2));


  //   Vector E_ca = Tb * E_cu; // strain vector [E11,E22,2*E12] referred to cartesian coor sys
  //   Vector K_ca = Tb * K_cu; //curvature vector [K11,K22,2*K12] referred to cartesian coor sys

  //   // opserr << "E_ca = " << E_ca << endln;
  //   // opserr << "K_ca = " << K_ca << endln;

  //   //
  //   // Compute the First variation of strain and curvature w.r.t the d.o.f.
  //   //

  //   // Local variables
  //   int kr = 0;
  //   int dirr = 0;
  //   Matrix dg(3, 2);
  //   Matrix dR = transpose(2, noFuncs, dr);
  //   ddR = transpose(3, noFuncs, ddR);
  //   Matrix dg3(3, 3 * noFuncs);
  //   Vector g3dg3(3 * noFuncs);
  //   Vector g3dg3lg3_3(3 * noFuncs);
  //   Matrix dn(3, 3 * noFuncs);

  //   // For "slicing" matrices
  //   ID ur_only(1);
  //   ID us_only(1);
  //   auto xyz_temp = arange<int>(0, 3); ID xyz(&(xyz_temp[0]), 3); // Takes value 0 1 2

  //   double lg3_3 = pow(lg3, 3);



  //   for (int ur = 0; ur < 3 * noFuncs; ++ur)
  //   {
  //     // Local node number kr and dof direction dirr
  //     kr = (ur + 1 + 2) / 3 - 1; // takes value 0 0 0, 1 1 1, 2 2 2, .... noFuncs-1 noFuncs-1 noFuncs-1
  //     dirr = ur - 3 * (kr); // Takes value 0 1 2 0 1 2 ...

  //     // opserr << "kr = " << kr << endln;
  //     // opserr << "dirr = " << dirr << endln;

  //     dg(dirr, 0) = dR(kr, 0);
  //     dg(dirr, 1) = dR(kr, 1);

  //     // Strain
  //     dE_cu(0) = dR(kr, 0) * g(dirr, 0);
  //     dE_cu(1) = dR(kr, 1) * g(dirr, 1);
  //     dE_cu(2) = 0.5 * (dR(kr, 0) * g(dirr, 1) + g(dirr, 0) * dR(kr, 1));

  //     dE_ca(0, ur) = (Tb * dE_cu)(0);
  //     dE_ca(1, ur) = (Tb * dE_cu)(1);
  //     dE_ca(2, ur) = (Tb * dE_cu)(2);


  //     // Curvature
  //     dg3(0, ur) = dg(1, 0) * g(2, 1) - dg(2, 0) * g(1, 1) + g(1, 0) * dg(2, 1) - g(2, 0) * dg(1, 1);
  //     dg3(1, ur) = dg(2, 0) * g(0, 1) - dg(0, 0) * g(2, 1) + g(2, 0) * dg(0, 1) - g(0, 0) * dg(2, 1);
  //     dg3(2, ur) = dg(0, 0) * g(1, 1) - dg(1, 0) * g(0, 1) + g(0, 0) * dg(1, 1) - g(1, 0) * dg(0, 1);

  //     g3dg3(ur) = g3(0) * dg3(0, ur) + g3(1) * dg3(1, ur) + g3(2) * dg3(2, ur);
  //     g3dg3lg3_3(ur) = g3dg3(ur) / lg3_3;


  //     dn(0, ur) = dg3(0, ur) / lg3 - g3(0) * g3dg3lg3_3(ur);
  //     dn(1, ur) = dg3(1, ur) / lg3 - g3(1) * g3dg3lg3_3(ur);
  //     dn(2, ur) = dg3(2, ur) / lg3 - g3(2) * g3dg3lg3_3(ur);


  //     dK_cu(0) = -(ddR(kr, 0) * n(dirr) + (h(0, 0) * dn(0, ur) + h(1, 0) * dn(1, ur) + h(2, 0) * dn(2, ur)));
  //     dK_cu(1) = -(ddR(kr, 1) * n(dirr) + (h(0, 1) * dn(0, ur) + h(1, 1) * dn(1, ur) + h(2, 1) * dn(2, ur)));
  //     dK_cu(2) = -(ddR(kr, 2) * n(dirr) + (h(0, 2) * dn(0, ur) + h(1, 2) * dn(1, ur) + h(2, 2) * dn(2, ur)));

  //     // dK_cu(0) = -(ddR(0, kr) * n(dirr) + (h(0, 0) * dn(0, ur) + h(1, 0) * dn(1, ur) + h(2, 0) * dn(2, ur)));
  //     // dK_cu(1) = -(ddR(1, kr) * n(dirr) + (h(0, 1) * dn(0, ur) + h(1, 1) * dn(1, ur) + h(2, 1) * dn(2, ur)));
  //     // dK_cu(2) = -(ddR(2, kr) * n(dirr) + (h(0, 2) * dn(0, ur) + h(1, 2) * dn(1, ur) + h(2, 2) * dn(2, ur)));

  //     // CHECK THIS OUT, I HAD TO SWAP THE INDICES FOR H

  //     // dK_cu(0) = -(ddR(0, kr) * n(dirr) + (h(0, 0) * dn(0, ur) + h(0, 1) * dn(1, ur) + h(0, 2) * dn(2, ur)));
  //     // dK_cu(1) = -(ddR(1, kr) * n(dirr) + (h(1, 0) * dn(0, ur) + h(1, 1) * dn(1, ur) + h(1, 2) * dn(2, ur)));
  //     // dK_cu(2) = -(ddR(2, kr) * n(dirr) + (h(2, 0) * dn(0, ur) + h(2, 1) * dn(1, ur) + h(2, 2) * dn(2, ur)));

  //     dK_ca(0, ur) = (Tb * dK_cu)(0);
  //     dK_ca(1, ur) = (Tb * dK_cu)(1);
  //     dK_ca(2, ur) = (Tb * dK_cu)(2);


  //   }


  //   // Second variation of strain and curvature w.r.t. dofs
  //   Vector ddE_cu(3);
  //   Matrix ddE_ca_0(3 * noFuncs, 3 * noFuncs); // Wanted to avoid declarin a 3D array, instead using 3 matrices
  //   Matrix ddE_ca_1(3 * noFuncs, 3 * noFuncs);
  //   Matrix ddE_ca_2(3 * noFuncs, 3 * noFuncs);

  //   Vector ddK_cu(3);
  //   Matrix ddK_ca_0(3 * noFuncs, 3 * noFuncs); // Wanted to avoid declarin a 3D array, instead using 3 matrices
  //   Matrix ddK_ca_1(3 * noFuncs, 3 * noFuncs);
  //   Matrix ddK_ca_2(3 * noFuncs, 3 * noFuncs);

  //   // Local variables
  //   int ks = 0;
  //   int dirs = 0;
  //   int dirt = 0;
  //   int ddir = 0;
  //   Vector ddg3(3);
  //   double tmp1 = 0.0;
  //   double tmp2 = 0.0;
  //   double lg3_5 = pow(lg3, 5);
  //   Vector ddn(3);



  //   if (nonLinearGeometry) // true if non-linear geometry
  //   {
  //     // opserr << "Using nonLinearGeometry!" << endln;
  //     for (int ur = 0; ur < 3 * noFuncs; ++ur)
  //     {
  //       // Local node number kr and dof direction dirr
  //       kr = (ur + 1 + 2) / 3 - 1; // takes value 0 0 0, 1 1 1, 2 2 2, .... noFuncs-1 noFuncs-1 noFuncs-1
  //       dirr = ur - 3 * (kr); // Takes value 0 1 2 0 1 2 ...

  //       for (int us = 0; us <= ur; ++us)
  //       {
  //         ks = (us + 1 + 2) / 3 - 1; // takes value 0 0 0, 1 1 1, 2 2 2, .... noFuncs-1 noFuncs-1 noFuncs-1
  //         dirs = ur - 3 * (ks); // Takes value 0 1 2 0 1 2 ...

  //         // Strain
  //         if (dirr == dirs)
  //         {
  //           ddE_cu(0) = dR(kr, 0) * dR(ks, 0);
  //           ddE_cu(1) = dR(kr, 1) * dR(ks, 1);
  //           ddE_cu(2) = 0.5 * (dR(kr, 0) * dR(ks, 1) + dR(kr, 1) * dR(ks, 0));
  //         }

  //         ddE_ca_0(ur, us) = (Tb * ddE_cu)(0);
  //         ddE_ca_1(ur, us) = (Tb * ddE_cu)(1);
  //         ddE_ca_2(ur, us) = (Tb * ddE_cu)(2);

  //         // Curvature

  //         dirt = 6 - dirr - dirs; // Check this
  //         ddir = dirr - dirs;

  //         if (ddir == -1 or ddir == 2)
  //         {
  //           ddg3(dirt) = dR(kr, 0) * dR(ks, 1) - dR(ks, 0) * dR(kr, 1);
  //         }
  //         else if (ddir == 1 or ddir == -2)
  //         {
  //           ddg3(dirt) = -dR(kr, 0) * dR(ks, 1) + dR(ks, 0) * dR(kr, 1);
  //         }

  //         ur_only.fill(ur);
  //         us_only.fill(us);

  //         tmp1 = -(ddg3 ^ g3) + (transpose(3, 1, dg3(xyz, ur_only)) * dg3(xyz, us_only) / lg3_3)(0, 0); // the (0,0) is to transform the 1x1 matrix into double
  //         tmp2 = 3 * g3dg3(ur) * g3dg3(us) / lg3_5;

  //         ddn(0) = ddg3(0) / lg3 - g3dg3lg3_3(us) * dg3(0, ur) - g3dg3lg3_3(ur) * dg3(0, us) + tmp1 * g3(0) + tmp2 * g3(0);
  //         ddn(1) = ddg3(1) / lg3 - g3dg3lg3_3(us) * dg3(1, ur) - g3dg3lg3_3(ur) * dg3(1, us) + tmp1 * g3(1) + tmp2 * g3(1);
  //         ddn(2) = ddg3(2) / lg3 - g3dg3lg3_3(us) * dg3(2, ur) - g3dg3lg3_3(ur) * dg3(2, us) + tmp1 * g3(2) + tmp2 * g3(2);

  //         // ddK_cu(0) = -(ddR(0, kr) * dn(dirr, us) + ddR(0, ks) * dn(dirs, ur) + (h(0, 0) * ddn(0) + h(0, 1) * ddn(1) + h(0, 2) * ddn(2)));
  //         // ddK_cu(1) = -(ddR(1, kr) * dn(dirr, us) + ddR(1, ks) * dn(dirs, ur) + (h(1, 0) * ddn(0) + h(1, 1) * ddn(1) + h(1, 2) * ddn(2)));
  //         // ddK_cu(2) = -(ddR(2, kr) * dn(dirr, us) + ddR(2, ks) * dn(dirs, ur) + (h(2, 0) * ddn(0) + h(2, 1) * ddn(1) + h(2, 2) * ddn(2)));

  //         // ddK_cu(0) = -(ddR(0, kr) * dn(dirr, us) + ddR(0, ks) * dn(dirs, ur) + (h(0, 0) * ddn(0) + h(1, 0) * ddn(1) + h(2, 0) * ddn(2)));
  //         // ddK_cu(1) = -(ddR(1, kr) * dn(dirr, us) + ddR(1, ks) * dn(dirs, ur) + (h(0, 1) * ddn(0) + h(1, 1) * ddn(1) + h(2, 1) * ddn(2)));
  //         // ddK_cu(2) = -(ddR(2, kr) * dn(dirr, us) + ddR(2, ks) * dn(dirs, ur) + (h(0, 2) * ddn(0) + h(1, 2) * ddn(1) + h(2, 2) * ddn(2)));

  //         ddK_cu(0) = -(ddR(kr, 0) * dn(dirr, us) + ddR(ks, 0) * dn(dirs, ur) + (h(0, 0) * ddn(0) + h(1, 0) * ddn(1) + h(2, 0) * ddn(2)));
  //         ddK_cu(1) = -(ddR(kr, 1) * dn(dirr, us) + ddR(ks, 1) * dn(dirs, ur) + (h(0, 1) * ddn(0) + h(1, 1) * ddn(1) + h(2, 1) * ddn(2)));
  //         ddK_cu(2) = -(ddR(kr, 2) * dn(dirr, us) + ddR(ks, 2) * dn(dirs, ur) + (h(0, 2) * ddn(0) + h(1, 2) * ddn(1) + h(2, 2) * ddn(2)));

  //         ddK_ca_0(ur, us) = (Tb * ddK_cu)(0);
  //         ddK_ca_1(ur, us) = (Tb * ddK_cu)(1);
  //         ddK_ca_2(ur, us) = (Tb * ddK_cu)(2);
  //       }
  //     }
  //   }

  //   // Leaving this here momentarily, further on it is done again, should delete and stick with this
  //   int nLayers = myPatch->getNLayers();

  //   Matrix A(3, 3);   // Membrane stiffness matrix
  //   Matrix B(3, 3);   // Coupling stifness matrix
  //   Matrix D(3, 3);   // Bending stiffness matrix
  //   Matrix T(3, 3);   // Rotating laminate matrix
  //   Matrix Cbar(3, 3); // Equivalent constitutive matrix
  //   for (int capa = 0; capa < nLayers; ++capa)
  //   {
  //     double iAngle     = myPatch->getAngle(capa);
  //     double iThickness = myPatch->getThickness(capa);
  //     double iZ         = myPatch->getZk(capa);// Mid center laminate position, 1 in the meantime

  //     const Matrix& Czig = materialPointers[gp][capa]->getTangent();

  //     T(0, 0) = pow(cos(iAngle), 2)                        ; T(0, 1) = pow(sin(iAngle), 2)           ; T(0, 2) = sin(iAngle) * cos(iAngle)               ;
  //     T(1, 0) = pow(sin(iAngle), 2)                        ; T(1, 1) = pow(cos(iAngle), 2)           ; T(1, 2) = -sin(iAngle) * cos(iAngle)              ;
  //     T(2, 0) = -2 * sin(iAngle) * cos(iAngle)             ; T(2, 1) = 2 * sin(iAngle) * cos(iAngle) ; T(2, 2) = pow(cos(iAngle), 2) - pow(sin(iAngle), 2) ;

  //     Cbar.addMatrixTripleProduct(0, T, Czig, 1.0);

  //     A.addMatrix(1.0, Cbar, iThickness);
  //     B.addMatrix(1.0, Cbar, iThickness * iZ);
  //     D.addMatrix(1.0, Cbar, (iThickness * pow(iZ, 2) + pow(iThickness, 3) / 12));
  //   }

  //   Vector N_ca = A * E_ca + B * K_ca; // Membrane forces
  //   Vector M_ca = B * E_ca + D * K_ca; // Bending moments

  //   // opserr << "A = " << A << endln;
  //   // opserr << "B = " << B << endln;
  //   // opserr << "D = " << D << endln;


  //   Matrix ke(3 * noFuncs, 3 * noFuncs);
  //   Matrix kem(3 * noFuncs, 3 * noFuncs);
  //   Matrix keb(3 * noFuncs, 3 * noFuncs);



  //   for (int ur = 0; ur < 3 * noFuncs ; ++ur)
  //   {
  //     ur_only.fill(ur);


  //     Matrix dN_ca = A * dE_ca(xyz, ur_only) + B * dK_ca(xyz, ur_only);
  //     Matrix dM_ca = B * dE_ca(xyz, ur_only) + D * dK_ca(xyz, ur_only);

  //     for (int us = 0; us <= ur; ++us)
  //     {
  //       us_only.fill(us);

  //       // Vectors for non-linear part
  //       Vector ddE_ca(3);
  //       ddE_ca(0) = ddE_ca_0(ur, us);
  //       ddE_ca(1) = ddE_ca_1(ur, us);
  //       ddE_ca(2) = ddE_ca_2(ur, us);

  //       Vector ddK_ca(3);
  //       ddK_ca(0) = ddK_ca_0(ur, us);
  //       ddK_ca(1) = ddK_ca_1(ur, us);
  //       ddK_ca(2) = ddK_ca_2(ur, us);

  //       // Membrane stiffness
  //       kem(ur, us) = (transpose(3, 1, dN_ca) * dE_ca(xyz, us_only))(0, 0); // Linear part

  //       // Bending stiffness
  //       keb(ur, us) = (transpose(3, 1, dM_ca) * dK_ca(xyz, us_only))(0, 0); // Linear part


  //       if (nonLinearGeometry) // 1 if nonLinear, 0 if Linear
  //       {
  //         kem(ur, us) += N_ca ^ ddE_ca; // Non-linear part
  //         keb(ur, us) += M_ca ^ ddK_ca; // Non-linear part
  //       }
  //       // Symmetric parts
  //       kem(us, ur) = kem(ur, us);
  //       keb(us, ur) = keb(ur, us);
  //     }
  //   }

  //   ke = (kem + keb) * (dA * J2 * wt); //Stiffness matrix on gauss point
  //   *stiff += ke; // Kiendl stiffness

  // }




  // // opserr << "K.noRows() = " << K.noRows() << endln;
  // // opserr << "K.noCols() = " << K.noCols() << endln;
  // // opserr << "ke.noRows() = " << ke.noRows() << endln;
  // // opserr << "ke.noCols() = " << ke.noCols() << endln;

  // // opserr << "Ke = " << memStiff * transpose(Bmem.noRows(), Bmem.noCols(), Bmem) * C * Bmem * J1 * J2 * wt + benStiff * transpose(Bben.noRows(), Bben.noCols(), Bben) * C * Bben * J1 * J2 * wt << endln;


  // // opserr << "Finished making Ke!!" << endln << endln;
  // // opserr << "K = " << K << endln;

  // // opserr << "K = " << K << endln;
  // // return K ;
  // // opserr << "*stiff = " << *stiff << endln;
  // return *stiff;
}

//return secant matrix
const Matrix&  IGAKLShell_BendingStrip::getInitialStiff( )
{
  opserr << "IGAKLShell_BendingStrip::getInitialStiff - eleTag" << this->getTag() << " called! " << endln;

  Matrix& K = *stiff;

  // K = getTangentStiff();

  // Rellenar K(i,j)  (inicial) aqui (no es necesario por ahora)

  return K ;
}


//return mass matrix
const Matrix&  IGAKLShell_BendingStrip::getMass( )
{
  // opserr << "IGAKLShell_BendingStrip::getMass - eleTag" << this->getTag() << " called! " << endln;

  // double t = 0.0; // thickness [m]
  // double rho = 0.0; //kg/m^3

  // double W = 0.0; //Total laminate weight
  // int nLayers = myPatch->getNLayers();

  // for (int capa = 0; capa < nLayers; ++capa)
  // {
  //   rho = (OPS_getNDMaterial(myPatch->getMatTag(capa)))->getRho(); // Density of material
  //   t = myPatch->getThickness(capa);
  //   W += rho * t;
  // }

  mass->Zero();

  // int noFuncs = myPatch->getNoFuncs();

  // Vector R(noFuncs);
  // Vector dRdxi(noFuncs);
  // Vector dRdeta(noFuncs);
  // Vector dR2dxi(noFuncs);
  // Vector dR2deta(noFuncs);
  // Vector dR2dxideta(noFuncs);

  // Matrix N(3, 3 * noFuncs);

  // float wt;
  // float ptU;
  // float ptV;
  // double xi;
  // double eta;

  // double J2;

  // Vector point(3);
  // Matrix pts(connectedExternalNodes.Size(), 3);

  // // Get pts span matrix for jacobian
  // for (int i = 0; i < connectedExternalNodes.Size(); ++i)
  // {
  //   point = nodePointers[i]->getCrds();
  //   for (int j = 0; j < 3; ++j)
  //   {
  //     pts(i, j) = point(j);
  //   }
  // }

  // // Loop over integrations points
  // for (int gp = 0; gp < ngauss; ++gp)
  // {
  //   wt = (*quadWeight)(gp);
  //   ptU = (*quadPoint)(gp, 0);
  //   ptV = (*quadPoint)(gp, 1);
  //   xi = myPatch->parent2ParametricSpace(xiE, ptU);
  //   eta = myPatch->parent2ParametricSpace(etaE, ptV);
  //   J2 = 0.5 * (xiE(1) - xiE(0)) * 0.5 * (etaE(1) - etaE(0));

  //   myPatch->Nurbs2DBasis2ndDers(xi, eta, R, dRdxi, dRdeta, dR2dxi, dR2deta, dR2dxideta);

  //   // Get first order jacobian
  //   Matrix dr(2, noFuncs);

  //   for (int i = 0; i < noFuncs; ++i)
  //   {
  //     dr(0, i) = dRdxi(i);
  //     dr(1, i) = dRdeta(i);
  //   }

  //   //compute the jacobian of physical and parameter domain mapping
  //   // then the derivative w.r.t spatial physical coordinates

  //   Matrix jacob = dr * pts;

  //   // a1, a2 and a3 vectors (surface basis vectors)
  //   // and its derivatives

  //   Vector a1(3);
  //   Vector a2(3);

  //   for (int i = 0; i < 3; ++i)
  //   {
  //     a1(i) = jacob(0, i);
  //     a2(i) = jacob(1, i);
  //   }
  //   Vector a3 = LovelyCrossProduct(a1, a2);
  //   double norma = a3.Norm();
  //   a3 /= norma;
  //   double J1 = norma;

  //   // Forming N matrix
  //   for (int i = 0; i < noFuncs; ++i)
  //   {
  //     N(0, 3 * i) = R(i);
  //     N(1, 3 * i + 1) = R(i);
  //     N(2, 3 * i + 2) = R(i);
  //   }

  //   // M += transpose(3, 3 * noFuncs, N) * N * rho * J1 * J2 * wt * t;
  //   // M.addMatrixTransposeProduct(1.0, N, N, rho * J1 * J2 * wt * t);
  //   (*mass).addMatrixTransposeProduct(1.0, N, N, W * J1 * J2 * wt );

  // }

  // opserr << "Finished making M!!! " << endln << endln;
  // opserr << "M = " << M << endln;

  // Rellenar M(i,j)   aqui

  return *mass ;
}

void IGAKLShell_BendingStrip::shellGeo(Matrix G, Matrix H, Vector& G3, double& dA, Vector& N, Matrix& Gab, Vector& Bv, Matrix& T_Gcon_E, Matrix& T_E_G, Matrix& T_G_E) // Get geometric quantities

{
  // opserr << "Starting shellGeo!" << endln;


  //G: Covariant base vectors, needed in the reference configuration
  //H: Hessian matrix, needed in the reference configuration


  // Basis vector G3 (covariant) on reference configuration

  // Vector G3(3); //Defined in input
  G3(0) = G(1, 0) * G(2, 1) - G(2, 0) * G(1, 1);
  G3(1) = G(2, 0) * G(0, 1) - G(0, 0) * G(2, 1);
  G3(2) = G(0, 0) * G(1, 1) - G(1, 0) * G(0, 1);

  // differential area dA = length of G3
  // double dA = G3.Norm(); //Defined in input
  dA = G3.Norm();

  // Normal vector N
  // Vector N = G3 / dA; //Defined in input
  N = G3 / dA;

  // Curvature coefficients in vector notation
  // Vector Bv(3); //Defined in input
  Bv.Zero();
  for (int i = 0; i < 3; ++i)
  {
    Bv(i) = H(0, i) * N(0) + H(1, i) * N(1) + H(2, i) * N(2) ;
    // Bv(i) = H(i, 0) * N(0) + H(i, 0) * N(1) + H(i, 2) * N(2) ; // Testing if the error is here
  }

  // Covariant metric coefficients Gab
  // Matrix Gab(2,2); //Defined in input
  Gab.addMatrixTransposeProduct(0, G, G, 1);

  // Contravariant metric coefficients Gab_con and contravariant base vectors G_con
  double invdetGab = 1.0 / (Gab(0, 0) * Gab(1, 1) - Gab(0, 1) * Gab(0, 1));
  // Testing this
  // double invdetGab = 1.0 / ( Gab(0, 0) * Gab(1, 1) - Gab(0, 1) * Gab(1, 0) );
  Matrix Gab_con(2, 2);
  Gab_con(0, 0) =  invdetGab * Gab(1, 1);
  Gab_con(0, 1) = -invdetGab * Gab(0, 1);
  Gab_con(1, 1) =  invdetGab * Gab(0, 0);
  Gab_con(1, 0) =  Gab_con(0, 1);
  // Gab_con(1, 0) =  -invdetGab*Gab(1, 0); // It should be this

  Matrix G_con(3, 2);
  G_con.Zero();
  G_con = G * transpose(2, 2, Gab_con);


  // Local cartesian coordinates
  Vector g1(3);
  Vector g2_con(3);
  for (int i = 0; i < 3; ++i)
  {
    g1(i) = G(i, 0);
    g2_con(i) = G_con(i, 1);
  }
  double lg1 = g1.Norm();
  double lg2_con = g2_con.Norm();

  Matrix E(3, 2);
  for (int i = 0; i < 3; ++i)
  {
    E(i, 0) = g1(i) / lg1;
    E(i, 1) = g2_con(i) / lg2_con;
  }

  // Quick Fix for bending strip
  int Q=myPatch->getOrders()(1);
  // if (Q==1)
  if (Q==2)
  {
    for (int i = 0; i < 3; ++i)
    {
      g1(i)=G(i,1);
      g2_con(i)=G_con(i,0);
    }
    lg1 = g1.Norm();
    lg2_con = g2_con.Norm();
    for (int i = 0; i < 3; ++i)
    {
      E(i, 0) = g1(i) / lg1;
      E(i, 1) = g2_con(i) / lg2_con;
    }
  }


  // Transformation matrix T_Gcon_E from G_con to E
  // with 2 in last row for strain in Voigt notation
  Matrix EG(2, 2);
  EG.addMatrixTransposeProduct(0, E, G_con, 1);

  // Matrix T_Gcon_E(3,3); //Defined in input

  T_Gcon_E(0, 0) = pow(EG(0, 0), 2)        ; T_Gcon_E(0, 1) = pow(EG(0, 1), 2)        ; T_Gcon_E(0, 2) = 2 * EG(0, 0) * EG(0, 1)                       ;
  T_Gcon_E(1, 0) = pow(EG(1, 0), 2)        ; T_Gcon_E(1, 1) = pow(EG(1, 1), 2)        ; T_Gcon_E(1, 2) = 2 * EG(1, 0) * EG(1, 1)                       ;
  T_Gcon_E(2, 0) = 2 * EG(0, 0) * EG(1, 0) ; T_Gcon_E(2, 1) = 2 * EG(0, 1) * EG(1, 1) ; T_Gcon_E(2, 2) = 2 * EG(0, 0) * EG(1, 1) + EG(0, 1) * EG(1, 0) ;


  // Transformation matrix T_E_G from E to G (for PK2 stress)

  T_E_G(0, 0) = pow(EG(0, 0), 2)  ; T_E_G(0, 1) = pow(EG(1, 0), 2) ;  T_E_G(0, 2) = 2 * EG(0, 0) * EG(1, 0)               ;
  T_E_G(1, 0) = pow(EG(0, 1), 2)  ; T_E_G(1, 1) = pow(EG(1, 1), 2) ;  T_E_G(1, 2) = 2 * EG(0, 1) * EG(1, 1)               ;
  T_E_G(2, 0) = EG(0, 0) * EG(0, 1) ; T_E_G(2, 1) = EG(1, 0) * EG(1, 1);  T_E_G(2, 2) = EG(0, 0) * EG(1, 1) + EG(1, 0) * EG(0, 1) ;


  // Transformation matrix T_g_e from g to e (for Cauchy stress)
  EG.addMatrixTransposeProduct(0, E, G, 1);

  T_G_E(0, 0) = pow(EG(0, 0), 2)  ; T_G_E(0, 1) = pow(EG(0, 1), 2) ;  T_G_E(0, 2) = 2 * EG(0, 0) * EG(0, 1)               ;
  T_G_E(1, 0) = pow(EG(1, 0), 2)  ; T_G_E(1, 1) = pow(EG(1, 1), 2) ;  T_G_E(1, 2) = 2 * EG(1, 0) * EG(1, 1)               ;
  T_G_E(2, 0) = EG(0, 0) * EG(1, 0) ; T_G_E(2, 1) = EG(0, 1) * EG(1, 1);  T_G_E(2, 2) = EG(0, 0) * EG(1, 1) + EG(0, 1) * EG(1, 0) ;


// opserr << "Finished shellGeo!" << endln;
  // End shell geo function
}


void  IGAKLShell_BendingStrip::zeroLoad( )
{
  // opserr << "IGAKLShell_BendingStrip::zeroLoad - eleTag" << this->getTag() << " called! " << endln;
  if (load != 0)
  {
    load->Zero(); // Uncomment this for twisted shell
  }

  applyLoad = 0;

  appliedB[0] = 0.0;
  appliedB[1] = 0.0;
  appliedB[2] = 0.0;

  return ;
}


int
IGAKLShell_BendingStrip::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  // opserr << "IGAKLShell_BendingStrip::addLoad - tag = " << this->getTag() << endln;
  int type;
  const Vector &data = theLoad->getData(type, loadFactor);

  if (type == LOAD_TAG_SelfWeight) {
    // opserr << "IGAKLShell_BendingStrip::addLoad - type = LOAD_TAG_SelfWeight - tag = " << this->getTag() << endln;
    // added compatability with selfWeight class implemented for all continuum elements, C.McGann, U.W.
    applyLoad = 1;
    appliedB[0] += loadFactor * data(0);
    appliedB[1] += loadFactor * data(1);
    appliedB[2] += loadFactor * data(2);
    return 0;
  } else if (type == LOAD_TAG_SurfaceLoader) {
    pressure = data(0);
    return 0;
  } else if (type == LOAD_TAG_IGAFollowerLoad) {
    // opserr << "IGAKLShell_BendingStrip::addLoad - type = LOAD_TAG_IGAFollowerLoad - tag = " << this->getTag() << endln;

    double xi = data(0);
    double eta = data(1);

    int Nnodes = connectedExternalNodes.Size();
    int NDOF = 3 * Nnodes;

    if (load == 0)
    {
      load = new Vector(NDOF);
    }


    if (pointInElement(xi, eta))
    {
      // opserr << "IGAKLShell_BendingStrip::addLoad - type = LOAD_TAG_IGAFollowerLoad -tag =  = " << this->getTag() << " called" << endln;
      Vector followerforce(3);

      followerforce(0) = data(2);
      followerforce(1) = data(3);
      followerforce(2) = data(4);

      // Obtain information from current configuration
      int noFuncs = myPatch->getNoFuncs();
      Vector R(noFuncs);
      Vector dRdxi(noFuncs);
      Vector dRdeta(noFuncs);
      Vector dR2dxi(noFuncs);
      Vector dR2deta(noFuncs);
      Vector dR2dxideta(noFuncs);
      R.Zero();
      dRdxi.Zero();
      dRdeta.Zero();
      dR2dxi.Zero();
      dR2deta.Zero();
      dR2dxideta.Zero();

      opserr << "xi = " << xi << endln;
      opserr << "eta = " << eta << endln;
      // if (eta == 0)
      // {
      //   eta = 1;
      // }
      // myPatch->Nurbs2DBasis2ndDers(xi, eta, R, dRdxi, dRdeta, dR2dxi, dR2deta, dR2dxideta);
      myPatch->Nurbs2DBasis2ndDers(xi, 1.0, R, dRdxi, dRdeta, dR2dxi, dR2deta, dR2dxideta); // I need this for twisted shell
      // myPatch->Nurbs2DBasis2ndDers(xi, eta, R, dRdxi, dRdeta, dR2dxi, dR2deta, dR2dxideta); // I need this for bending plate circle

      // dRdxi /= 16; //I need this for twisted shell
      // dRdeta /= 2; //I need this for twisted shell

      // opserr << "R = " << R << endln;
      // opserr << "dRdxi = " << dRdxi << endln;
      // opserr << "dRdeta = " << dRdeta << endln;

      // myPatch->Nurbs2DBasis2ndDers(1.0, 1.0, R, dRdxi, dRdeta, dR2dxi, dR2deta, dR2dxideta); // Hardcoding this to get the normal on the midplane

      opserr << "R = " << R << endln;
      opserr << "dRdxi = " << dRdxi << endln;
      opserr << "dRdeta = " << dRdeta << endln;

      // myPatch->Nurbs2DBasis2ndDers(xi, eta, R, dRdxi, dRdeta, dR2dxi, dR2deta, dR2dxideta);

      // double wt;
      // double ptU;
      // double ptV;
      // double xi;
      // double eta;
      // double J2;

      // opserr << "R = " << R << endln;
      // opserr << "dRdxi = " << dRdxi << endln;

      // for (int gp = 0; gp < ngauss; ++gp)
      // {
      //   wt = (*quadWeight)(gp);
      //   ptU = (*quadPoint)(gp, 0);
      //   ptV = (*quadPoint)(gp, 1);
      //   xi = myPatch->parent2ParametricSpace(xiE, ptU);
      //   eta = myPatch->parent2ParametricSpace(etaE, ptV);
      //   J2 = 0.5 * (xiE(1) - xiE(0)) * 0.5 * (etaE(1) - etaE(0));

      //   opserr << "wt = " << wt << endln;
      //   opserr << "ptU = " << ptU << endln;
      //   opserr << "ptV = " << ptV << endln;

      //   R.Zero();
      //   dRdxi.Zero();
      //   dRdeta.Zero();
      //   dR2dxi.Zero();
      //   dR2deta.Zero();
      //   dR2dxideta.Zero();
      //   myPatch->Nurbs2DBasis2ndDers(xi, eta, R, dRdxi, dRdeta, dR2dxi, dR2deta, dR2dxideta);
      // }

      Vector point(3);
      Vector point_disp(3);
      Matrix pts(connectedExternalNodes.Size(), 3);
      Matrix pts_d(connectedExternalNodes.Size(), 3);

      bool nonLinearGeometry = myPatch->getAnalysisType();

      // Get pts span matrix for jacobian
      for (int i = 0; i < connectedExternalNodes.Size(); ++i)
      {
        point = nodePointers[i]->getCrds();
        point_disp = nodePointers[i]->getTrialDisp(); //getTrialDisp
        for (int j = 0; j < 3; ++j)
        {
          pts(i, j) = point(j);
          if (isnan(point_disp(j)))
          {
            opserr << "Nan found on addLoad = " << endln;
          }
          if (nonLinearGeometry)
          {
            pts_d(i, j) = pts(i, j) + point_disp(j);
          }
          else
          {
            pts_d(i, j) = pts(i, j);
          }
        }
      }

      // Get first and second order jacobians
      Matrix dr(2, noFuncs);
      Matrix dr2(3, noFuncs);


      for (int i = 0; i < noFuncs; ++i)
      {
        dr(0, i) = dRdxi(i);
        dr(1, i) = dRdeta(i);

        dr2(0, i) = dR2dxi(i);
        dr2(1, i) = dR2deta(i);
        dr2(2, i) = dR2dxideta(i);
      }

      Matrix ddR(3, noFuncs);
      for (int i = 0; i < noFuncs; ++i)
      {
        ddR(0, i) = dR2dxi(i);
        ddR(1, i) = dR2deta(i);
        ddR(2, i) = dR2dxideta(i);
      }


      Matrix g = dr  * pts_d; //G: Covariant base vectors, needed in the current configuration
      Matrix h = ddR * pts_d; //H: Hessian matrix, needed in the current configuration

      // Transposing because it is needed in cols drdxi, drdeta
      g = transpose(2, 3, g);
      h = transpose(3, 3, h); // Shouldn't need this


      // Shell geo deformed parameters and temporal
      Vector g3(3); // = g3
      Vector n(3); // = n
      double lg3_tmp; // = dA
      Matrix gab_tmp(2, 2); // = Gab
      Vector bv_tmp(3); // = Bv
      Matrix T_Gcon_E(3, 3);
      Matrix T_E_G_tmp(3, 3);
      Matrix T_G_E(3, 3);


      // Get T_G_E for base transformation and g3
      shellGeo(g , h , g3     , lg3_tmp , n     , gab_tmp   , bv_tmp   , T_Gcon_E, T_E_G_tmp, T_G_E);

      // Local covariant base vectors
      Vector g1(3);
      Vector g2(3);
      for (int i = 0; i < 3; ++i)
      {
        g1(i) = g(i, 0);
        g2(i) = g(i, 1);
      }

      // opserr << "g3 = " << g3 << endln;
      // To transform from local covariant g to local cartesian e
      Vector e1 = T_G_E * g1;
      Vector e2 = T_G_E * g2;

      Vector e3 = n;

      // for (int i = 0; i < 3; ++i)
      // {
      //   if (e3(i)<0)
      //   {
      //     e3(i)*=-1;
      //   }
      // }

      // Vector e1 = T_Gcon_E * g1; // Trying
      // Vector e2 = T_Gcon_E * g2;

      // Normalizing vectors
      e1 /= e1.Norm();
      e2 /= e2.Norm();
      e3 /= e3.Norm();

      opserr << "e1 = " << e1 << endln;
      opserr << "e2 = " << e2 << endln;
      opserr << "e3 = " << e3 << endln;
      opserr << "data(2) = " << data(2) << endln;
      opserr << "data(3) = " << data(3) << endln;
      opserr << "data(4) = " << data(4) << endln;

      // followerforce = loadFactor * (data(2) * e1 + data(3) * e2 + data(4) * e3);
      followerforce = (data(2) * e1 + data(3) * e2 + data(4) * e3);


      // followerforce(0) = 0;

      // opserr << "e3 = " << e3 << endln;

      // opserr << "followerforce = " << followerforce << endln;

      // N debe ser evaluado en xi, eta
      myPatch->Nurbs2DBasis2ndDers(xi, eta, R, dRdxi, dRdeta, dR2dxi, dR2deta, dR2dxideta);
      Matrix N(3, 3 * noFuncs);
      for (int i = 0; i < noFuncs; ++i)
      {
        N(0, 3 * i) = R(i);
        N(1, 3 * i + 1) = R(i);
        N(2, 3 * i + 2) = R(i);
      }

      // opserr << "old load = " << *load << endln;

      // Error found here
      load->addMatrixTransposeVector(1.0, N, followerforce, 1.0);

      opserr << "load = " << *load << endln;

      opserr << "e3 = " << e3 << endln;
      opserr << "followerforce = " << followerforce << endln;
      opserr << "loadFactor = " << loadFactor << endln;


      // resid->addMatrixTransposeVector(1.0, N, followerforce, -1.0);

      // opserr << "*load = " << *load << endln;

      // opserr << "*load = " << *load << endln;
      // The error is that the shape functions are not interpolatory, so the above multiplication does not apply the force in the desired DOF, as the shape functions are greater elsewhere.

      // opserr << "xi = " << xi << endln;
      // opserr << "eta = " << eta << endln;
      // opserr << "xiE = " << xiE << endln;
      // opserr << "etaE = " << etaE << endln;
      // opserr << "R = " << R << endln;
      // opserr << "-1*followerforce = " << -1 * followerforce << endln;
      // opserr << "*load = " << *load << endln;
    }

    return 0;

  } else {
    opserr << "ShellMITC4::addLoad() - ele with tag: " << this->getTag() << " does not deal with load type: " << type << "\n";
    return -1;
  }
}



int
IGAKLShell_BendingStrip::addInertiaLoadToUnbalance(const Vector &accel)
{
  // opserr << "IGAKLShell_BendingStrip::addInertiaLoadToUnbalance - eleTag" << this->getTag() << " called! " << endln;
  return 0;
}



//get residual
const Vector&  IGAKLShell_BendingStrip::getResistingForce() 
{
  // opserr << "IGAKLShell::getResistingForce - eleTag" << this->getTag() << " called! " << endln;

  int tang_flag = 1 ; //get the tangent
  bool nonLinearGeometry = myPatch->getAnalysisType();

  formResidAndTangent( tang_flag ) ;

  // Adding K*NodalDisplacements to resid

  int Nnodes = connectedExternalNodes.Size();
  int NDOF = 3 * Nnodes;
  static Vector res(NDOF);   /// ngdl
  res.resize(NDOF);
  res.Zero();

  // res = this->getResistingForce();

  // Calculo desplazamientos en nodos
  static Vector NodalDisplacements(NDOF);
  NodalDisplacements.resize(NDOF);
  NodalDisplacements.Zero();
  static Vector disp_i(3);
  disp_i.Zero();


  if (nonLinearGeometry == false)
  {
    // Loop through nodes
    int k = 0;
    for (int i = 0; i < Nnodes; ++i)
    {
      Node *node_i = nodePointers[i];
      disp_i = node_i->getTrialDisp();
      for (int j = 0; j < 3; ++j)
      {
        NodalDisplacements(k) = disp_i(j);
        k += 1;
      }
    }

    // Adding K*nodalDisplacements;
    // opserr << "Adding K*nodalDisplacements" << endln;
    // *resid += 1*(*stiff) * NodalDisplacements;
    // *resid += 1*(this->getTangentStiff()) * NodalDisplacements;
    *resid += (*stiff) * NodalDisplacements;

  }

  // *resid += this->getMass() * NodalAccelerations;


  // subtract external loads
  if (load != 0)
  {
    *resid -= *load;
    // opserr << "this->getTag() = " << this->getTag() << endln;
    // opserr << "*resid = " << *resid << endln;
  }

  // opserr << "*resid = " << *resid << endln;

  return *resid ;

  // resid->Zero();
  // // opserr << "connectedExternalNodes = " << connectedExternalNodes << endln;

  // // opserr << "Element number:  = " << this->getTag() << endln;

  // // opserr << "xiE = " << xiE << endln;
  // // opserr << "etaE = " << etaE << endln;

  // Vector gFact = myPatch->getGravityFactors();
  // // opserr << "gFact = " << gFact << endln;


  // double t = 0.0; // thickness [m]
  // double rho = 0.0; //kg/m^3

  // double W = 0.0; //Total laminate weight (kg/m^2)
  // int nLayers = myPatch->getNLayers();

  // for (int capa = 0; capa < nLayers; ++capa)
  // {
  //   rho = (OPS_getNDMaterial(myPatch->getMatTag(capa)))->getRho(); // Density of material
  //   t = myPatch->getThickness(capa); // Thickness
  //   W += rho * t;
  // }



  // bool nonLinearGeometry = myPatch->getAnalysisType();


  // float wt;
  // float ptU;
  // float ptV;
  // double xi;
  // double eta;

  // int noFuncs = myPatch->getNoFuncs();


  // Vector R(noFuncs);
  // Vector dRdxi(noFuncs);
  // Vector dRdeta(noFuncs);
  // Vector dR2dxi(noFuncs);
  // Vector dR2deta(noFuncs);
  // Vector dR2dxideta(noFuncs);

  // double J2; //Jacobian for the parametric mapping

  // Vector point(3);
  // Vector point_disp(3);
  // Matrix pts(connectedExternalNodes.Size(), 3);
  // Matrix pts_d(connectedExternalNodes.Size(), 3);
  // // Get pts span matrix for jacobian
  // for (int i = 0; i < connectedExternalNodes.Size(); ++i)
  // {
  //   point = nodePointers[i]->getCrds();
  //   point_disp = nodePointers[i]->getTrialDisp(); //getTrialDisp
  //   for (int j = 0; j < 3; ++j)
  //   {
  //     pts(i, j) = point(j);
  //     if (isnan(point_disp(j)))
  //     {
  //       opserr << "Nan found on formResidAndTangent = " << endln;
  //     }
  //     if (nonLinearGeometry)
  //     {
  //       pts_d(i, j) = pts(i, j) + point_disp(j);
  //     }
  //     else
  //     {
  //       pts_d(i, j) = pts(i, j);
  //     }
  //   }

  // }

  // // opserr << "pts = " << pts << endln;

  // // Loop over integrations points
  // for (int gp = 0; gp < ngauss; ++gp)
  // {
  //   wt = (*quadWeight)(gp);
  //   ptU = (*quadPoint)(gp, 0);
  //   ptV = (*quadPoint)(gp, 1);
  //   xi = myPatch->parent2ParametricSpace(xiE, ptU);
  //   eta = myPatch->parent2ParametricSpace(etaE, ptV);
  //   J2 = 0.5 * (xiE(1) - xiE(0)) * 0.5 * (etaE(1) - etaE(0));
  //   // opserr << "xi = " << xi << endln;
  //   // opserr << "eta = " << eta << endln;
  //   R.Zero();
  //   dRdxi.Zero();
  //   dRdeta.Zero();
  //   dR2dxi.Zero();
  //   dR2deta.Zero();
  //   dR2dxideta.Zero();
  //   myPatch->Nurbs2DBasis2ndDers(xi, eta, R, dRdxi, dRdeta, dR2dxi, dR2deta, dR2dxideta);
  //   // opserr << "Got derivatives! " <<  endln;

  //   // Get first and second order jacobians
  //   Matrix dr(2, noFuncs);
  //   Matrix dr2(3, noFuncs);


  //   for (int i = 0; i < noFuncs; ++i)
  //   {
  //     dr(0, i) = dRdxi(i);
  //     dr(1, i) = dRdeta(i);

  //     dr2(0, i) = dR2dxi(i);
  //     dr2(1, i) = dR2deta(i);
  //     dr2(2, i) = dR2dxideta(i);
  //   }

  //   //compute the jacobian of physical and parameter domain mapping
  //   // then the derivative w.r.t spatial physical coordinates, current configuration


  //   Matrix jacob = dr * pts;  //G: Covariant base vectors, needed in the reference configuration
  //   Matrix jacob2 = dr2 * pts; //H: Hessian matrix, needed in the reference configuration


  //   // !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
  //   // !    DESCRIPTION OF THE VARIABLES !
  //   // !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

  //   // ! g: matrix with the 3 components of the in plane vectors g1, g2
  //   // ! dR: first derivatives of the shape functions
  //   // ! ddR: second derivatives of the shape functions
  //   // ! dn: first variation of the normal vector g3
  //   // ! lg3: length of the normal vector g3

  //   // ! E_cu: strain vector [E11,E22,E12] referred to the curvilinear coord. sys
  //   // ! E_ca: strain vector [E11,E22,2*E12] referred to cartesian coord. sys
  //   // ! K_cu: curvature vector [K11,K22,K12] referred to the curvilinear coord. sys
  //   // ! K_ca: curvature vector [K11,K22,2*K12] referred to the cartesian coord. sys

  //   // ! dE_cu: first variation of the membrane strain in local coord. system
  //   // ! dE_ca: first variation of the membrane strain in global coord. system
  //   // ! dK_cu: first variation of the curvature in local coord. system
  //   // ! dK_ca: first variation of the curvature in global coord. system

  //   // ! ddE_cu: second variation of the membrane strain in local coord. system
  //   // ! ddE_ca: second variation of the membrane strain in global coord. system
  //   // ! ddK_cu: second variation of the curvature in local coord. system
  //   // ! ddK_ca: second variation of the curvature in global coord. system

  //   Vector E_cu(3);
  //   Vector K_cu(3);

  //   Vector dE_cu(3);
  //   Matrix dE_ca(3, 3 * noFuncs);
  //   Vector dK_cu(3);
  //   Matrix dK_ca(3, 3 * noFuncs);

  //   Vector fiem(3 * noFuncs);
  //   Vector fieb(3 * noFuncs);
  //   Vector fie(3 * noFuncs);

  //   Matrix ddR(3, noFuncs);
  //   for (int i = 0; i < noFuncs; ++i)
  //   {
  //     ddR(0, i) = dR2dxi(i);
  //     ddR(1, i) = dR2deta(i);
  //     ddR(2, i) = dR2dxideta(i);
  //   }


  //   Matrix G = dr  * pts;
  //   Matrix H = ddR * pts;

  //   Matrix g = dr  * pts_d;
  //   Matrix h = ddR * pts_d;


  //   G = transpose(2, 3, G); // Transposing because it is needed in cols drdxi, drdeta
  //   g = transpose(2, 3, g);
  //   H = transpose(3, 3, H);
  //   h = transpose(3, 3, h);


  //   // Shell geo reference parameters and temporal
  //   Vector g3_tmp(3);
  //   Vector n_tmp(3);
  //   double dA; // = lG3
  //   Matrix Gab_r(2, 2); // = Gab
  //   Vector Bv_r(3); // = Bv
  //   Matrix Tb(3, 3); // = T_Gcon_E
  //   Matrix T_E_G_tmp(3, 3); // = T_Gcon_E
  //   Matrix T_G_E_tmp(3, 3); // = T_Gcon_E

  //   // Shell geo deformed parameters and temporal
  //   Vector g3(3); // = g3
  //   Vector n(3); // = n
  //   double lg3; // = dA
  //   Matrix gab(2, 2); // = Gab
  //   Vector bv(3); // = Bv
  //   Matrix Tb_tmp(3, 3);

  //   // Call shell geometry functions
  //   shellGeo(G , H , g3_tmp , dA  , n_tmp , Gab_r , Bv_r , Tb, T_E_G_tmp, T_G_E_tmp);
  //   shellGeo(g , h , g3     , lg3 , n     , gab   , bv   , Tb_tmp, T_E_G_tmp, T_G_E_tmp);

  //   //
  //   // Compute the membrane strains and the curvatures according to the eqns 1.161 - 1.162 of Computational FSI (Yuri's book)
  //   //

  //   //strain vector [E11,E22,E12] referred to curvilinear coor sys
  //   E_cu(0) = 0.5 * (gab(0, 0) - Gab_r(0, 0));
  //   E_cu(1) = 0.5 * (gab(1, 1) - Gab_r(1, 1));
  //   E_cu(2) = 0.5 * (gab(0, 1) - Gab_r(0, 1));

  //   //curvature vector [K11,K22,K12] referred to curvilinear coor sys
  //   K_cu(0) = -(bv(0) - Bv_r(0));
  //   K_cu(1) = -(bv(1) - Bv_r(1));
  //   K_cu(2) = -(bv(2) - Bv_r(2));


  //   Vector E_ca = Tb * E_cu; // strain vector [E11,E22,2*E12] referred to cartesian coor sys
  //   Vector K_ca = Tb * K_cu; //curvature vector [K11,K22,2*K12] referred to cartesian coor sys

  //   //
  //   // Compute the First variation of strain and curvature w.r.t the d.o.f.
  //   //

  //   // Local variables
  //   int kr = 0;
  //   int dirr = 0;
  //   Matrix dg(3, 2);
  //   Matrix dR = transpose(2, noFuncs, dr);
  //   ddR = transpose(3, noFuncs, ddR);
  //   Matrix dg3(3, 3 * noFuncs);
  //   Vector g3dg3(3 * noFuncs);
  //   Vector g3dg3lg3_3(3 * noFuncs);
  //   Matrix dn(3, 3 * noFuncs);

  //   // For "slicing" matrices
  //   ID ur_only(1);
  //   ID us_only(1);
  //   auto xyz_temp = arange<int>(0, 3); ID xyz(&(xyz_temp[0]), 3); // Takes value 0 1 2

  //   double lg3_3 = pow(lg3, 3);



  //   for (int ur = 0; ur < 3 * noFuncs; ++ur)
  //   {
  //     // Local node number kr and dof direction dirr
  //     kr = (ur + 1 + 2) / 3 - 1; // takes value 0 0 0, 1 1 1, 2 2 2, .... noFuncs-1 noFuncs-1 noFuncs-1
  //     dirr = ur - 3 * (kr); // Takes value 0 1 2 0 1 2 ...

  //     // opserr << "kr = " << kr << endln;
  //     // opserr << "dirr = " << dirr << endln;

  //     dg(dirr, 0) = dR(kr, 0);
  //     dg(dirr, 1) = dR(kr, 1);

  //     // Strain
  //     dE_cu(0) = dR(kr, 0) * g(dirr, 0);
  //     dE_cu(1) = dR(kr, 1) * g(dirr, 1);
  //     dE_cu(2) = 0.5 * (dR(kr, 0) * g(dirr, 1) + g(dirr, 0) * dR(kr, 1));

  //     dE_ca(0, ur) = (Tb * dE_cu)(0);
  //     dE_ca(1, ur) = (Tb * dE_cu)(1);
  //     dE_ca(2, ur) = (Tb * dE_cu)(2);


  //     // Curvature
  //     dg3(0, ur) = dg(1, 0) * g(2, 1) - dg(2, 0) * g(1, 1) + g(1, 0) * dg(2, 1) - g(2, 0) * dg(1, 1);
  //     dg3(1, ur) = dg(2, 0) * g(0, 1) - dg(0, 0) * g(2, 1) + g(2, 0) * dg(0, 1) - g(0, 0) * dg(2, 1);
  //     dg3(2, ur) = dg(0, 0) * g(1, 1) - dg(1, 0) * g(0, 1) + g(0, 0) * dg(1, 1) - g(1, 0) * dg(0, 1);

  //     g3dg3(ur) = g3(0) * dg3(0, ur) + g3(1) * dg3(1, ur) + g3(2) * dg3(2, ur);
  //     g3dg3lg3_3(ur) = g3dg3(ur) / lg3_3;


  //     dn(0, ur) = dg3(0, ur) / lg3 - g3(0) * g3dg3lg3_3(ur);
  //     dn(1, ur) = dg3(1, ur) / lg3 - g3(1) * g3dg3lg3_3(ur);
  //     dn(2, ur) = dg3(2, ur) / lg3 - g3(2) * g3dg3lg3_3(ur);

  //     dK_cu(0) = -(ddR(kr, 0) * n(dirr) + (h(0, 0) * dn(0, ur) + h(1, 0) * dn(1, ur) + h(2, 0) * dn(2, ur)));
  //     dK_cu(1) = -(ddR(kr, 1) * n(dirr) + (h(0, 1) * dn(0, ur) + h(1, 1) * dn(1, ur) + h(2, 1) * dn(2, ur)));
  //     dK_cu(2) = -(ddR(kr, 2) * n(dirr) + (h(0, 2) * dn(0, ur) + h(1, 2) * dn(1, ur) + h(2, 2) * dn(2, ur)));

  //     // CHECK THIS OUT, I HAD TO SWAP THE INDICES FOR H

  //     // dK_cu(0) = -(ddR(0, kr) * n(dirr) + (h(0, 0) * dn(0, ur) + h(0, 1) * dn(1, ur) + h(0, 2) * dn(2, ur)));
  //     // dK_cu(1) = -(ddR(1, kr) * n(dirr) + (h(1, 0) * dn(0, ur) + h(1, 1) * dn(1, ur) + h(1, 2) * dn(2, ur)));
  //     // dK_cu(2) = -(ddR(2, kr) * n(dirr) + (h(2, 0) * dn(0, ur) + h(2, 1) * dn(1, ur) + h(2, 2) * dn(2, ur)));

  //     dK_ca(0, ur) = (Tb * dK_cu)(0);
  //     dK_ca(1, ur) = (Tb * dK_cu)(1);
  //     dK_ca(2, ur) = (Tb * dK_cu)(2);
  //   }

  //   int nLayers = myPatch->getNLayers();

  //   Matrix A(3, 3);   // Membrane stiffness matrix
  //   Matrix B(3, 3);   // Coupling stifness matrix
  //   Matrix D(3, 3);   // Bending stiffness matrix
  //   Matrix T(3, 3);   // Rotating laminate matrix
  //   Matrix Cbar(3, 3); // Equivalent constitutive matrix
  //   for (int capa = 0; capa < nLayers; ++capa)
  //   {
  //     double iAngle     = myPatch->getAngle(capa);
  //     double iThickness = myPatch->getThickness(capa);
  //     double iZ         = myPatch->getZk(capa);// Mid center laminate position, 1 in the meantime

  //     const Matrix& Czig = materialPointers[gp][capa]->getTangent();

  //     T(0, 0) = pow(cos(iAngle), 2)                        ; T(0, 1) = pow(sin(iAngle), 2)           ; T(0, 2) = sin(iAngle) * cos(iAngle)               ;
  //     T(1, 0) = pow(sin(iAngle), 2)                        ; T(1, 1) = pow(cos(iAngle), 2)           ; T(1, 2) = -sin(iAngle) * cos(iAngle)              ;
  //     T(2, 0) = -2 * sin(iAngle) * cos(iAngle)             ; T(2, 1) = 2 * sin(iAngle) * cos(iAngle) ; T(2, 2) = pow(cos(iAngle), 2) - pow(sin(iAngle), 2) ;

  //     Cbar.addMatrixTripleProduct(0, T, Czig, 1.0);

  //     A.addMatrix(1.0, Cbar, iThickness);
  //     B.addMatrix(1.0, Cbar, iThickness * iZ);
  //     D.addMatrix(1.0, Cbar, (iThickness * pow(iZ, 2) + pow(iThickness, 3) / 12));
  //   }

  //   Vector N_ca = A * E_ca + B * K_ca; // Membrane forces
  //   Vector M_ca = B * E_ca + D * K_ca; // Bending moments


  //   for (int ur = 0; ur < 3 * noFuncs ; ++ur)
  //   {
  //     // Local node number kr and dof direction dirr
  //     kr = (ur + 1 + 2) / 3 - 1; // takes value 0 0 0, 1 1 1, 2 2 2, .... noFuncs-1 noFuncs-1 noFuncs-1
  //     dirr = ur - 3 * (kr); // Takes value 0 1 2 0 1 2 ...

  //     ur_only.fill(ur);

  //     Matrix dN_ca = A * dE_ca(xyz, ur_only) + B * dK_ca(xyz, ur_only);
  //     Matrix dM_ca = B * dE_ca(xyz, ur_only) + D * dK_ca(xyz, ur_only);

  //     Vector dE_ca_here(3);
  //     dE_ca_here(0) = dE_ca(xyz, ur_only)(0, 0);
  //     dE_ca_here(1) = dE_ca(xyz, ur_only)(1, 0);
  //     dE_ca_here(2) = dE_ca(xyz, ur_only)(2, 0);

  //     Vector dK_ca_here(3);
  //     dK_ca_here(0) = dK_ca(xyz, ur_only)(0, 0);
  //     dK_ca_here(1) = dK_ca(xyz, ur_only)(1, 0);
  //     dK_ca_here(2) = dK_ca(xyz, ur_only)(2, 0);

  //     // Force vectors
  //     fiem(ur) = (N_ca ^ dE_ca_here);
  //     fieb(ur) = (M_ca ^ dK_ca_here);

  //   }

  //   Matrix N(3, 3 * noFuncs);
  //   Vector pressure_vector = pressure * g3;

  //   for (int i = 0; i < noFuncs; ++i)
  //   {
  //     N(0, 3 * i) = R(i);
  //     N(1, 3 * i + 1) = R(i);
  //     N(2, 3 * i + 2) = R(i);
  //   }


  //   fie = (fiem + fieb) * (dA * J2 * wt); //Residual vector on gauss point

  //   fie.addMatrixTransposeVector(1.0, N, pressure_vector, dA * J2 * wt);

  //   // Body forces
  //   Vector b(appliedB, 3);

  //   if (applyLoad == 1)
  //   {
  //     fie.addMatrixTransposeVector(1.0, N, b, W * dA * J2 * wt);
  //   } else {
  //     fie.addMatrixTransposeVector(1.0, N, gFact, W * dA * J2 * wt);
  //   }

  //   *resid += fie; // Kiendl stiffness


  // }

  // // Aqui se rellena el vector de fuerzas residuales...

  // // // subtract external loads
  // // if (load != 0)
  // //   *resid -= *load;

  // // Adding follower load vector to resid (substracting actually to make it positive)
  // // opserr << "*resid = " << *resid << endln;
  // // opserr << "*load = " << -1*(*load) << endln;
  // *resid-=(*load);
  // // load->Zero();

  // return *resid ;


}

void IGAKLShell_BendingStrip::formResidAndTangent( int tang_flag )
{
  //zero stiffness and residual
  stiff->Zero();
  resid->Zero();

  bool nonLinearGeometry = myPatch->getAnalysisType();

  Vector gFact = myPatch->getGravityFactors();

  // Gauss integration parameters
  float wt;
  float ptU;
  float ptV;
  double xi;
  double eta;

  // Number of shape functions and layers
  int noFuncs = myPatch->getNoFuncs();
  int nLayers = myPatch->getNLayers();




  // Start of static declarations
  static Matrix ke(3 * noFuncs, 3 * noFuncs);
  static Matrix kem(3 * noFuncs, 3 * noFuncs);
  static Matrix keb(3 * noFuncs, 3 * noFuncs);

  static Vector E_cu(3);
  static Vector E_ca(3);
  static Vector K_cu(3);
  static Vector K_ca(3);
  static Vector dE_cu(3);
  static Matrix dE_ca(3, 3 * noFuncs);
  static Vector dK_cu(3);
  static Matrix dK_ca(3, 3 * noFuncs);
  static Vector N_ca(3);
  static Vector M_ca(3);
  static Vector fiem(3 * noFuncs);
  static Vector fieb(3 * noFuncs);
  static Vector fie(3 * noFuncs);

  static Matrix pts(noFuncs, 3); // Matrix for storing element nodes coordinates, noFuncs,3
  static Matrix pts_d(noFuncs, 3); // Matrix for storing element nodes displaced coordinates

  static Matrix G_t(2, 3);
  static Matrix H_t(3, 3);
  static Matrix g_t(2, 3);
  static Matrix h_t(3, 3);
  // Shell geo reference parameters and temporal
  static Vector g3_tmp(3);
  static Vector n_tmp(3);
  static double dA = 0; // = lG3
  static Matrix Gab_r(2, 2);// = Gab
  static Vector Bv_r(3); // = Bv
  static Matrix Tb(3, 3); // = T_Gcon_E
  static Matrix T_E_G_tmp(3, 3); // = T_E_G
  static Matrix T_G_E_tmp(3, 3); // = T_G_E
  static Matrix T_G_E(3, 3); // = T_G_E for testing
  // Shell geo deformed parameters and temporal
  static Vector g3(3); // = g3
  static Vector n(3); // = n
  static double lg3 = 0; // = dA
  static Matrix gab(2, 2); // = Gab
  static Vector bv(3); // = Bv
  static Matrix Tb_tmp(3, 3);

  static Matrix dR(noFuncs, 3);
  static Matrix ddR(noFuncs, 3);

  static Matrix dr(3, noFuncs);
  static Matrix ddr(3, noFuncs);

    // Getting ply matrices

  static Matrix A(3, 3);   // Membrane stiffness matrix
  static Matrix B(3, 3);   // Coupling stifness matrix
  static Matrix D(3, 3);   // Bending stiffness matrix
  static Matrix T(3, 3);   // Rotating laminate matrix
  static Matrix C(3, 3); // Equivalent constitutive matrix
  static double W; //Total laminate weight (kg/m^2)
  static Vector Ezeta(3);
  static double iAngle;    
  static double iThickness;
  static double iZ;        
  static double iRho;      

  // Variables for stiffness assembly
  static double lg3_3;
  static double lg3_5;
  static int kr;
  static int dirr;
  static Matrix dg(3, 2);
  static Matrix dg3(3, 3 * noFuncs);
  static Vector g3dg3(3 * noFuncs);
  static Vector g3dg3lg3_3(3 * noFuncs);
  static Matrix dn(3, 3 * noFuncs);

  static Vector ddE_cu(3);
  static Matrix ddE_ca_0(3 * noFuncs, 3 * noFuncs); // Wanted to avoid declaring a 3D array, instead using 3 matrices
  static Matrix ddE_ca_1(3 * noFuncs, 3 * noFuncs);
  static Matrix ddE_ca_2(3 * noFuncs, 3 * noFuncs);

  static Vector ddK_cu(3);
  static Matrix ddK_ca_0(3 * noFuncs, 3 * noFuncs); // Wanted to avoid declaring a 3D array, instead using 3 matrices
  static Matrix ddK_ca_1(3 * noFuncs, 3 * noFuncs);
  static Matrix ddK_ca_2(3 * noFuncs, 3 * noFuncs);

  static int ks;
  static int dirs;
  static int dirt;
  static int ddir;
  static double tmp1;
  static double tmp2;
  static Vector ddn(3);
  static Vector ddg3(3);

  static Vector Szeta(3);



  static Vector dE_ca_ur(3);
  static Vector dK_ca_ur(3);
  static Vector dN_ca(3);
  static Vector dM_ca(3);
  static Vector ddE_ca(3);
  static Vector ddK_ca(3);
  static Vector dK_ca_us(3);
  static Vector dE_ca_us(3);


  static Vector R(noFuncs);  // 
  static Vector dRdxi(noFuncs);
  static Vector dRdeta(noFuncs);
  static Vector dR2dxi(noFuncs);
  static Vector dR2deta(noFuncs);
  static Vector dR2dxideta(noFuncs);

  //Jacobian for the parametric mapping
  static double J2;

  // Start of resizing static Vectors and Matrices
  ke.resize(3 * noFuncs, 3 * noFuncs);
  kem.resize(3 * noFuncs, 3 * noFuncs);
  keb.resize(3 * noFuncs, 3 * noFuncs);
  dE_ca.resize(3, 3 * noFuncs);
  dK_ca.resize(3, 3 * noFuncs);
  fiem.resize(3 * noFuncs);
  fieb.resize(3 * noFuncs);
  fie.resize(3 * noFuncs);
  pts.resize(noFuncs, 3); // Matrix for storing element nodes coordinates, noFuncs,3
  pts_d.resize(noFuncs, 3); // Matrix for storing element nodes displaced coordinates
  dR.resize(noFuncs, 3);
  ddR.resize(noFuncs, 3);
  dr.resize(3, noFuncs);
  ddr.resize(3, noFuncs);
  dg3.resize(3, 3 * noFuncs);
  g3dg3.resize(3 * noFuncs);
  g3dg3lg3_3.resize(3 * noFuncs);
  dn.resize(3, 3 * noFuncs);
  ddE_ca_0.resize(3 * noFuncs, 3 * noFuncs); // Wanted to avoid declaring a 3D array, instead using 3 matrices
  ddE_ca_1.resize(3 * noFuncs, 3 * noFuncs);
  ddE_ca_2.resize(3 * noFuncs, 3 * noFuncs);
  ddK_ca_0.resize(3 * noFuncs, 3 * noFuncs); // Wanted to avoid declaring a 3D array, instead using 3 matrices
  ddK_ca_1.resize(3 * noFuncs, 3 * noFuncs);
  ddK_ca_2.resize(3 * noFuncs, 3 * noFuncs);
  R.resize(noFuncs);  // 
  dRdxi.resize(noFuncs);
  dRdeta.resize(noFuncs);
  dR2dxi.resize(noFuncs);
  dR2deta.resize(noFuncs);
  dR2dxideta.resize(noFuncs);




  // START OF CALCULATION
  pts.Zero();
  pts_d.Zero();


  // Get pts span matrix for jacobian
  for (int i = 0; i < connectedExternalNodes.Size(); ++i)
  {
    const Vector& point  = nodePointers[i]->getCrds() ; // Vector for storing node coordinates
    const Vector& point_disp = nodePointers[i]->getTrialDisp(); // Vector for storing node displacement
    for (int j = 0; j < 3; ++j)
    {
      pts(i, j) = point(j);
      if (isnan(point_disp(j)))
      {
        opserr << "Nan found on formResidAndTangent = " << endln;
      }
      if (nonLinearGeometry)
      {
        pts_d(i, j) = pts(i, j) + point_disp(j);
      }
      else
      {
        pts_d(i, j) = pts(i, j);
      }
    }
  }


  // Loop over integrations points
  for (int gp = 0; gp < ngauss; ++gp)
  {
    wt = (*quadWeight)(gp);
    ptU = (*quadPoint)(gp, 0);
    ptV = (*quadPoint)(gp, 1);

    xi = myPatch->parent2ParametricSpace(xiE, ptU);
    eta = myPatch->parent2ParametricSpace(etaE, ptV);
    J2 = 0.5 * (xiE(1) - xiE(0)) * 0.5 * (etaE(1) - etaE(0));

    R.Zero();
    dRdxi.Zero();
    dRdeta.Zero();
    dR2dxi.Zero();
    dR2deta.Zero();
    dR2dxideta.Zero();
    myPatch->Nurbs2DBasis2ndDers(xi, eta, R, dRdxi, dRdeta, dR2dxi, dR2deta, dR2dxideta);




    // !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!
    // !    DESCRIPTION OF THE VARIABLES !
    // !~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~!

    // ! g: matrix with the 3 components of the in plane vectors g1, g2
    // ! dR: first derivatives of the shape functions
    // ! ddR: second derivatives of the shape functions
    // ! dn: first variation of the normal vector g3
    // ! lg3: length of the normal vector g3

    // ! E_cu: strain vector [E11,E22,E12] referred to the curvilinear coord. sys
    // ! E_ca: strain vector [E11,E22,2*E12] referred to cartesian coord. sys
    // ! K_cu: curvature vector [K11,K22,K12] referred to the curvilinear coord. sys
    // ! K_ca: curvature vector [K11,K22,2*K12] referred to the cartesian coord. sys

    // ! dE_cu: first variation of the membrane strain in local coord. system
    // ! dE_ca: first variation of the membrane strain in global coord. system
    // ! dK_cu: first variation of the curvature in local coord. system
    // ! dK_ca: first variation of the curvature in global coord. system

    // ! ddE_cu: second variation of the membrane strain in local coord. system
    // ! ddE_ca: second variation of the membrane strain in global coord. system
    // ! ddK_cu: second variation of the curvature in local coord. system
    // ! ddK_ca: second variation of the curvature in global coord. system

    

    ke.Zero(); 

    kem.Zero(); E_cu.Zero(); dE_cu.Zero(); E_ca.Zero(); dE_ca.Zero();
    N_ca.Zero();

    keb.Zero(); K_cu.Zero(); dK_cu.Zero(); K_ca.Zero();dK_ca.Zero();
    M_ca.Zero();


    // Get the 1st and the 2nd derivatives of the shape functions and call them dR, ddR
    for (int i = 0; i < noFuncs; ++i)
    {
      dr(0, i) = dRdxi(i);
      dr(1, i) = dRdeta(i);
      dr(2, i) = 0;
      ddr(0, i) = dR2dxi(i);
      ddr(1, i) = dR2deta(i);
      ddr(2, i) = dR2dxideta(i);
    }

    G_t = dr  * pts;
    H_t = ddr * pts;

    g_t = dr  * pts_d;
    h_t = ddr * pts_d;

    // Transposing because it is needed in cols drdxi, drdeta, etc
    // G is finally in shape (3,2) and H(3,3)
    Matrix G(3,2); G.Zero();
    Matrix g(3,2); g.Zero();
    Matrix H(3,3); H.Zero();
    Matrix h(3,3); h.Zero();

    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 2; ++j)
      {
        G(i,j) = G_t(j,i);
        g(i,j) = g_t(j,i);
      }
    }

    for (int i = 0; i < 3; ++i)
    {
      for (int j = 0; j < 3; ++j)
      {
        H(i,j) = H_t(j,i);
        h(i,j) = h_t(j,i);
      }
    }

    for (int i = 0; i < noFuncs; ++i)
    {
      dR(i, 0) = dRdxi(i);
      dR(i, 1) = dRdeta(i);
      dR(i, 2) = 0;

      ddR(i, 0) = dR2dxi(i);
      ddR(i, 1) = dR2deta(i);
      ddR(i, 2) = dR2dxideta(i);
    }


    // Call shell geometry functions
    Gab_r.Zero(); Bv_r.Zero(); Tb.Zero(); T_G_E.Zero(); g3.Zero(); n.Zero(); gab.Zero(); bv.Zero();
    shellGeo(G , H , g3_tmp , dA  , n_tmp , Gab_r , Bv_r , Tb, T_E_G_tmp, T_G_E);  // Reference configuration
    shellGeo(g , h , g3     , lg3 , n     , gab   , bv   , Tb_tmp, T_E_G_tmp, T_G_E_tmp); // Deformed or actual configuration

    // Compute the in-plane rotation of the plies
    double theta1 = 0.0;
    Vector bG1(3);
    bG1(0) = G(0, 0);
    bG1(1) = G(1, 0);
    bG1(2) = G(2, 0);

    if (bG1(0) > 0 && bG1(1) > 0) {
      theta1 = atan(fabs(bG1(1) / bG1(0)));

    } else if (bG1(0) > 0 && bG1(1) < 0) {
      theta1 = 2 * M_PI - atan(fabs(bG1(1) / bG1(0)));

    } else if (bG1(0)<0 && bG1(1)>0) {
      theta1 = M_PI - atan(fabs(bG1(1) / bG1(0)));

    } else if (bG1(0) < 0 && bG1(1) < 0) {
      theta1 = M_PI + atan(fabs(bG1(1) / bG1(0)));

    } else if (bG1(0) > 0 && bG1(1) == 0) {
      theta1 = 0;

    } else if (bG1(0) < 0 && bG1(1) == 0) {
      theta1 = M_PI;

    } else if (bG1(0) == 0 && bG1(1) > 0) {
      theta1 = 90.0 / 180.0 * M_PI;

    } else if (bG1(0) == 0 && bG1(1) < 0) {
      theta1 = 270.0 / 180.0 * M_PI;
    }

    theta1 = 0; // Not working, making zero in-plane rotation
    //

    //
    // Compute the membrane strains and the curvatures according to the eqns 1.161 - 1.162 of Computational FSI (Yuri's book)
    //

    E_cu.Zero();
    K_cu.Zero();

    //strain vector [E11,E22,E12] referred to curvilinear coor sys
    E_cu(0) = 0.5 * (gab(0, 0) - Gab_r(0, 0));
    E_cu(1) = 0.5 * (gab(1, 1) - Gab_r(1, 1));
    E_cu(2) = 0.5 * (gab(0, 1) - Gab_r(0, 1));

    //curvature vector [K11,K22,K12] referred to curvilinear coor sys
    K_cu(0) = -(bv(0) - Bv_r(0));
    K_cu(1) = -(bv(1) - Bv_r(1));
    K_cu(2) = -(bv(2) - Bv_r(2));


    E_ca = Tb * E_cu; // strain vector [E11,E22,2*E12] referred to cartesian coor sys, Tb is T_Gcon_E
    K_ca = Tb * K_cu; //curvature vector [K11,K22,2*K12] referred to cartesian coor sys, , Tb is T_Gcon_E


    A.Zero();         // Membrane stiffness matrix
    B.Zero();         // Coupling stifness matrix
    D.Zero();         // Bending stiffness matrix
    W = 0.0; //Total laminate weight (kg/m^2)

    // opserr << "Getting Tangent" << endln;

    for (int capa = 0; capa < nLayers; ++capa)
    {
      iAngle     = myPatch -> getAngle(capa) - theta1;
      iThickness = myPatch -> getThickness(capa);
      iZ         = myPatch -> getZk(capa);// Mid center laminate position
      iRho       = (OPS_getNDMaterial(myPatch -> getMatTag(capa))) -> getRho(); // Density of material kg/m^3;
      // opserr << "iRho = " << iRho << endln;

      W += iRho * iThickness;

      T(0, 0) = pow(cos(iAngle), 2)                        ; T(0, 1) = pow(sin(iAngle), 2)           ; T(0, 2) = sin(iAngle) * cos(iAngle)               ;
      T(1, 0) = pow(sin(iAngle), 2)                        ; T(1, 1) = pow(cos(iAngle), 2)           ; T(1, 2) = -sin(iAngle) * cos(iAngle)              ;
      T(2, 0) = -2 * sin(iAngle) * cos(iAngle)             ; T(2, 1) = 2 * sin(iAngle) * cos(iAngle) ; T(2, 2) = pow(cos(iAngle), 2) - pow(sin(iAngle), 2) ;

      // Aqui le pedimos al material su C

      // calcular Ezeta
      // Vector Ezeta = E_cu + iZ * K_cu; // Curvilinear
      Ezeta.Zero();
      Ezeta = T * (E_ca + iZ * K_ca); // Cartesian

      // Am i rotating twice?
      materialPointers[gp][capa] -> setTrialStrain(Ezeta); // Giving the trial strain to the material, rotated and in cartesian coordinates

      Matrix Czig = materialPointers[gp][capa]->getTangent();

      // int Q=myPatch->getOrders()(1);
      // if (Q==2) // THIS WAS Q==1, CHANGED IT BECAUSE OF THE BENCHMARKS, CHECK OUT
      // {
      //   // Czig(0,1) = Czig(1,0) = Czig(1,1) = Czig(2,2) = 0;
      //   Czig(0,1) = Czig(1,0) = Czig(0,0) = Czig(2,2) = 0;
      // }
      // else{
      //   Czig(0,1) = Czig(1,0) = Czig(1,1) = Czig(2,2) = 0;
      // }
      // Czig(0,1) = Czig(1,0) = Czig(2,2) = 0;

      C.addMatrixTripleProduct(0, T, Czig, 1.0);

      A.addMatrix(0, C, iThickness);
      B.addMatrix(0, C, iThickness * iZ);
      D.addMatrix(1.0, C, (iThickness * pow(iZ, 2) + pow(iThickness, 3) / 12) );
    }

    // Only bending stiffness
    A.Zero();
    B.Zero(); 

    // int Q=myPatch->getOrders()(1);
    // if (Q==2) // THIS WAS Q==1, CHANGED IT BECVAUSE OF THE BENCHAMRKS, CHECK OUT
    //   {
    //     D(1,1) = max(D(0,0),D(1,1));
    //     D(0,0) = 0;
    //   }
    // else{
    //     D(0,0) = max(D(0,0) , D(1,1));
    //     D(1,1) = 0;
    //   }

    D(0,0) = max(D(0,0),D(1,1));
    // D(1,1) = D(0,0);
    D(2,2) = 0;
    // opserr << "D = " << D << endln;


    double lg3_3 = pow(lg3, 3);
    double lg3_5 = pow(lg3, 5);

    //
    // Compute the First variation of strain and curvature w.r.t the d.o.f.
    //

    // Local variables
    int kr = 0;
    int dirr = 0;

    dg.Zero();

    // For "slicing" matrices
    ID ur_only(1);
    ID us_only(1);
    auto xyz_temp = arange<int>(0, 3); ID xyz(&(xyz_temp[0]), 3); // Takes value 0 1 2

    dg.Zero();

    for (int ur = 0; ur < 3 * noFuncs; ++ur)
    {
      // Local node number kr and dof direction dirr
      kr = (ur + 1 + 2) / 3 - 1; // takes value 0 0 0, 1 1 1, 2 2 2, .... noFuncs-1 noFuncs-1 noFuncs-1
      dirr = ur - 3 * (kr); // Takes value 0 1 2 0 1 2 ...

      dg(dirr, 0) = dR(kr, 0);
      dg(dirr, 1) = dR(kr, 1);

      // ADDING THIS, WEIRD
      // if (dirr==2)
      // {
      //   dg.Zero();
      // }
      // ADDING THIS, WEIRD

      // Strain
      dE_cu(0) = dR(kr, 0) * g(dirr, 0);
      dE_cu(1) = dR(kr, 1) * g(dirr, 1);
      dE_cu(2) = 0.5 * ( dR(kr, 0) * g(dirr, 1) + g(dirr, 0) * dR(kr, 1) );

      dE_ca(0, ur) = (Tb * dE_cu)(0);
      dE_ca(1, ur) = (Tb * dE_cu)(1);
      dE_ca(2, ur) = (Tb * dE_cu)(2);


      // Curvature
      dg3(0, ur) = dg(1, 0) * g(2, 1) - dg(2, 0) * g(1, 1) + g(1, 0) * dg(2, 1) - g(2, 0) * dg(1, 1);
      dg3(1, ur) = dg(2, 0) * g(0, 1) - dg(0, 0) * g(2, 1) + g(2, 0) * dg(0, 1) - g(0, 0) * dg(2, 1);
      dg3(2, ur) = dg(0, 0) * g(1, 1) - dg(1, 0) * g(0, 1) + g(0, 0) * dg(1, 1) - g(1, 0) * dg(0, 1);
      // dg3(2, ur) = 0;



      // Kratos code THIS GETS CLOSE!
      if (dirr == 0)
      {
        dg3(0, ur) = 0;
        dg3(1, ur) = -dR(kr,0) * g(2,1) + dR(kr,1) * g(2,0);
        dg3(2, ur) = dR(kr,0) * g(1,1) - dR(kr,1) * g(1,0);
      }
      else if (dirr == 1)
      {
        dg3(0, ur) = dR(kr,0) * g(2,1) - dR(kr,1) * g(2,0); 
        dg3(1, ur) = 0;
        dg3(2, ur) = -dR(kr,0) * g(0,1) + dR(kr,1) * g(0,0);
      }
      else if (dirr == 2)
      {
        dg3(0, ur) = -dR(kr,0) * g(1,1) + dR(kr,1) * g(1,0); 
        dg3(1, ur) = dR(kr,0) * g(0,1) - dR(kr,1) * g(0,0);
        dg3(2, ur) = 0;
      }

      g3dg3(ur) = g3(0) * dg3(0, ur) + g3(1) * dg3(1, ur) + g3(2) * dg3(2, ur);
      g3dg3lg3_3(ur) = g3dg3(ur) / lg3_3;

      dn(0, ur) = dg3(0, ur) / lg3 - g3(0) * g3dg3lg3_3(ur);
      dn(1, ur) = dg3(1, ur) / lg3 - g3(1) * g3dg3lg3_3(ur);
      dn(2, ur) = dg3(2, ur) / lg3 - g3(2) * g3dg3lg3_3(ur);
      // dn(2, ur) = 0;


      // // Original
      // dK_cu(0) = -( ddR(kr, 0) * n(dirr) + (h(0, 0) * dn(0, ur) + h(1, 0) * dn(1, ur) + h(2, 0) * dn(2, ur)) );
      // dK_cu(1) = -( ddR(kr, 1) * n(dirr) + (h(0, 1) * dn(0, ur) + h(1, 1) * dn(1, ur) + h(2, 1) * dn(2, ur)) );
      // dK_cu(2) = -( ddR(kr, 2) * n(dirr) + (h(0, 2) * dn(0, ur) + h(1, 2) * dn(1, ur) + h(2, 2) * dn(2, ur)) ); 


      // testing
      dK_cu(0) = 0 - (ddR(kr, 0) * n(dirr) + (h(0, 0) * dn(0, ur) + h(1, 0) * dn(1, ur) + h(2, 0) * dn(2, ur))) ;
      dK_cu(1) = 0 - (ddR(kr, 1) * n(dirr) + (h(0, 1) * dn(0, ur) + h(1, 1) * dn(1, ur) + h(2, 1) * dn(2, ur))) ;
      dK_cu(2) = 0 - (ddR(kr, 2) * n(dirr) + (h(0, 2) * dn(0, ur) + h(1, 2) * dn(1, ur) + h(2, 2) * dn(2, ur))) ;

      // // Kratos code
      if (dirr == 0)
      {
        dK_cu(0) = 0 - (ddR(kr,0)*n(0) + h(0,0)*dn(0,ur) + h(1,0)*dn(1,ur)+h(2,0)*dn(2,ur));
        dK_cu(1) = 0 - (ddR(kr,1)*n(0) + h(0,1)*dn(0,ur) + h(1,1)*dn(1,ur)+h(2,1)*dn(2,ur));
        dK_cu(2) = 0 - (ddR(kr,2)*n(0) + h(0,2)*dn(0,ur) + h(1,2)*dn(1,ur)+h(2,2)*dn(2,ur));
      }
      else if (dirr == 1)
      {
        dK_cu(0) = 0 - (ddR(kr,0)*n(1) + h(0,0)*dn(0,ur) + h(1,0)*dn(1,ur)+h(2,0)*dn(2,ur));
        dK_cu(1) = 0 - (ddR(kr,1)*n(1) + h(0,1)*dn(0,ur) + h(1,1)*dn(1,ur)+h(2,1)*dn(2,ur));
        dK_cu(2) = 0 - (ddR(kr,2)*n(1) + h(0,2)*dn(0,ur) + h(1,2)*dn(1,ur)+h(2,2)*dn(2,ur));
      }
      else if (dirr == 2)
      {
        dK_cu(0) = 0 - (ddR(kr,0)*n(2) + h(0,0)*dn(0,ur) + h(1,0)*dn(1,ur)+h(2,0)*dn(2,ur));
        dK_cu(1) = 0 - (ddR(kr,1)*n(2) + h(0,1)*dn(0,ur) + h(1,1)*dn(1,ur)+h(2,1)*dn(2,ur));
        dK_cu(2) = 0 - (ddR(kr,2)*n(2) + h(0,2)*dn(0,ur) + h(1,2)*dn(1,ur)+h(2,2)*dn(2,ur));
      }

      // dK_cu(0) = -( h(dirr, 0) * n(dirr) + (ddR(0, 0) * dn(0, ur) + ddR(1, 0) * dn(1, ur) + ddR(2, 0) * dn(2, ur)) );
      // dK_cu(1) = -( h(dirr, 1) * n(dirr) + (ddR(0, 1) * dn(0, ur) + ddR(1, 1) * dn(1, ur) + ddR(2, 1) * dn(2, ur)) );
      // dK_cu(2) = -( h(dirr, 2) * n(dirr) + (ddR(0, 2) * dn(0, ur) + ddR(1, 2) * dn(1, ur) + ddR(2, 2) * dn(2, ur)) );

      // Maybe it is this
      // dK_cu(0) = -( h(dirr, 0) * n(dirr) + (ddR(0, 0) * dn(0, ur) + ddR(1, 0) * dn(1, ur) + ddR(2, 0) * dn(2, ur)) );
      // dK_cu(1) = -( h(dirr, 1) * n(dirr) + (ddR(0, 1) * dn(0, ur) + ddR(1, 1) * dn(1, ur) + ddR(2, 1) * dn(2, ur)) );
      // dK_cu(2) = -( h(dirr, 2) * n(dirr) + (ddR(0, 2) * dn(0, ur) + ddR(1, 2) * dn(1, ur) + ddR(2, 2) * dn(2, ur)) );

      // dK_cu(0) = -(ddR(0, kr) * n(dirr) + (h(0, 0) * dn(0, ur) + h(1, 0) * dn(1, ur) + h(2, 0) * dn(2, ur)));
      // dK_cu(1) = -(ddR(1, kr) * n(dirr) + (h(0, 1) * dn(0, ur) + h(1, 1) * dn(1, ur) + h(2, 1) * dn(2, ur)));
      // dK_cu(2) = -(ddR(2, kr) * n(dirr) + (h(0, 2) * dn(0, ur) + h(1, 2) * dn(1, ur) + h(2, 2) * dn(2, ur)));

      // CHECK THIS OUT, I HAD TO SWAP THE INDICES FOR H

      // dK_cu(0) = -(ddR(0, kr) * n(dirr) + (h(0, 0) * dn(0, ur) + h(0, 1) * dn(1, ur) + h(0, 2) * dn(2, ur)));
      // dK_cu(1) = -(ddR(1, kr) * n(dirr) + (h(1, 0) * dn(0, ur) + h(1, 1) * dn(1, ur) + h(1, 2) * dn(2, ur)));
      // dK_cu(2) = -(ddR(2, kr) * n(dirr) + (h(2, 0) * dn(0, ur) + h(2, 1) * dn(1, ur) + h(2, 2) * dn(2, ur)));

      dK_ca(0, ur) = (Tb * dK_cu)(0);
      dK_ca(1, ur) = (Tb * dK_cu)(1);
      dK_ca(2, ur) = (Tb * dK_cu)(2);

    }


    // Second variation of strain and curvature w.r.t. dofs

    // Local variables
    int ks = 0;
    int dirs = 0;
    int dirt = 0;
    int ddir = 0;

    double tmp1 = 0.0;
    double tmp2 = 0.0;



    if (nonLinearGeometry) // true if non-linear geometry
    {
      // opserr << "Using nonLinearGeometry!" << endln;
      for (int ur = 0; ur < 3 * noFuncs; ++ur)
      {
        // Local node number kr and dof direction dirr
        kr = (ur + 1 + 2) / 3 - 1; // takes value 0 0 0, 1 1 1, 2 2 2, .... noFuncs-1 noFuncs-1 noFuncs-1
        dirr = ur - 3 * (kr); // Takes value 0 1 2 0 1 2 ...

        for (int us = 0; us <= ur; ++us)
        {
          ks = (us + 1 + 2) / 3 - 1; // takes value 0 0 0, 1 1 1, 2 2 2, .... noFuncs-1 noFuncs-1 noFuncs-1
          // ks = (us + 2) / 3 ; // takes value 0 0 0, 1 1 1, 2 2 2, .... noFuncs-1 noFuncs-1 noFuncs-1
          dirs = us - 3 * (ks); // Takes value 0 1 2 0 1 2 ...

          // Strain
          ddE_cu.Zero();
          if (dirr == dirs)
          {
            ddE_cu(0) = dR(kr, 0) * dR(ks, 0);
            ddE_cu(1) = dR(kr, 1) * dR(ks, 1);
            ddE_cu(2) = 0.5 * ( dR(kr, 0) * dR(ks, 1) + dR(kr, 1) * dR(ks, 0) );
          }

          ddE_ca_0(ur, us) = (Tb * ddE_cu)(0);
          ddE_ca_1(ur, us) = (Tb * ddE_cu)(1);
          ddE_ca_2(ur, us) = (Tb * ddE_cu)(2);

          // Curvature
          ddg3.Zero();

          dirt = 6 - dirr - dirs; // Check this
          ddir = dirr - dirs;

          if (ddir == -1 || ddir == 2)
          {
            // opserr << "dirt - 3 = " << dirt-3 << endln;
            ddg3(dirt - 3) = dR(kr, 0) * dR(ks, 1) - dR(ks, 0) * dR(kr, 1);
          }
          else if (ddir == 1 || ddir == -2)
          {
            // opserr << "dirt - 3 = " << dirt-3 << endln;
            ddg3(dirt - 3) = -dR(kr, 0) * dR(ks, 1) + dR(ks, 0) * dR(kr, 1);
          }

          ur_only.fill(ur);
          us_only.fill(us);

          tmp1 = -( (ddg3 ^ g3) + (transpose(3, 1, dg3(xyz, ur_only)) * dg3(xyz, us_only)) (0, 0) ) / lg3_3 ; // the (0,0) is to transform the 1x1 matrix into double
          tmp2 = 3 * g3dg3(ur) * g3dg3(us) / lg3_5;



          ddn(0) = ddg3(0) / lg3 - g3dg3lg3_3(us) * dg3(0, ur) - g3dg3lg3_3(ur) * dg3(0, us) + tmp1 * g3(0) + tmp2 * g3(0);
          ddn(1) = ddg3(1) / lg3 - g3dg3lg3_3(us) * dg3(1, ur) - g3dg3lg3_3(ur) * dg3(1, us) + tmp1 * g3(1) + tmp2 * g3(1);
          ddn(2) = ddg3(2) / lg3 - g3dg3lg3_3(us) * dg3(2, ur) - g3dg3lg3_3(ur) * dg3(2, us) + tmp1 * g3(2) + tmp2 * g3(2);

          // ddK_cu(0) = -(ddR(0, kr) * dn(dirr, us) + ddR(0, ks) * dn(dirs, ur) + (h(0, 0) * ddn(0) + h(0, 1) * ddn(1) + h(0, 2) * ddn(2)));
          // ddK_cu(1) = -(ddR(1, kr) * dn(dirr, us) + ddR(1, ks) * dn(dirs, ur) + (h(1, 0) * ddn(0) + h(1, 1) * ddn(1) + h(1, 2) * ddn(2)));
          // ddK_cu(2) = -(ddR(2, kr) * dn(dirr, us) + ddR(2, ks) * dn(dirs, ur) + (h(2, 0) * ddn(0) + h(2, 1) * ddn(1) + h(2, 2) * ddn(2)));

          // ddK_cu(0) = -(ddR(0, kr) * dn(dirr, us) + ddR(0, ks) * dn(dirs, ur) + (h(0, 0) * ddn(0) + h(1, 0) * ddn(1) + h(2, 0) * ddn(2)));
          // ddK_cu(1) = -(ddR(1, kr) * dn(dirr, us) + ddR(1, ks) * dn(dirs, ur) + (h(0, 1) * ddn(0) + h(1, 1) * ddn(1) + h(2, 1) * ddn(2)));
          // ddK_cu(2) = -(ddR(2, kr) * dn(dirr, us) + ddR(2, ks) * dn(dirs, ur) + (h(0, 2) * ddn(0) + h(1, 2) * ddn(1) + h(2, 2) * ddn(2)));

          ddK_cu(0) = -( ddR(kr, 0) * dn(dirr, us) + ddR(ks, 0) * dn(dirs, ur) + h(0, 0) * ddn(0) + h(1, 0) * ddn(1) + h(2, 0) * ddn(2) );
          ddK_cu(1) = -( ddR(kr, 1) * dn(dirr, us) + ddR(ks, 1) * dn(dirs, ur) + h(0, 1) * ddn(0) + h(1, 1) * ddn(1) + h(2, 1) * ddn(2) );
          ddK_cu(2) = -( ddR(kr, 2) * dn(dirr, us) + ddR(ks, 2) * dn(dirs, ur) + h(0, 2) * ddn(0) + h(1, 2) * ddn(1) + h(2, 2) * ddn(2) );

          ddK_ca_0(ur, us) = (Tb * ddK_cu)(0);
          ddK_ca_1(ur, us) = (Tb * ddK_cu)(1);
          ddK_ca_2(ur, us) = (Tb * ddK_cu)(2);
        }
      }
    }

    N_ca.Zero();
    M_ca.Zero();
 

    //alternativamente // necesitas para plasticidad
    for (int capa = 0; capa < nLayers; ++capa)
    {
      Matrix Czig = materialPointers[gp][capa]->getTangent();

      double iThickness = myPatch->getThickness(capa);
      double iZ         = myPatch->getZk(capa);// Mid center laminate position
      double iAngle     = myPatch->getAngle(capa) - theta1;

      T(0, 0) = pow(cos(iAngle), 2)                        ; T(0, 1) = pow(sin(iAngle), 2)           ; T(0, 2) = sin(iAngle) * cos(iAngle)               ;
      T(1, 0) = pow(sin(iAngle), 2)                        ; T(1, 1) = pow(cos(iAngle), 2)           ; T(1, 2) = -sin(iAngle) * cos(iAngle)              ;
      T(2, 0) = -2 * sin(iAngle) * cos(iAngle)             ; T(2, 1) = 2 * sin(iAngle) * cos(iAngle) ; T(2, 2) = pow(cos(iAngle), 2) - pow(sin(iAngle), 2) ;

      // int Q=myPatch->getOrders()(1);
      // if (Q==2) // THIS WAS Q==1, CHANGED IT BECVAUSE OF THE BENCHAMRKS, CHECK OUT
      // {
      //   Czig(0,1) = Czig(1,0) = Czig(0,0) = Czig(2,2) = 0;
      //   // Czig(0,1) = Czig(1,0) = Czig(1,1) = Czig(2,2) = 0;
      // }
      // else{
      //   Czig(0,1) = Czig(1,0) = Czig(1,1) = Czig(2,2) = 0;
      // }
      // Czig(0,1) = Czig(1,0) = Czig(2,2) = 0;


      Matrix C(3, 3);
      C.addMatrixTripleProduct(0, T, Czig, 1.0);


      const Vector Szeta = transpose(3, 3, T) * materialPointers[gp][capa]->getStress();  // Rotating stresses
      // Szeta = transpose(3, 3, T) * Szeta; // Rotating stresses
      // Szeta = T * Szeta; // Rotating stresses

      // Integrating stress in every ply
      // N_ca += iThickness * Szeta; Only bending stiffness

      // MAYBE WRONG HERE, 
      // M_ca += iThickness * iZ * Szeta + (C * K_ca) * pow(iThickness, 3) / 12.0 ;
      M_ca += (C * K_ca) * pow(iThickness, 3) / 12.0 ;

    }

    N_ca = A * E_ca + B * K_ca;
    M_ca = B * E_ca + D * K_ca;

    kem.Zero();
    keb.Zero();

    for (int ur = 0; ur < 3 * noFuncs ; ++ur)
    {

      dE_ca_ur(0) = dE_ca(0, ur);
      dE_ca_ur(1) = dE_ca(1, ur);
      dE_ca_ur(2) = dE_ca(2, ur);

      dK_ca_ur(0) = dK_ca(0, ur);
      dK_ca_ur(1) = dK_ca(1, ur);
      dK_ca_ur(2) = dK_ca(2, ur);

      dN_ca = A * dE_ca_ur + B * dK_ca_ur;
      dM_ca = B * dE_ca_ur + D * dK_ca_ur;


      for (int us = 0; us <= ur; ++us)  
      {
    
        // Vectors for non-linear part
        ddE_ca(0) = ddE_ca_0(ur, us);
        ddE_ca(1) = ddE_ca_1(ur, us);
        ddE_ca(2) = ddE_ca_2(ur, us);

        ddK_ca(0) = ddK_ca_0(ur, us);
        ddK_ca(1) = ddK_ca_1(ur, us);
        ddK_ca(2) = ddK_ca_2(ur, us);

        // Membrane stiffness
        dE_ca_us(0) = dE_ca(0, us);
        dE_ca_us(1) = dE_ca(1, us);
        dE_ca_us(2) = dE_ca(2, us);

        kem(ur, us) = dN_ca ^ dE_ca_us; // Linear part
        // No membrane stiffness in bending strip
        

        // Bending stiffness
        dK_ca_us(0) = dK_ca(0, us);
        dK_ca_us(1) = dK_ca(1, us);
        dK_ca_us(2) = dK_ca(2, us);

        keb(ur, us) = dM_ca ^ dK_ca_us; // Linear part




        if (nonLinearGeometry) // 1 if nonLinear, 0 if Linear
        {
          // SHOULD THIS BE HERE?
          // kem(ur, us) += N_ca ^ ddE_ca; // Non-linear part
          keb(ur, us) += M_ca ^ ddK_ca; // Non-linear part
        }

        // Symmetric parts
        kem(us, ur) = kem(ur, us);
        keb(us, ur) = keb(ur, us);
      }

      // Force vectors (Residual)
      // fiem(ur) = (N_ca ^ dE_ca_ur);
      fieb(ur) = (M_ca ^ dK_ca_ur);

    }

    // for (int ur = 0; ur < 3*noFuncs; ++ur)
    // {
    //   // Local node number kr and dof direction dirr
    //   kr = (ur + 1 + 2) / 3 - 1; // takes value 0 0 0, 1 1 1, 2 2 2, .... noFuncs-1 noFuncs-1 noFuncs-1
    //   dirr = ur - 3 * (kr); // Takes value 0 1 2 0 1 2 ...

    //   opserr << "kr, dirr = " << kr <<' ' << dirr << endln;
    //   for (int us = 0; us < 3*noFuncs; ++us)
    //   {
    //     if (us>=ur) { ke(ur,us) = ke(us,ur); }
    //     ks = (us + 1 + 2) / 3 - 1; // takes value 0 0 0, 1 1 1, 2 2 2, .... noFuncs-1 noFuncs-1 noFuncs-1
    //     dirs = us - 3 * (ks); // Takes value 0 1 2 0 1 2 ...
    //     opserr << "ks,dirs = " << ks <<' '<< dirs << endln;
    //   }
    // }

    // opserr << "kem*(dA * J2 * wt) = " << kem*(dA * J2 * wt) << endln;
    // opserr << "keb*(dA * J2 * wt) = " << keb*(dA * J2 * wt) << endln;

    Matrix N(3, 3 * noFuncs);
    Vector pressure_vector = pressure * g3;

    for (int i = 0; i < noFuncs; ++i)
    {
      N(0, 3 * i) = R(i);
      N(1, 3 * i + 1) = R(i);
      N(2, 3 * i + 2) = R(i);
    }

    // keb.Zero(); // Just debugging, delete after
    // kem.Zero(); // Just debugging, delete after
    fie = (fiem + fieb) * (dA * J2 * wt); //Residual vector on gauss point
    // fie.addMatrixTransposeVector(1.0, N, pressure_vector, dA * J2 * wt);
    ke = (kem + keb) * (dA * J2 * wt); //Stiffness matrix on gauss point
    // opserr << "Adding just bending stiffness!"<< endln;

    // opserr << "kem*(dA * J2 * wt) = " << kem*(dA * J2 * wt) << endln;
    
    // opserr << "ke = " << ke << endln;

    // Body forces
    // Vector b(appliedB, 3);

    // if (applyLoad == 1)
    // {
    //   fie.addMatrixTransposeVector(1.0, N, b, -1 * W * dA * J2 * wt);
    // } else {
    //   fie.addMatrixTransposeVector(1.0, N, gFact, -1 * W * dA * J2 * wt);
    // }

    // Add contribution to residual vector
    *resid += fie;

    // Add contribution to stiffness matrix if tang_flag == 1
    if (tang_flag == 1 ) {
      *stiff += ke;
    }
  }
  // opserr << "*stiff = " << *stiff << endln;
}


//get residual with inertia terms
const Vector&  IGAKLShell_BendingStrip::getResistingForceIncInertia( )
{
  // opserr << "IGAKLShell::getResistingForceIncInertia - eleTag" << this->getTag() << " called! " << endln;
  int Nnodes = connectedExternalNodes.Size();
  int NDOF = 3 * Nnodes;
  static Vector res(NDOF);   /// ngdl
  res.resize(NDOF);
  res.Zero();
  // Calculo aceleraciones en nodos y los guardo en res o accel
  static Vector NodalAccelerations(NDOF);
  static Vector NodalDisplacements(NDOF);
  NodalAccelerations.resize(NDOF);
  NodalDisplacements.resize(NDOF);

  NodalAccelerations.Zero();
  NodalDisplacements.Zero();
  static Vector accel_i(3);
  static Vector disp_i(3);
  accel_i.Zero();
  disp_i.Zero();


  // Loop through nodes
  int k = 0;
  for (int i = 0; i < Nnodes; ++i)
  {
    Node *node_i = nodePointers[i];
    accel_i = node_i->getTrialAccel();
    disp_i = node_i->getTrialDisp();
    for (int j = 0; j < 3; ++j)
    {
      NodalAccelerations(k) = accel_i(j);
      NodalDisplacements(k) = disp_i(j);
      k += 1;
      // opserr << "NodalAccelerations(Nnodes * j + i) = " << NodalAccelerations(Nnodes * j + i) << endln;
    }
  }

  int tang_flag = 1 ; //get the tangent
  bool nonLinearGeometry = myPatch->getAnalysisType();

  // formResidAndTangent( tang_flag ) ;

  // res = *resid;
  res = this->getResistingForce(); // + this->getMass() * NodalAccelerations;
  // Calculo masa M
  // res += this->getMass() * NodalAccelerations;

  // if (not nonLinearGeometry)
  // {
  //   // opserr << "IGAKLShell::getResistingForceIncInertia - eleTag" << this->getTag() << " called! " << endln;
  //   // res += 1*(*stiff) * NodalDisplacements;
  //   res += 1*(this->getTangentStiff()) * NodalDisplacements;
  //   // *resid += 1*(*stiff) * NodalDisplacements;
  // }

  // // add the damping forces if rayleigh damping
  // if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
  //   res += this->getRayleighDampingForces();
    // *resid += this->getRayleighDampingForces();

  // opserr << "M = " << this->getMass() << endln;

  // opserr << "NodalAccelerations = " << NodalAccelerations << endln;

  // opserr << "this->getMass() * NodalAccelerations = " << this->getMass() * NodalAccelerations << endln;

  // *resid += this->getMass() * NodalAccelerations;
  // res = getResistingForce() + this->getMass() * NodalAccelerations;

  // res = this->resid

  // opserr << "Finished IGAKLShell::getResistingForceIncInertia - eleTag" << this->getTag() << " called! " << endln;
  *resid=res;
  // return res;
  return *resid;
}

//**********************************************************************

Matrix
IGAKLShell_BendingStrip::transpose( int dim1,
                       int dim2,
                       const Matrix &M )
{
  int i ;
  int j ;

  Matrix Mtran( dim2, dim1 ) ;

  for ( i = 0; i < dim1; i++ ) {
    for ( j = 0; j < dim2; j++ )
      Mtran(j, i) = M(i, j) ;
  } // end for i

  return Mtran ;
}

//**********************************************************************

int  IGAKLShell_BendingStrip::sendSelf (int commitTag, Channel &theChannel)
{
  int res = 0;

  return res;
}

int  IGAKLShell_BendingStrip::recvSelf (int commitTag,
                           Channel &theChannel,
                           FEM_ObjectBroker &theBroker)
{
  int res = 0;



  return res;
}
//**************************************************************************

int
IGAKLShell_BendingStrip::displaySelf(Renderer &theViewer, int displayMode, float fact, const char **modes, int numMode)
{

  int error = 0;

  return error;
}


bool IGAKLShell_BendingStrip::pointInElement(double xi, double eta) const
{
  bool in_xi = xiE(0) <= xi && xi <= xiE(1);
  bool in_eta = etaE(0) <= eta && eta <= etaE(1);
  return in_xi && in_eta;
}
