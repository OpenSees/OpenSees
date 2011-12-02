///////////////////////////////////////////////////////////////////////////////
//
// COPYRIGHT (C):     :-))
// PROJECT:           Object Oriented Finite Element Program
// FILE:              TwentyNodeBrick_u_p_U.cpp
// CLASS:             TwentyNodeBrick_u_p_U
// MEMBER FUNCTIONS:
//
// MEMBER VARIABLES
//
// PURPOSE:           Finite Element Class for coupled system
//  "Coupled system": Solid and fluid coexist.
//                    u-- Solid displacement
//                    p-- Pore pressure
//                    U-- Absolute fluid displacement
//
// RETURN:
// VERSION:
// LANGUAGE:          C++
// TARGET OS:         DOS || UNIX || . . .
// DESIGNER:          Boris Jeremic, Zhao Cheng
// PROGRAMMER:        Boris Jeremic, Zhaohui Yang, Xiaoyan Wu, Zhao Cheng
// DATE:              Aug. 2001
// UPDATE HISTORY:    Modified from EightNodeBrick.cpp  reorganized a lot by Xiaoyan
//                    01/24/2002    Xiaoyan
//                    Add the permeability tensor and ks, kf to the constructor  Xiaoyan
//
//                    Clean-up and re-write by Zhao Cheng, 10/20/2004
//
//                    Fixed a bug, and some small modification, ZC 05/2006
///////////////////////////////////////////////////////////////////////////////
//
#ifndef TWENTYNODEBRICK_U_P_U_CPP
#define TWENTYNODEBRICK_U_P_U_CPP

#include <TwentyNodeBrick_u_p_U.h>

const int TwentyNodeBrick_u_p_U::Num_IntegrationPts = 3;
const int TwentyNodeBrick_u_p_U::Num_TotalGaussPts = 27;
const int TwentyNodeBrick_u_p_U::Num_Nodes = 20;
const int TwentyNodeBrick_u_p_U::Num_Dim = 3;
const int TwentyNodeBrick_u_p_U::Num_Dof = 7;
const int TwentyNodeBrick_u_p_U::Num_ElemDof = 140;
const double TwentyNodeBrick_u_p_U::pts[3] = {-0.774596669241483, 0.0, +0.774596669241483};
const double TwentyNodeBrick_u_p_U::wts[3] = {5.0/9.0, 8.0/9.0, 5.0/9.0};
Matrix TwentyNodeBrick_u_p_U::MCK(Num_ElemDof, Num_ElemDof);
Vector TwentyNodeBrick_u_p_U::P(Num_ElemDof);

//======================================================================
TwentyNodeBrick_u_p_U::TwentyNodeBrick_u_p_U(int element_number,
                                             int node_numb_1,
                                             int node_numb_2,
                                             int node_numb_3,
                                             int node_numb_4,
                                             int node_numb_5,
                                             int node_numb_6,
                                             int node_numb_7,
                                             int node_numb_8,
                                             int node_numb_9,
                                             int node_numb_10,
                                             int node_numb_11,
                                             int node_numb_12,
                                             int node_numb_13,
                                             int node_numb_14,
                                             int node_numb_15,
                                             int node_numb_16,
                                             int node_numb_17,
                                             int node_numb_18,
                                             int node_numb_19,
                                             int node_numb_20,
                                           NDMaterial *Globalmmodel,
                                           double b1,
                                           double b2,
                                           double b3,
                                           double nn,
                                           double alf,
                                           double rs,
                                           double rf,
                                           double permb_x,
                                           double permb_y,
                                           double permb_z,
                                           double kks,
                                           double kkf,
                                           double pp)
 : Element(element_number,
           ELE_TAG_TwentyNodeBrick_u_p_U ),
           connectedExternalNodes(Num_Nodes),
           perm(Num_Dim),
           bf(Num_Dim),
           poro(nn),
           alpha(alf),
           rho_s(rs),
           rho_f(rf),
           ks(kks),
           kf(kkf),
           pressure(pp),
           Q(0),
           Ki(0)
{
    // permeability
    perm(0) = permb_x;
    perm(1) = permb_y;
    perm(2) = permb_z;

    if (perm(0)==0.0 || perm(1)==0.0 || perm(2)==0.0) {
       opserr<<" Error, TwentyNodeBrick_u_p_U:: permeability (kx/ky/kz) is zero! \n";
       exit(-1);
    }

    // body forces
    bf(0) = b1;
    bf(1) = b2;
    bf(2) = b3;

    connectedExternalNodes( 0) = node_numb_1;
    connectedExternalNodes( 1) = node_numb_2;
    connectedExternalNodes( 2) = node_numb_3;
    connectedExternalNodes( 3) = node_numb_4;
    connectedExternalNodes( 4) = node_numb_5;
    connectedExternalNodes( 5) = node_numb_6;
    connectedExternalNodes( 6) = node_numb_7;
    connectedExternalNodes( 7) = node_numb_8;
    connectedExternalNodes( 8) = node_numb_9;
    connectedExternalNodes( 9) = node_numb_10;
    connectedExternalNodes(10) = node_numb_11;
    connectedExternalNodes(11) = node_numb_12;
    connectedExternalNodes(12) = node_numb_13;
    connectedExternalNodes(13) = node_numb_14;
    connectedExternalNodes(14) = node_numb_15;
    connectedExternalNodes(15) = node_numb_16;
    connectedExternalNodes(16) = node_numb_17;
    connectedExternalNodes(17) = node_numb_18;
    connectedExternalNodes(18) = node_numb_19;
    connectedExternalNodes(19) = node_numb_20;

    theMaterial = new NDMaterial *[Num_TotalGaussPts];
    if (theMaterial == 0) {
       opserr<<" TwentyNodeBrick_u_p_U::TwentyNodeBrick_u_p_U -- failed allocate material model pointer\n";
       exit(-1);
    }
    for (int i=0; i<Num_TotalGaussPts; i++) {
       theMaterial[i] = Globalmmodel->getCopy();
       if (theMaterial[i] == 0) {
        opserr<<" TwentyNodeBrick_u_p_U::TwentyNodeBrick_u_p_U -- failed allocate material model pointer\n";
        exit(-1);
       }
    }

}


//======================================================================
TwentyNodeBrick_u_p_U::TwentyNodeBrick_u_p_U ()
 : Element(0, ELE_TAG_TwentyNodeBrick_u_p_U ),
   connectedExternalNodes(Num_Nodes), perm(Num_Dim), bf(Num_Dim),
   poro(0.0), alpha(1.0), rho_s(0.0),rho_f(0.0), ks(0.0), kf(0.0), pressure(0.0), Q(0), Ki(0)
{
   theMaterial = 0;

   for (int j=0; j<Num_Nodes; j++)
     theNodes[j] = 0;
}

//======================================================================
TwentyNodeBrick_u_p_U::~TwentyNodeBrick_u_p_U ()
{
   if (theMaterial)
     delete [] theMaterial;

   for (int j=0; j<Num_Nodes; j++)
     theNodes[j] = 0;

   if (Q != 0)
     delete Q;

   if (Ki != 0)
     delete Ki;

}

//======================================================================
int TwentyNodeBrick_u_p_U::getNumExternalNodes (void) const
{
    return Num_Nodes;
}

//======================================================================
const ID& TwentyNodeBrick_u_p_U::getExternalNodes (void)
{
    return connectedExternalNodes;
}

//======================================================================
Node ** TwentyNodeBrick_u_p_U::getNodePtrs (void)
{
        return theNodes;
}

//======================================================================
int TwentyNodeBrick_u_p_U::getNumDOF (void)
{
    return Num_ElemDof;
}

//======================================================================
void TwentyNodeBrick_u_p_U::setDomain (Domain *theDomain)
{
  int i, Ndof;

  if (theDomain == 0) {
    for (i=0; i<Num_Nodes; i++) {
      theNodes[i] = 0;
    }
  }

  for (i=0; i<Num_Nodes; i++) {
    theNodes[i] = theDomain->getNode(connectedExternalNodes(i));
    if (theNodes[i] == 0) {
      opserr << "Error TwentyNodeBrick_u_p_U : node not found in the domain" << "\n";
      return ;
    }
    Ndof = theNodes[i]->getNumberDOF();
    if( Ndof != Num_Dof) {
      opserr << "Error TwentyNodeBrick_u_p_U : has wrong number of DOFs at its nodes" << "\n";
      return ;
    }
  }

  this->DomainComponent::setDomain(theDomain);

}

//======================================================================
int TwentyNodeBrick_u_p_U::commitState (void)
{
    int retVal = 0;
    int i;

    if ((retVal = this->Element::commitState()) != 0) {
      opserr << "EightNodeBrick-u_p_U::commitState () - failed in base class";
      return (-1);
    }

    for (i=0; i<Num_TotalGaussPts; i++ )
      retVal += theMaterial[i]->commitState();

    return retVal;
}

//======================================================================
int TwentyNodeBrick_u_p_U::revertToLastCommit (void)
{
    int retVal = 0;
    int i;

    for (i=0; i<Num_TotalGaussPts; i++ )
      retVal += theMaterial[i]->revertToLastCommit() ;

    return retVal;
}

//======================================================================
int TwentyNodeBrick_u_p_U::revertToStart (void)
{
    int retVal = 0;
    int i;

    for (i=0; i<Num_TotalGaussPts; i++ )
      retVal += theMaterial[i]->revertToStart() ;

    return retVal;
}

//======================================================================
const Matrix& TwentyNodeBrick_u_p_U::getTangentStiff (void)
{
    return this->getStiff(1);
}

//======================================================================
const Matrix& TwentyNodeBrick_u_p_U::getInitialStiff (void)
{
    return this->getStiff(0);
}

//======================================================================
const Matrix& TwentyNodeBrick_u_p_U::getDamp (void)
{
    MCK.Zero();  // necessary
        
	tensor tC = getDampTensorC123();
    
    int i, j, m, n;

    double Ctemp = 0.0;
    tensor CRm;
    tensor CRk;    
    
    if (alphaM != 0.0) 
      CRm = getMassTensorMsf();
    
    if (betaK != 0.0) 
      CRk = getStiffnessTensorKep();
    
    if (betaK0 != 0.0 || betaKc != 0.0) {
	  opserr << "Warning: EightNodeBrick-u_p_U:: betaK0 or betaKc are not used" << "\n";  
    }

    for ( i=0 ; i<Num_Nodes; i++ ) {
      for ( j=0; j<Num_Nodes; j++ ) {
        for( m=0; m<Num_Dim; m++) {
          for( n=0; n<Num_Dim; n++) 
            {
              Ctemp = tC.cval(i+1, m+1, n+1, j+1) *(poro*poro);
	          
              //C1
              MCK(i*Num_Dof+m, j*Num_Dof+n) = Ctemp;
              
              if (alphaM != 0.0)
	            MCK(i*Num_Dof+m, j*Num_Dof+n) += CRm.cval(i+1, j+1) *((1.0-poro)*rho_s *alphaM);
              
              if (betaK != 0.0)
	            MCK(i*Num_Dof+m, j*Num_Dof+n) += CRk.cval(i+1, m+1, n+1, j+1) * betaK;
              
              //C3
              MCK(i*Num_Dof+m+4, j*Num_Dof+n+4) = Ctemp;
              
              if (alphaM != 0.0)
                MCK(i*Num_Dof+m+4, j*Num_Dof+n+4) += CRm.cval(i+1, j+1) *(poro*rho_f *alphaM);
              
              //C2 and C2^T
              MCK(i*Num_Dof+m, j*Num_Dof+n+4) = - Ctemp;
              MCK(j*Num_Dof+n+4, i*Num_Dof+m) = - Ctemp;
            }
        }
      }
    }

    return MCK;
}

//======================================================================
const Matrix& TwentyNodeBrick_u_p_U::getMass (void)
{
    MCK.Zero();  // necessary
    
	tensor tM = getMassTensorMsf();

    int i, j;
    double Mtemp1 = 0.0;
    double Mtemp2 = 0.0;

    for ( i=0 ; i<Num_Nodes; i++ ) {
      for ( j=0; j<Num_Nodes; j++ ) {
                
        //Ms, Note *(1.0-poro)*rho_s here!
        Mtemp1 = tM.cval(i+1, j+1) *(1.0-poro)*rho_s;

        MCK(i*Num_Dof+0, j*Num_Dof+0) = Mtemp1;
        MCK(i*Num_Dof+1, j*Num_Dof+1) = Mtemp1;
        MCK(i*Num_Dof+2, j*Num_Dof+2) = Mtemp1;
        
        //Mf, Note *poro*rho_f here!
        Mtemp2 = tM.cval(i+1, j+1) *poro*rho_f;

        MCK(i*Num_Dof+4, j*Num_Dof+4) = Mtemp2;
        MCK(i*Num_Dof+5, j*Num_Dof+5) = Mtemp2;
        MCK(i*Num_Dof+6, j*Num_Dof+6) = Mtemp2;
      }
    }

    return MCK;
}

//======================================================================
void TwentyNodeBrick_u_p_U::zeroLoad()
{
   if ( Q != 0 )
     Q->Zero();
}

//======================================================================
int TwentyNodeBrick_u_p_U::addLoad(ElementalLoad *theLoad, double loadFactor)
{
   int type;
   const Vector& data = theLoad->getData(type, loadFactor);

   if ( type == LOAD_TAG_BrickSelfWeight ) {    
     Vector Fbody = this->getBodyForce();     
     if ( Q == 0 )
       Q = new Vector(Num_ElemDof);     
     *Q = Fbody * loadFactor;
   }
   
   else {
     opserr << "TwentyNodeBrick_u_p_U::addLoad() " << this->getTag() << ", load type unknown\n";
     return -1;
   }

   return 0;
}

//======================================================================
int TwentyNodeBrick_u_p_U::addInertiaLoadToUnbalance(const Vector &accel)
{
    static Vector avu(Num_ElemDof);

    int i, ik;

    for (i=0; i<Num_Nodes; i++) {
      const Vector &RA = theNodes[i]->getRV(accel);

      if ( RA.Size() != Num_Dof ) {
        opserr << "TwentyNodeBrick_u_p_U::addInertiaLoadToUnbalance matrix and vector sizes are incompatable\n";
        return (-1);
      }

      ik = i*Num_Dof;

      avu(ik +0) = RA(0);
      avu(ik +1) = RA(1);
      avu(ik +2) = RA(2);
      avu(ik +3) = 0.0;
      avu(ik +4) = RA(4);
      avu(ik +5) = RA(5);
      avu(ik +6) = RA(6);
    }

    if (Q == 0)
      Q = new Vector(Num_ElemDof);
    
    this->getMass();
    Q->addMatrixVector(1.0, MCK, avu, -1.0);

    return 0;
}

//========================================================================
const Vector& TwentyNodeBrick_u_p_U::getResistingForce ()
{            
    static Vector avu(Num_ElemDof);

    P.Zero();
       
    // Internal
    P = this->getInternalForce();

    // + K0*u
    int i, j;
    for (i=0; i<Num_Nodes; i++) {      
      const Vector &disp = theNodes[i]->getTrialDisp();
      if ( disp.Size() != Num_Dof ) {
        opserr << "EightNode_Brick_u_p_U::getResistingForce(): matrix and vector sizes are incompatable \n";
        exit(-1);
      }
      for (j=0; j<Num_Dof; j++) {
        avu(i*Num_Dof +j) = disp(j);
      }
    }

    this->getStiff00();
    P.addMatrixVector(1.0, MCK, avu, 1.0);

    if (Q != 0)
      P.addVector(1.0, *Q, -1.0);

    return P;
}

//========================================================================
const Vector& TwentyNodeBrick_u_p_U::getResistingForceIncInertia ()
{
    static Vector avu(Num_ElemDof);
    
    this->getResistingForce();    

    // + M*a
    int i, j;
    for (i=0; i<Num_Nodes; i++) {
      const Vector &acc = theNodes[i]->getTrialAccel();
      if ( acc.Size() != Num_Dof ) {
        opserr << "EightNode_Brick_u_p_U::getResistingForceIncInertia matrix and vector sizes are incompatable \n";
        exit(-1);
      }
      for (j=0; j<Num_Dof; j++) {
        avu(i*Num_Dof +j) = acc(j);
      }
    }

    this->getMass();
    P.addMatrixVector(1.0, MCK, avu, 1.0);

    // + C*v
    for (i=0; i<Num_Nodes; i++) {
      const Vector &vel = theNodes[i]->getTrialVel();
      if ( vel.Size() != Num_Dof ) {
        opserr << "EightNode_Brick_u_p_U::getResistingForceIncInertia matrix and vector sizes are incompatable \n";
        exit(-1);											       
      }
      for (j=0; j<Num_Dof; j++) {
        avu(i*Num_Dof +j) = vel(j);
      }
    }

    this->getDamp();
    P.addMatrixVector(1.0, MCK, avu, 1.0);

    return P;
}

//=============================================================================
int TwentyNodeBrick_u_p_U::sendSelf (int commitTag, Channel &theChannel)
{
     // Not implemented yet
     return 0;
}

//=============================================================================
int TwentyNodeBrick_u_p_U::recvSelf (int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
     // Not implemented yet
     return 0;
}

//=============================================================================
int TwentyNodeBrick_u_p_U::displaySelf (Renderer &theViewer, int displayMode, float fact)
{
     // Not implemented yet
     return 0;
}


//=============================================================================
Response* TwentyNodeBrick_u_p_U::setResponse(const char **argv, int argc, Information &eleInfo, OPS_Stream &output)
{
  Response *theResponse = 0;

  char outputData[32];
  output.tag("ElementOutput");
  output.attr("eleType","TwentyNodeBrick_u_p_U");
  output.attr("eleTag",this->getTag());
  for (int i=1; i<=Num_Nodes; i++) {
    sprintf(outputData,"node%d",i);
    output.attr(outputData,connectedExternalNodes[i-1]);
  }

  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0) {

    for (int i=1; i<=Num_Nodes; i++)
      for (int j=1; j<=Num_Dof; j++) {
	sprintf(outputData,"P%d_%d",j,i);
	output.tag("ResponseType",outputData);
      }

    theResponse = new ElementResponse(this, 1, this->getResistingForce());
    
  } else if (strcmp(argv[0],"material") == 0 || strcmp(argv[0],"integrPoint") == 0) { 
    int pointNum = atoi(argv[1]);
    if (pointNum > 0 && pointNum <= Num_TotalGaussPts) {

      output.tag("GaussPoint");
      output.attr("number",pointNum);

      theResponse =  theMaterial[pointNum-1]->setResponse(&argv[2], argc-2, eleInfo, output);

      output.endTag(); // GaussPoint
    }
  }

  else if (strcmp(argv[0],"stresses") ==0) {

    for (int i=0; i<Num_TotalGaussPts; i++) {
      output.tag("GaussPoint");
      output.attr("number",i+1);
      output.tag("NdMaterialOutput");
      output.attr("classType", theMaterial[i]->getClassTag());
      output.attr("tag", theMaterial[i]->getTag());

      output.tag("ResponseType","sigma11");
      output.tag("ResponseType","sigma22");
      output.tag("ResponseType","sigma33");
      output.tag("ResponseType","sigma23");
      output.tag("ResponseType","sigma31");
      output.tag("ResponseType","sigma12");      

      output.endTag(); // NdMaterialOutput
      output.endTag(); // GaussPoint
    }
     
    theResponse = new ElementResponse(this, 5, Vector(Num_TotalGaussPts*6) );

  } else if (strcmp(argv[0],"gausspoint") == 0 || strcmp(argv[0],"GaussPoint") == 0) {
     theResponse = new ElementResponse(this, 6, Vector(Num_TotalGaussPts*Num_Dim) );
  }
  output.endTag(); // ElementOutput

  return theResponse; 

}

//=============================================================================
int TwentyNodeBrick_u_p_U::getResponse(int responseID, Information &eleInfo)
{
  if (responseID == 1)
    return eleInfo.setVector(this->getResistingForce());

  else if (responseID == 5) {
    static Vector stresses(Num_TotalGaussPts*6);
    stresstensor sigma;
    int cnt = 0;
    int i;
    for (i=0; i<Num_TotalGaussPts; i++) {
      sigma = theMaterial[i]->getStressTensor();
      stresses(cnt++) = sigma.cval(1,1);  //xx
      stresses(cnt++) = sigma.cval(2,2);  //yy
      stresses(cnt++) = sigma.cval(3,3);  //zz
      stresses(cnt++) = sigma.cval(2,3);  //yz
      stresses(cnt++) = sigma.cval(3,1);  //zx
      stresses(cnt++) = sigma.cval(1,2);  //xy
    }
    return eleInfo.setVector(stresses);
  }

  else if (responseID == 6) {
    static Vector Gpts(Num_TotalGaussPts*Num_Dim);
    tensor GCoord;
    int cnt = 0;
    int i,j;
    GCoord = getGaussPts();
    for (i=0; i<Num_TotalGaussPts; i++) {
      for (j=0; j<Num_Dim; j++) {
        Gpts(cnt++) = GCoord.cval(i+1,j+1);
      }
    }
    return eleInfo.setVector(Gpts);
  }

  else if (responseID == 7) {
    static Vector Gpts(Num_TotalGaussPts*2);
    stresstensor sigma;
    int i;
    for (i=0; i<Num_TotalGaussPts; i++) {
        sigma = theMaterial[i]->getStressTensor();
        Gpts(i*2   ) = sigma.p_hydrostatic();
        Gpts(i*2 +1) = sigma.q_deviatoric(); 
    }
    return eleInfo.setVector(Gpts);
  }

  else
    return (-1);
}

//=============================================================================
void TwentyNodeBrick_u_p_U::Print(OPS_Stream &s, int flag)
{
    s << "TwentyNodeBrick_u_p_U, element id:  " << this->getTag() << "\n";
    s << "Connected external nodes:  " << connectedExternalNodes << "\n";

    s << "Node 1: " << connectedExternalNodes(0) << "\n";
    s << "Node 2: " << connectedExternalNodes(1) << "\n";
    s << "Node 3: " << connectedExternalNodes(2) << "\n";
    s << "Node 4: " << connectedExternalNodes(3) << "\n";
    s << "Node 5: " << connectedExternalNodes(4) << "\n";
    s << "Node 6: " << connectedExternalNodes(5) << "\n";
    s << "Node 7: " << connectedExternalNodes(6) << "\n";
    s << "Node 8: " << connectedExternalNodes(7) << "\n";
    s << "Node 9: " << connectedExternalNodes(8) << "\n";
    s << "Node 10: " << connectedExternalNodes(9) << "\n";
    s << "Node 11: " << connectedExternalNodes(10) << "\n";
    s << "Node 12: " << connectedExternalNodes(11) << "\n";
    s << "Node 13: " << connectedExternalNodes(12) << "\n";
    s << "Node 14: " << connectedExternalNodes(13) << "\n";
    s << "Node 15: " << connectedExternalNodes(14) << "\n";
    s << "Node 16: " << connectedExternalNodes(15) << "\n";
    s << "Node 17: " << connectedExternalNodes(16) << "\n";
    s << "Node 18: " << connectedExternalNodes(17) << "\n";
    s << "Node 19: " << connectedExternalNodes(18) << "\n";
    s << "Node 20: " << connectedExternalNodes(19) << "\n";

    s << "Material model:  " << "\n";

    int GP_c_r, GP_c_s, GP_c_t, where;

    for( GP_c_r = 0 ; GP_c_r < Num_IntegrationPts ; GP_c_r++ ) {
      for( GP_c_s = 0 ; GP_c_s < Num_IntegrationPts ; GP_c_s++ ) {
        for( GP_c_t = 0 ; GP_c_t < Num_IntegrationPts ; GP_c_t++ ) {
          where = (GP_c_r*Num_IntegrationPts+GP_c_s)*Num_IntegrationPts+GP_c_t;
          s << "\n where = " << where+1 << "\n";
          s << " r= " << GP_c_r << " s= " << GP_c_s << " t= " << GP_c_t << "\n";
          theMaterial[where]->Print(s);
        }
      }
    }

}

//======================================================================
int TwentyNodeBrick_u_p_U::update()
{
    int ret = 0;

    double r  = 0.0;
    double s  = 0.0;
    double t  = 0.0;

    int Tdisp_dim[] = {Num_Nodes, Num_Dof};
    tensor total_displacements(2, Tdisp_dim, 0.0);
    int tdisp_dim[] = {Num_Nodes, Num_Dim};
    tensor total_disp(2, tdisp_dim, 0.0);

    int dh_dim[] = {Num_Nodes, Num_Dim};
    tensor dh(2, dh_dim, 0.0);

    straintensor eps;

    tensor dhGlobal;

    total_displacements = getNodesDisp();
    int i;
    for (i=1; i<=Num_Nodes; i++) {
      total_disp.val(i,1) = total_displacements.cval(i,1);
      total_disp.val(i,2) = total_displacements.cval(i,2);
      total_disp.val(i,3) = total_displacements.cval(i,3);
    }

    int GP_c_r, GP_c_s, GP_c_t, where;

    for( GP_c_r = 0 ; GP_c_r < Num_IntegrationPts; GP_c_r++ ) {
      r = pts[GP_c_r];
      for( GP_c_s = 0 ; GP_c_s < Num_IntegrationPts; GP_c_s++ ) {
        s = pts[GP_c_s];
        for( GP_c_t = 0 ; GP_c_t < Num_IntegrationPts; GP_c_t++ ) {
          t = pts[GP_c_t];
          where = (GP_c_r*Num_IntegrationPts+GP_c_s)*Num_IntegrationPts+GP_c_t;
          dhGlobal = dh_Global(r,s,t);
          eps = total_disp("ia") * dhGlobal("ib");
          eps.null_indices();
          eps.symmetrize11();
          if ( (theMaterial[where]->setTrialStrain(eps) ) )
            opserr << "TwentyNodeBrick_u_p_U::update (tag: " << this->getTag() << "), not converged\n";
        }
      }
    }

  return ret;
}

//======================================================================
const Vector& TwentyNodeBrick_u_p_U::getInternalForce ()
{
    static Vector Pforce(Num_ElemDof);

    double r  = 0.0;
    double rw = 0.0;
    double s  = 0.0;
    double sw = 0.0;
    double t  = 0.0;
    double tw = 0.0;
    double weight = 0.0;
    double det_of_Jacobian = 0.0;
    int where = 0;

    int dh_dim[] = {Num_Nodes,Num_Dim};
    tensor dh(2, dh_dim, 0.0);
    tensor Pins(2, dh_dim, 0.0);

    tensor Jacobian;
    stresstensor Mstress;

    int GP_c_r, GP_c_s, GP_c_t;
    for( GP_c_r = 0 ; GP_c_r < Num_IntegrationPts; GP_c_r++ ) {
      r = pts[GP_c_r];
      rw = wts[GP_c_r];
      for( GP_c_s = 0 ; GP_c_s < Num_IntegrationPts; GP_c_s++ ) {
        s = pts[GP_c_s];
        sw = wts[GP_c_s];
        for( GP_c_t = 0 ; GP_c_t < Num_IntegrationPts; GP_c_t++ ) {
          t = pts[GP_c_t];
          tw = wts[GP_c_t];
          where = (GP_c_r*Num_IntegrationPts+GP_c_s)*Num_IntegrationPts+GP_c_t;
          Jacobian = Jacobian_3D(r,s,t);
          det_of_Jacobian = Jacobian.determinant();
          dh = dh_Global(r,s,t);
          weight = rw * sw * tw * det_of_Jacobian;
          Mstress = this->theMaterial[where]->getStressTensor();
          Pins += (dh("Kj")*Mstress("ij"))*weight;
        }
      }
    }
    
    Pforce.Zero(); // necessary
    
    int i, j;
    for (i=0; i<Num_Nodes; i++) {
      for (j=0; j<Num_Dim; j++) {
        Pforce(i*Num_Dof+j) = Pins.cval(i+1, j+1);
      }
    }

    return Pforce;
}


//======================================================================
const Vector& TwentyNodeBrick_u_p_U::getBodyForce ()
{
    static Vector Pforce(Num_ElemDof);

    double r  = 0.0;
    double rw = 0.0;
    double s  = 0.0;
    double sw = 0.0;
    double t  = 0.0;
    double tw = 0.0;
    double weight = 0.0;
    double det_of_Jacobian = 0.0;

    int hp_dim[] = {Num_Nodes};
    tensor hp(1, hp_dim, 0.0);
    tensor Pexf(1, hp_dim, 0.0);

    tensor Jacobian;

    int GP_c_r, GP_c_s, GP_c_t;

    for( GP_c_r = 0 ; GP_c_r < Num_IntegrationPts; GP_c_r++ ) {
      r = pts[GP_c_r];
      rw = wts[GP_c_r];
      for( GP_c_s = 0 ; GP_c_s < Num_IntegrationPts; GP_c_s++ ) {
        s = pts[GP_c_s];
        sw = wts[GP_c_s];
        for( GP_c_t = 0 ; GP_c_t < Num_IntegrationPts; GP_c_t++ ) {
          t = pts[GP_c_t];
          tw = wts[GP_c_t];
          hp = shapeFunction(r,s,t);
          Jacobian = Jacobian_3D(r,s,t);
          det_of_Jacobian = Jacobian.determinant();
          weight = rw * sw * tw * det_of_Jacobian;
          Pexf += hp *weight;
        }
      }
    }

    Pforce.Zero(); // necessary
    
    int i, j;
    for (i=0; i<Num_Nodes; i++) {
      for (j=0; j<Num_Dim; j++) {
        Pforce(i*Num_Dof +j) = Pexf.cval(i+1) * bf(j) * (1.0-poro) * rho_s;
        Pforce(i*Num_Dof +j +4) = Pexf.cval(i+1) * bf(j) * poro * rho_f;
      }
    }

    return Pforce;
}


//======================================================================
const Matrix& TwentyNodeBrick_u_p_U::getStiff00 (void)
{
    MCK.Zero();  // necessary
	
    int i, j, m;

	tensor tG   = getStiffnessTensorG12();

    //G1 and G1^T, Note *(alpha-poro) here!
    for ( i=0 ; i<Num_Nodes; i++ ) {
      for ( j=0; j<Num_Nodes; j++ ) {
        for( m=0; m<Num_Dim; m++)  {
            MCK(i*Num_Dof+m, j*Num_Dof+3) = -tG.cval(i+1, m+1, j+1) *(alpha-poro);
            MCK(j*Num_Dof+3, i*Num_Dof+m) = -tG.cval(i+1, m+1, j+1) *(alpha-poro);
        }
      }
    }

    //G2 and G2^T, Note *poro here!
    for ( i=0 ; i<Num_Nodes; i++ ) {
      for ( j=0; j<Num_Nodes; j++ ) {
        for( m=0; m<Num_Dim; m++) {
            MCK(i*Num_Dof+m+4, j*Num_Dof+3) = -tG.cval(i+1, m+1, j+1) *poro;
            MCK(j*Num_Dof+3, i*Num_Dof+m+4) = -tG.cval(i+1, m+1, j+1) *poro;
        }
      }
    }


    //P
    if (ks == 0.0 || kf == 0.0) {
       opserr<<" Error, TwentyNodeBrick_u_p_U::getStiffnessTensorP -- solid and/or fluid bulk modulus is zero\n";
       exit(-1);
    }

    double  oneOverQ = poro/kf + (alpha-poro)/ks;

    tensor tP = getMassTensorMsf();

    for ( i=0 ; i<Num_Nodes; i++ ) {
      for ( j=0; j<Num_Nodes; j++ ) {
        MCK(i*Num_Dof+3, j*Num_Dof+3) = tP.cval(i+1, j+1) * (-oneOverQ);
      }
    }

    return MCK;
}

//======================================================================
const Matrix& TwentyNodeBrick_u_p_U::getStiff (int Ki_flag)
{
    if (Ki_flag != 0 && Ki_flag != 1) {
      opserr << "Error TwentyNodeBrick_u_p_U::getStiff() - illegal use\n";
      exit(-1);
    }

    if (Ki_flag == 0 && Ki != 0)
      return *Ki;

    tensor tKep = getStiffnessTensorKep();

    int i, j, m, n;
    
    this->getStiff00();

    // + Kep
    for ( i=0 ; i<Num_Nodes; i++ ) {
      for ( j=0; j<Num_Nodes; j++ ) {
        for( m=0; m<Num_Dim; m++) {
          for( n=0; n<Num_Dim; n++) 
            MCK(i*Num_Dof+m, j*Num_Dof+n) = tKep.cval(i+1, m+1, n+1, j+1);
        }
      }
    }

    if( Ki_flag == 1)
      return MCK;

    Ki = new Matrix(MCK);

    if (Ki == 0) {
      opserr << "Error TwentyNodeBrick_u_p_U::getStiff() -";
      opserr << "ran out of memory\n";
      exit(-1);
    }

    return *Ki;
}

//======================================================================
tensor TwentyNodeBrick_u_p_U::getStiffnessTensorKep( )
{
    int K_dim[] = {Num_Nodes, Num_Dim, Num_Dim, Num_Nodes};
    tensor Kep(4, K_dim, 0.0);
    tensor Kkt(4, K_dim, 0.0);

    double r  = 0.0;
    double rw = 0.0;
    double s  = 0.0;
    double sw = 0.0;
    double t  = 0.0;
    double tw = 0.0;
    int where = 0;
    double weight = 0.0;
    double det_of_Jacobian = 0.0;

    int dh_dim[] = {Num_Nodes,Num_Dim};
    tensor dh(2, dh_dim, 0.0);

    tensor Constitutive;

    tensor Jacobian;
    tensor dhGlobal;

    int GP_c_r, GP_c_s, GP_c_t;

    for( GP_c_r = 0 ; GP_c_r < Num_IntegrationPts; GP_c_r++ ) {
      r = pts[GP_c_r];
      rw = wts[GP_c_r];
      for( GP_c_s = 0 ; GP_c_s < Num_IntegrationPts; GP_c_s++ ) {
        s = pts[GP_c_s];
        sw = wts[GP_c_s];
        for( GP_c_t = 0 ; GP_c_t < Num_IntegrationPts; GP_c_t++ ) {
          t = pts[GP_c_t];
          tw = wts[GP_c_t];
          where = (GP_c_r*Num_IntegrationPts+GP_c_s)*Num_IntegrationPts+GP_c_t;
          Jacobian = Jacobian_3D(r,s,t);
          det_of_Jacobian = Jacobian.determinant();
          dhGlobal = dh_Global(r,s,t);
          weight = rw * sw * tw * det_of_Jacobian;
          Constitutive = theMaterial[where]->getTangentTensor();
          Kkt = dhGlobal("kj")*Constitutive("ijml");
          Kkt = Kkt("kiml")*dhGlobal("pl")*weight;
          Kep = Kep + Kkt;
        }
      }
    }

    return Kep;
}

//======================================================================
tensor TwentyNodeBrick_u_p_U::getStiffnessTensorG12()
{
    // This is for G1 and G2
    // G1 = (alpha-poro) *G;
    // G2 = poro *G;

    int G_dim[] = {Num_Nodes, Num_Dim, Num_Nodes};
    tensor G(3, G_dim, 0.0);

    double r  = 0.0;
    double rw = 0.0;
    double s  = 0.0;
    double sw = 0.0;
    double t  = 0.0;
    double tw = 0.0;
    double weight = 0.0;

    int dh_dim[] = {Num_Nodes,Num_Dim};
    tensor dh(2, dh_dim, 0.0);
    int hp_dim[] = {Num_Nodes};
    tensor hp(1, hp_dim, 0.0);

    double det_of_Jacobian = 0.0;

    tensor Jacobian;
    tensor dhGlobal;

    int GP_c_r, GP_c_s, GP_c_t;

    for( GP_c_r = 0 ; GP_c_r < Num_IntegrationPts; GP_c_r++ ) {
      r = pts[GP_c_r];
      rw = wts[GP_c_r];
      for( GP_c_s = 0 ; GP_c_s < Num_IntegrationPts; GP_c_s++ ) {
        s = pts[GP_c_s];
        sw = wts[GP_c_s];
        for( GP_c_t = 0 ; GP_c_t < Num_IntegrationPts; GP_c_t++ ) {
          t = pts[GP_c_t];
          tw = wts[GP_c_t];
          hp= shapeFunction(r,s,t);
          Jacobian = Jacobian_3D(r,s,t);
          dhGlobal = dh_Global(r,s,t);
          det_of_Jacobian = Jacobian.determinant();
          weight = rw * sw * tw * det_of_Jacobian;
          G += dhGlobal("ki")*hp("m") * weight;
        }
      }
    }

    return G;
}

//======================================================================
tensor TwentyNodeBrick_u_p_U::getDampTensorC123()
{
    // This is for C1, C2 and C3, C1 = C2 = c3
    // Since solid and fluid shape function the same

    int perm_dim[] = {Num_Dim, Num_Dim};
    tensor perm_inv(2, perm_dim, 0.0);
    
    perm_inv.val(1,1) = 1.0/perm(0);
    perm_inv.val(2,2) = 1.0/perm(1);
    perm_inv.val(3,3) = 1.0/perm(2);

    int C_dim[] = {Num_Nodes,Num_Dim,Num_Dim,Num_Nodes};
    tensor C123(4,C_dim,0.0);
    int c_dim[] = {Num_Nodes,Num_Dim,Num_Dim};
    tensor c123(3,c_dim,0.0);

    double r  = 0.0;
    double rw = 0.0;
    double s  = 0.0;
    double sw = 0.0;
    double t  = 0.0;
    double tw = 0.0;
    double weight = 0.0;
    double det_of_Jacobian = 0.0;

    int hp_dim[] = {Num_Nodes};
    tensor hp(1, hp_dim,0.0);

    tensor Jacobian;

    int GP_c_r, GP_c_s, GP_c_t;

    for( GP_c_r = 0 ; GP_c_r < Num_IntegrationPts; GP_c_r++ ) {
      r = pts[GP_c_r];
      rw = wts[GP_c_r];
      for( GP_c_s = 0 ; GP_c_s < Num_IntegrationPts; GP_c_s++ ) {
        s = pts[GP_c_s];
        sw = wts[GP_c_s];
        for( GP_c_t = 0 ; GP_c_t < Num_IntegrationPts; GP_c_t++ ) {
          t = pts[GP_c_t];
          tw = wts[GP_c_t];
          hp = shapeFunction(r,s,t);
          Jacobian = Jacobian_3D(r,s,t);
          det_of_Jacobian = Jacobian.determinant();
          weight = rw * sw * tw * det_of_Jacobian;
          c123 = hp("k")*perm_inv("ij");
          C123 += c123("kij")*hp("m") *weight;
        }
      }
    }

    return C123;
}

//======================================================================
tensor TwentyNodeBrick_u_p_U::getMassTensorMsf()
{
    // This is for Ms and Mf -> M_kl    
    // Ms = Msf * (1.0-poro)*rho_s
    // Mf = Msf * poro*rho_f

    // Also this is for Compression term Pc
    // Pc = Msf * oneOverQ
    //    = Msf * (poro/kf + (alpha-poro)/ks)

    int M_dim[] = {Num_Nodes, Num_Nodes};
    tensor Msf(2, M_dim, 0.0);

    double r  = 0.0;
    double rw = 0.0;
    double s  = 0.0;
    double sw = 0.0;
    double t  = 0.0;
    double tw = 0.0;
    double weight = 0.0;
    double det_of_Jacobian = 0.0;

    int dh_dim[] = {Num_Nodes,Num_Dim};
    tensor dh(2, dh_dim, 0.0);
    int hp_dim[] = {Num_Nodes};
    tensor hp(1, hp_dim,0.0);

    tensor Jacobian;

    int GP_c_r, GP_c_s, GP_c_t;

    for( GP_c_r = 0 ; GP_c_r < Num_IntegrationPts; GP_c_r++ ) {
      r = pts[GP_c_r];
      rw = wts[GP_c_r];
      for( GP_c_s = 0 ; GP_c_s < Num_IntegrationPts; GP_c_s++ ) {
        s = pts[GP_c_s];
        sw = wts[GP_c_s];
        for( GP_c_t = 0 ; GP_c_t < Num_IntegrationPts; GP_c_t++ ) {
          t = pts[GP_c_t];
          tw = wts[GP_c_t];
          hp = shapeFunction(r,s,t);
          Jacobian = Jacobian_3D(r,s,t); 
          det_of_Jacobian = Jacobian.determinant();
          weight = rw * sw * tw * det_of_Jacobian;
          Msf += hp("m")*hp("n")*weight;
        }
      }
    }

    return Msf;
}

//======================================================================
tensor TwentyNodeBrick_u_p_U::Jacobian_3D(double x, double y, double z)
{
     tensor N_C = this->getNodesCrds();
     tensor dh = this->shapeFunctionDerivative(x, y, z);
     
     tensor J3D = N_C("ki") * dh("kj");
       J3D.null_indices();
     return J3D;
}

//======================================================================
tensor TwentyNodeBrick_u_p_U::Jacobian_3Dinv(double x, double y, double z)
{
     tensor J = this->Jacobian_3D(x,y,z);
     return J.inverse();
}

//======================================================================
tensor TwentyNodeBrick_u_p_U::dh_Global(double x, double y, double z)
{
     tensor JacobianINV0 = this->Jacobian_3Dinv(x, y, z);
     tensor dh = this->shapeFunctionDerivative(x, y, z);
     tensor  dhGlobal = dh("ik") * JacobianINV0("kj");
       dhGlobal.null_indices();
     return dhGlobal;
}

//======================================================================
tensor TwentyNodeBrick_u_p_U::getNodesCrds(void)
{
  int i,j;
  int dimX[] = {Num_Nodes,Num_Dim};
  tensor N_coord(2, dimX, 0.0);

  for (i=0; i<Num_Nodes; i++) {
    const Vector&TNodesCrds = theNodes[i]->getCrds();
    for (j=0; j<Num_Dim; j++) {
      N_coord.val(i+1,j+1) = TNodesCrds(j);
    }
  }

  return N_coord;
}

//======================================================================
tensor TwentyNodeBrick_u_p_U::getNodesDisp(void)
{
  int i,j;
  int dimU[] = {Num_Nodes, Num_Dof};
  tensor total_disp(2, dimU, 0.0);

  for (i=0; i<Num_Nodes; i++) {
    const Vector&TNodesDisp = theNodes[i]->getTrialDisp();
    for (j=0; j<Num_Dof; j++) {
      total_disp.val(i+1,j+1) = TNodesDisp(j);
    }
  }

  return total_disp;
}

//======================================================================
double TwentyNodeBrick_u_p_U::getPorePressure(double x1, double x2, double x3)
{
  double pp = 0.0;
  int i;

  for (i=0; i<Num_Nodes; i++) {
    const Vector& T_disp = theNodes[i]->getTrialDisp();
    pp += shapeFunction(x1,x2,x3).cval(i+1) * T_disp(3);
  }

  return pp;
}

//======================================================================
tensor TwentyNodeBrick_u_p_U::shapeFunction(double r1, double r2, double r3)
{
    int Hfun[] = {Num_Nodes};
    tensor h(1, Hfun, 0.0);

      // influence of the node number 20
    h.val(20)=(1.0+r1)*(1.0-r2)*(1.0-r3*r3)/4.0;
      // influence of the node number 19
    h.val(19)=(1.0-r1)*(1.0-r2)*(1.0-r3*r3)/4.0;
      // influence of the node number 18
    h.val(18)=(1.0-r1)*(1.0+r2)*(1.0-r3*r3)/4.0;
      // influence of the node number 17
    h.val(17)=(1.0+r1)*(1.0+r2)*(1.0-r3*r3)/4.0;

      // influence of the node number 16
    h.val(16)=(1.0+r1)*(1.0-r2*r2)*(1.0-r3)/4.0;
      // influence of the node number 15
    h.val(15)=(1.0-r1*r1)*(1.0-r2)*(1.0-r3)/4.0;
      // influence of the node number 14
    h.val(14)=(1.0-r1)*(1.0-r2*r2)*(1.0-r3)/4.0;
      // influence of the node number 13
    h.val(13)=(1.0-r1*r1)*(1.0+r2)*(1.0-r3)/4.0;

      // influence of the node number 12
    h.val(12)=(1.0+r1)*(1.0-r2*r2)*(1.0+r3)/4.0;
      // influence of the node number 11
    h.val(11)=(1.0-r1*r1)*(1.0-r2)*(1.0+r3)/4.0;
      // influence of the node number 10
    h.val(10)=(1.0-r1)*(1.0-r2*r2)*(1.0+r3)/4.0;
      // influence of the node number 9
    h.val( 9)=(1.0-r1*r1)*(1.0+r2)*(1.0+r3)/4.0;

      // influence of the node number 8
    h.val(8)=(1.0+r1)*(1.0-r2)*(1.0-r3)/8.0 - (h.val(15)+h.val(16)+h.val(20))/2.0;
      // influence of the node number 7
    h.val(7)=(1.0-r1)*(1.0-r2)*(1.0-r3)/8.0 - (h.val(14)+h.val(15)+h.val(19))/2.0;
      // influence of the node number 6
    h.val(6)=(1.0-r1)*(1.0+r2)*(1.0-r3)/8.0 - (h.val(13)+h.val(14)+h.val(18))/2.0;
      // influence of the node number 5
    h.val(5)=(1.0+r1)*(1.0+r2)*(1.0-r3)/8.0 - (h.val(13)+h.val(16)+h.val(17))/2.0;

      // influence of the node number 4
    h.val(4)=(1.0+r1)*(1.0-r2)*(1.0+r3)/8.0 - (h.val(11)+h.val(12)+h.val(20))/2.0;
      // influence of the node number 3
    h.val(3)=(1.0-r1)*(1.0-r2)*(1.0+r3)/8.0 - (h.val(10)+h.val(11)+h.val(19))/2.0;
      // influence of the node number 2
    h.val(2)=(1.0-r1)*(1.0+r2)*(1.0+r3)/8.0 - (h.val(10)+h.val(18)+h.val(9))/2.0;
      // influence of the node number 1
    h.val(1)=(1.0+r1)*(1.0+r2)*(1.0+r3)/8.0 - (h.val(12)+h.val(17)+h.val(9))/2.0;

    return h;
}


//==============================================================
tensor TwentyNodeBrick_u_p_U::shapeFunctionDerivative(double r1, double r2, double r3)
{
    int DHfun[] = {Num_Nodes,Num_Dim};
    tensor dh(2, DHfun, 0.0);

      // influence of the node number 20
    dh.val(20,1) =   (1.0-r2)*(1.0-r3*r3)/4.0;
    dh.val(20,2) = - (1.0+r1)*(1.0-r3*r3)/4.0;
    dh.val(20,3) = - (1.0+r1)*(1.0-r2)*r3/2.0;
      // influence of the node number 19
    dh.val(19,1) = - (1.0-r2)*(1.0-r3*r3)/4.0;
    dh.val(19,2) = - (1.0-r1)*(1.0-r3*r3)/4.0;
    dh.val(19,3) = - (1.0-r1)*(1.0-r2)*r3/2.0;
      // influence of the node number 18
    dh.val(18,1) = - (1.0+r2)*(1.0-r3*r3)/4.0;
    dh.val(18,2) =   (1.0-r1)*(1.0-r3*r3)/4.0;
    dh.val(18,3) = - (1.0-r1)*(1.0+r2)*r3/2.0;
      // influence of the node number 17
    dh.val(17,1) =   (1.0+r2)*(1.0-r3*r3)/4.0;
    dh.val(17,2) =   (1.0+r1)*(1.0-r3*r3)/4.0;
    dh.val(17,3) = - (1.0+r1)*(1.0+r2)*r3/2.0;

      // influence of the node number 16
    dh.val(16,1) =   (1.0-r2*r2)*(1.0-r3)/4.0;
    dh.val(16,2) = - (1.0+r1)*r2*(1.0-r3)/2.0;
    dh.val(16,3) = - (1.0+r1)*(1.0-r2*r2)/4.0;
      // influnce of the node number 15
    dh.val(15,1) = - r1*(1.0-r2)*(1.0-r3)/2.0;
    dh.val(15,2) = - (1.0-r1*r1)*(1.0-r3)/4.0;
    dh.val(15,3) = - (1.0-r1*r1)*(1.0-r2)/4.0;
      // influence of the node number 14
    dh.val(14,1) = - (1.0-r2*r2)*(1.0-r3)/4.0;
    dh.val(14,2) = - (1.0-r1)*r2*(1.0-r3)/2.0;
    dh.val(14,3) = - (1.0-r1)*(1.0-r2*r2)/4.0;
      // influence of the node number 13
    dh.val(13,1) = - r1*(1.0+r2)*(1.0-r3)/2.0;
    dh.val(13,2) =   (1.0-r1*r1)*(1.0-r3)/4.0;
    dh.val(13,3) = - (1.0-r1*r1)*(1.0+r2)/4.0;

      // influence of the node number 12
    dh.val(12,1) =   (1.0-r2*r2)*(1.0+r3)/4.0;
    dh.val(12,2) = - (1.0+r1)*r2*(1.0+r3)/2.0;
    dh.val(12,3) =   (1.0+r1)*(1.0-r2*r2)/4.0;
      // influence of the node number 11
    dh.val(11,1) = - r1*(1.0-r2)*(1.0+r3)/2.0;
    dh.val(11,2) = - (1.0-r1*r1)*(1.0+r3)/4.0;
    dh.val(11,3) =   (1.0-r1*r1)*(1.0-r2)/4.0;
      // influence of the node number 10
    dh.val(10,1) = - (1.0-r2*r2)*(1.0+r3)/4.0;
    dh.val(10,2) = - (1.0-r1)*r2*(1.0+r3)/2.0;
    dh.val(10,3) =   (1.0-r1)*(1.0-r2*r2)/4.0;
      // influence of the node number 9
    dh.val(9,1)  = - r1*(1.0+r2)*(1.0+r3)/2.0;
    dh.val(9,2)  =   (1.0-r1*r1)*(1.0+r3)/4.0;
    dh.val(9,3)  =   (1.0-r1*r1)*(1.0+r2)/4.0;

      // influence of the node number 8
    dh.val(8,1)= (1.0-r2)*(1.0-r3)/8.0 - (dh.val(15,1)+dh.val(16,1)+dh.val(20,1))/2.0;
    dh.val(8,2)=-(1.0+r1)*(1.0-r3)/8.0 - (dh.val(15,2)+dh.val(16,2)+dh.val(20,2))/2.0;
    dh.val(8,3)=-(1.0+r1)*(1.0-r2)/8.0 - (dh.val(15,3)+dh.val(16,3)+dh.val(20,3))/2.0;
      // influence of the node number 7
    dh.val(7,1)=-(1.0-r2)*(1.0-r3)/8.0 - (dh.val(14,1)+dh.val(15,1)+dh.val(19,1))/2.0;
    dh.val(7,2)=-(1.0-r1)*(1.0-r3)/8.0 - (dh.val(14,2)+dh.val(15,2)+dh.val(19,2))/2.0;
    dh.val(7,3)=-(1.0-r1)*(1.0-r2)/8.0 - (dh.val(14,3)+dh.val(15,3)+dh.val(19,3))/2.0;
      // influence of the node number 6
    dh.val(6,1)=-(1.0+r2)*(1.0-r3)/8.0 - (dh.val(13,1)+dh.val(14,1)+dh.val(18,1))/2.0;
    dh.val(6,2)= (1.0-r1)*(1.0-r3)/8.0 - (dh.val(13,2)+dh.val(14,2)+dh.val(18,2))/2.0;
    dh.val(6,3)=-(1.0-r1)*(1.0+r2)/8.0 - (dh.val(13,3)+dh.val(14,3)+dh.val(18,3))/2.0;
      // influence of the node number 5
    dh.val(5,1)= (1.0+r2)*(1.0-r3)/8.0 - (dh.val(13,1)+dh.val(16,1)+dh.val(17,1))/2.0;
    dh.val(5,2)= (1.0+r1)*(1.0-r3)/8.0 - (dh.val(13,2)+dh.val(16,2)+dh.val(17,2))/2.0;
    dh.val(5,3)=-(1.0+r1)*(1.0+r2)/8.0 - (dh.val(13,3)+dh.val(16,3)+dh.val(17,3))/2.0;

      // influence of the node number 4
    dh.val(4,1)= (1.0-r2)*(1.0+r3)/8.0 - (dh.val(11,1)+dh.val(12,1)+dh.val(20,1))/2.0;
    dh.val(4,2)=-(1.0+r1)*(1.0+r3)/8.0 - (dh.val(11,2)+dh.val(12,2)+dh.val(20,2))/2.0;
    dh.val(4,3)= (1.0+r1)*(1.0-r2)/8.0 - (dh.val(11,3)+dh.val(12,3)+dh.val(20,3))/2.0;
      // influence of the node number 3
    dh.val(3,1)=-(1.0-r2)*(1.0+r3)/8.0 - (dh.val(10,1)+dh.val(11,1)+dh.val(19,1))/2.0;
    dh.val(3,2)=-(1.0-r1)*(1.0+r3)/8.0 - (dh.val(10,2)+dh.val(11,2)+dh.val(19,2))/2.0;
    dh.val(3,3)= (1.0-r1)*(1.0-r2)/8.0 - (dh.val(10,3)+dh.val(11,3)+dh.val(19,3))/2.0;
      // influence of the node number 2
    dh.val(2,1)=-(1.0+r2)*(1.0+r3)/8.0 - (dh.val(10,1)+dh.val(18,1)+dh.val(9,1))/2.0;
    dh.val(2,2)= (1.0-r1)*(1.0+r3)/8.0 - (dh.val(10,2)+dh.val(18,2)+dh.val(9,2))/2.0;
    dh.val(2,3)= (1.0-r1)*(1.0+r2)/8.0 - (dh.val(10,3)+dh.val(18,3)+dh.val(9,3))/2.0;
      // influence of the node number 1
    dh.val(1,1)= (1.0+r2)*(1.0+r3)/8.0 - (dh.val(12,1)+dh.val(17,1)+dh.val(9,1))/2.0;
    dh.val(1,2)= (1.0+r1)*(1.0+r3)/8.0 - (dh.val(12,2)+dh.val(17,2)+dh.val(9,2))/2.0;
    dh.val(1,3)= (1.0+r1)*(1.0+r2)/8.0 - (dh.val(12,3)+dh.val(17,3)+dh.val(9,3))/2.0;

    return dh;
}

//==============================================================
tensor TwentyNodeBrick_u_p_U::getGaussPts(void)
{
    int dimensions1[] = {Num_TotalGaussPts, Num_Dim};
    tensor Gs(2, dimensions1, 0.0);
    int dimensions2[] = {Num_Nodes};
    tensor shp(1, dimensions2, 0.0);

    double r = 0.0;
    double s = 0.0;
    double t = 0.0;
    int i, j, where;

    int GP_c_r, GP_c_s, GP_c_t;

    for( GP_c_r = 0 ; GP_c_r < Num_IntegrationPts; GP_c_r++ ) {
      r = pts[GP_c_r];
      for( GP_c_s = 0 ; GP_c_s < Num_IntegrationPts; GP_c_s++ ) {
        s = pts[GP_c_s];
        for( GP_c_t = 0 ; GP_c_t < Num_IntegrationPts; GP_c_t++ ) {
          t = pts[GP_c_t];
          where = (GP_c_r*Num_IntegrationPts+GP_c_s)*Num_IntegrationPts+GP_c_t;
          shp = shapeFunction(r,s,t);
          for (i=0; i<Num_Nodes; i++) {
            const Vector& T_Crds = theNodes[i]->getCrds();
            for (j=0; j<Num_Dim; j++) {
              Gs.val(where+1, j+1) += shp.cval(i+1) * T_Crds(j);
            }
          }
        }
      }
    }

    return Gs;
}


#endif



