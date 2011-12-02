//===============================================================================

//# COPYRIGHT (C): Woody's license (by BJ):

//                 ``This    source  code is Copyrighted in

//                 U.S.,  for  an  indefinite  period,  and anybody

//                 caught  using it without our permission, will be

//                 mighty good friends of ourn, cause we don't give

//                 a  darn.  Hack it. Compile it. Debug it. Run it.

//                 Yodel  it.  Enjoy it. We wrote it, that's all we

//                 wanted to do.''

//

//# PROJECT:           Object Oriented Finite Element Program

//# PURPOSE:           Finite Deformation Hyper-Elastic classes

//# CLASS:

//#

//# VERSION:           0.6_(1803398874989) (golden section)

//# LANGUAGE:          C++

//# TARGET OS:         all...

//# DESIGN:            Zhao Cheng, Boris Jeremic (jeremic@ucdavis.edu)

//# PROGRAMMER(S):     Zhao Cheng, Boris Jeremic

//#

//#

//# DATE:              Sept2005              

//# UPDATE HISTORY:   

//#                     

//#                    

//#

//===============================================================================

#ifndef TOTALLAGRANGIANFD8NODEBRICK_CPP

#define TOTALLAGRANGIANFD8NODEBRICK_CPP



#include <TotalLagrangianFD8NodeBrick.h>



const int TotalLagrangianFD8NodeBrick::NumIntegrationPts = 2;

const int TotalLagrangianFD8NodeBrick::NumTotalGaussPts = 8;

const int TotalLagrangianFD8NodeBrick::NumNodes = 8;

const int TotalLagrangianFD8NodeBrick::NumDof = 3;

const int TotalLagrangianFD8NodeBrick::NumElemDof = NumNodes*NumDof;



Matrix TotalLagrangianFD8NodeBrick::K(NumElemDof, NumElemDof);

Matrix TotalLagrangianFD8NodeBrick::M(NumElemDof, NumElemDof);

Vector TotalLagrangianFD8NodeBrick::P(NumElemDof);

const double TotalLagrangianFD8NodeBrick::pts[2] = {-0.577350269189626, +0.577350269189626};

const double TotalLagrangianFD8NodeBrick::wts[2] = {1.0, 1.0};



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



TotalLagrangianFD8NodeBrick::TotalLagrangianFD8NodeBrick(int tag,

int node_numb_1,  int node_numb_2,  int node_numb_3,  int node_numb_4,

int node_numb_5,  int node_numb_6,  int node_numb_7,  int node_numb_8,

NDMaterial &m, double b1, double b2, double b3)

:Element(tag, ELE_TAG_TotalLagrangianFD8NodeBrick ),

 theMaterial(0), connectedExternalNodes(NumNodes), Q(0), bf(NumDof), Ki(0)

{

      connectedExternalNodes( 0) = node_numb_1;

      connectedExternalNodes( 1) = node_numb_2;

      connectedExternalNodes( 2) = node_numb_3;

      connectedExternalNodes( 3) = node_numb_4;

      connectedExternalNodes( 4) = node_numb_5;

      connectedExternalNodes( 5) = node_numb_6;

      connectedExternalNodes( 6) = node_numb_7;

      connectedExternalNodes( 7) = node_numb_8;



      bf(0) = b1;

      bf(1) = b2;

      bf(2) = b3;



      theMaterial = new NDMaterial *[NumTotalGaussPts];



      if (theMaterial == 0) {

       opserr<<"FiniteDeformationElastic3D::FiniteDeformationElastic3D -- failed allocate material model pointer\n";

       exit(-1);

      }

      

      int i;

      for (i=0; i<NumTotalGaussPts; i++) {

       theMaterial[i] = m.getCopy();

       if (theMaterial[i] == 0) {

        opserr<<"FiniteDeformationElastic3D::FiniteDeformationElastic3D -- failed allocate material model pointer\n";

        exit(-1);

       }

      }



      rho = m.getRho();



      for (i=0; i<NumNodes; i++) theNodes[i] = 0;



}



//-------------------------------------------------------------------------------------------

TotalLagrangianFD8NodeBrick::TotalLagrangianFD8NodeBrick ()

:Element(0, ELE_TAG_TotalLagrangianFD8NodeBrick ),

 theMaterial(0), connectedExternalNodes(NumNodes), Q(0), bf(NumDof), Ki(0)

{    

	 int i;

     for (i=0; i<NumNodes; i++) {  

       theNodes[i] = 0;

     }



     bf(0) = 0.0;

     bf(1) = 0.0;

     bf(2) = 0.0;



     rho = 0.0;

}



//-------------------------------------------------------------------------------------------------

TotalLagrangianFD8NodeBrick::~TotalLagrangianFD8NodeBrick ()

{   

	int i;

	for (i=0; i<NumTotalGaussPts; i++) {

      if (theMaterial[i]) 

        delete theMaterial[i];

    }



    if(theMaterial) 

      delete [] theMaterial;



    if(Ki) 

      delete Ki;



}

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++





//=============================================================================

int TotalLagrangianFD8NodeBrick::getNumExternalNodes () const

{

    return NumNodes;

}



//=============================================================================

const ID& TotalLagrangianFD8NodeBrick::getExternalNodes ()

{

    return connectedExternalNodes;

}



//=============================================================================

Node **TotalLagrangianFD8NodeBrick::getNodePtrs(void)

{

    return theNodes;

}



//=============================================================================

int TotalLagrangianFD8NodeBrick::getNumDOF ()

{

    return NumElemDof;

}



//=============================================================================

void TotalLagrangianFD8NodeBrick::setDomain (Domain *theDomain)

{

    int i;



    // Check Domain is not null - invoked when object removed from a domain

    if (theDomain == 0) {

      for (i=0; i<NumNodes; i++) {

	     theNodes[i] = 0; 

      }

      return;

    }

    

    for (i=0; i<NumNodes; i++) {

	  theNodes[i] = theDomain->getNode(connectedExternalNodes(i));

	  if ( theNodes[i]==0 ) {

		  opserr << "FATAL ERROR TotalLagrangianFD8NodeBrick (tag: " << this->getTag() <<

          " ), node not found in domain\n";

          exit(-1);		  

	  }	   

    }



    this->DomainComponent::setDomain(theDomain);  // Very Important!!



    for (i=0; i<NumNodes; i++) {

	  if ( theNodes[i]->getNumberDOF() != NumDof ) {

        opserr << "FATAL ERROR TotalLagrangianFD8NodeBrick (tag: " << this->getTag() <<

        "), has differing number of DOFs at its nodes\n";

        exit(-1);		    		    

	  }	      

    }



}



//=============================================================================

int TotalLagrangianFD8NodeBrick::commitState ()

{

  int retVal = 0;



  if ((retVal = this->Element::commitState()) != 0)

    opserr << "TotalLagrangianFD8NodeBrick::commitState () - failed in base class";



  int i;

  for (i=0; i<NumTotalGaussPts; i++)

    retVal += theMaterial[i]->commitState();



  return retVal;

}



//=============================================================================

int TotalLagrangianFD8NodeBrick::revertToLastCommit ()

{

    int retVal = 0;

    

    int i;

    for (i=0; i<NumTotalGaussPts; i++)

       retVal += theMaterial[i]->revertToLastCommit();



    return retVal;

}



//=============================================================================

int TotalLagrangianFD8NodeBrick::revertToStart ()

{

    int retVal = 0;

    

    int i;

    for (i=0; i<NumTotalGaussPts; i++)

       retVal += theMaterial[i]->revertToStart();



    return retVal;



}



//=============================================================================

int TotalLagrangianFD8NodeBrick::update ()

{

    int ret = 0;

    tensor dh;

    tensor dH_dX;

    int where = 0;

    int GP_c_r, GP_c_s, GP_c_t;

    double r = 0.0;

    double s = 0.0;

    double t = 0.0;

    tensor I_ij("I", 2, def_dim_2);

    tensor currentF;

    tensor updatedF;



    tensor InitialNodesCrds = this->getNodesCrds();

    tensor CurrentNodesDisp = this->getNodesDisp();



    for( GP_c_r = 0 ; GP_c_r < NumIntegrationPts ; GP_c_r++ ) {

      r = pts[GP_c_r ];

      for( GP_c_s = 0 ; GP_c_s < NumIntegrationPts ; GP_c_s++ ) {

        s = pts[GP_c_s ];

        for( GP_c_t = 0 ; GP_c_t < NumIntegrationPts ; GP_c_t++ ) {

          t = pts[GP_c_t ];

          where =(GP_c_r * NumIntegrationPts + GP_c_s) * NumIntegrationPts + GP_c_t;

          //dh = shapeFunctionDerivative(r,s,t);

          dH_dX = this->dh_Global(r,s,t);

          currentF = CurrentNodesDisp("Ia") * dH_dX("Ib");

            currentF.null_indices();

          updatedF = currentF + I_ij;

          ret += theMaterial[where]->setTrialF(updatedF);

        }

      }

    }

    return ret;

}



//======================================================================

tensor TotalLagrangianFD8NodeBrick::Jacobian_3D(double x, double y, double z)

  {

     tensor N_C = this->getNodesCrds();

     tensor dh = this->shapeFunctionDerivative(x, y, z);

     

     tensor J3D = N_C("ki") * dh("kj");

       J3D.null_indices();

     return J3D;

  }



//======================================================================

tensor TotalLagrangianFD8NodeBrick::Jacobian_3Dinv(double x, double y, double z)

  {

     tensor J = this->Jacobian_3D(x,y,z);

     return J.inverse();

  }



//======================================================================

tensor TotalLagrangianFD8NodeBrick::dh_Global(double x, double y, double z)

  {

     tensor JacobianINV0 = this->Jacobian_3Dinv(x, y, z);

     tensor dh = this->shapeFunctionDerivative(x, y, z);

     tensor  dhGlobal = dh("ik") * JacobianINV0("kj");

       dhGlobal.null_indices();

     return dhGlobal;

  }



//======================================================================

tensor TotalLagrangianFD8NodeBrick::getStiffnessTensor(void)

  {

    tensor tI2("I", 2, def_dim_2);

	  

	int K_dim[] = {NumNodes,NumDof,NumDof,NumNodes};

    tensor Kk(4,K_dim,0.0);



    double r  = 0.0;

    double rw = 0.0;

    double s  = 0.0;

    double sw = 0.0;

    double t  = 0.0;

    double tw = 0.0;

   

    int where = 0;

    int GP_c_r, GP_c_s, GP_c_t;

    double weight = 0.0;



    int dh_dim[] = {NumNodes,NumDof};

    tensor dh(2, dh_dim, 0.0);

    stresstensor PK2Stress;

    tensor L2;



    double det_of_Jacobian = 0.0;



    tensor Jacobian;

    tensor dhGlobal;

    tensor nodesDisp;

    tensor F;

    //tensor temp01;

    tensor temp02;

    tensor temp03;

    tensor temp04; 

    tensor temp05;

    tensor temp06;



    nodesDisp = this->getNodesDisp( );



    for( GP_c_r = 0 ; GP_c_r < NumIntegrationPts ; GP_c_r++ ) {

      r = pts[GP_c_r ];

      rw = wts[GP_c_r ];

      for( GP_c_s = 0 ; GP_c_s < NumIntegrationPts ; GP_c_s++ ) {

        s = pts[GP_c_s ];

        sw = wts[GP_c_s ];

        for( GP_c_t = 0 ; GP_c_t < NumIntegrationPts ; GP_c_t++ ) {

          t = pts[GP_c_t ];

          tw = wts[GP_c_t ];

          where =(GP_c_r * NumIntegrationPts + GP_c_s) * NumIntegrationPts + GP_c_t;

          //dh = shapeFunctionDerivative(r,s,t);

          Jacobian = this->Jacobian_3D(r,s,t);

          det_of_Jacobian  = Jacobian.determinant();

          dhGlobal = this->dh_Global(r,s,t);

          weight = rw * sw * tw * det_of_Jacobian;

          PK2Stress = theMaterial[where]->getStressTensor();

          L2 = theMaterial[where]->getTangentTensor();

          F = theMaterial[where]->getF();

                        

          //K1

          temp04 = dhGlobal("Pb") * tI2("mn");

            temp04.null_indices(); 

          temp02 = PK2Stress("bd") * dhGlobal("Qd");   

            temp02.null_indices();

          temp06 = temp04("Pbmn") * temp02("bQ") * weight;

            temp06.null_indices(); 

          Kk += temp06;

                        

          //K2

          temp03 =  dhGlobal("Pb") * F("ma");

            temp03.null_indices(); 

          temp04 = F("nc") * L2("abcd");

            temp04.null_indices(); 

          temp05 = temp04("nabd") * dhGlobal("Qd"); 

            temp05.null_indices(); 

          temp06 = temp03("Pbma") * temp05("nabQ") * weight;

            temp06.null_indices(); 

          Kk += temp06;

        }

      }

    }



    return Kk;

  }



//======================================================================

tensor TotalLagrangianFD8NodeBrick::getRtensor(void)

  {

    int R_dim[] = {NumNodes,NumDof};

    tensor Rr(2,R_dim,0.0);



    double r  = 0.0;

    double rw = 0.0;

    double s  = 0.0;

    double sw = 0.0;

    double t  = 0.0;

    double tw = 0.0;



    int where = 0;

    int GP_c_r, GP_c_s, GP_c_t;

    double weight = 0.0;



    int dh_dim[] = {NumNodes,NumDof};

    tensor dh(2, dh_dim, 0.0);



    double det_of_Jacobian = 0.0;



    tensor Jacobian;

    tensor JacobianINV;

    tensor dhGlobal;

    tensor currentF;

    tensor nodesDisp;

    //stresstensor PK2Stress;

    tensor temp01;

    tensor temp02;

    tensor F;



    nodesDisp = this->getNodesDisp( );



    for( GP_c_r = 0 ; GP_c_r < NumIntegrationPts ; GP_c_r++ ) {

	  r = pts[GP_c_r ];

      rw = wts[GP_c_r ];

      for( GP_c_s = 0 ; GP_c_s < NumIntegrationPts ; GP_c_s++ ) {

        s = pts[GP_c_s ];

        sw = wts[GP_c_s ];

        for( GP_c_t = 0 ; GP_c_t < NumIntegrationPts ; GP_c_t++ ) {

          t = pts[GP_c_t ];

          tw = wts[GP_c_t ];

          where =(GP_c_r * NumIntegrationPts + GP_c_s) * NumIntegrationPts + GP_c_t;

          //dh = shapeFunctionDerivative(r,s,t);

          Jacobian = this->Jacobian_3D(r,s,t);

          det_of_Jacobian  = Jacobian.determinant();

          dhGlobal = this->dh_Global(r,s,t);

          weight = rw * sw * tw * det_of_Jacobian;

          //PK2Stress = theMaterial[where]->getStressTensor();

          //F = theMaterial[where]->getF();

          //temp01 = PK2Stress("ik") * F("jk");

          //  temp01.null_indices();

          temp01 = theMaterial[where]->getPK1StressTensor();

          temp02 = dhGlobal("PJ") * temp01("iJ") * weight;

            temp02.null_indices(); 

          Rr += temp02;

        }

      }

    }



    return Rr;

  }



//======================================================================

tensor TotalLagrangianFD8NodeBrick::getBodyForce(void)

  {

    int B_dim[] = {NumNodes,NumDof};

    tensor Bb(2,B_dim,0.0);



    double r  = 0.0;

    double rw = 0.0;

    double s  = 0.0;

    double sw = 0.0;

    double t  = 0.0;

    double tw = 0.0;



    int where = 0;

    int GP_c_r, GP_c_s, GP_c_t;

    double weight = 0.0;



    int h_dim[] = {20};

    tensor h(1, h_dim, 0.0);

    int dh_dim[] = {NumNodes,NumDof};

    tensor dh(2, dh_dim, 0.0);

    int bodyforce_dim[] = {3};

    tensor bodyforce(1, bodyforce_dim, 0.0);



    double det_of_Jacobian = 0.0;



    tensor Jacobian;

    tensor JacobianINV;



    bodyforce.val(1) = bf(0);

    bodyforce.val(2) = bf(1);

    bodyforce.val(3) = bf(2);



    for( GP_c_r = 0 ; GP_c_r < NumIntegrationPts ; GP_c_r++ ) {

      r = pts[GP_c_r ];

      rw = wts[GP_c_r ];

      for( GP_c_s = 0 ; GP_c_s < NumIntegrationPts ; GP_c_s++ ) {

        s = pts[GP_c_s ];

        sw = wts[GP_c_s ];

        for( GP_c_t = 0 ; GP_c_t < NumIntegrationPts ; GP_c_t++ ) {

          t = pts[GP_c_t ];

          tw = wts[GP_c_t ];

          where =(GP_c_r * NumIntegrationPts + GP_c_s) * NumIntegrationPts + GP_c_t;

          h = shapeFunction(r,s,t);

          dh = shapeFunctionDerivative(r,s,t);

          Jacobian = this->Jacobian_3D(r,s,t);

          det_of_Jacobian  = Jacobian.determinant();

          weight = rw * sw * tw * det_of_Jacobian;

          Bb = Bb +  h("P") * bodyforce("i") * rho *weight;

             Bb.null_indices();

        }

      }

    }

    return Bb;

  }



//======================================================================

tensor TotalLagrangianFD8NodeBrick::getSurfaceForce(void)

  {

    int S_dim[] = {NumNodes, NumDof};

    tensor Ss(2,S_dim,0.0);

    // Need Work Here!



    return Ss;

  }



//============================================================================

tensor TotalLagrangianFD8NodeBrick::getForces(void)

  {

    int F_dim[] = {NumNodes,NumDof};

    tensor Ff(2,F_dim,0.0);



    Ff = this->getBodyForce( ) + this->getSurfaceForce( );



    return Ff;

  }



//=============================================================================

const Matrix &TotalLagrangianFD8NodeBrick::getTangentStiff ()

{

     K.Zero();



     tensor stifftensor = this->getStiffnessTensor();



     int kki=0;

     int kkj=0;

     

     int i, j, k, l;

     for (i=1 ; i<=NumNodes ; i++ ) {

        for (j=1 ; j<=NumNodes ; j++ ) {

           for (k=1 ; k<=NumDof ; k++ ) {

              for (l=1 ; l<=NumDof ; l++ ) {

                 kki = k + NumDof*(i-1);

                 kkj = l + NumDof*(j-1);

                 K(kki-1 , kkj-1) = stifftensor.cval(i,k,l,j); 

              }

           }

        }

     }



     return K;

}



//=============================================================================

const Matrix &TotalLagrangianFD8NodeBrick::getInitialStiff ()

{

     if (Ki != 0) return *Ki;



     K.Zero();

     K = this->getTangentStiff ();



     Ki = new Matrix(K);



     return K;

}



//=============================================================================

const Matrix &TotalLagrangianFD8NodeBrick::getMass ()

{

    // Need Work Here

    M.Zero();

    return M;

}



//======================================================================

tensor TotalLagrangianFD8NodeBrick::getNodesCrds(void) 

{

    const int dimensions[] = {NumNodes, NumDof};

    tensor N_coord(2, dimensions, 0.0);



    int i, j;

    for (i=0; i<NumNodes; i++) {

	  const Vector &TNodesCrds = theNodes[i]->getCrds();

      for (j=0; j<NumDof; j++) {

        N_coord.val(i+1, j+1) = TNodesCrds(j);

	  }		    

    }

    

    return N_coord;

}



//=============================================================================================

tensor TotalLagrangianFD8NodeBrick::getNodesDisp(void)

  {

    int i, j;

    int dimU[] = {NumNodes, NumDof};

    tensor total_disp(2, dimU, 0.0);



    for (i=0; i<NumNodes; i++) {

      const Vector &TNodesDisp = theNodes[i]->getTrialDisp();

      for (j=0; j<NumDof; j++) {

        total_disp.val(i+1, j+1) = TNodesDisp(j);

      }

    }



    return total_disp;

  }



//=============================================================================

void TotalLagrangianFD8NodeBrick::zeroLoad(void)

{



	if ( Q != 0 )

	  Q->Zero();

    

	return;

}





//=============================================================================

int

TotalLagrangianFD8NodeBrick::addLoad(ElementalLoad *theLoad, double loadFactor)

{

    opserr<<"TotalLagrangianFD8NodeBrick::addLoad - load type unknown for ele with tag: "<<this->getTag();          

    return -1;

}



//=============================================================================

int TotalLagrangianFD8NodeBrick::addInertiaLoadToUnbalance(const Vector &accel)

{

    // Check for a quick return

    if (rho == 0.0) return 0;



    static Vector ra(NumElemDof);

    int i, j;



    for (i=0; i<NumNodes; i++) {

      const Vector &RA = theNodes[i]->getRV(accel);

      if ( RA.Size() != NumDof ) {

        opserr << "TotalLagrangianFD8NodeBrick::addInertiaLoadToUnbalance(): matrix and vector sizes are incompatable \n";

        return (-1);

      }

      

      for (j=0; j<NumDof; j++) {

	    ra(i*NumDof +j) = RA(j);

      }



    }



    this->getMass();



    if (Q == 0)  

      Q = new Vector(NumElemDof);



    Q->addMatrixVector(1.0, M, ra, -1.0);



    return 0;  

    

}





//=============================================================================

const Vector &TotalLagrangianFD8NodeBrick::getResistingForce ()

{   

	int i, j;

    int f_dim[] = {NumNodes, NumDof};

    tensor NodalForces_in(2, f_dim, 0.0);

    NodalForces_in = this->getRtensor() - this->getForces();

    

    for (i=0; i<NumNodes; i++) {

      for (j=0; j<NumDof; j++) {

         P(i*NumDof +j) = NodalForces_in.cval(i+1, j+1);

      }

    }

        

    if ( Q != 0 )

      P.addVector(1.0, *Q, -1.0);

    

    return P;

}



//=============================================================================

const Vector &TotalLagrangianFD8NodeBrick::getResistingForceIncInertia ()

{

    int i, j;

    Vector a(NumElemDof);

    

    this->getResistingForce();



    if (rho != 0.0)

    {

      for (i=0; i<NumNodes; i++) {

        const Vector &acc = theNodes[i]->getTrialAccel();

        if ( acc.Size() != NumDof ) {

          opserr << "TotalLagrangianFD8NodeBrick::getResistingForceIncInertia matrix and vector sizes are incompatable \n";

          exit(-1);

        }

      for (j=0; j<NumDof; j++) {

        a(i*NumDof +j) = acc(j);

      }

    }



    this->getMass();

    P.addMatrixVector(1.0, M, a, 1.0);



  }



  return P;

}



//=============================================================================

int TotalLagrangianFD8NodeBrick::sendSelf (int commitTag, Channel &theChannel)

{

     // Not implemtented yet

     return 0;

}



//=============================================================================

int TotalLagrangianFD8NodeBrick::recvSelf (int commitTag, Channel &theChannel,

FEM_ObjectBroker &theBroker)

{

     // Not implemtented yet

     return 0;

}





//=============================================================================

int TotalLagrangianFD8NodeBrick::displaySelf (Renderer &theViewer, int displayMode, float fact)

{

     // Not implemtented yet

     return 0;

}



//=============================================================================

void TotalLagrangianFD8NodeBrick::Print(OPS_Stream &s, int flag)

{

    s << "\nTotalLagrangianFD8NodeBrick, element id:  " << this->getTag() << endln;

    s << "\nConnected external nodes:  " << connectedExternalNodes;

    s << "\nBody forces:  " << bf(0) << " " << bf(1) << " " << bf(2) << endln;



    theMaterial[0]->Print(s,flag);



    tensor sigma;

    Vector P00(6);

    

    int i;

    for (i=0; i<NumTotalGaussPts; i++)

    {

      sigma = theMaterial[i]->getCauchyStressTensor();

      P00(0) = sigma.val(1,1);

      P00(1) = sigma.val(2,2);

      P00(2) = sigma.val(3,3);

      P00(3) = sigma.val(2,3);

      P00(4) = sigma.val(3,1);

      P00(5) = sigma.val(1,2);



      s << "\n where = " << i << endln;

      s << " Stress (Cauchy): xx yy zz yz zx xy) " << P00 << endln;

    }



}



//=============================================================================

Response * TotalLagrangianFD8NodeBrick::setResponse (const char **argv, int argc, Information &eleInfo, OPS_Stream &output)

{

  Response *theResponse = 0;



  char outputData[32];



  output.tag("ElementOutput");

  output.attr("eleType","TotalLagrangianFD8NodeBrick");

  output.attr("eleTag",this->getTag());

  for (int i=1; i<=NumNodes; i++) {

    sprintf(outputData,"node%d",i);

    output.attr(outputData,connectedExternalNodes[i-1]);

  }



  if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0) {

    for (int i=1; i<=NumNodes; i++)

      for (int j=1; j<=NumDof; j++) {

	sprintf(outputData,"P%d_%d",j,i);

	output.tag("ResponseType",outputData);

      }

    

    theResponse = new ElementResponse(this, 1, this->getResistingForce());

    

  } else if (strcmp(argv[0],"CauchyStress") == 0 || strcmp(argv[0],"stress") == 0)

      theResponse = new ElementResponse(this, 3, Vector(NumTotalGaussPts*6));

    

    else if (strcmp(argv[0],"PK2Stress") == 0 || strcmp(argv[0],"PK2stress") == 0)

      theResponse = new ElementResponse(this, 4, Vector(NumTotalGaussPts*6));



    else if (strcmp(argv[0],"EulerianStrain") == 0 || strcmp(argv[0],"strain") == 0)

      theResponse = new ElementResponse(this, 5, Vector(NumTotalGaussPts*6));

    

    else if (strcmp(argv[0],"LagrangianStrain") == 0 || strcmp(argv[0],"iniStrain") == 0)

      theResponse = new ElementResponse(this, 6, Vector(NumTotalGaussPts*6));



    output.endTag(); // ElementOutput



  return theResponse; 

}



//=============================================================================

int TotalLagrangianFD8NodeBrick::getResponse (int responseID, Information &eleInfo)

{    

	 int i;

     static Vector P0(NumTotalGaussPts*6);

     

     switch (responseID) {

     

     case 1:

          return eleInfo.setVector(this->getResistingForce() );



     case 3: { 

        Vector P0(NumTotalGaussPts*6);

        tensor sigma; 

        for (i=0; i<NumTotalGaussPts; i++) {

          sigma = theMaterial[i]->getCauchyStressTensor();

          P0(i*6 +0 ) = sigma.val(1,1);

          P0(i*6 +1 ) = sigma.val(2,2);

          P0(i*6 +2 ) = sigma.val(3,3);

          P0(i*6 +3 ) = sigma.val(2,3);

          P0(i*6 +4 ) = sigma.val(3,1);

          P0(i*6 +5 ) = sigma.val(1,2);

        }

        return eleInfo.setVector(P0);

     }



     case 4: { 

        Vector P0(NumTotalGaussPts*6);

        tensor sigma; 

        for (i=0; i<NumTotalGaussPts; i++) {

          sigma = theMaterial[i]->getStressTensor();

          P0(i*6 +0 ) = sigma.val(1,1);

          P0(i*6 +1 ) = sigma.val(2,2);

          P0(i*6 +2 ) = sigma.val(3,3);

          P0(i*6 +3 ) = sigma.val(2,3);

          P0(i*6 +4 ) = sigma.val(3,1);

          P0(i*6 +5 ) = sigma.val(1,2);

        }

        return eleInfo.setVector(P0);

     }



     case 5: { 

        Vector P0(NumTotalGaussPts*6);

        tensor e;

	tensor E;

	tensor F;

	tensor tI2("I", 2, def_dim_2); 

        for (i=0; i<NumTotalGaussPts; i++) {

          E = theMaterial[i]->getStrainTensor();

	  F = theMaterial[i]->getF();

	  F = F.inverse();

	  e = F("ki")*F("kj"); e.null_indices();

	  e = (tI2-e) *0.5;

          P0(i*6 +0 ) = e.val(1,1);

          P0(i*6 +1 ) = e.val(2,2);

          P0(i*6 +2 ) = e.val(3,3);

          P0(i*6 +3 ) = e.val(2,3);

          P0(i*6 +4 ) = e.val(3,1);

          P0(i*6 +5 ) = e.val(1,2);

        }

        return eleInfo.setVector(P0);

     }



     case 6: { 

        Vector P0(NumTotalGaussPts*6);

        tensor E; 

        for (i=0; i<NumTotalGaussPts; i++) {

          E = theMaterial[i]->getStrainTensor();

          P0(i*6 +0 ) = E.val(1,1);

          P0(i*6 +1 ) = E.val(2,2);

          P0(i*6 +2 ) = E.val(3,3);

          P0(i*6 +3 ) = E.val(2,3);

          P0(i*6 +4 ) = E.val(3,1);

          P0(i*6 +5 ) = E.val(1,2);

        }

        return eleInfo.setVector(P0);

     }

    

     default:

     return -1;



   }    

}





//#############################################################################

//===================================================================

tensor TotalLagrangianFD8NodeBrick::shapeFunction(double r1, double r2, double r3)

  {



    int dimension[] = {NumNodes};

    tensor h(1, dimension, 0.0);



      // influence of the node number 8

    h.val(8)=(1.0+r1)*(1.0-r2)*(1.0-r3)*0.125;

      // influence of the node number 7

    h.val(7)=(1.0-r1)*(1.0-r2)*(1.0-r3)*0.125;

      // influence of the node number 6

    h.val(6)=(1.0-r1)*(1.0+r2)*(1.0-r3)*0.125;

      // influence of the node number 5

    h.val(5)=(1.0+r1)*(1.0+r2)*(1.0-r3)*0.125;



      // influence of the node number 4

    h.val(4)=(1.0+r1)*(1.0-r2)*(1.0+r3)*0.125;

      // influence of the node number 3

    h.val(3)=(1.0-r1)*(1.0-r2)*(1.0+r3)*0.125;

      // influence of the node number 2

    h.val(2)=(1.0-r1)*(1.0+r2)*(1.0+r3)*0.125;

      // influence of the node number 1

    h.val(1)=(1.0+r1)*(1.0+r2)*(1.0+r3)*0.125;



    return h;

  }





//==============================================================

tensor TotalLagrangianFD8NodeBrick::shapeFunctionDerivative(double r1, double r2, double r3)

  {



    int dimensions[] = {NumNodes, NumDof};

    tensor dh(2, dimensions, 0.0);



      // influence of the node number 8

    dh.val(8,1)= (1.0-r2)*(1.0-r3)*0.125;

    dh.val(8,2)=-(1.0+r1)*(1.0-r3)*0.125;

    dh.val(8,3)=-(1.0+r1)*(1.0-r2)*0.125;

      // influence of the node number 7

    dh.val(7,1)=-(1.0-r2)*(1.0-r3)*0.125;

    dh.val(7,2)=-(1.0-r1)*(1.0-r3)*0.125;

    dh.val(7,3)=-(1.0-r1)*(1.0-r2)*0.125;

      // influence of the node number 6

    dh.val(6,1)=-(1.0+r2)*(1.0-r3)*0.125;

    dh.val(6,2)= (1.0-r1)*(1.0-r3)*0.125;

    dh.val(6,3)=-(1.0-r1)*(1.0+r2)*0.125;

      // influence of the node number 5

    dh.val(5,1)= (1.0+r2)*(1.0-r3)*0.125;

    dh.val(5,2)= (1.0+r1)*(1.0-r3)*0.125;

    dh.val(5,3)=-(1.0+r1)*(1.0+r2)*0.125;



      // influence of the node number 4

    dh.val(4,1)= (1.0-r2)*(1.0+r3)*0.125;

    dh.val(4,2)=-(1.0+r1)*(1.0+r3)*0.125;

    dh.val(4,3)= (1.0+r1)*(1.0-r2)*0.125;

      // influence of the node number 3

    dh.val(3,1)=-(1.0-r2)*(1.0+r3)*0.125;

    dh.val(3,2)=-(1.0-r1)*(1.0+r3)*0.125;

    dh.val(3,3)= (1.0-r1)*(1.0-r2)*0.125;

      // influence of the node number 2

    dh.val(2,1)=-(1.0+r2)*(1.0+r3)*0.125;

    dh.val(2,2)= (1.0-r1)*(1.0+r3)*0.125;

    dh.val(2,3)= (1.0-r1)*(1.0+r2)*0.125;

      // influence of the node number 1

    dh.val(1,1)= (1.0+r2)*(1.0+r3)*0.125;

    dh.val(1,2)= (1.0+r1)*(1.0+r3)*0.125;

    dh.val(1,3)= (1.0+r1)*(1.0+r2)*0.125;



    return dh;

  }





#endif



