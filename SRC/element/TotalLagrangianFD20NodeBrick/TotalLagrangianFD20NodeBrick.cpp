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
//# DATE:              Sept2003
//# UPDATE HISTORY:    28May2004, Zhao put all Ks & Rs in the integration cycle and
//#                          and use minor symmetries to make concise and efficient
//#                    April2005 Zhao adds new output options
//#
//===============================================================================
#ifndef TOTALLAGRANGIANFD20NODEBRICK_CPP
#define TOTALLAGRANGIANFD20NODEBRICK_CPP

#include <TotalLagrangianFD20NodeBrick.h>

#define NumIntegrationPts 3
#define NumTotalGaussPts 27
#define NumNodes 20
#define NumDof 3
#define NumElemDof 60

double TotalLagrangianFD20NodeBrick::matrixData[NumElemDof*NumElemDof];
Matrix TotalLagrangianFD20NodeBrick::K(matrixData, NumElemDof, NumElemDof);
Matrix TotalLagrangianFD20NodeBrick::M(matrixData, NumElemDof, NumElemDof);
Vector TotalLagrangianFD20NodeBrick::P(NumElemDof);
double TotalLagrangianFD20NodeBrick::pts[NumDof] = {-0.774596669241483, 0.0, +0.774596669241483};
double TotalLagrangianFD20NodeBrick::wts[NumDof] = {+0.55555555555555556, +0.88888888888888889,+0.55555555555555556};

//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

TotalLagrangianFD20NodeBrick::TotalLagrangianFD20NodeBrick(int tag,
int node_numb_1,  int node_numb_2,  int node_numb_3,  int node_numb_4,
int node_numb_5,  int node_numb_6,  int node_numb_7,  int node_numb_8,
int node_numb_9,  int node_numb_10, int node_numb_11, int node_numb_12,
int node_numb_13, int node_numb_14, int node_numb_15, int node_numb_16,
int node_numb_17, int node_numb_18, int node_numb_19, int node_numb_20,
NDMaterial &m, double b1, double b2, double b3)
:Element(tag, ELE_TAG_TotalLagrangianFD20NodeBrick ),
 theMaterial(0), connectedExternalNodes(NumNodes), Q(NumElemDof), bf(NumDof), Ki(0)
{
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

      bf(0) = b1;
      bf(1) = b2;
      bf(2) = b3;

      theMaterial = new NDMaterial *[NumTotalGaussPts];

      if (theMaterial == 0)
      {
       opserr<<"FiniteDeformationElastic3D::FiniteDeformationElastic3D -- failed allocate material model pointer\n";
       exit(-1);
      }

      for (int i=0; i<NumTotalGaussPts; i++)
      {
       theMaterial[i] = m.getCopy();
       if (theMaterial[i] == 0)
       {
        opserr<<"FiniteDeformationElastic3D::FiniteDeformationElastic3D -- failed allocate material model pointer\n";
        exit(-1);
       }
      }

      rho = m.getRho();

      for (int j=0; j<NumNodes; j++) theNodes[j] = 0;

}

//-------------------------------------------------------------------------------------------
TotalLagrangianFD20NodeBrick::TotalLagrangianFD20NodeBrick ()
:Element(0, ELE_TAG_TotalLagrangianFD20NodeBrick ),
 theMaterial(0), connectedExternalNodes(NumNodes), Q(NumElemDof), bf(NumDof), Ki(0)
{
     for (int i=0; i<NumNodes; i++)  theNodes[i] = 0;

     bf(0) = 0.0;
     bf(1) = 0.0;
     bf(2) = 0.0;

     rho = 0.0;
}

//-------------------------------------------------------------------------------------------------
TotalLagrangianFD20NodeBrick::~TotalLagrangianFD20NodeBrick ()
{
    for (int i=0; i<NumTotalGaussPts; i++)
    {
     if (theMaterial[i]) delete theMaterial[i];
    }

    if(theMaterial) delete [] theMaterial;

    if(Ki) delete Ki;

}
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


//=============================================================================
int TotalLagrangianFD20NodeBrick::getNumExternalNodes () const
{
    return NumNodes;
}

//=============================================================================
const ID& TotalLagrangianFD20NodeBrick::getExternalNodes ()
{
    return connectedExternalNodes;
}

//=============================================================================
Node **TotalLagrangianFD20NodeBrick::getNodePtrs(void)
{
  return theNodes;
}

//=============================================================================
int TotalLagrangianFD20NodeBrick::getNumDOF ()
{
    return NumElemDof;
}

//=============================================================================
void TotalLagrangianFD20NodeBrick::setDomain (Domain *theDomain)
{
    // Check Domain is not null - invoked when object removed from a domain
    if (theDomain == 0)
    {
        theNodes[0] = 0;
        theNodes[1] = 0;
        theNodes[2] = 0;
        theNodes[3] = 0;
        theNodes[4] = 0;
        theNodes[5] = 0;
        theNodes[6] = 0;
        theNodes[7] = 0;
        theNodes[8] = 0;
        theNodes[9] = 0;
        theNodes[10] = 0;
        theNodes[11] = 0;
        theNodes[12] = 0;
        theNodes[13] = 0;
        theNodes[14] = 0;
        theNodes[15] = 0;
        theNodes[16] = 0;
        theNodes[17] = 0;
        theNodes[18] = 0;
        theNodes[19] = 0;
  return;
    }


      int Nd1 = connectedExternalNodes(0);
      int Nd2 = connectedExternalNodes(1);
      int Nd3 = connectedExternalNodes(2);
      int Nd4 = connectedExternalNodes(3);
      int Nd5 = connectedExternalNodes(4);
      int Nd6 = connectedExternalNodes(5);
      int Nd7 = connectedExternalNodes(6);
      int Nd8 = connectedExternalNodes(7);
      int Nd9  = connectedExternalNodes( 8);
      int Nd10 = connectedExternalNodes( 9);
      int Nd11 = connectedExternalNodes(10);
      int Nd12 = connectedExternalNodes(11);
      int Nd13 = connectedExternalNodes(12);
      int Nd14 = connectedExternalNodes(13);
      int Nd15 = connectedExternalNodes(14);
      int Nd16 = connectedExternalNodes(15);
      int Nd17 = connectedExternalNodes(16);
      int Nd18 = connectedExternalNodes(17);
      int Nd19 = connectedExternalNodes(18);
      int Nd20 = connectedExternalNodes(19);


      theNodes[0] = theDomain->getNode(Nd1);
      theNodes[1] = theDomain->getNode(Nd2);
      theNodes[2] = theDomain->getNode(Nd3);
      theNodes[3] = theDomain->getNode(Nd4);
      theNodes[4] = theDomain->getNode(Nd5);
      theNodes[5] = theDomain->getNode(Nd6);
      theNodes[6] = theDomain->getNode(Nd7);
      theNodes[7] = theDomain->getNode(Nd8);
      theNodes[8] = theDomain->getNode(Nd9);
      theNodes[9] = theDomain->getNode(Nd10);
      theNodes[10] = theDomain->getNode(Nd11);
      theNodes[11] = theDomain->getNode(Nd12);
      theNodes[12] = theDomain->getNode(Nd13);
      theNodes[13] = theDomain->getNode(Nd14);
      theNodes[14] = theDomain->getNode(Nd15);
      theNodes[15] = theDomain->getNode(Nd16);
      theNodes[16] = theDomain->getNode(Nd17);
      theNodes[17] = theDomain->getNode(Nd18);
      theNodes[18] = theDomain->getNode(Nd19);
      theNodes[19] = theDomain->getNode(Nd20);

      if (theNodes[0]==0  || theNodes[1]==0  || theNodes[2]==0  || theNodes[3]==0  ||
          theNodes[4]==0  || theNodes[5]==0  || theNodes[6]==0  || theNodes[7]==0  ||
          theNodes[8]==0  || theNodes[9]==0  || theNodes[10]==0 || theNodes[11]==0 ||
          theNodes[12]==0 || theNodes[13]==0 || theNodes[14]==0 || theNodes[15]==0 ||
          theNodes[16]==0 || theNodes[17]==0 || theNodes[18]==0 || theNodes[19]==0 )
      {
        opserr << "FATAL ERROR TotalLagrangianFD20NodeBrick (tag: " << this->getTag() <<
        " ), node not found in domain\n";
        exit(-1);
      }

      this->DomainComponent::setDomain(theDomain);  // Very Important!!

      int dofNd1 = theNodes[0]->getNumberDOF();
      int dofNd2 = theNodes[1]->getNumberDOF();
      int dofNd3 = theNodes[2]->getNumberDOF();
      int dofNd4 = theNodes[3]->getNumberDOF();
      int dofNd5 = theNodes[4]->getNumberDOF();
      int dofNd6 = theNodes[5]->getNumberDOF();
      int dofNd7 = theNodes[6]->getNumberDOF();
      int dofNd8 = theNodes[7]->getNumberDOF();
      int dofNd9 = theNodes[8]->getNumberDOF();
      int dofNd10 = theNodes[9]->getNumberDOF();
      int dofNd11 = theNodes[10]->getNumberDOF();
      int dofNd12 = theNodes[11]->getNumberDOF();
      int dofNd13 = theNodes[12]->getNumberDOF();
      int dofNd14 = theNodes[13]->getNumberDOF();
      int dofNd15 = theNodes[14]->getNumberDOF();
      int dofNd16 = theNodes[15]->getNumberDOF();
      int dofNd17 = theNodes[16]->getNumberDOF();
      int dofNd18 = theNodes[17]->getNumberDOF();
      int dofNd19 = theNodes[18]->getNumberDOF();
      int dofNd20 = theNodes[19]->getNumberDOF();

      if (dofNd1  != NumDof || dofNd2  != NumDof || dofNd3  != NumDof || dofNd4  != NumDof ||
          dofNd5  != NumDof || dofNd6  != NumDof || dofNd7  != NumDof || dofNd8  != NumDof ||
          dofNd9  != NumDof || dofNd10 != NumDof || dofNd11 != NumDof || dofNd12 != NumDof ||
          dofNd13 != NumDof || dofNd14 != NumDof || dofNd15 != NumDof || dofNd16 != NumDof ||
          dofNd17 != NumDof || dofNd18 != NumDof || dofNd19 != NumDof || dofNd20 != NumDof )
      {
        opserr << "FATAL ERROR TotalLagrangianFD20NodeBrick (tag: " << this->getTag() <<
        "), has differing number of DOFs at its nodes\n";
        exit(-1);
      }


}

//=============================================================================
int TotalLagrangianFD20NodeBrick::commitState ()
{
  int retVal = 0;

  if ((retVal = this->Element::commitState()) != 0)
  {
    opserr << "TotalLagrangianFD20NodeBrick::commitState () - failed in base class";
  }


  for (int i = 0; i < NumTotalGaussPts ; i++)
         retVal += theMaterial[i]->commitState();

  return retVal;
}

//=============================================================================
int TotalLagrangianFD20NodeBrick::revertToLastCommit ()
{
    int retVal = 0;

    for (int i = 0; i < NumTotalGaussPts; i++)
       retVal += theMaterial[i]->revertToLastCommit();

    return retVal;
}

//=============================================================================
int TotalLagrangianFD20NodeBrick::revertToStart ()
{
    int retVal = 0;

    for (int i = 0; i < NumTotalGaussPts; i++)
       retVal += theMaterial[i]->revertToStart();

    return retVal;

}

//=============================================================================
int TotalLagrangianFD20NodeBrick::update ()
{
    int ret = 0;
    tensor dh;
    tensor dH_dX;
    int where;
    double r = 0.0;
    double s = 0.0;
    double t = 0.0;
    tensor I_ij("I", 2, def_dim_2);
    tensor currentF;
    tensor updatedF;

    tensor InitialNodesCrds = this->getNodesCrds();
    tensor CurrentNodesDisp = this->getNodesDisp();

    for( short GP_c_r = 0 ; GP_c_r < NumIntegrationPts ; GP_c_r++ )
      {
        r = pts[GP_c_r ];
        for( short GP_c_s = 0 ; GP_c_s < NumIntegrationPts ; GP_c_s++ )
          {
            s = pts[GP_c_s ];
            for( short GP_c_t = 0 ; GP_c_t < NumIntegrationPts ; GP_c_t++ )
              {
                t = pts[GP_c_t ];
                        where =(GP_c_r * NumIntegrationPts + GP_c_s) * NumIntegrationPts + GP_c_t;
                        dh = shapeFunctionDerivative(r,s,t);
                        dH_dX = this->dh_Global(dh);
                        currentF = CurrentNodesDisp("Ia") * dH_dX("Ib");
                          currentF.null_indices(); //CurrentNodesDisp.null_indices(); dH_dX.null_indices();
                        updatedF = currentF + I_ij;
                        ret += theMaterial[where]->setTrialF(updatedF);
        }
    }
       }
    return ret;
}

//======================================================================
tensor TotalLagrangianFD20NodeBrick::Jacobian_3D(tensor dh)
  {
     tensor N_C = this->getNodesCrds();
     tensor J3D = dh("ij") * N_C("ik");
       J3D.null_indices(); //dh.null_indices(); N_C.null_indices();
     return J3D;
  }

//======================================================================
tensor TotalLagrangianFD20NodeBrick::Jacobian_3Dinv(tensor dh)
  {
     tensor N_C = this->getNodesCrds();
     tensor J3D = dh("ij") * N_C("ik");
       J3D.null_indices(); //dh.null_indices(); N_C.null_indices();
     tensor J3Dinv = J3D.inverse();
     return J3Dinv;
  }

//======================================================================
tensor TotalLagrangianFD20NodeBrick::dh_Global(tensor dh)
  {
     tensor  JacobianINV0 = this->Jacobian_3Dinv(dh);
     tensor  dhGlobal_0 = dh("ij") * JacobianINV0("kj");
       dhGlobal_0.null_indices(); //dh.null_indices(); JacobianINV0.null_indices();
     return dhGlobal_0;
  }

//======================================================================
tensor TotalLagrangianFD20NodeBrick::getStiffnessTensor(void)
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

    short where = 0;
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

    for( short GP_c_r = 0 ; GP_c_r < NumIntegrationPts ; GP_c_r++ )
      {
        r = pts[GP_c_r ];
        rw = wts[GP_c_r ];
        for( short GP_c_s = 0 ; GP_c_s < NumIntegrationPts ; GP_c_s++ )
          {
            s = pts[GP_c_s ];
            sw = wts[GP_c_s ];
            for( short GP_c_t = 0 ; GP_c_t < NumIntegrationPts ; GP_c_t++ )
              {
                t = pts[GP_c_t ];
                tw = wts[GP_c_t ];
                        where =(GP_c_r * NumIntegrationPts + GP_c_s) * NumIntegrationPts + GP_c_t;
                        dh = shapeFunctionDerivative(r,s,t);
                        Jacobian = this->Jacobian_3D(dh);
                        det_of_Jacobian  = Jacobian.determinant();
                        dhGlobal = this->dh_Global(dh);
                        weight = rw * sw * tw * det_of_Jacobian;
                        PK2Stress = theMaterial[where]->getStressTensor();
                        L2 = theMaterial[where]->getTangentTensor();
                        F = theMaterial[where]->getF();
                        
           //K1
                        temp04 = dhGlobal("Pb") * tI2("mn");
                          temp04.null_indices(); //dhGlobal.null_indices(); tI2.null_indices();
                        temp02 = PK2Stress("bd") * dhGlobal("Qd");   
                          temp02.null_indices(); //dhGlobal.null_indices(); PK2Stress.null_indices();
                        temp06 = temp04("Pbmn") * temp02("bQ") * weight;
                          temp06.null_indices(); //temp04.null_indices(); temp02.null_indices();
                        Kk += temp06;
                        
           //K2
                        temp03 =  dhGlobal("Pb") * F("ma");
                          temp03.null_indices(); //dhGlobal.null_indices(); F.null_indices();
                        temp04 = F("nc") * L2("abcd");
                          temp04.null_indices(); //F.null_indices(); L2.null_indices();
                        temp05 = temp04("nabd") * dhGlobal("Qd"); 
                          temp05.null_indices(); //temp04.null_indices(); dhGlobal.null_indices();
                        temp06 = temp03("Pbma") * temp05("nabQ") * weight;
                          temp06.null_indices(); //temp03.null_indices(); temp05.null_indices();
                        Kk += temp06;
                                                                    
                        
///*
//			// K1
//                        temp01 = dhGlobal("Pj") * L2("ijkl");
//                        temp01.null_indices(); dhGlobal.null_indices(); L2.null_indices();
//                        temp05 = temp01("Pikl") * dhGlobal("Qk") * weight;
//                        temp05.null_indices(); dhGlobal.null_indices(); temp01.null_indices();
//                        Kk += temp05;
//
//			// K2
//                        temp01 = dhGlobal("nk") * nodesDisp("nm");
//                        temp01.null_indices(); dhGlobal.null_indices(); nodesDisp.null_indices();
//                        temp02 = temp01("km") * dhGlobal("Qm");
//                        temp02.null_indices(); dhGlobal.null_indices(); temp01.null_indices();
//                        temp05 = dhGlobal("Pj") * L2("ijkl");
//                        temp05.null_indices(); dhGlobal.null_indices(); L2.null_indices();
//                        temp06 = temp05("Pikl") * temp02("kQ") * weight;
//                        temp06.null_indices(); temp05.null_indices(); temp02.null_indices();
//                        Kk += temp06;
//
//
//			// K3
//                        temp01 = dhGlobal("mk") * nodesDisp("mx");
//                        temp01.null_indices( ); dhGlobal.null_indices( ); nodesDisp.null_indices( );
//                        temp02 = temp01("kx") * dhGlobal("Qx");
//                        temp02.null_indices( ); dhGlobal.null_indices( ); temp01.null_indices( );
//                        temp03 = dhGlobal("nk") * nodesDisp("ny");
//                        temp03.null_indices( ); dhGlobal.null_indices( ); nodesDisp.null_indices( );
//                        temp04 = dhGlobal("Py") * temp03("jy");
//                        temp04.null_indices( ); temp03.null_indices( ); dhGlobal.null_indices( );
//                        temp05 = temp04("Pj") * L2("ijkl");
//                        temp05.null_indices( ); temp04.null_indices( ); L2.null_indices( );
//                        temp06 = temp05("Pikl") * temp02("kQ") * weight;
//                        temp06.null_indices( ); temp05.null_indices( ); temp05.null_indices( );
//                        Kk += temp06;
//
//			// K4
//                        temp01 = dhGlobal("mj") * nodesDisp("mx");
//                        temp01.null_indices(); dhGlobal.null_indices(); nodesDisp.null_indices();
//                        temp02 =  dhGlobal("Px") * temp01("jx");
//                        temp02.null_indices(); dhGlobal.null_indices(); temp01.null_indices();
//                        temp05 = temp02("Pj") * L2("ijkl");
//                        temp05.null_indices(); temp02.null_indices(); L2.null_indices();
//                        temp06 = temp05("Pikl") * dhGlobal("Qk") * weight;
//                        temp06.null_indices(); dhGlobal.null_indices(); temp05.null_indices();
//                        Kk += temp06;
//
//			// K5
//                        temp01 =  dhGlobal("Qj") * PK2Stress("ij");
//                        temp01.null_indices(); dhGlobal.null_indices(); PK2Stress.null_indices();
//                        temp05 = dhGlobal("Pl") * temp01("Qi") * weight;
//                        temp05.null_indices( ); dhGlobal.null_indices(); temp01.null_indices();
//                        for ( int i=1 ; i<=NumNodes ; i++ ) {
//                           for ( int j=1 ; j<=NumNodes ; j++ ) {
//                              for ( int k=1 ; k<=NumDof ; k++ ) {
//                                 for ( int l=1 ; l<=NumDof ; l++ ) {
//                                    temp06.val(i,k,l,j) = temp05.cval(i,l,j,k); 
//                                 }
//                              }
//                           }
//                        }
//			Kk += temp06;
//*/			            

              }
          }
      }


    return Kk;
  }

//======================================================================
tensor TotalLagrangianFD20NodeBrick::getRtensor(void)
  {
    int R_dim[] = {NumNodes,NumDof};
    tensor Rr(2,R_dim,0.0);

    double r  = 0.0;
    double rw = 0.0;
    double s  = 0.0;
    double sw = 0.0;
    double t  = 0.0;
    double tw = 0.0;

    short where = 0;
    double weight = 0.0;

    int dh_dim[] = {NumNodes,NumDof};
    tensor dh(2, dh_dim, 0.0);

    double det_of_Jacobian = 0.0;

    tensor Jacobian;
    tensor JacobianINV;
    tensor dhGlobal;
    tensor currentF;
    tensor nodesDisp;
    stresstensor PK2Stress;
    tensor temp01;
    tensor temp02;
    //tensor temp03;
    tensor F;

    nodesDisp = this->getNodesDisp( );

    for( short GP_c_r = 0 ; GP_c_r < NumIntegrationPts ; GP_c_r++ )
      {
        r = pts[GP_c_r ];
        rw = wts[GP_c_r ];
        for( short GP_c_s = 0 ; GP_c_s < NumIntegrationPts ; GP_c_s++ )
          {
            s = pts[GP_c_s ];
            sw = wts[GP_c_s ];
            for( short GP_c_t = 0 ; GP_c_t < NumIntegrationPts ; GP_c_t++ )
              {
                t = pts[GP_c_t ];
                tw = wts[GP_c_t ];
                        where =(GP_c_r * NumIntegrationPts + GP_c_s) * NumIntegrationPts + GP_c_t;
                        dh = shapeFunctionDerivative(r,s,t);
                        Jacobian = this->Jacobian_3D(dh);
                        det_of_Jacobian  = Jacobian.determinant();
                        dhGlobal = this->dh_Global(dh);
                        weight = rw * sw * tw * det_of_Jacobian;
                        PK2Stress = theMaterial[where]->getStressTensor();
                        F = theMaterial[where]->getF();                        
                        
                        temp01 = PK2Stress("ik") * F("jk");
                          temp01.null_indices(); //PK2Stress.null_indices(); F.null_indices();
                        temp02 = dhGlobal("Pj") * temp01("ij") * weight;
                          temp02.null_indices(); //dhGlobal.null_indices(); temp01.null_indices();
                        Rr += temp02; 
                        
///*                                              
//		// R1
//			temp01 = dhGlobal("Pj") * PK2Stress("ij") * weight;
//                        temp01.null_indices(); dhGlobal.null_indices(); PK2Stress.null_indices();
//                        Rr += temp01;			
//			
//		// R2
//                        temp01 = dhGlobal("Bm") * nodesDisp("Bn");
//                        temp01.null_indices(); dhGlobal.null_indices(); nodesDisp.null_indices();
//                        temp02 = dhGlobal("Pn") * temp01("mn");
//                        temp02.null_indices(); dhGlobal.null_indices(); temp01.null_indices();
//                        temp03 = temp02("Pm") * PK2Stress("im") * weight;
//                        temp03.null_indices(); temp02.null_indices(); PK2Stress.null_indices();
//                        Rr += temp03;
//*/                        

              }
          }
      }

    return Rr;
  }

//======================================================================
tensor TotalLagrangianFD20NodeBrick::getBodyForce(void)
  {
    int B_dim[] = {NumNodes,NumDof};
    tensor Bb(2,B_dim,0.0);

    double r  = 0.0;
    double rw = 0.0;
    double s  = 0.0;
    double sw = 0.0;
    double t  = 0.0;
    double tw = 0.0;

    short where = 0;
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

    for( short GP_c_r = 0 ; GP_c_r < NumIntegrationPts ; GP_c_r++ )
      {
        r = pts[GP_c_r ];
        rw = wts[GP_c_r ];
        for( short GP_c_s = 0 ; GP_c_s < NumIntegrationPts ; GP_c_s++ )
          {
            s = pts[GP_c_s ];
            sw = wts[GP_c_s ];
            for( short GP_c_t = 0 ; GP_c_t < NumIntegrationPts ; GP_c_t++ )
              {
                t = pts[GP_c_t ];
                tw = wts[GP_c_t ];
                        where =(GP_c_r * NumIntegrationPts + GP_c_s) * NumIntegrationPts + GP_c_t;
                        h = shapeFunction(r,s,t);
                        dh = shapeFunctionDerivative(r,s,t);
                        Jacobian = this->Jacobian_3D(dh);
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
tensor TotalLagrangianFD20NodeBrick::getSurfaceForce(void)
  {
    int S_dim[] = {NumNodes,NumDof};
    tensor Ss(2,S_dim,0.0);
    // Need Work Here!

    return Ss;
  }

//============================================================================
tensor TotalLagrangianFD20NodeBrick::getForces(void)
  {
    int F_dim[] = {NumNodes,NumDof};
    tensor Ff(2,F_dim,0.0);

    Ff = this->getBodyForce( ) + this->getSurfaceForce( );

    return Ff;
  }

//=============================================================================
const Matrix &TotalLagrangianFD20NodeBrick::getTangentStiff ()
{
     K.Zero();

     tensor stifftensor = this->getStiffnessTensor();

     int kki=0;
     int kkj=0;

     for ( int i=1 ; i<=NumNodes ; i++ ) {
        for ( int j=1 ; j<=NumNodes ; j++ ) {
           for ( int k=1 ; k<=NumDof ; k++ ) {
              for ( int l=1 ; l<=NumDof ; l++ ) {
                 kki = k+NumDof*(i-1);
                 kkj = l+NumDof*(j-1);
                 K( kki-1 , kkj-1 ) = stifftensor.cval(i,k,l,j); 
              }
           }
        }
     }

     return K;
}

//=============================================================================
const Matrix &TotalLagrangianFD20NodeBrick::getInitialStiff ()
{
     if (Ki != 0) return *Ki;

     K.Zero();
     K = this->getTangentStiff ();

     Ki = new Matrix(K);

     return K;

}

//=============================================================================
const Matrix &TotalLagrangianFD20NodeBrick::getMass ()
{
    // Need Work Here
    M.Zero();
    return M;
}

//======================================================================
tensor TotalLagrangianFD20NodeBrick::getNodesCrds(void)
  {
    const int dimensions[] = {NumNodes,NumDof};
    tensor N_coord(2, dimensions, 0.0);

// Using node pointers, which come from the Domain
    const Vector &nd1Crds = theNodes[0]->getCrds();
    const Vector &nd2Crds = theNodes[1]->getCrds();
    const Vector &nd3Crds = theNodes[2]->getCrds();
    const Vector &nd4Crds = theNodes[3]->getCrds();
    const Vector &nd5Crds = theNodes[4]->getCrds();
    const Vector &nd6Crds = theNodes[5]->getCrds();
    const Vector &nd7Crds = theNodes[6]->getCrds();
    const Vector &nd8Crds = theNodes[7]->getCrds();

    const Vector &nd9Crds  = theNodes[8]->getCrds();
    const Vector &nd10Crds = theNodes[9]->getCrds();
    const Vector &nd11Crds = theNodes[10]->getCrds();
    const Vector &nd12Crds = theNodes[11]->getCrds();

    const Vector &nd13Crds = theNodes[12]->getCrds();
    const Vector &nd14Crds = theNodes[13]->getCrds();
    const Vector &nd15Crds = theNodes[14]->getCrds();
    const Vector &nd16Crds = theNodes[15]->getCrds();


    const Vector &nd17Crds = theNodes[16]->getCrds();
    const Vector &nd18Crds = theNodes[17]->getCrds();
    const Vector &nd19Crds = theNodes[18]->getCrds();
    const Vector &nd20Crds = theNodes[19]->getCrds();

    N_coord.val(1,1)=nd1Crds(0); N_coord.val(1,2)=nd1Crds(1); N_coord.val(1,3)=nd1Crds(2);
    N_coord.val(2,1)=nd2Crds(0); N_coord.val(2,2)=nd2Crds(1); N_coord.val(2,3)=nd2Crds(2);
    N_coord.val(3,1)=nd3Crds(0); N_coord.val(3,2)=nd3Crds(1); N_coord.val(3,3)=nd3Crds(2);
    N_coord.val(4,1)=nd4Crds(0); N_coord.val(4,2)=nd4Crds(1); N_coord.val(4,3)=nd4Crds(2);
    N_coord.val(5,1)=nd5Crds(0); N_coord.val(5,2)=nd5Crds(1); N_coord.val(5,3)=nd5Crds(2);
    N_coord.val(6,1)=nd6Crds(0); N_coord.val(6,2)=nd6Crds(1); N_coord.val(6,3)=nd6Crds(2);
    N_coord.val(7,1)=nd7Crds(0); N_coord.val(7,2)=nd7Crds(1); N_coord.val(7,3)=nd7Crds(2);
    N_coord.val(8,1)=nd8Crds(0); N_coord.val(8,2)=nd8Crds(1); N_coord.val(8,3)=nd8Crds(2);

    N_coord.val(9 ,1)=nd9Crds(0);  N_coord.val(9 ,2)=nd9Crds(1);  N_coord.val(9 ,3)=nd9Crds(2);
    N_coord.val(10,1)=nd10Crds(0); N_coord.val(10,2)=nd10Crds(1); N_coord.val(10,3)=nd10Crds(2);
    N_coord.val(11,1)=nd11Crds(0); N_coord.val(11,2)=nd11Crds(1); N_coord.val(11,3)=nd11Crds(2);
    N_coord.val(12,1)=nd12Crds(0); N_coord.val(12,2)=nd12Crds(1); N_coord.val(12,3)=nd12Crds(2);

    N_coord.val(13,1)=nd13Crds(0); N_coord.val(13,2)=nd13Crds(1); N_coord.val(13,3)=nd13Crds(2);
    N_coord.val(14,1)=nd14Crds(0); N_coord.val(14,2)=nd14Crds(1); N_coord.val(14,3)=nd14Crds(2);
    N_coord.val(15,1)=nd15Crds(0); N_coord.val(15,2)=nd15Crds(1); N_coord.val(15,3)=nd15Crds(2);
    N_coord.val(16,1)=nd16Crds(0); N_coord.val(16,2)=nd16Crds(1); N_coord.val(16,3)=nd16Crds(2);

    N_coord.val(17,1)=nd17Crds(0); N_coord.val(17,2)=nd17Crds(1); N_coord.val(17,3)=nd17Crds(2);
    N_coord.val(18,1)=nd18Crds(0); N_coord.val(18,2)=nd18Crds(1); N_coord.val(18,3)=nd18Crds(2);
    N_coord.val(19,1)=nd19Crds(0); N_coord.val(19,2)=nd19Crds(1); N_coord.val(19,3)=nd19Crds(2);
    N_coord.val(20,1)=nd20Crds(0); N_coord.val(20,2)=nd20Crds(1); N_coord.val(20,3)=nd20Crds(2);

    return N_coord;

  }

//=============================================================================================
tensor TotalLagrangianFD20NodeBrick::getNodesDisp(void)
  {
    const int dimensions[] = {NumNodes,NumDof};
    tensor total_disp(2, dimensions, 0.0);

    const Vector &TotDis1 = theNodes[0]->getTrialDisp();
    const Vector &TotDis2 = theNodes[1]->getTrialDisp();
    const Vector &TotDis3 = theNodes[2]->getTrialDisp();
    const Vector &TotDis4 = theNodes[3]->getTrialDisp();
    const Vector &TotDis5 = theNodes[4]->getTrialDisp();
    const Vector &TotDis6 = theNodes[5]->getTrialDisp();
    const Vector &TotDis7 = theNodes[6]->getTrialDisp();
    const Vector &TotDis8 = theNodes[7]->getTrialDisp();
    const Vector &TotDis9 = theNodes[8]->getTrialDisp();
    const Vector &TotDis10 = theNodes[9]->getTrialDisp();
    const Vector &TotDis11 = theNodes[10]->getTrialDisp();
    const Vector &TotDis12 = theNodes[11]->getTrialDisp();
    const Vector &TotDis13 = theNodes[12]->getTrialDisp();
    const Vector &TotDis14 = theNodes[13]->getTrialDisp();
    const Vector &TotDis15 = theNodes[14]->getTrialDisp();
    const Vector &TotDis16 = theNodes[15]->getTrialDisp();
    const Vector &TotDis17 = theNodes[16]->getTrialDisp();
    const Vector &TotDis18 = theNodes[17]->getTrialDisp();
    const Vector &TotDis19 = theNodes[18]->getTrialDisp();
    const Vector &TotDis20 = theNodes[19]->getTrialDisp();

    total_disp.val(1,1)=TotDis1(0); total_disp.val(1,2)=TotDis1(1);total_disp.val(1,3)=TotDis1(2);
    total_disp.val(2,1)=TotDis2(0); total_disp.val(2,2)=TotDis2(1);total_disp.val(2,3)=TotDis2(2);
    total_disp.val(3,1)=TotDis3(0); total_disp.val(3,2)=TotDis3(1);total_disp.val(3,3)=TotDis3(2);
    total_disp.val(4,1)=TotDis4(0); total_disp.val(4,2)=TotDis4(1);total_disp.val(4,3)=TotDis4(2);
    total_disp.val(5,1)=TotDis5(0); total_disp.val(5,2)=TotDis5(1);total_disp.val(5,3)=TotDis5(2);
    total_disp.val(6,1)=TotDis6(0); total_disp.val(6,2)=TotDis6(1);total_disp.val(6,3)=TotDis6(2);
    total_disp.val(7,1)=TotDis7(0); total_disp.val(7,2)=TotDis7(1);total_disp.val(7,3)=TotDis7(2);
    total_disp.val(8,1)=TotDis8(0); total_disp.val(8,2)=TotDis8(1);total_disp.val(8,3)=TotDis8(2);
    total_disp.val(9,1)=TotDis9(0); total_disp.val(9,2)=TotDis9(1);total_disp.val(9,3)=TotDis9(2);
    total_disp.val(10,1)=TotDis10(0); total_disp.val(10,2)=TotDis10(1);total_disp.val(10,3)=TotDis10(2);
    total_disp.val(11,1)=TotDis11(0); total_disp.val(11,2)=TotDis11(1);total_disp.val(11,3)=TotDis11(2);
    total_disp.val(12,1)=TotDis12(0); total_disp.val(12,2)=TotDis12(1);total_disp.val(12,3)=TotDis12(2);
    total_disp.val(13,1)=TotDis13(0); total_disp.val(13,2)=TotDis13(1);total_disp.val(13,3)=TotDis13(2);
    total_disp.val(14,1)=TotDis14(0); total_disp.val(14,2)=TotDis14(1);total_disp.val(14,3)=TotDis14(2);
    total_disp.val(15,1)=TotDis15(0); total_disp.val(15,2)=TotDis15(1);total_disp.val(15,3)=TotDis15(2);
    total_disp.val(16,1)=TotDis16(0); total_disp.val(16,2)=TotDis16(1);total_disp.val(16,3)=TotDis16(2);
    total_disp.val(17,1)=TotDis17(0); total_disp.val(17,2)=TotDis17(1);total_disp.val(17,3)=TotDis17(2);
    total_disp.val(18,1)=TotDis18(0); total_disp.val(18,2)=TotDis18(1);total_disp.val(18,3)=TotDis18(2);
    total_disp.val(19,1)=TotDis19(0); total_disp.val(19,2)=TotDis19(1);total_disp.val(19,3)=TotDis19(2);
    total_disp.val(20,1)=TotDis20(0); total_disp.val(20,2)=TotDis20(1);total_disp.val(20,3)=TotDis20(2);

    return total_disp;
  }

//=============================================================================
void TotalLagrangianFD20NodeBrick::zeroLoad(void)
{
    Q.Zero();
    return;
}


//=============================================================================
int
TotalLagrangianFD20NodeBrick::addLoad(ElementalLoad *theLoad, double loadFactor)
{
    opserr<<"TotalLagrangianFD20NodeBrick::addLoad - load type unknown for ele with tag: "<<this->getTag();          
    return -1;
}

//=============================================================================
int TotalLagrangianFD20NodeBrick::addInertiaLoadToUnbalance(const Vector &accel)
{
  // Check for a quick return
  if (rho == 0.0) return 0;

  // Get R * accel from the nodes
  const Vector &Raccel1  = theNodes[0]->getRV(accel);
  const Vector &Raccel2  = theNodes[1]->getRV(accel);
  const Vector &Raccel3  = theNodes[2]->getRV(accel);
  const Vector &Raccel4  = theNodes[3]->getRV(accel);
  const Vector &Raccel5  = theNodes[4]->getRV(accel);
  const Vector &Raccel6  = theNodes[5]->getRV(accel);
  const Vector &Raccel7  = theNodes[6]->getRV(accel);
  const Vector &Raccel8  = theNodes[7]->getRV(accel);
  const Vector &Raccel9  = theNodes[8]->getRV(accel);
  const Vector &Raccel10 = theNodes[9]->getRV(accel);
  const Vector &Raccel11 = theNodes[10]->getRV(accel);
  const Vector &Raccel12 = theNodes[11]->getRV(accel);
  const Vector &Raccel13 = theNodes[12]->getRV(accel);
  const Vector &Raccel14 = theNodes[13]->getRV(accel);
  const Vector &Raccel15 = theNodes[14]->getRV(accel);
  const Vector &Raccel16 = theNodes[15]->getRV(accel);
  const Vector &Raccel17 = theNodes[16]->getRV(accel);
  const Vector &Raccel18 = theNodes[17]->getRV(accel);
  const Vector &Raccel19 = theNodes[18]->getRV(accel);
  const Vector &Raccel20 = theNodes[19]->getRV(accel);

    if (NumDof != Raccel1.Size()  || NumDof != Raccel2.Size()  || NumDof != Raccel3.Size()  || NumDof != Raccel4.Size() ||
        NumDof != Raccel5.Size()  || NumDof != Raccel6.Size()  || NumDof != Raccel7.Size()  || NumDof != Raccel8.Size() ||
        NumDof != Raccel9.Size()  || NumDof != Raccel10.Size() || NumDof != Raccel11.Size() || NumDof != Raccel12.Size()||
        NumDof != Raccel13.Size() || NumDof != Raccel14.Size() || NumDof != Raccel15.Size() || NumDof != Raccel16.Size()||
        NumDof != Raccel17.Size() || NumDof != Raccel18.Size() || NumDof != Raccel19.Size() || NumDof != Raccel20.Size() )
  {
  opserr << "TotalLagrangianFD20NodeBrick::addInertiaLoadToUnbalance " << "matrix and vector sizes are incompatable\n";
  return -1;
  }

  static Vector a(NumElemDof);

    a( 0) = Raccel1(0);    a( 1) = Raccel1(1);    a( 2) = Raccel1(2);
    a( 3) = Raccel2(0);    a( 4) = Raccel2(1);    a( 5) = Raccel2(2);
    a( 6) = Raccel3(0);    a( 7) = Raccel3(1);    a( 8) = Raccel3(2);
    a( 9) = Raccel4(0);    a(10) = Raccel4(1);    a(11) = Raccel4(2);
    a(12) = Raccel5(0);    a(13) = Raccel5(1);    a(14) = Raccel5(2);
    a(15) = Raccel6(0);    a(16) = Raccel6(1);    a(17) = Raccel6(2);
    a(18) = Raccel7(0);    a(19) = Raccel7(1);    a(20) = Raccel7(2);
    a(21) = Raccel8(0);    a(22) = Raccel8(1);    a(23) = Raccel8(2);
    a(24) = Raccel9(0);    a(25) = Raccel9(1);    a(26) = Raccel9(2);
    a(27) = Raccel10(0);   a(28) = Raccel10(1);   a(29) = Raccel10(2);
    a(30) = Raccel11(0);   a(31) = Raccel11(1);   a(32) = Raccel11(2);
    a(33) = Raccel12(0);   a(34) = Raccel12(1);   a(35) = Raccel12(2);
    a(36) = Raccel13(0);   a(37) = Raccel13(1);   a(38) = Raccel13(2);
    a(39) = Raccel14(0);   a(40) = Raccel14(1);   a(41) = Raccel14(2);
    a(42) = Raccel15(0);   a(43) = Raccel15(1);   a(44) = Raccel15(2);
    a(45) = Raccel16(0);   a(46) = Raccel16(1);   a(47) = Raccel16(2);
    a(48) = Raccel17(0);   a(49) = Raccel17(1);   a(50) = Raccel17(2);
    a(51) = Raccel18(0);   a(52) = Raccel18(1);   a(53) = Raccel18(2);
    a(54) = Raccel19(0);   a(55) = Raccel19(1);   a(56) = Raccel19(2);
    a(57) = Raccel20(0);   a(58) = Raccel20(1);   a(59) = Raccel20(2);

    Q.addMatrixVector(1.0, M, a, -1.0);

    return 0;
}


//=============================================================================
const Vector &TotalLagrangianFD20NodeBrick::getResistingForce ()
{
    int f_dim[] = {NumNodes,NumDof};
    tensor NodalForces_in(2, f_dim, 0.0);
    NodalForces_in = this->getRtensor() - this->getForces();
    
    for (int i = 0; i < NumNodes; i++) {
      for (int j = 0; j < NumDof; j++) {
         P(i *3 + j) = NodalForces_in.cval(i+1, j+1);
      }
    }
    
    P.addVector(1.0, Q,-1.0);
    
    return P;
}

//=============================================================================
const Vector &TotalLagrangianFD20NodeBrick::getResistingForceIncInertia ()
{

  this->getResistingForce();

  if (rho != 0.0)
  {

    // form the mass matrix
    this->getMass();

    const Vector &accel1 = theNodes[0]->getTrialAccel();
    const Vector &accel2 = theNodes[1]->getTrialAccel();
    const Vector &accel3 = theNodes[2]->getTrialAccel();
    const Vector &accel4 = theNodes[3]->getTrialAccel();
    const Vector &accel5 = theNodes[4]->getTrialAccel();
    const Vector &accel6 = theNodes[5]->getTrialAccel();
    const Vector &accel7 = theNodes[6]->getTrialAccel();
    const Vector &accel8 = theNodes[7]->getTrialAccel();
    const Vector &accel9 = theNodes[8]->getTrialAccel();
    const Vector &accel10 = theNodes[9]->getTrialAccel();
    const Vector &accel11 = theNodes[10]->getTrialAccel();
    const Vector &accel12 = theNodes[11]->getTrialAccel();
    const Vector &accel13 = theNodes[12]->getTrialAccel();
    const Vector &accel14 = theNodes[13]->getTrialAccel();
    const Vector &accel15 = theNodes[14]->getTrialAccel();
    const Vector &accel16 = theNodes[15]->getTrialAccel();
    const Vector &accel17 = theNodes[16]->getTrialAccel();
    const Vector &accel18 = theNodes[17]->getTrialAccel();
    const Vector &accel19 = theNodes[18]->getTrialAccel();
    const Vector &accel20 = theNodes[19]->getTrialAccel();

    static Vector a(60);

    a( 0) = accel1(0);    a( 1) = accel1(1);    a( 2) = accel1(2);
    a( 3) = accel2(0);    a( 4) = accel2(1);    a( 5) = accel2(2);
    a( 6) = accel3(0);    a( 7) = accel3(1);    a( 8) = accel3(2);
    a( 9) = accel4(0);    a(10) = accel4(1);    a(11) = accel4(2);
    a(12) = accel5(0);    a(13) = accel5(1);    a(14) = accel5(2);
    a(15) = accel6(0);    a(16) = accel6(1);    a(17) = accel6(2);
    a(18) = accel7(0);    a(19) = accel7(1);    a(20) = accel7(2);
    a(21) = accel8(0);    a(22) = accel8(1);    a(23) = accel8(2);
    a(24) = accel9(0);    a(25) = accel9(1);    a(26) = accel9(2);
    a(27) = accel10(0);   a(28) = accel10(1);   a(29) = accel10(2);
    a(30) = accel11(0);   a(31) = accel11(1);   a(32) = accel11(2);
    a(33) = accel12(0);   a(34) = accel12(1);   a(35) = accel12(2);
    a(36) = accel13(0);   a(37) = accel13(1);   a(38) = accel13(2);
    a(39) = accel14(0);   a(40) = accel14(1);   a(41) = accel14(2);
    a(42) = accel15(0);   a(43) = accel15(1);   a(44) = accel15(2);
    a(45) = accel16(0);   a(46) = accel16(1);   a(47) = accel16(2);
    a(48) = accel17(0);   a(49) = accel17(1);   a(50) = accel17(2);
    a(51) = accel18(0);   a(52) = accel18(1);   a(53) = accel18(2);
    a(54) = accel19(0);   a(55) = accel19(1);   a(56) = accel19(2);
    a(57) = accel20(0);   a(58) = accel20(1);   a(59) = accel20(2);

    M = this->getMass();

        // P += M * a
    P.addMatrixVector(1.0, M, a, 1.0);

  }

  return P;
}

//=============================================================================
int TotalLagrangianFD20NodeBrick::sendSelf (int commitTag, Channel &theChannel)
{
     // Not implemtented yet
     return 0;
}

//=============================================================================
int TotalLagrangianFD20NodeBrick::recvSelf (int commitTag, Channel &theChannel,
FEM_ObjectBroker &theBroker)
{
     // Not implemtented yet
     return 0;
}


//=============================================================================
int TotalLagrangianFD20NodeBrick::displaySelf (Renderer &theViewer, int displayMode, float fact)
{
    //Needs more work...
    // first determine the end points of the quad based on
    // the display factor (a measure of the distorted image)
    // store this information in 4 3d vectors v1 through v4
    const Vector &end1Crd = theNodes[0]->getCrds();
    const Vector &end2Crd = theNodes[1]->getCrds();
    const Vector &end3Crd = theNodes[2]->getCrds();
    const Vector &end4Crd = theNodes[3]->getCrds();
    const Vector &end5Crd = theNodes[4]->getCrds();
    const Vector &end6Crd = theNodes[5]->getCrds();
    const Vector &end7Crd = theNodes[6]->getCrds();
    const Vector &end8Crd = theNodes[7]->getCrds();

    const Vector &end1Disp = theNodes[0]->getDisp();
    const Vector &end2Disp = theNodes[1]->getDisp();
    const Vector &end3Disp = theNodes[2]->getDisp();
    const Vector &end4Disp = theNodes[3]->getDisp();
    const Vector &end5Disp = theNodes[4]->getDisp();
    const Vector &end6Disp = theNodes[5]->getDisp();
    const Vector &end7Disp = theNodes[6]->getDisp();
    const Vector &end8Disp = theNodes[7]->getDisp();

    static Vector v1(3);
    static Vector v2(3);
    static Vector v3(3);
    static Vector v4(3);
    static Vector v5(3);
    static Vector v6(3);
    static Vector v7(3);
    static Vector v8(3);

    for (int i = 0; i < 2; i++)
    {
      v1(i) = end1Crd(i) + end1Disp(i)*fact;
      v2(i) = end2Crd(i) + end2Disp(i)*fact;
      v3(i) = end3Crd(i) + end3Disp(i)*fact;
      v4(i) = end4Crd(i) + end4Disp(i)*fact;
      v5(i) = end5Crd(i) + end5Disp(i)*fact;
      v6(i) = end6Crd(i) + end6Disp(i)*fact;
      v7(i) = end7Crd(i) + end7Disp(i)*fact;
      v8(i) = end8Crd(i) + end8Disp(i)*fact;
    }

    int error = 0;

    error += theViewer.drawLine (v1, v2, 1.0, 1.0);
    error += theViewer.drawLine (v2, v3, 1.0, 1.0);
    error += theViewer.drawLine (v3, v4, 1.0, 1.0);
    error += theViewer.drawLine (v4, v1, 1.0, 1.0);

    error += theViewer.drawLine (v5, v6, 1.0, 1.0);
    error += theViewer.drawLine (v6, v7, 1.0, 1.0);
    error += theViewer.drawLine (v7, v8, 1.0, 1.0);
    error += theViewer.drawLine (v8, v5, 1.0, 1.0);

    error += theViewer.drawLine (v1, v5, 1.0, 1.0);
    error += theViewer.drawLine (v2, v6, 1.0, 1.0);
    error += theViewer.drawLine (v3, v7, 1.0, 1.0);
    error += theViewer.drawLine (v4, v8, 1.0, 1.0);

    return error;

}

//=============================================================================
void TotalLagrangianFD20NodeBrick::Print(OPS_Stream &s, int flag)
{
    s << "\nTotalLagrangianFD20NodeBrick, element id:  " << this->getTag() << endln;
    s << "\nConnected external nodes:  " << connectedExternalNodes;
    s << "\nBody forces:  " << bf(0) << " " << bf(1) << " " << bf(2) << endln;

    theMaterial[0]->Print(s,flag);

    tensor sigma;
    Vector P00(6);
    
    for (int i=0; i<NumTotalGaussPts; i++)
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
Response * TotalLagrangianFD20NodeBrick::setResponse (const char **argv, int argc, Information &eleInfo)
{
    if (strcmp(argv[0],"force") == 0 || strcmp(argv[0],"forces") == 0)
      return new ElementResponse(this, 1, Vector(NumElemDof));
    
    else if (strcmp(argv[0],"stiff") == 0 || strcmp(argv[0],"stiffness") == 0)
      return new ElementResponse(this, 2, Matrix(matrixData, NumElemDof, NumElemDof));
    
    else if (strcmp(argv[0],"CauchyStress") == 0 || strcmp(argv[0],"stress") == 0)
      return new ElementResponse(this, 3, Vector(NumTotalGaussPts*6));
    
    else if (strcmp(argv[0],"PK2Stress") == 0 || strcmp(argv[0],"PK2stress") == 0)
      return new ElementResponse(this, 4, Vector(NumTotalGaussPts*6));

    // Added ZC 01/18/2005 to output strains
    else if (strcmp(argv[0],"EulerianStrain") == 0 || strcmp(argv[0],"strain") == 0)
      return new ElementResponse(this, 5, Vector(NumTotalGaussPts*6));
    
    else if (strcmp(argv[0],"LagrangianStrain") == 0 || strcmp(argv[0],"iniStrain") == 0)
      return new ElementResponse(this, 6, Vector(NumTotalGaussPts*6));

    else
      return 0;
}

//=============================================================================
int TotalLagrangianFD20NodeBrick::getResponse (int responseID, Information &eleInfo)
{
     static Vector P0(NumTotalGaussPts*6);
     
     switch (responseID) {
     
     case 1:
          return eleInfo.setVector(this->getResistingForce() );

     case 2:
          return eleInfo.setMatrix(this->getTangentStiff() );

     case 3: { 
        Vector P0(NumTotalGaussPts*6);
        tensor sigma; 
        for (int i=0; i<NumTotalGaussPts; i++) {
          sigma = theMaterial[i]->getCauchyStressTensor();
          P0( i*6+0 ) = sigma.val(1,1);
          P0( i*6+1 ) = sigma.val(2,2);
          P0( i*6+2 ) = sigma.val(3,3);
          P0( i*6+3 ) = sigma.val(2,3);
          P0( i*6+4 ) = sigma.val(3,1);
          P0( i*6+5 ) = sigma.val(1,2);
        }
        return eleInfo.setVector(P0);
     }

     case 4: { 
        Vector P0(NumTotalGaussPts*6);
        tensor sigma; 
        for (int i=0; i<NumTotalGaussPts; i++) {
          sigma = theMaterial[i]->getStressTensor();
          P0( i*6+0 ) = sigma.val(1,1);
          P0( i*6+1 ) = sigma.val(2,2);
          P0( i*6+2 ) = sigma.val(3,3);
          P0( i*6+3 ) = sigma.val(2,3);
          P0( i*6+4 ) = sigma.val(3,1);
          P0( i*6+5 ) = sigma.val(1,2);
        }
        return eleInfo.setVector(P0);
     }

    // Added ZC 01/18/2005 to output strains
     case 5: { 
        Vector P0(NumTotalGaussPts*6);
        tensor e;
	tensor E;
	tensor F;
	tensor tI2("I", 2, def_dim_2); 
        for (int i=0; i<NumTotalGaussPts; i++) {
          E = theMaterial[i]->getStrainTensor();
	  F = theMaterial[i]->getF();
	  F = F.inverse();
	  e = F("ki")*F("kj"); e.null_indices();
	  e = (tI2-e) *0.5;
          P0( i*6+0 ) = e.val(1,1);
          P0( i*6+1 ) = e.val(2,2);
          P0( i*6+2 ) = e.val(3,3);
          P0( i*6+3 ) = e.val(2,3);
          P0( i*6+4 ) = e.val(3,1);
          P0( i*6+5 ) = e.val(1,2);
        }
        return eleInfo.setVector(P0);
     }

     case 6: { 
        Vector P0(NumTotalGaussPts*6);
        tensor E; 
        for (int i=0; i<NumTotalGaussPts; i++) {
          E = theMaterial[i]->getStrainTensor();
          P0( i*6+0 ) = E.val(1,1);
          P0( i*6+1 ) = E.val(2,2);
          P0( i*6+2 ) = E.val(3,3);
          P0( i*6+3 ) = E.val(2,3);
          P0( i*6+4 ) = E.val(3,1);
          P0( i*6+5 ) = E.val(1,2);
        }
        return eleInfo.setVector(P0);
     }
    
     default:
        return -1;

     }    
}


//#############################################################################
//===================================================================
tensor TotalLagrangianFD20NodeBrick::shapeFunction(double r1, double r2, double r3)
  {

    int dimension[] = {NumNodes};
    tensor h(1, dimension, 0.0);

    // influence of the node number 20
        h.val(20)=(1.0+r1)*(1.0-r2)*(1.0-r3*r3)*0.25;
    // influence of the node number 19
        h.val(19)=(1.0-r1)*(1.0-r2)*(1.0-r3*r3)*0.25;
    // influence of the node number 18
        h.val(18)=(1.0-r1)*(1.0+r2)*(1.0-r3*r3)*0.25;
    // influence of the node number 17
        h.val(17)=(1.0+r1)*(1.0+r2)*(1.0-r3*r3)*0.25;

    // influence of the node number 16
        h.val(16)=(1.0+r1)*(1.0-r2*r2)*(1.0-r3)*0.25;
    // influence of the node number 15
        h.val(15)=(1.0-r1*r1)*(1.0-r2)*(1.0-r3)*0.25;
    // influence of the node number 14
        h.val(14)=(1.0-r1)*(1.0-r2*r2)*(1.0-r3)*0.25;
    // influence of the node number 13
        h.val(13)=(1.0-r1*r1)*(1.0+r2)*(1.0-r3)*0.25;

    // influence of the node number 12
        h.val(12)=(1.0+r1)*(1.0-r2*r2)*(1.0+r3)*0.25;
    // influence of the node number 11
        h.val(11)=(1.0-r1*r1)*(1.0-r2)*(1.0+r3)*0.25;
    // influence of the node number 10
        h.val(10)=(1.0-r1)*(1.0-r2*r2)*(1.0+r3)*0.25;
    // influence of the node number 9
        h.val( 9)=(1.0-r1*r1)*(1.0+r2)*(1.0+r3)*0.25;

      // influence of the node number 8
    h.val(8)=(1.0+r1)*(1.0-r2)*(1.0-r3)*0.125 - (h.val(15)+h.val(16)+h.val(20))*0.5;
      // influence of the node number 7
    h.val(7)=(1.0-r1)*(1.0-r2)*(1.0-r3)*0.125 - (h.val(14)+h.val(15)+h.val(19))*0.5;
      // influence of the node number 6
    h.val(6)=(1.0-r1)*(1.0+r2)*(1.0-r3)*0.125 - (h.val(13)+h.val(14)+h.val(18))*0.5;
      // influence of the node number 5
    h.val(5)=(1.0+r1)*(1.0+r2)*(1.0-r3)*0.125 - (h.val(13)+h.val(16)+h.val(17))*0.5;

      // influence of the node number 4
    h.val(4)=(1.0+r1)*(1.0-r2)*(1.0+r3)*0.125 - (h.val(11)+h.val(12)+h.val(20))*0.5;
      // influence of the node number 3
    h.val(3)=(1.0-r1)*(1.0-r2)*(1.0+r3)*0.125 - (h.val(10)+h.val(11)+h.val(19))*0.5;
      // influence of the node number 2
    h.val(2)=(1.0-r1)*(1.0+r2)*(1.0+r3)*0.125 - (h.val(10)+h.val(18)+h.val( 9))*0.5;
      // influence of the node number 1
    h.val(1)=(1.0+r1)*(1.0+r2)*(1.0+r3)*0.125 - (h.val(12)+h.val(17)+h.val( 9))*0.5;

    return h;
  }


//==============================================================
tensor TotalLagrangianFD20NodeBrick::shapeFunctionDerivative(double r1, double r2, double r3)
  {

    int dimensions[] = {NumNodes,NumDof};
    tensor dh(2, dimensions, 0.0);

    // influence of the node number 20
        dh.val(20,1) =   (1.0-r2)*(1.0-r3*r3)*0.25;
        dh.val(20,2) = - (1.0+r1)*(1.0-r3*r3)*0.25;
        dh.val(20,3) = - (1.0+r1)*(1.0-r2)*r3*0.50;
    // influence of the node number 19
        dh.val(19,1) = - (1.0-r2)*(1.0-r3*r3)*0.25;
        dh.val(19,2) = - (1.0-r1)*(1.0-r3*r3)*0.25;
        dh.val(19,3) = - (1.0-r1)*(1.0-r2)*r3*0.50;
    // influence of the node number 18
        dh.val(18,1) = - (1.0+r2)*(1.0-r3*r3)*0.25;
        dh.val(18,2) =   (1.0-r1)*(1.0-r3*r3)*0.25;
        dh.val(18,3) = - (1.0-r1)*(1.0+r2)*r3*0.50;
    // influence of the node number 17
        dh.val(17,1) =   (1.0+r2)*(1.0-r3*r3)*0.25;
        dh.val(17,2) =   (1.0+r1)*(1.0-r3*r3)*0.25;
        dh.val(17,3) = - (1.0+r1)*(1.0+r2)*r3*0.50;

    // influence of the node number 16
        dh.val(16,1) =   (1.0-r2*r2)*(1.0-r3)*0.25;
        dh.val(16,2) = - (1.0+r1)*r2*(1.0-r3)*0.50;
        dh.val(16,3) = - (1.0+r1)*(1.0-r2*r2)*0.25;
    // influnce of the node number 15
        dh.val(15,1) = - r1*(1.0-r2)*(1.0-r3)*0.50;
        dh.val(15,2) = - (1.0-r1*r1)*(1.0-r3)*0.25;
        dh.val(15,3) = - (1.0-r1*r1)*(1.0-r2)*0.25;
    // influence of the node number 14
        dh.val(14,1) = - (1.0-r2*r2)*(1.0-r3)*0.25;
        dh.val(14,2) = - (1.0-r1)*r2*(1.0-r3)*0.50;
        dh.val(14,3) = - (1.0-r1)*(1.0-r2*r2)*0.25;
    // influence of the node number 13
        dh.val(13,1) = - r1*(1.0+r2)*(1.0-r3)*0.50;
        dh.val(13,2) =   (1.0-r1*r1)*(1.0-r3)*0.25;
        dh.val(13,3) = - (1.0-r1*r1)*(1.0+r2)*0.25;

    // influence of the node number 12
        dh.val(12,1) =   (1.0-r2*r2)*(1.0+r3)*0.25;
        dh.val(12,2) = - (1.0+r1)*r2*(1.0+r3)*0.50;
        dh.val(12,3) =   (1.0+r1)*(1.0-r2*r2)*0.25;
    // influence of the node number 11
        dh.val(11,1) = - r1*(1.0-r2)*(1.0+r3)*0.50;
        dh.val(11,2) = - (1.0-r1*r1)*(1.0+r3)*0.25;
        dh.val(11,3) =   (1.0-r1*r1)*(1.0-r2)*0.25;
    // influence of the node number 10
        dh.val(10,1) = - (1.0-r2*r2)*(1.0+r3)*0.25;
        dh.val(10,2) = - (1.0-r1)*r2*(1.0+r3)*0.50;
        dh.val(10,3) =   (1.0-r1)*(1.0-r2*r2)*0.25;
    // influence of the node number 9
        dh.val(9,1)  = - r1*(1.0+r2)*(1.0+r3)*0.50;
        dh.val(9,2)  =   (1.0-r1*r1)*(1.0+r3)*0.25;
        dh.val(9,3)  =   (1.0-r1*r1)*(1.0+r2)*0.25;

      // influence of the node number 8
    dh.val(8,1)= (1.0-r2)*(1.0-r3)*0.125 - (dh.val(15,1)+dh.val(16,1)+dh.val(20,1))*0.5;
    dh.val(8,2)=-(1.0+r1)*(1.0-r3)*0.125 - (dh.val(15,2)+dh.val(16,2)+dh.val(20,2))*0.5;
    dh.val(8,3)=-(1.0+r1)*(1.0-r2)*0.125 - (dh.val(15,3)+dh.val(16,3)+dh.val(20,3))*0.5;
      // influence of the node number 7
    dh.val(7,1)=-(1.0-r2)*(1.0-r3)*0.125 - (dh.val(14,1)+dh.val(15,1)+dh.val(19,1))*0.5;
    dh.val(7,2)=-(1.0-r1)*(1.0-r3)*0.125 - (dh.val(14,2)+dh.val(15,2)+dh.val(19,2))*0.5;
    dh.val(7,3)=-(1.0-r1)*(1.0-r2)*0.125 - (dh.val(14,3)+dh.val(15,3)+dh.val(19,3))*0.5;
      // influence of the node number 6
    dh.val(6,1)=-(1.0+r2)*(1.0-r3)*0.125 - (dh.val(13,1)+dh.val(14,1)+dh.val(18,1))*0.5;
    dh.val(6,2)= (1.0-r1)*(1.0-r3)*0.125 - (dh.val(13,2)+dh.val(14,2)+dh.val(18,2))*0.5;
    dh.val(6,3)=-(1.0-r1)*(1.0+r2)*0.125 - (dh.val(13,3)+dh.val(14,3)+dh.val(18,3))*0.5;
      // influence of the node number 5
    dh.val(5,1)= (1.0+r2)*(1.0-r3)*0.125 - (dh.val(13,1)+dh.val(16,1)+dh.val(17,1))*0.5;
    dh.val(5,2)= (1.0+r1)*(1.0-r3)*0.125 - (dh.val(13,2)+dh.val(16,2)+dh.val(17,2))*0.5;
    dh.val(5,3)=-(1.0+r1)*(1.0+r2)*0.125 - (dh.val(13,3)+dh.val(16,3)+dh.val(17,3))*0.5;

      // influence of the node number 4
    dh.val(4,1)= (1.0-r2)*(1.0+r3)*0.125 - (dh.val(11,1)+dh.val(12,1)+dh.val(20,1))*0.5;
    dh.val(4,2)=-(1.0+r1)*(1.0+r3)*0.125 - (dh.val(11,2)+dh.val(12,2)+dh.val(20,2))*0.5;
    dh.val(4,3)= (1.0+r1)*(1.0-r2)*0.125 - (dh.val(11,3)+dh.val(12,3)+dh.val(20,3))*0.5;
      // influence of the node number 3
    dh.val(3,1)=-(1.0-r2)*(1.0+r3)*0.125 - (dh.val(10,1)+dh.val(11,1)+dh.val(19,1))*0.5;
    dh.val(3,2)=-(1.0-r1)*(1.0+r3)*0.125 - (dh.val(10,2)+dh.val(11,2)+dh.val(19,2))*0.5;
    dh.val(3,3)= (1.0-r1)*(1.0-r2)*0.125 - (dh.val(10,3)+dh.val(11,3)+dh.val(19,3))*0.5;
      // influence of the node number 2
    dh.val(2,1)=-(1.0+r2)*(1.0+r3)*0.125 - (dh.val(10,1)+dh.val(18,1)+dh.val( 9,1))*0.5;
    dh.val(2,2)= (1.0-r1)*(1.0+r3)*0.125 - (dh.val(10,2)+dh.val(18,2)+dh.val( 9,2))*0.5;
    dh.val(2,3)= (1.0-r1)*(1.0+r2)*0.125 - (dh.val(10,3)+dh.val(18,3)+dh.val( 9,3))*0.5;
      // influence of the node number 1
    dh.val(1,1)= (1.0+r2)*(1.0+r3)*0.125 - (dh.val(12,1)+dh.val(17,1)+dh.val( 9,1))*0.5;
    dh.val(1,2)= (1.0+r1)*(1.0+r3)*0.125 - (dh.val(12,2)+dh.val(17,2)+dh.val( 9,2))*0.5;
    dh.val(1,3)= (1.0+r1)*(1.0+r2)*0.125 - (dh.val(12,3)+dh.val(17,3)+dh.val( 9,3))*0.5;

    return dh;
  }


#endif

