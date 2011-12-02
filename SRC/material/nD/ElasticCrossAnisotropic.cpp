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
//# PURPOSE:           Elastic Cross Anisotropic Material implementation:
//# CLASS:             ElasticIsotropic3D
//#
//# VERSION:           0.61803398874989 (golden section)
//# LANGUAGE:          C++
//# TARGET OS:         all...
//# DESIGN:            Zhaohui Yang, Boris Jeremic (jeremic@ucdavis.edu)
//# PROGRAMMER(S):     Zhaohui Yang, Yi Bian, Boris Jeremic
//#
//#
//# DATE:              10Oct2002 Yi Bian
//# UPDATE HISTORY:    March 20, 2003 Revised by Joey Yang & Boris Jeremic, UC Davis
//#                    Aug2006   Z.Cheng
//#
//===============================================================================
//


#include <ElasticCrossAnisotropic.h>

Tensor ElasticCrossAnisotropic::Dt(4, def_dim_4, 0.0 );
stresstensor ElasticCrossAnisotropic::Stress;


///////////////////////////////////////////////////////////////////////////////
ElasticCrossAnisotropic::ElasticCrossAnisotropic(int tag, 
                                                 double Ehp, 
                                                 double Evp, 
                                                 double nuhvp, 
                                                 double nuhhp, 
                                                 double Ghvp, 
                                                 double rhop):
NDMaterial(tag, ND_TAG_ElasticCrossAnisotropic3D),
Eh(Ehp), 
Ev(Evp), 
nuhv(nuhvp),
nuhh(nuhhp),
Ghv(Ghvp),
rho(rhop)
{

}

///////////////////////////////////////////////////////////////////////////////
ElasticCrossAnisotropic::ElasticCrossAnisotropic()
{

}

///////////////////////////////////////////////////////////////////////////////
ElasticCrossAnisotropic::~ElasticCrossAnisotropic ()
{

}

///////////////////////////////////////////////////////////////////////////////
double ElasticCrossAnisotropic::getrho()
{
	return rho;
}

///////////////////////////////////////////////////////////////////////////////
double  ElasticCrossAnisotropic::getMatParameter(int MatParameterID)
{
	switch (MatParameterID) {
	    case (1):
		    return Eh;
	    case (2):
		    return Ev;
	    case (3):
		    return nuhv;
	    case (4):
		    return nuhh;
	    case (5):
		return Ghv;
	    case (6):
		    return rho;
	    default: {
		    opserr << "Warning! ElasticIsotropic3D:: incorrect Materal parameter ID" << "\n";
		    return 0.0;
        }
    }
}

///////////////////////////////////////////////////////////////////////////////
NDMaterial* ElasticCrossAnisotropic::getCopy (const char *type)
{
    if (strcmp(type,"ThreeDimensional") == 0) {
		ElasticCrossAnisotropic *theModel;
		theModel = new ElasticCrossAnisotropic (this->getTag(), Eh, Ev, nuhv, nuhh, Ghv, rho);
		return theModel;
    }

    else {
		opserr <<"ElasticCrossAnisotropic::getModel failed to get model " << type << "\n";
		return 0;
    }
}

int
ElasticCrossAnisotropic::setTrialStrain (const Tensor &v)
{
    Strain = v;
    return 0;
}

int
ElasticCrossAnisotropic::setTrialStrain (const Tensor &v, const Tensor &r)
{
    Strain = v;
    return 0;
}

int
ElasticCrossAnisotropic::setTrialStrainIncr (const Tensor &v)
{
    Strain += v;
    return 0;
}

int
ElasticCrossAnisotropic::setTrialStrainIncr (const Tensor &v, const Tensor &r)
{
    Strain += v;
    return 0;
}

const Tensor&
ElasticCrossAnisotropic::getTangentTensor (void)
{
   //Old codes from Yi Bian
   //double A = 1/((1 + nuhv) * (1 - 2 * nuhv));
   //D(0, 0) = A * (1 - nuhv) * Ev;
   //D(1, 1) = D(2, 2) = A * (1 - nuhv) * Eh;
   //D(0, 1) = D(0, 2) = D(1, 0) = D(2, 0) = A * Eh * nuhh;
   //D(3, 3) = Eh/(1 + nuhv);
   //D(4, 4) = D(5, 5) = 2 * Ghv;

   //  Compliance matrix C Refer to Thor C. Heimdahl and Andrew Drescher
   //  "Elastic Anisotropy of Tire Shreds", {ASCE} Journal of Geotechnical
   //  and Geoenvironmental Engineering, 1999, pp. 383-389
   //  |Sxx|    | 1/Eh     -nuhh/Eh -nuhv/Ev     0       0     0    |
   //  |Syy|    |-nuhh/Eh    1/Eh   -nuhv/Ev     0       0     0    |
   //  |Szz|    |-nuhv/Ev  -nuhv/Ev   1/Ev       0       0     0    |
   //  |Sxy| C= |   0         0       0  2(1+nuhh)/Eh     0     0    |
   //  |Sxz|    |   0         0       0         0     1/(Ghv)   0    |
   //  |Syz|    |   0         0       0         0         0   1/(Ghv)|
   
   // Form compliance matrix D
   Matrix D(6,6);
   double A = 1.0/Eh;
   double B = 1.0/Ev;
   D(0,0) = D(1,1) = A;
   D(2,2) = B;
   D(0,1) = D(1,0) = -nuhh*A;
   D(0,2) = D(2,0) = D(1,2) = D(2,1) = -nuhv*B;
   D(3,3) = 2*(1.0+nuhh)*A;
   D(4,4) = D(5,5) = 1.0/Ghv;

   D.Invert( D );

   Dt.val(1,1,1,1) = D(0,0); 
   Dt.val(1,1,2,2) = D(0,1); 
   Dt.val(1,1,3,3) = D(0,2); // --> Sigma_xx

   Dt.val(1,2,1,2) = D(3,3); 
   Dt.val(1,2,2,1) = D(3,3); // --> Sigma_xy

   Dt.val(1,3,1,3) = D(4,4); 
   Dt.val(1,3,3,1) = D(4,4); // --> Sigma_xz

   Dt.val(2,1,1,2) = D(3,3); 
   Dt.val(2,1,2,1) = D(3,3); // --> Sigma_yx

   Dt.val(2,2,1,1) = D(1,0); 
   Dt.val(2,2,2,2) = D(1,1); 
   Dt.val(2,2,3,3) = D(1,2); // --> Sigma_yy

   Dt.val(2,3,2,3) = D(5,5); 
   Dt.val(2,3,3,2) = D(5,5); // --> Sigma_yz

   Dt.val(3,1,1,3) = D(4,4); 
   Dt.val(3,1,3,1) = D(4,4); // --> Sigma_zx

   Dt.val(3,2,2,3) = D(5,5); 
   Dt.val(3,2,3,2) = D(5,5); // --> Sigma_zy

   Dt.val(3,3,1,1) = D(2,0); 
   Dt.val(3,3,2,2) = D(2,1); 
   Dt.val(3,3,3,3) = D(2,2); // --> Sigma_zz
   
   return Dt;
}

const stresstensor& ElasticCrossAnisotropic::getStressTensor (void)
{
    Tensor Dt0 = getTangentTensor();
    Stress = Dt0("ijkl") * Strain("kl");
    
    return Stress;
}

const straintensor& ElasticCrossAnisotropic::getStrainTensor (void)
{
	return Strain;
}

int
ElasticCrossAnisotropic::commitState (void)
{
	return 0;
}

int
ElasticCrossAnisotropic::revertToLastCommit (void)
{
	return 0;
}

int
ElasticCrossAnisotropic::revertToStart (void)
{
	// added: C.McGann, U.Washington for InitialStateAnalysis
	if (ops_InitialStateAnalysis) {
		// do nothing, keep state variables from last step
	} else {
		// normal call for revertToStart (not initialStateAnalysis)
		Strain = Strain*0.0;
	}
	
	return 0;
}

NDMaterial*
ElasticCrossAnisotropic::getCopy (void)
{
	ElasticCrossAnisotropic *theCopy =
	  new ElasticCrossAnisotropic (this->getTag(), Eh, Ev, nuhv, nuhh, Ghv, rho);
	
    return theCopy;
}

const char*
ElasticCrossAnisotropic::getType (void) const
{
	return "ThreeDimensional";
}

int
ElasticCrossAnisotropic::sendSelf(int commitTag, Channel &theChannel)
{
	// Need work here.
	return 0;
}

int
ElasticCrossAnisotropic::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
	// Need work here.
	return 0;
}

void
ElasticCrossAnisotropic::Print (OPS_Stream &s, int flag)
{
	s << "Elastic Cross-Anisotropic Material Model\n";
	s << "\tEh:  " << Eh << "\tEv:  " << Ev << "\n";
	s << "\tnuhv:  " << nuhv << "\tnuhh:  " << nuhh << "\n";
	s << "\tGhv:  " << Ghv << "\trho:  " << rho << "\n";
	
	return;
}

