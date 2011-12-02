//#############################################################################
//# COPYRIGHT (C):     :-))                                                   #
//# PROJECT:           Object Oriented Finite Element Program                 #
//# PURPOSE:                                                                  #
//#                                                                           #
//# CLASS:                                                                    #
//#                                                                           #
//# VERSION:                                                                  #
//# LANGUAGE:          C++																																																				#
//# TARGET OS:         DOS || UNIX || . . .                                   #
//# DESIGNER(S):       Boris Jeremic, Zhaohui Yang                            #
//# PROGRAMMER(S):     Boris Jeremic, Zhaohui Yang                            #
//# CONTACT:           jeremic@ucdavis.edu                                    #
//#                                                                           #
//#                                                                           #
//# DATE:              Aug, Sept, Oct 2000                                    #
//# UPDATE HISTORY:                                                           #
//#                                                                           #
//#                                                                           #
//#                                                                           #
//#                                                                           #
//# SHORT EXPLANATION: 																																																							#
//#                                                                           #
//#                                                                           #
//#                                                                           #
//#                                                                           #
//#############################################################################

//$Source: /usr/local/cvs/OpenSees/SRC/material/nD/ElasticIsotropic3D.cpp,v $
//$Date: 2000-12-18 10:50:41 $
//$Revision: 1.2 $
                                                                        
#include <ElasticIsotropic3D.h>
#include <Channel.h>
#include <Tensor.h>

ElasticIsotropic3D::ElasticIsotropic3D
(int tag, double E, double nu) :
 ElasticIsotropicMaterial (tag, ND_TAG_ElasticIsotropic3D, E, nu),
 sigma(6), D(6,6), epsilon(6)
{
	// Set up the elastic constant matrix for 3D elastic isotropic 
	D.Zero();
        Dt = tensor( 4, def_dim_4, 0.0 ); 
	setElasticStiffness();

}

ElasticIsotropic3D::ElasticIsotropic3D():
 ElasticIsotropicMaterial (0, ND_TAG_ElasticIsotropic3D, 0.0, 0.0),
 sigma(6), D(6,6), epsilon(6)
{
       Dt = tensor( 4, def_dim_4, 0.0 );
}

ElasticIsotropic3D::~ElasticIsotropic3D ()
{

}

int
ElasticIsotropic3D::setTrialStrain (const Vector &v)
{
	epsilon = v;

	return 0;
}

int
ElasticIsotropic3D::setTrialStrain (const Vector &v, const Vector &r)
{
	epsilon = v;

	return 0;
}

int
ElasticIsotropic3D::setTrialStrainIncr (const Vector &v)
{
	epsilon += v;

	return 0;
}

int
ElasticIsotropic3D::setTrialStrainIncr (const Vector &v, const Vector &r)
{
	epsilon += v;

	return 0;
}

const Matrix&
ElasticIsotropic3D::getTangent (void)
{
	return D;
}

const Vector&
ElasticIsotropic3D::getStress (void)
{
	sigma = D*epsilon;
	return sigma;
}

const Vector&
ElasticIsotropic3D::getStrain (void)
{
	return epsilon;
}

int
ElasticIsotropic3D::setTrialStrain (const Tensor &v)
{
    Strain = v;
    return 0;
}

int
ElasticIsotropic3D::setTrialStrain (const Tensor &v, const Tensor &r)
{
    Strain = v;
    return 0;
}

int
ElasticIsotropic3D::setTrialStrainIncr (const Tensor &v)
{
    Strain = Strain + v;
    return 0;
}

int
ElasticIsotropic3D::setTrialStrainIncr (const Tensor &v, const Tensor &r)
{
    Strain = Strain + v;
    return 0;
}

const Tensor&
ElasticIsotropic3D::getTangentTensor (void)
{
    return Dt;
}

const Tensor&
ElasticIsotropic3D::getStressTensor (void)
{
    Stress = Dt("ijkl") * Strain("kl");
    return Stress;
}

const Tensor&
ElasticIsotropic3D::getStrainTensor (void)
{
    return Strain;
}

int
ElasticIsotropic3D::commitState (void)
{
	// Nothing to commit
	return 0;
}

int
ElasticIsotropic3D::revertToLastCommit (void)
{
	// Nothing was previously committed
	return 0;
}

int
ElasticIsotropic3D::revertToStart (void)
{
	// Nothing was previously committed
	return 0;
}

NDMaterial*
ElasticIsotropic3D::getCopy (void)
{
	ElasticIsotropic3D *theCopy =
		new ElasticIsotropic3D (this->getTag(), E, v);
	theCopy->epsilon = this->epsilon;
	theCopy->sigma = this->sigma;
	theCopy->Strain = this->Strain;
	theCopy->Stress = this->Stress;
	// D and Dt are created in the constructor call

	return theCopy;
}

const char*
ElasticIsotropic3D::getType (void) const
{
	return "ElasticIsotropic3D";
}

int
ElasticIsotropic3D::getOrder (void) const
{
	return 6;
}

int 
ElasticIsotropic3D::sendSelf(int commitTag, Channel &theChannel)
{
	int res = 0;

	static Vector data(3);

	data(0) = this->getTag();
	data(1) = E;
	data(2) = v;

    	res += theChannel.sendVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
		g3ErrorHandler->warning("%s -- could not send Vector",
			"ElasticIsotropic3D::sendSelf");
		return res;
	}

	return res;
}

int
ElasticIsotropic3D::recvSelf(int commitTag, Channel &theChannel, 
		 FEM_ObjectBroker &theBroker)
{
	int res = 0;

	static Vector data(6);

	res += theChannel.recvVector(this->getDbTag(), commitTag, data);
	if (res < 0) {
		g3ErrorHandler->warning("%s -- could not receive Vector",
			"ElasticIsotropic3D::recvSelf");
		return res;
	}
    
	this->setTag((int)data(0));
    	E = data(1);
	v = data(2);

	// Set up the elastic constant matrix for 3D elastic isotropic
	D.Zero();
	setElasticStiffness();
	
	return res;
}
    
void
ElasticIsotropic3D::Print(ostream &s, int flag)
{
	s << "ElasticIsotropic3D" << endl;
	s << "\ttag: " << this->getTag() << endl;
	s << "\tE: " << E << endl;
	s << "\tv: " << v << endl;
	//s << "\tD: " << D << endl;
}


//================================================================================
void ElasticIsotropic3D::setElasticStiffness(void)
{    
    tensor ret( 4, def_dim_4, 0.0 );
    				       
    // Kronecker delta tensor
    tensor I2("I", 2, def_dim_2);

    tensor I_ijkl = I2("ij")*I2("kl");


    //I_ijkl.null_indices();
    tensor I_ikjl = I_ijkl.transpose0110();
    tensor I_iljk = I_ijkl.transpose0111();
    tensor I4s = (I_ikjl+I_iljk)*0.5;

    // Building elasticity tensor
    ret = I_ijkl*( E*v / ( (1.0+v)*(1.0 - 2.0*v) ) ) + I4s*( E / (1.0 + v) );
    
    //ret.print();
    Dt = ret;
    //D = Dt;

    return;

}


