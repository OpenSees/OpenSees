/***************************************************************************
                          NullEvolution.cpp  -  description
                             -------------------
    begin                : Thu Aug 1 2002
    email                : rkaul@ce-blume215-pent-2.stanford.edu
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *                                                                         *
 ***************************************************************************/

#include "NullEvolution.h"
#define NULL_EVOL_CLASS_TAG -1

Vector NullEvolution::vec_dim_1(1);
Vector NullEvolution::vec_dim_2(2);
Vector NullEvolution::vec_dim_3(3);

NullEvolution::NullEvolution(int tag, double isox)
:YS_Evolution(tag, NULL_EVOL_CLASS_TAG, 0.0, 0.0, 1, 0.0, 0.0)
{
	isotropicFactor(0) = isox;
	isotropicFactor_hist(0) = isox;
}


NullEvolution::NullEvolution(int tag, double isox, double isoy)
:YS_Evolution(tag, NULL_EVOL_CLASS_TAG, 0.0, 0.0, 2, 0.0, 0.0)
{
	isotropicFactor(0) = isox;
	isotropicFactor(1) = isoy;
	isotropicFactor_hist(0) = isox;
	isotropicFactor_hist(1) = isoy;

}


NullEvolution::NullEvolution(int tag, double isox, double isoy, double isoz)
:YS_Evolution(tag, NULL_EVOL_CLASS_TAG, 0.0, 0.0, 3, 0.0, 0.0)
{
	isotropicFactor(0) = isox;
	isotropicFactor(1) = isoy;
	isotropicFactor(2) = isoz;
	isotropicFactor_hist(0) = isox;
	isotropicFactor_hist(1) = isoy;
	isotropicFactor_hist(2) = isoz;
}

NullEvolution::~NullEvolution()
{
}



/** No descriptions */
int NullEvolution::evolveSurface(YieldSurface_BC * ys, double magPlasticDefo,
                                Vector & G, Vector & F_Surface, int flag)
{
	// does nothing
	return 0;
}

/** No descriptions */
YS_Evolution * NullEvolution::getCopy()
{
NullEvolution *copy=0;

	if(dimension == 1)
		copy = new NullEvolution(getTag(), isotropicFactor(0));
	else if(dimension == 2)
		copy = new NullEvolution(getTag(), isotropicFactor(0), isotropicFactor(1));
	else if(dimension == 3)
		copy = new NullEvolution(getTag(), isotropicFactor(0), isotropicFactor(1), isotropicFactor(2));
	else
		copy = 0;
	return copy;
}

/** No descriptions */
const Vector & NullEvolution::getEquiPlasticStiffness()
{
	if(dimension == 1)
		return vec_dim_1;
	else if (dimension == 2)
		return vec_dim_2;
	else if (dimension == 3)
		return vec_dim_3;
	else
		opserr << "NullEvolution::getEquiPlasticStiffness() - error dimension\n";

	return vec_dim_3;
}

double NullEvolution::getTrialPlasticStrains(int dof)
{
	return 0;
}

double NullEvolution::getCommitPlasticStrains(int dof)
{
	return 0;
}


/** No descriptions */
int NullEvolution::getResponse(int responseID, Information & info)
{
	return 0;
}

/** No descriptions */
Response * NullEvolution::setResponse(char * * argv, int argc, Information & info)
{
	return 0;
}

/** No descriptions */
int NullEvolution::commitState(int status)
{
	return 0;
}
/** No descriptions */
int NullEvolution::revertToLastCommit()
{
	return 0;
}

/*
// implement these or else the base class will complain
void	NullEvolution::toDeformedCoord(double &x){
	}
void	NullEvolution::toDeformedCoord(double &x, double &y){
	}
void	NullEvolution::toDeformedCoord(double &x, double &y, double &z){
	}

void	NullEvolution::toOriginalCoord(double &x){
	}
void	NullEvolution::toOriginalCoord(double &x, double &y){
	}
void	NullEvolution::toOriginalCoord(double &x, double &y, double &z){
	}
*/




