// MultiLinearKp.cpp: implementation of the MultiLinearKp class.
//
//////////////////////////////////////////////////////////////////////

#include "MultiLinearKp.h"
#include <stdlib.h>

#define MAT_TAG_MULTILINEAR -1
#define DEBG 0
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

MultiLinearKp::MultiLinearKp(int tag, Vector &sum_plas_defo, Vector &kp)
:PlasticHardeningMaterial(tag,MAT_TAG_MULTILINEAR),
sumPlasDefo(sum_plas_defo.Size()+1), Kp(kp.Size()+1)
{
 	if(sumPlasDefo.Size() != Kp.Size())
 	{
 		opserr << "ERROR MultiLinear() - incompatible vector sizes\n";
 	}

 	numPoints = sum_plas_defo.Size();

 	for(int i=0; i < numPoints; i++)
 	{
 		sumPlasDefo(i) =  sum_plas_defo(i);
 		if(sumPlasDefo(i) < 0)
 			opserr << "ERROR MultiLinear() -  sumPlasDefo < 0\n";

 		Kp(i) = kp(i);
 	}
 	if(sumPlasDefo(0) != 0)
 		opserr << "WARNING MultiLinear() -  sumPlasDefo(0) != 0\n";

 	Kp(numPoints)          = Kp(numPoints -1);
 	sumPlasDefo(numPoints) = sumPlasDefo(numPoints -1)*1000;

 	//cout << *this;
}

MultiLinearKp::~MultiLinearKp()
{

}


double MultiLinearKp::getTrialPlasticStiffness()
{
double K = 0;
double sumDisp = val_trial;

	if( sumDisp > sumPlasDefo(numPoints-1))
	{
		K = residual*Kp(numPoints-1);
		
		if(sFactor != 1.0)
			K = Kp(0)*sFactor;
		return K;
	}
		
	for(int i=0; i< numPoints; i++)
	{
		double x1 = sumPlasDefo(i);
		double y1 = Kp(i);
		double x2 = sumPlasDefo(i+1);
		double y2 = Kp(i+1);
		
		if(sumDisp < x2 && sumDisp >= x1)
		{
			if (sumDisp == x1)
				return y1;
				
			if(x2 == x1)
			{
				opserr << "WARNING - MultiLinear::getTangent() x2 -x1 = 0 \n";
				return 0;
			}
			
			double m = (y2 - y1)/(x2 -x1);
		    double b = y1 - m*x1;
		    K = m*sumDisp + b;
		 	break;
		}
	}

	if(sFactor != 1.0)
		K = Kp(0)*sFactor;
	else
		K = residual*K;

	// opserr << "K = " << K << ", sFactor = " << sFactor << endln; // opserr << "\a";
	return K;
}


void MultiLinearKp::Print(OPS_Stream &s, int flag)
{
	this->PlasticHardeningMaterial::Print(s, flag);
	
	s << "+-MultiLinear" << endln;
	s << "    Kp = " << this->getTrialPlasticStiffness();
	s << "    SumPlasDefo Vector = " <<  sumPlasDefo;
	s << "    Kp Vector          = " <<  Kp << endln;
}

PlasticHardeningMaterial *MultiLinearKp::getCopy(void)
{
	Vector spd(numPoints);
    Vector kp(numPoints);

    for(int i =0; i < numPoints; i++)
    {
    	spd(i) =  sumPlasDefo(i);
    	kp(i)  =  Kp(i);
    }

    // Don't want to pass the actual vectors or else the size will
    // keep on increasing by 1
 	PlasticHardeningMaterial *theMat = new MultiLinearKp(getTag(), spd, kp);
    return theMat;
}

