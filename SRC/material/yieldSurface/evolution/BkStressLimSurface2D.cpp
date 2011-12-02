//BkStressLimSurface2D
//
//////////////////////////////////////////////////////////////////////

#include "BkStressLimSurface2D.h"
#include <YieldSurface_BC.h>
#include <math.h>

#define evolDebug 0
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

BkStressLimSurface2D::BkStressLimSurface2D(int tag, int classTag, double min_iso_factor,
						double iso_ratio, double kin_ratio, YieldSurface_BC &lim_surface,
						PlasticHardeningMaterial &kinX,
						PlasticHardeningMaterial &kinY,
						PlasticHardeningMaterial &isoXPos,
						PlasticHardeningMaterial &isoXNeg,
						PlasticHardeningMaterial &isoYPos,
						PlasticHardeningMaterial &isoYNeg,
						int restype, double res_Fact, double app_Fact, double dir
						)
:YS_Evolution2D(tag, classTag, min_iso_factor, iso_ratio, kin_ratio),
	defPosX(true), defPosY(true), resAlgo(restype),
	resFactor(res_Fact), appFactor(app_Fact), direction(dir), direction_orig(dir)
{
	if (dir < -1.0 )
	{
		opserr << "WARNING: BkStressLimSurface2D() - Dir should be between -1 and +1\n";
		opserr << "Set to variable \n";
		direction_orig = 10;
	}

	if(direction_orig > 1)
		direction = 0.0;
		
	kinMatX = kinX.getCopy();
	kinMatY = kinY.getCopy();
	isoMatXPos = isoXPos.getCopy();
	isoMatXNeg = isoXNeg.getCopy();
	isoMatYPos = isoYPos.getCopy();
	isoMatYNeg = isoYNeg.getCopy();
	
	// don't really need a copy for kinMats /LimSurface because only getDrift() is used
	// even if the user reuses.. a separate copy will be created
	// that said ..
	limSurface = lim_surface.getCopy();
	limSurface->setTransformation(0, 1, 1, 1);
}

BkStressLimSurface2D::~BkStressLimSurface2D()
{
  if (kinMatX != 0)
    delete kinMatX;

  if (kinMatY != 0)
    delete kinMatY;

  if (isoMatXPos != 0)
    delete isoMatXPos;

  if (isoMatXNeg != 0)
    delete isoMatXNeg;

  if (isoMatYPos != 0)
    delete isoMatYPos;

  if (isoMatYNeg != 0)
    delete isoMatYNeg;
	
  if (limSurface != 0)
    delete limSurface;
}

void BkStressLimSurface2D::setResidual(double res)
{
	kinMatX->setResidual(res);
	kinMatY->setResidual(res);
}

	        
int BkStressLimSurface2D::commitState()
{
	this->YS_Evolution2D::commitState();

    int res  = kinMatX->commitState();
	res += kinMatY->commitState();
	res += isoMatXPos->commitState();
	res += isoMatXNeg->commitState();
	res += isoMatYPos->commitState();
	res += isoMatYNeg->commitState();
	
	return res;
}

int BkStressLimSurface2D::revertToLastCommit(void)
{
	this->YS_Evolution2D::revertToLastCommit();

	kinMatX->revertToLastCommit();
	kinMatY->revertToLastCommit();
	isoMatXPos->revertToLastCommit();
	isoMatXNeg->revertToLastCommit();
	isoMatYPos->revertToLastCommit();
	isoMatYNeg->revertToLastCommit();

	return 0;
}


void BkStressLimSurface2D::setTrialPlasticStrains(double lamda, const Vector &f, const Vector &g)
{
//	double epx = isotropicRatio*lamda*g(0);
//	double epy = isotropicRatio*lamda*g(1);

	// set wrt absolute for easier calibration
	double epx = lamda*g(0);
	double epy = lamda*g(1);

//	opserr << "epx = " << epx << ", epy = " << epy << endln;
//	opserr << "gx  = " << g(0)  << ", gy  = " << g(1)  << endln;
//	opserr << "\a";

	if(epx > 0)
		defPosX = true;
	else
		defPosX = false;

	if(epy > 0)
		defPosY = true;
	else
		defPosY = false;
		
	isoMatXPos->setTrialIncrValue(epx);
	isoMatXNeg->setTrialIncrValue(-1*epx);
	isoMatYPos->setTrialIncrValue(epy);
	isoMatYNeg->setTrialIncrValue(-1*epy);
	//!! when should sumIsoEp be reset? - using same as plastic hardening

	double x0 = translate_hist(0);
	double y0 = translate_hist(1);
	double fx = f(0);
	double fy = f(1);

	limSurface->hModel->toOriginalCoord(x0, y0);
	double drift = limSurface->getDrift(x0, y0);

	if(direction_orig > 1)
		direction = fabs(y0);

	if(fabs(y0) >= 0.80)
		direction = 1.0; // constP
	
int resType = resAlgo;
	
	double dR = fabs(drift); // in-case outside    pinching starts late

	switch (resType)
	{
	case 1:
	{
		if(drift >= 0.0)
		{
			if(sign(g(0)) != sign(translate_hist(0)))
				dR = 1.5 + drift; // approx value = 2   metal case
			 else
				dR = 0.0;
		}
		else     // no pinching
		{
			//old limSurface->hModel->toOriginalCoord(fx, fy);
			//old  dR = limSurface->interpolate(x0, y0, fx, fy);


			// y0 range -1 to +1
			if(sign(g(0)) != sign(translate_hist(0)))
			{
				  // dR = 2.0 + drift; // drift < 0
				dR = fabs(limSurface->getDrift(0.0, y0))*2 - fabs(drift);
//          not required
//			else
//				dR = fabs(drift);

//			opserr << "!!drift 0, y0 = " << limSurface->getDrift(0.0, y0)
//			     << ", drift = " << drift << endln;
			}
		}

		break;
	} //case 1 - Metals

	case 2:
	{
		if(drift >= 0.0)
			dR = 0.0;
		else     // pinching starts early
		{
//			limSurface->hModel->toOriginalCoord(fx, fy);
//			dR = limSurface->interpolate(x0, y0, fx, fy);
			if(sign(g(0)) != sign(translate_hist(0)))
				dR = fabs(limSurface->getDrift(0.0, y0))*2 - fabs(drift);
		}

		break;
	}//case 2 - Pinching,  Kp =  Kp0 -> 0

	case 3:
	{
		if(drift >= 0.0)
			dR = 0.0;
		break;
	}//case 3 - Pinching, Kp = 0 -> Kp0 -> 0

	case 4:
	{
		if(drift >= 0.0)
		{
			if(sign(g(0)) == sign(translate_hist(0)))
				dR = 0.0;
		}
/*		else
		{
			if(sign(g(0)) != sign(translate_hist(0)))
				dR = 0.0;
		}
*/		
		break;
	}
	
	default:
	{
		opserr << "WARNING - Unknown residual algo\n";
		opserr << *this;
		if(drift >= 0.0)
			dR = 0.0;
	}
	
	} // switch - algo

double sfactor = 1.0;
	
	resHardening = false;
	resApproach  = false;
    if(drift >= 0.0)
    {
		if(sign(g(0)) == sign(translate_hist(0)))
		{
			resHardening = true;
			if(resType > 1)
				sfactor = resFactor; 
		}
		else
		{
			resApproach  = true;
			if(resType > 1)
				sfactor = appFactor;

		}

//		opserr << "----- Drift > 0 --- ( " << sfactor << ")\n";
    }
    
	// absolute values - no need to have history 
	kinMatX->setTrialValue(dR, sfactor);
	kinMatY->setTrialValue(dR, sfactor);
}

const Vector &BkStressLimSurface2D::getEquiPlasticStiffness()
{	
	double kp_kin_x =  kinMatX->getTrialPlasticStiffness();
	double kp_kin_y =  kinMatY->getTrialPlasticStiffness();
	double kp_iso_x =  isoMatXPos->getTrialPlasticStiffness();
	double kp_iso_y =  isoMatYPos->getTrialPlasticStiffness();

	if(!defPosX)
		kp_iso_x =  isoMatXNeg->getTrialPlasticStiffness();
	if(!defPosY)
		kp_iso_y =  isoMatYNeg->getTrialPlasticStiffness();

//	opserr << *isoMatYPos;
//	opserr << *isoMatXPos;
		
	v2(0) =isotropicRatio*kp_iso_x + kinematicRatio*kp_kin_x;
	v2(1) =isotropicRatio*kp_iso_y + kinematicRatio*kp_kin_y;
	
	if(isotropicFactor(0) <=minIsoFactor)  
		v2(0) = 0;

	if(isotropicFactor(1) <=minIsoFactor)
		v2(1) = 0;

	  
	return v2;
}

double BkStressLimSurface2D::getTrialPlasticStrains(int dir)
{
	if(dir == 0 && defPosX)
		return isoMatXPos->getTrialValue();
	else if(dir == 0 && !defPosX)
		return isoMatXNeg->getTrialValue();
	else if (dir == 1 && defPosY)
		return isoMatYPos->getTrialValue();
	else if (dir == 1 && !defPosY)
		return isoMatYNeg->getTrialValue();
	else
		opserr << "BkStressLimSurface2D::getTrialPlasticStrains(double dir) - incorrect dir||condition \n";
	return 0;
}

double BkStressLimSurface2D::getCommitPlasticStrains(int dir)
{
	opserr << "WARNING: BkStressLimSurface2D::getCommitPlasticStrains(.) "
		  << " not yet implemented" << endln;
	return this->getTrialPlasticStrains(dir);
}


double BkStressLimSurface2D::getIsoPlasticStiffness(int dir)
{
	if(dir == 0 && defPosX)
		return isoMatXPos->getTrialPlasticStiffness();
	else if(dir == 0 && !defPosX)
		return isoMatXNeg->getTrialPlasticStiffness();
	else if (dir == 1 && defPosY)
		return isoMatYPos->getTrialPlasticStiffness();
	else if (dir == 1 && !defPosY)
		return isoMatYNeg->getTrialPlasticStiffness();
	else
		opserr << "BkStressLimSurface2D::getIsoPlasticStiffness(double dir) - incorrect dir/condition \n";
	return 0;	
}

double BkStressLimSurface2D::getKinPlasticStiffness(int dir)
{
	if(dir == 0)
		return kinMatX->getTrialPlasticStiffness();
	else if (dir == 1)
		return kinMatY->getTrialPlasticStiffness();
	else
		opserr << "BkStressLimSurface2D::getKinPlasticStiffness(double dir) - incorrect dir\n";
	return 0;

}

Vector& BkStressLimSurface2D::getEvolDirection(Vector &f_new)
{
	// -1 => Radial Evolution
	//  0 => From geometric center (~ normal)
	//  1 => Constant-P
	
	v2(0) = 0.0;
	if(direction >= 0)
		v2(1) = direction*f_new(1);
	else
		v2(1) = direction*translate_init(1);
	return v2;
}

int BkStressLimSurface2D::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
	
	limSurface->displaySelf(theViewer, limSurface->SurfOnly, fact);
	return  0;
}

void BkStressLimSurface2D::Print(OPS_Stream &s, int flag)
{
	s << "BkStressLimSurface2D \n";
	s << "iso_Ratio = " << isotropicRatio << "\n";
	s << "isotropicFactor_hist = " << isotropicFactor_hist;
	s << "translateX       = " << translate(0) << ",\ttranslateY = " << translate(1) << "\n";
	s << "\n";

}
	
