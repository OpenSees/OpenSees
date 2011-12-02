//
//
//////////////////////////////////////////////////////////////////////

#include "PlasticHardening2D.h"
#include <YieldSurface_BC.h>
#include <math.h>

#define strnDebug 0
#define stifDebug 0
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

PlasticHardening2D::PlasticHardening2D(int tag, int classTag, double min_iso_factor,
				double iso_ratio, double kin_ratio,
				 PlasticHardeningMaterial &kpx_pos,
				 PlasticHardeningMaterial &kpx_neg,
				 PlasticHardeningMaterial &kpy_pos,
				 PlasticHardeningMaterial &kpy_neg, double dir)
:YS_Evolution2D(tag, classTag, min_iso_factor, iso_ratio, kin_ratio),
direction(dir), defPosX(true), defPosY(true)
{
	if (dir < -1.0 || dir > 1.0)
	{
		opserr << "WARNING: PlasticHardening2D() - Dir should be between -1 and +1\n";
		opserr << "Set to 0 \n";
		direction = 0.0;
	}

	kpMatXPos = kpx_pos.getCopy();
	kpMatXNeg = kpx_neg.getCopy();
	
	kpMatYPos = kpy_pos.getCopy();
	kpMatYNeg = kpy_neg.getCopy();
}

PlasticHardening2D::~PlasticHardening2D()
{
  if (kpMatXPos != 0)
    delete kpMatXPos;

  if (kpMatXNeg != 0)
    delete kpMatXNeg;

  if (kpMatYPos != 0)
    delete kpMatYPos;

  if (kpMatYNeg != 0)
    delete kpMatYNeg;
}

int PlasticHardening2D::commitState()
{
	this->YS_Evolution2D::commitState();
	
    int res  = kpMatXPos->commitState();
		res += kpMatXNeg->commitState();
		res += kpMatYPos->commitState();
		res += kpMatYNeg->commitState();

	/*if(getTag() ==1)
	{
	if(stifDebug)
	{
		double v0 = getIsoPlasticStiffness(0);
		double v1 = getIsoPlasticStiffness(1);
		opserr << v0 << "\t " << v1 << endln;
	}
    }*/
    
	return res;
}

int PlasticHardening2D::revertToLastCommit(void)
{
	this->YS_Evolution2D::revertToLastCommit();
	
	kpMatXPos->revertToLastCommit();
	kpMatXNeg->revertToLastCommit();
	kpMatYPos->revertToLastCommit();
	kpMatYNeg->revertToLastCommit();

	return 0;
}


// In a plastic hardening material, Kp_iso == Kp_kin == Kp_equivalent
void PlasticHardening2D::setTrialPlasticStrains(double lamda, const Vector &f, const Vector &g)
{

//	opserr << *tmpYSPtr;
	
	double epx = lamda*g(0);
	double epy = lamda*g(1);
//	double val = f(0);
//	double chk = 0.8;
//	val = translate_hist(0);

	defPosX = true;
	if(epx < 0)
		defPosX = false;

//	kpMatXPos->setTrialIncrValue(epx);
//	kpMatXNeg->setTrialIncrValue(-1*epx);
// no need of if.. else - just for remembering the condition

//	pinchX = false;
	if(defPosX)
	{
//		if(val < -1*chk)
//		{
//			opserr << "+Pinch [";
//			opserr << tmpYSPtr->ele_Location << "]\n";
//			pinchX = true;
//		}
//		else
			kpMatXPos->setTrialIncrValue(epx);

		kpMatXNeg->setTrialIncrValue(-1*epx);
	}
	else
	{
//		if(val >  chk)
//		{
//			opserr << "-Pinch [";
//			opserr << tmpYSPtr->ele_Location << "]\n";
//			pinchX = true;
//		}
//		else
			kpMatXNeg->setTrialIncrValue(-1*epx);

		kpMatXPos->setTrialIncrValue(epx);
	}

	defPosY = true;
	if(epy < 0)
		defPosY = false;

    if(defPosY)
    {
	//	 if(translate_hist(1) >= 0)
			kpMatYPos->setTrialIncrValue(epy);

		kpMatYNeg->setTrialIncrValue(-1*epy);
	}
	else
	{
	//	 if(translate_hist(1) <= 0)
			kpMatYNeg->setTrialIncrValue(-1*epy);

		kpMatYPos->setTrialIncrValue(epy);
	}

	if(strnDebug)
	{
		opserr << "epx = " << epx << ", epy = " << epy << endln;
		opserr << "bool defPosX = " << defPosX << ", bool defPosY = " << defPosY << endln;
	}

}

const Vector &PlasticHardening2D::getEquiPlasticStiffness()
{
	if(freezeEvolution)
	{
		v2(0) = 0.0;
		v2(1) = 0.0;
		return v2;
	}

	if(defPosX == true)
	  v2(0) =  kpMatXPos->getTrialPlasticStiffness();
	 else
	  v2(0) =  kpMatXNeg->getTrialPlasticStiffness();

//	if(pinchX)
//		v2(0) = 0.123*v2(0);

	if(defPosY == true)
	  v2(1) =  kpMatYPos->getTrialPlasticStiffness();
	else
	  v2(1) =  kpMatYNeg->getTrialPlasticStiffness();


	if(strnDebug)
		opserr << "Kp " << v2;

	return v2;
}

double PlasticHardening2D::getTrialPlasticStrains(int dir)
{
	if(dir == 0 && defPosX)
		return kpMatXPos->getTrialValue();
	else if(dir == 0 && !defPosX)
		return kpMatXNeg->getTrialValue();
	else if (dir == 1 && defPosY)
		return kpMatYPos->getTrialValue();
	else if (dir == 1 && !defPosY)
		return kpMatYNeg->getTrialValue();
	else
		opserr << "PlasticHardening2D::getTrialPlasticStrains(double dir) - incorrect dir||condition \n";
	return 0;
}

double PlasticHardening2D::getCommitPlasticStrains(int dir)
{
	opserr << "PlasticHardening2D::getCommitPlasticStrains(double dir) - not yet implemented \n";
	this->getTrialPlasticStrains(dir);
	return 0;
}

double PlasticHardening2D::getIsoPlasticStiffness(int dir)
{
double kp =0;

	if(dir == 0)
	{
		if(defPosX)
		  kp = kpMatXPos->getTrialPlasticStiffness();
		else
		  kp = kpMatXNeg->getTrialPlasticStiffness();

//		if(pinchX)
//			kp = 0.123*kp;
	}
	else if (dir == 1)
	{
		if(defPosY)
		  kp =  kpMatYPos->getTrialPlasticStiffness();
		 else
		  kp = kpMatYNeg->getTrialPlasticStiffness();
	}
	else
		opserr << "WARNING: PlasticHardening2D::getPlasticStiffness(int dir) - incorrect dir\n";
	return kp;
}

double PlasticHardening2D::getKinPlasticStiffness(int dir)
{
	return this->getIsoPlasticStiffness(dir);
}

Vector& PlasticHardening2D::getEvolDirection(Vector &f_new)
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

void PlasticHardening2D::Print(OPS_Stream &s, int flag)
{
	s << "PlasticHardening2D \n";
	s << "iso_Ratio = " << isotropicRatio << "\n";
	s << "isotropicFactor_hist = " << isotropicFactor_hist;
	s << "translateX       = " << translate(0) << ",\ttranslateY = " << translate(1) << "\n";
	s << "\n";
}

