// YS_Evolution.cpp: implementation of the YS_Evolution class.
//
//////////////////////////////////////////////////////////////////////

#include "YS_Evolution.h"

Vector YS_Evolution::crd1(1);
Vector YS_Evolution::crd2(2);
Vector YS_Evolution::crd3(3);

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

YS_Evolution::YS_Evolution(int tag, int classtag,
                                     double iso_ratio, double kin_ratio,
                                     int _dimension,
                                     double shr_iso, double shr_kin)
                  :TaggedObject(tag), MovableObject(classtag),
           isotropicRatio(iso_ratio), kinematicRatio(kin_ratio),
           isotropicRatio_orig(iso_ratio), kinematicRatio_orig(kin_ratio),
           dimension(_dimension), deformable(false),
			translate_hist(_dimension), translate(_dimension), translate_init(_dimension),
			isotropicFactor(_dimension), isotropicFactor_hist(_dimension),
			isotropicRatio_shrink(shr_iso), kinematicRatio_shrink(shr_kin),
			freezeEvolution(false)
{
	/*if( (isotropicRatio + kinematicRatio) != 1.0 || (isotropicRatio + kinematicRatio) != 0.0)
	{
	 	opserr << "WARNING YS_Evolution(..) -  isotropicRatio + kinematicRatio != 1 or 0 -> "
		  << isotropicRatio + kinematicRatio << endln;
		// this always gives warning message
	}*/
	translate.Zero();
	translate_hist.Zero();
	translate_init.Zero();

//	opserr << "DIM = " << dimension << endln;
	for(int i=0; i<dimension;i++)
	{
		isotropicFactor(i)=1;
		isotropicFactor_hist(i)=1;
	}
}

YS_Evolution::~YS_Evolution()
{
}

void YS_Evolution::setInitTranslation(Vector &initTranslate)
{
	if(initTranslate.Size() > dimension)
	{
	 	opserr << "WARNING -  newTranslate" << initTranslate << " outside the dimensions\n";
	}

	translate = initTranslate;
	translate_hist = initTranslate;
	translate_init = initTranslate;
}

const Vector &YS_Evolution::getInitTranslation(void)
{
	return translate_init;
}


int YS_Evolution::update(int flag)
{
	// does nothing
	return 0;
}

void YS_Evolution::setResidual(double res)
{
	// does nothing
}

	
int YS_Evolution::commitState()
{
	isotropicFactor_hist   = isotropicFactor;
	translate_hist         = translate;
	return 0;
}

int YS_Evolution::revertToLastCommit(void)
{
	isotropicFactor   =  isotropicFactor_hist;
	translate         =  translate_hist;
	return 0;
}

void  YS_Evolution::toDeformedCoord(Vector &coord)
{
//	if(getTag() == 10)
//	{
//		opserr << "DIM: " << dimension << endln;
//		opserr << "ISO: " << isotropicFactor;
//		opserr << "\a";
//		}
	
    // no vector multipication in opensees
    for(int i=0; i< coord.Size(); i++)
    {
		coord(i) = coord(i)*isotropicFactor(i) + translate(i);
	}

}
	
void  YS_Evolution::toOriginalCoord(Vector &coord)
{
    for(int i=0; i< coord.Size(); i++)
    {
		double notrans = coord(i) - translate(i);
		coord(i) = notrans/isotropicFactor(i);
	}

}

void YS_Evolution::toDeformedCoord(double &x)
{
//	opserr << "WARNING YS_Evolution::toDeformedCoord(.) - This method should not be called\n";
	crd1(0) = x;
	this->toDeformedCoord(crd1);
	x = crd1(0);
}

void YS_Evolution::toDeformedCoord(double &x, double &y)
{
//  opserr << "WARNING YS_Evolution::toDeformedCoord(..) - This method should not be called\n";
	crd2(0) = x;
	crd2(1) = y;
	this->toDeformedCoord(crd2);
	x = crd2(0);
	y = crd2(1);
}

void YS_Evolution::toDeformedCoord(double &x, double &y, double &z)
{
//	opserr << "WARNING YS_Evolution::toDeformedCoord(...) - This method should not be called\n";
	crd3(0) = x;
	crd3(1) = y;
	crd3(2) = z;
	this->toDeformedCoord(crd3);
	x = crd3(0);
	y = crd3(1);
	z = crd3(2);

}

void YS_Evolution::toOriginalCoord(double &x)
{
//	opserr << "WARNING YS_Evolution::toOriginalCoord(.) - This method should not be called\n";
	crd1(0) = x;
	this->toOriginalCoord(crd1);
	x = crd1(0);

}

void YS_Evolution::toOriginalCoord(double &x, double &y)
{
//	opserr << "WARNING YS_Evolution::toOriginalCoord(..) - This method should not be called\n";
	crd2(0) = x;
	crd2(1) = y;
	this->toOriginalCoord(crd2);
	x = crd2(0);
	y = crd2(1);
}

void YS_Evolution::toOriginalCoord(double &x, double &y, double &z)
{
//	opserr << "WARNING YS_Evolution::toOriginalCoord(...) - This method should not be called\n";
	crd3(0) = x;
	crd3(1) = y;
	crd3(2) = z;
	this->toOriginalCoord(crd3);
	x = crd3(0);
	y = crd3(1);
	z = crd3(2);

}

double YS_Evolution::getCommitTranslation(int dir)
{
#ifdef _G3DEBUG
	checkDimension(dir);
#endif

	return translate_hist(dir);
}

double YS_Evolution::getTrialTranslation(int dir)
{
#ifdef _G3DEBUG
	checkDimension(dir);
#endif

	return translate(dir);
}

double YS_Evolution::getCommitIsotropicFactor(int dir)
{
	return isotropicFactor_hist(dir);
}

double YS_Evolution::getTrialIsotropicFactor(int dir)
{
	return isotropicFactor(dir);
}

void YS_Evolution::checkDimension(int dir)
{
	if(dir < 0 || dir >= dimension)
	{
	 	opserr << "WARNING - Direction " << dir << " outside the dimensions\n";
	}

}

void  YS_Evolution::setDeformable(bool defo)
{
	this->deformable = defo;
}

void YS_Evolution::Print(OPS_Stream &s, int flag)
{
	opserr << " YS_Evolution - tag = " << getTag() << endln;
}

/*OPS_Stream &operator<<(OPS_Stream &s, const YS_Evolution &hModel)
{
	hModel.Print(s);
  return s;
}*/

	
	
	
