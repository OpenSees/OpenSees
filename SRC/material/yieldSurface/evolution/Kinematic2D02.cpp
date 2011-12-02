//Kinematic2D02
//
//////////////////////////////////////////////////////////////////////

#include "Kinematic2D02.h"
#include <math.h>
//#include "YieldSurface_BC.h"
 
#define evolDebug 0
#define KINEMATIC2D02_CLASSTAG -1

NullPlasticMaterial Kinematic2D02::nullMat(-1);

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Kinematic2D02::Kinematic2D02(int tag, double min_iso_factor, 
									YieldSurface_BC &lim_surface,
									PlasticHardeningMaterial &kpx,
									PlasticHardeningMaterial &kpy,
									int algo, double resfact, double appfact, double dir)
:BkStressLimSurface2D(tag, KINEMATIC2D02_CLASSTAG, min_iso_factor,0.0, 1.0,
                        lim_surface, kpx, kpy,
                        nullMat, nullMat, nullMat, nullMat, algo, resfact, appfact, dir)
{
	
}

Kinematic2D02::~Kinematic2D02()
{

}

YS_Evolution *Kinematic2D02::getCopy()
{
YS_Evolution *theCopy = new Kinematic2D02(getTag(), minIsoFactor,
                                          *limSurface, *kinMatX, *kinMatY, resAlgo,
                                          resFactor, appFactor, direction);

	return theCopy;
}

/*int Kinematic2D02::displaySelf(Renderer &theViewer, int displayMode, float fact)
{

 limSurface->displaySelf(theViewer, limSurface->SurfOnly, fact);
//	limSurface->displaySelf(theViewer, displayMode, fact);
	return  0;
}
*/
void Kinematic2D02::Print(OPS_Stream &s, int flag)
{
	s << "Kinematic2D02 \n";
	s << "iso_Ratio = " << isotropicRatio << "\n";
	s << "isotropicFactor_hist = " << isotropicFactor_hist;
	s << "translateX       = " << translate(0) << ",\ttranslateY = " << translate(1) << "\n";
	s << "\n";

}
	
