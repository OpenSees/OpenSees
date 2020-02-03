// YS_Evolution2D.cpp: implementation of the YS_HardeningModel class.
//
//////////////////////////////////////////////////////////////////////

#include "YS_Evolution2D.h"
#include <YieldSurface_BC.h>
#include <math.h>
#include <MaterialResponse.h>

#define evolDebug  0
#define modifDebug 0
#define transDebug 0

Vector YS_Evolution2D::v2(2);

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

YS_Evolution2D::YS_Evolution2D(int tag, int classtag,
							   double min_iso_factor,
                               double iso_ratio, double kin_ratio)
  :YS_Evolution(tag, classtag, iso_ratio, kin_ratio, 2),
   minIsoFactor(min_iso_factor), softening(false)
{
//	sumPlasticDeformX = 0;
//	sumPlasticDeformX_hist = 0;
//	sumPlasticDeformY = 0;
//	sumPlasticDeformY_hist = 0;
	isotropicFactor(0)   = 1;
	isotropicFactor(1)   = 1;
    
	isotropicFactor_hist(0) = 1;
	isotropicFactor_hist(1) = 1;

}

YS_Evolution2D::~YS_Evolution2D()
{

}

int YS_Evolution2D::commitState()
{
	this->YS_Evolution::commitState();
	
//	sumPlasticDeformX_hist = sumPlasticDeformX;
//	sumPlasticDeformY_hist = sumPlasticDeformY;
	freezeEvolution = false;
	return 0;
}

int YS_Evolution2D::update(int flag)
{	
	isotropicFactor   =  isotropicFactor_hist;
	translate     =  translate_hist;
	
	return 0;
}

int YS_Evolution2D::revertToLastCommit(void)
{
	this->YS_Evolution::revertToLastCommit();
//	sumPlasticDeformX =  sumPlasticDeformX_hist;
//	sumPlasticDeformY =  sumPlasticDeformY_hist;
	return 0;
}

/*
void YS_Evolution2D::toDeformedCoord(double &x, double &y)
{
double x1 = x, y1 = y;
double  translateX = translate(0);
double  translateY = translate(1);

	if(transDebug)
	{
		opserr << "\nYieldSurfaceBC2D::toDeformedCoord(double &x, double &y)\n";
		opserr << "Original - fx = " << x1 << ",\tfy = " << y1 << "\n";
	}

	// isotropic
    x1 = x1*isotropicFactor(0);
    y1 = y1*isotropicFactor(1);

	// kinmatic
    x1 += translateX;
    y1 += translateY;

	x = x1;
	y = y1;

	if(transDebug)
	{
		opserr << "Deformed - fx = " << x << ",\tfy = " << y << "\n";
		opserr << "isotropicFactor =  " <<  isotropicFactor << "\n";
		opserr << "translateX       = " << translateX << ",\ttranslateY = " << translateY << "\n";
	}
}
*/

/*
void YS_Evolution2D::toOriginalCoord(double &x, double &y)
{
double  translateX = translate(0);
double  translateY = translate(1);

  	if(transDebug)
	{
		opserr << "\nYieldSurfaceBC2D::toOriginalCoord(double &x, double &y)\n";
		opserr << "Deformed - fx = " << x << ",\tfy = " << y << "\n";
		opserr << "isotropicFactor =  " <<  isotropicFactor << "\n";
		opserr << "translateX       = " << translateX << ",\ttranslateY = " << translateY << "\n";
	}

	// kinematic
    x -= translateX;
    y -= translateY;

	// isotropic
	x = x/isotropicFactor(0);
    y = y/isotropicFactor(1);

	if(transDebug)
	{
	 	opserr << "Original - fx = " << x << ",\tfy = " << y << "\n";
		opserr << "\a";
	}
}
*/

int YS_Evolution2D::evolveSurface(YieldSurface_BC *ys, double lamda, 
                                  Vector &G, Vector &F_Surface, int flag)
{
	// first and fore-most:
	tmpYSPtr = ys;
//	int loc = ys->ele_Location;
//	opserr << " evolve surface [" << opserr << loc << "]\n";
//	opserr << *ys;

	//freezeEvolution = false; -> have set this in commitState -> don't change
	// first save the values on stack
	// static vectors could get reallocated elsewhere
	Vector f_sur(2);
		f_sur(0) = F_Surface(0);
		f_sur(1) = F_Surface(1);
	Vector gl(2);
		gl(0) = G(0);
		gl(1) = G(1);
	
	setTrialPlasticStrains(lamda, f_sur, gl);
	if(freezeEvolution)
		return 0;

	//deformable = true;
	
	double kinX = gl(0)*getKinPlasticStiffness(0)/ys->getCap(0);
	double kinY = gl(1)*getKinPlasticStiffness(1)/ys->getCap(1);
	double isoX = gl(0)*getIsoPlasticStiffness(0)/ys->getCap(0);
	double isoY = gl(1)*getIsoPlasticStiffness(1)/ys->getCap(1);

	// opserr << "isoX = " << isoX << ", isoY = " << isoY << endln;
	
	// kinematic hardening
	double lamda_kin = kinematicRatio*lamda;	
	double dfx_kin = lamda_kin*kinX;
	double dfy_kin = lamda_kin*kinY;

	// isotropic hardening
	double lamda_iso = isotropicRatio*lamda;
	double dfx_iso = lamda_iso*isoX;
	double dfy_iso = lamda_iso*isoY;

	double dfx_total =  dfx_kin + dfx_iso;
	double dfy_total =  dfy_kin + dfy_iso;
	
	double fx_new = f_sur(0) + dfx_total;
	double fy_new = f_sur(1) + dfy_total;

	double fx_iso_new = f_sur(0) + dfx_iso;
	double fy_iso_new = f_sur(1) + dfy_iso;

    // check 1: for cross-over -> same as: |F_sur + df| > |F_sur| && Kp < 0
	// sign change

	//cout << "f_sur = " << f_sur(0) << ", f new = " << fx_new << endln;
	bool ys_harden = true;
	toOriginalCoord(fx_new, fy_new);
	if(ys->getDrift(fx_new, fy_new) < 0)
		ys_harden = false;

	bool iso_harden = true;
	toOriginalCoord(fx_iso_new, fy_iso_new);
	if(ys->getDrift(fx_iso_new, fy_iso_new) < 0)
		iso_harden = false;

	
 	if(!ys_harden && (sign(f_sur(0)) != sign(fx_new))) // risky to use for fy -> P
 	{
	  // need to fix this
 		dfx_iso = 0.0;
 		dfx_kin = 0.0;
	    opserr << "Condition happened..\n";
	    opserr << *ys;	    
		freezeEvolution = true;
	    
	    //cerr << "F_Surface(0) = " << f_sur(0) << ", F_New = " << fx_new << endln;
	    //cin.get();

		// if(!deformable) //nothing to do
		// return anyway for now
		{

			return 0;
//			dfy_iso = 0.0;
//			dfy_kin = 0.0;
		}
	}
		
	if(!ys_harden && (kinematicRatio != kinematicRatio_shrink)
				  && (isotropicRatio != isotropicRatio_shrink)  )
	{
		double kinRatio = kinematicRatio_shrink;
		double isoRatio = isotropicRatio_shrink;
		// here it might be a good idea to redo the above step
		// will not make difference for peak-oriented but for
		// others may cause convergence problems  (Kp_iso != Kp_kin)
		// what if its not softening anymore?
		lamda_iso = isoRatio*lamda;
		dfx_iso = lamda_iso*isoX;
		dfy_iso = lamda_iso*isoY;
		lamda_kin = kinRatio*lamda;
		dfx_kin = lamda_kin*kinX;
		dfy_kin = lamda_kin*kinY;

		dfx_total =  dfx_kin + dfx_iso;
		dfy_total =  dfy_kin + dfy_iso;
		fx_new = f_sur(0) + dfx_total;
		fy_new = f_sur(1) + dfy_total;

		toOriginalCoord(fx_new, fy_new);
		if(ys->getDrift(fx_new, fy_new) > 0)
		{
			opserr << "oops: YS_Evolution2D::evolveSurface() - softens->hardens\n";
			ys_harden = true;
		}
	}


   	// Update the isotropicFactor vector
	/*
	// This way does not work!
	int x_grow = 1;
	if(fabs(f_sur(0) + dfx_iso) < fabs(f_sur(0)))
	{
		// opserr << "Softening!\n";
		x_grow = -1;
	}
	
	int y_grow = 1;
	if(fabs(f_sur(1) + dfy_iso) < fabs(f_sur(1)))
		y_grow = -1;
	*/

	int x_grow = 1, y_grow = 1;
	if(getIsoPlasticStiffness(0) < 0)
		x_grow = -1;
	if(getIsoPlasticStiffness(1) < 0)
		y_grow = -1;
	
	if(evolDebug)
	{
		opserr << "F_Surface Vector: " << f_sur;
		opserr << "Fx_new = " << fx_new << ", Fy_new = " << fy_new << endln;
		opserr << "Gradient: " << gl;
	    opserr << "KpxIso = " << getIsoPlasticStiffness(0) << ", KpyIso = " <<  getIsoPlasticStiffness(1) << endln;
	    opserr << "F_surf(1) = " << f_sur(1) << ", dfy_iso = " << dfy_iso << endln;
	    opserr << "x_grow = " <<  x_grow << ", y_grow = " <<  y_grow << endln;
	    opserr << "---------------------------------------------------------" << endln;
     }

	Vector mgnf(2);
	mgnf = isotropicFactor_hist;
	if(flag==1)
		mgnf = isotropicFactor;
	Vector delMag(2);
	
	if(deformable)
	{

		delMag(0) = x_grow*fabs(dfx_iso);
		delMag(1) = y_grow*fabs(dfy_iso);
	}
	else
	{
		double dR = sqrt(dfx_iso*dfx_iso + dfy_iso*dfy_iso);
		if(!iso_harden)
		dR = -1*dR;
		delMag(0) = dR;
		delMag(1) = dR;
	}

	Vector isoFact = mgnf + delMag;

	//check 2: For min isotropic factor
	 if( (isotropicFactor(0) + delMag(0)) <= minIsoFactor)
	{
		delMag(0) = 0.0;
		dfx_kin = 0.0;
        freezeEvolution = true;
		if(!deformable)// nothing to do
			return 0;
	}
	if( (isotropicFactor(1) + delMag(1)) <= minIsoFactor)
	{
		delMag(1) = 0.0;
		dfy_kin = 0.0;
		freezeEvolution = true;

		if(!deformable)
			return 0;
	}

	// update the translation vector
	double fx_aim = f_sur(0) + dfx_kin;
	double fy_aim = f_sur(1) + dfy_kin;

    //cout << "YS_Evolution2D - F_Surface = " << F_Surface;

	toOriginalCoord(fx_aim, fy_aim);
	Vector f_aim(2);
	f_aim(0) = fx_aim;
	f_aim(1) = fy_aim;
	v2 = getEvolDirection(f_aim);

	Vector df_kin = ys->translationTo(f_aim, v2);
	// correct for isotropic factor
	Vector trans(2);
	trans = translate_hist;
	if(flag==1)
		trans = translate;

    // Update the quantities
	translate(0) = trans(0) + df_kin(0)*isotropicFactor(0);
	translate(1) = trans(1) + df_kin(1)*isotropicFactor(1);
	isotropicFactor = mgnf + delMag;

	return 0;
}

int YS_Evolution2D::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
	return -1;
}
	
Response* YS_Evolution2D::setResponse(char **argv, int argc, Information &matInfo)
{
	return 0;
/*
UniaxialMaterial *mat;
int id0 = 0;

	if(strcmp(argv[1],"materialX") == 0)
	{
		mat = kpMatX;
		id0 = 0;
	}
	else if (strcmp(argv[1],"materialY") == 0)
	{
		mat = kpMatY;
		id0 = 3;
	}
	else
		return 0; // for now


	// responseIds from 1 - 3 for matX, 4 - 6 for matY
    // stress
    if (strcmp(argv[0],"stress") == 0)
		return new MaterialResponse(mat, id0+1, mat->getStress());

    // tangent
    else if (strcmp(argv[0],"tangent") == 0)
		return new MaterialResponse(mat, id0+2, mat->getTangent());

    // strain
	else if (strcmp(argv[0],"strain") == 0)
		return new MaterialResponse(mat, id0+3, mat->getStrain());

    // otherwise unknown
    else
		return 0;
	*/	
}

int YS_Evolution2D::getResponse(int responseID, Information &matInfo)
{
	return 0;
  /*
  // each subclass must implement its own stuff
  switch (responseID) {
    case 1:
      matInfo.setDouble(kpMatX->getStress());
      break;

    case 2:
      matInfo.setDouble(kpMatX->getTangent());
      break;

    case 3:
      matInfo.setDouble(kpMatX->getStrain());
      break;

    case 4:
      matInfo.setDouble(kpMatY->getStress());
      break;

    case 5:
      matInfo.setDouble(kpMatY->getTangent());
      break;

    case 6:
      matInfo.setDouble(kpMatY->getStrain());
      break;

    default:
      return -1;
  }

  return 0;
  */
}


