// YieldSurface_BC2D.cpp: implementation of the YieldSurfaceBC class.
//
//////////////////////////////////////////////////////////////////////

#include "YieldSurface_BC2D.h"
#define modifDebug  0
#define returnDebug 0
#define transDebug  0
#define driftDebug  0

#include <MaterialResponse.h>

double YieldSurface_BC2D::error(1.0e-6);
Vector YieldSurface_BC2D::v6(6);
Vector YieldSurface_BC2D::T2(2);
Vector YieldSurface_BC2D::F2(2);
Vector YieldSurface_BC2D::g2(2);
Vector YieldSurface_BC2D::v2(2);
Vector YieldSurface_BC2D::v4(4);


//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

YieldSurface_BC2D::YieldSurface_BC2D(int tag, int classTag, double xmax, double ymax,
                  YS_Evolution &model)
:YieldSurface_BC(tag, classTag, model, xmax, ymax)
{
	status_hist = -1; //elastic

	fx_hist = 0;
	fy_hist = 0;

	offset     = 0.01;
	increment  = 0.022;

	state = 0;

}

YieldSurface_BC2D::~YieldSurface_BC2D()
{

}

const Vector &YieldSurface_BC2D::getExtent(void)
{
	v4(0) = xPos;
	v4(1) = xNeg;
	v4(2) = yPos;
	v4(3) = yNeg;

	return v4;
}


//////////////////////////////////////////////////////////////////////
// Transformation
//////////////////////////////////////////////////////////////////////

// Use this function to initialize other stuff
void YieldSurface_BC2D::setTransformation(int xDof, int yDof, int xFact, int yFact)
{
	this->YieldSurface_BC::setTransformation(xDof, yDof, xFact, yFact);

	this->setExtent();
	if(xPos == 0 && yPos == 0 && xNeg ==0 && yNeg == 0)
	{
	 	opserr << "WARNING - YieldSurface_BC2D - surface extent not set correctly\n";
	}

	if(xPos == 0 || xNeg == 0)
		opserr << "Error - YieldSurface_BC2D no X extent\n";

	////////////////////////////////////////////////////////
	// Next set the 'a' and 'b' for the internal quad
	////////////////////////////////////////////////////////

double x1, y1, x2, y2;

	// QUAD 1
	x1 = 0; y1 = yPos - offset ; x2 = xPos - offset; y2 = 0;
	a1 = (y1 - y2)/(x1 - x2);
	b1 = y1 - a1*x1;

	// QUAD 2
	x1 = 0; y1 = yPos - offset; x2 = xNeg + offset; y2 = 0;
	a2 = (y1 - y2)/(x1 - x2);
	b2 = y1 - a2*x1;

	// QUAD 3
	x1 = 0; y1 = yNeg + offset; x2 = xNeg + offset; y2 = 0;
	a3 = (y1 - y2)/(x1 - x2);
	b3 = y1 - a3*x1;

	// QUAD 4
	x1 = 0; y1 = yNeg + offset; x2 = xPos - offset; y2 = 0;
	a4 = (y1 - y2)/(x1 - x2);
	b4 = y1 - a4*x1;

}


//////////////////////////////////////////////////////////////////////
// Overridden methods
//////////////////////////////////////////////////////////////////////

int YieldSurface_BC2D::getState(int stateInfo)
{
	if(stateInfo == this->StateLoading)
		return isLoading;
	else
		return -1;
}

int YieldSurface_BC2D::commitState(Vector &force)
{
	this->YieldSurface_BC::commitState(force);
	
	status_hist = this->getTrialForceLocation(force);
//    opserr << "YieldSurface_BC2D::commitState(..), status = " << status_hist << endln;

//	opserr << "YieldSurface_BC2D::commitState " << getTag() <<" - force location: " << status_hist << endln;
//	opserr << " FORCE = " << force;
	if(status_hist > 0)
	{
		opserr << "WARNING - YieldSurface_BC2D::commitState(..) [" << getTag()<<"]\n";
		opserr << "Can't commit with force outside the surface\n";
		opserr << "\a";
	}

	double driftOld = this->getDrift(fx_hist, fy_hist);
	double driftNew = this->getTrialDrift(force);
	isLoading = 0;
	if(status_hist >= 0)
	{
		isLoading = 1;
	}
	else
	{
		if(driftOld < driftNew) // drift will be negative
			isLoading = 1;
	}
		
	//hModel->commitState(status_hist);
	hModel->commitState();

	toLocalSystem(force, fx_hist, fy_hist, true);
	hModel->toOriginalCoord(fx_hist, fy_hist);
	
	// opserr << "yPos = " << yPos << ", fy = " << fy_hist << endln;
	if(fy_hist/yPos > 0.85)
		hModel->setDeformable(true);
	else
		hModel->setDeformable(false); 

	gx_hist = 0;
	gy_hist = 0;

	if(status_hist==0)
	{
		// double fx = fx_hist;
		// double fy = fy_hist;

		// hModel->toOriginalCoord(fx, fy);
		getGradient(gx_hist, gy_hist, fx_hist, fy_hist);
	}

	return 0;
}


int YieldSurface_BC2D::update(int flag)
{
	return hModel->update(flag);
}


int YieldSurface_BC2D::revertToLastCommit(void)
{
	hModel->revertToLastCommit();
	return 0;
}

// assuming in original coordinates, return value does not account
// for any isotropic factor
Vector& YieldSurface_BC2D::translationTo(Vector &f_new, Vector &f_dir)
{
	double x1 = f_dir(0);
	double y1 = f_dir(1);

	double x2 = f_new(0);
	double y2 = f_new(1);

	// first find the point to interpolate to
	bool is_hardening = true;
	state = 1;
	
	double hi = getDrift(x2, y2);
		if(hi < 0)
		{
			is_hardening = false;
			state = -1;
		}
	if(fabs(hi) < 1e-12) state = 0;

//	opserr << "Drift New = " << hi <<endl;
//	opserr << "F new = " << f_new;
	
	hi = 5*fabs(hi); //approx half to interpolate - that didn't work for all cases
	// changed from factor of 2 to 5
	// radial/normal evolution -> no problem
	// for bounding surface, dir. of evolution with hardening may cause the
	// interpolation end point (xi, yi) to be out side

	double h = sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1));
	double c = hi/h;
	double sign= -1.0;

    if(c > 1.0)
    {
		opserr << "oops - YieldSurface_BC2D::translationTo - c > 1.0 \n";
		c = 1.0;
	}

	if(!is_hardening)
    	sign = 1.0;

//	double xi = x2 + sign*c*fabs(x2 - x1);
//	double yi = y2 + sign*c*fabs(y2 - y1); 
	double xi = x2 + sign*c*(x2 - x1);
	double yi = y2 + sign*c*(y2 - y1);
        
    if(is_hardening) // Fi inside, Fnew outside
    {
        double t = interpolate(xi, yi, x2, y2);
        T2(0) = (1 - t)*(x2 - xi);
        T2(1) = (1 - t)*(y2 - yi);
    }
    else // Fnew inside, Fi outside
    {
        double t = interpolate(x2, y2, xi, yi);
        T2(0) = t*(x2 - xi);
        T2(1) = t*(y2 - yi);
    }

//	opserr << "F New = " << f_new;
//	opserr << "F dir = " << f_dir;
//	opserr << "Translation vector = " << v2;
	
	return T2;
}

// magPlasticDefo is based on incremental force
// we always modify the surface only if the last conv force was on the
// surface, otherwise we'll just set to surface if it shoots through

// return value >  0 --> surface is expanding
// return value == 0 --> surface unchanged
// return value <  0 --> surface is shrinking
int YieldSurface_BC2D::modifySurface(double magPlasticDefo, Vector &Fsurface, Matrix &G, int flag)
{

// check 1
	if( this->getTrialForceLocation(Fsurface) !=0)
	{
		opserr << "Can't modify surface with Force Location = " << getTrialForceLocation(Fsurface) << endln;
		return 0;
	}


// check 2
	if(magPlasticDefo < 0)
	{
		opserr << "\nYieldSurface_BC2D::modifySurface(..) \n";
		opserr << "Warning -   magPlasticDefo < 0 " << magPlasticDefo << "\n";
		// magPlasticDefo = 0;
		return 0;
	}


	double fx, fy, fx_def, fy_def, gx, gy;

	toLocalSystem(Fsurface, fx_def, fy_def, true);
	// set to ys coordinates
	toLocalSystem(G, gx, gy, false, true);

	F2(0) = fx_def;
	F2(1) = fy_def;
	g2(0) = gx;
	g2(1) = gy;

	hModel->evolveSurface(this, magPlasticDefo, g2, F2, flag);
	
	/*double  isotropicFactor      = hModel->getTrialIsotropicFactor(0);
	double  isotropicFactor_hist = hModel->getCommitIsotropicFactor(0);

	int res;
	if( fabs (fabs(isotropicFactor) - fabs(isotropicFactor_hist)) < 1e-13) // remains same
	{
		res = 0;
	}
	else if(isotropicFactor > isotropicFactor_hist) // surface is expanding
	{
	 	res = 1;
	}
	else if(isotropicFactor < isotropicFactor_hist)
	{
		res = -1;
	}
	else
		opserr << "ModifySurface - this condition should not happen\n";
    */
	return state;

}


void YieldSurface_BC2D::getCommitGradient(Matrix &G)
{
double gx = gx_hist;
double gy = gy_hist;

	toElementSystem(G, gx, gy, false, true);
}

void YieldSurface_BC2D::getTrialGradient(Matrix &G, Vector &force)
{
double gx, gy, fx, fy;



	toLocalSystem(force, fx, fy, true);
	hModel->toOriginalCoord(fx, fy);

//	double drift = getDrift(fx, fy);
//	opserr << "YieldSurface_BC2D::getTrialGradient  trial_drift: " << this->getTrialDrift(force)
//		 << "  drift: " << drift << endln;
//	opserr << "----> calling getGradient " << endln;

    getGradient(gx, gy, fx, fy);
//	opserr << "<---- done calling getGradient " << endln;

	toElementSystem(G, gx, gy, false, true);
}


double YieldSurface_BC2D::getTrialDrift(Vector &force)
{
double fx, fy;
	toLocalSystem(force, fx, fy, true);
	hModel->toOriginalCoord(fx, fy);
double drift = getDrift(fx, fy);

    return drift;
}


// used by the element
int YieldSurface_BC2D::getTrialForceLocation(Vector &force)
{
double drift = this->getTrialDrift(force);
		return this->forceLocation(drift);
}

int YieldSurface_BC2D::getCommitForceLocation()
{
	return status_hist;
}

// used by the yield surface
// here we decide on what we define as outside, inside
// and on the surface
// -1 -> inside
//  0 -> within
// +1 -> outside
int  YieldSurface_BC2D::forceLocation(double drift)
{
double tolNeg = 0.00;
double tolPos = 1e-5;

int	   status = -2;

	// first get rid of -1e-13, etc
	if(fabs(drift) < 1e-7)
    	drift = 0;

	if(drift <  -tolNeg)
	{
		status = -1; //inside
	}
	else if(drift >= -tolNeg && drift <= tolPos) // on the surface
    {
     	status = 0;
    }
	else if(drift > tolPos) // outside the surface
	{
	 	status = 1;
	}
	else
	{
	 	opserr << "YieldSurface_BC2D::forceLocation(double drift) - this condition not possible\n";
	 	opserr << "\a";
	}

	return status;
}

void YieldSurface_BC2D::addPlasticStiffness(Matrix &K)
{
Vector v2 = hModel->getEquiPlasticStiffness();

       v6.Zero();
double kpX =  v2(0);
double kpY =  v2(1);

//void toElementSystem(Vector &eleVector, double &x, double &y);
      toElementSystem(v6, kpX, kpY, false, false);

      for(int i=0; i<6; i++)
      {
        K(i,i) += v6(i);
      }

}


double YieldSurface_BC2D::getDrift(double x, double y)
{

double sdrift = getSurfaceDrift(x, y);

// if sdrift > 0, trust it to be outside but if
// sdrift < 0, the point may be either inside
// or outside.  Reason is that surfaces close
// on themselves nicely, but may deform or
// open up while expanding.

//	if(sdrift > 0)
//		return sdrift;

// Now sdrift < 0, if it is inside the internal
// quad, sdrift is correct

double R_xy = sqrt(x*x + y*y);
//double x_in, y_in, R_in; //internal
double x0, y0, R0; // intersection

// find the point where the line joining 0,0 and x,y
// intersects the internal quad -> x_in, y_in (R_in)

	if(x!=0)
	{
	double x1 =0, y1 = 0;
	double x2 = x, y2 = y;
	double a_ext = (y1 - y2)/(x1 - x2);
	double b_ext =  y1 - a_ext*x1;
	double a_int, b_int;

	// y = a_ext*x +b_ext
		if(x > 0 && y >=0)
		{
			a_int = a1;
			b_int = b1;
		}
		else if(x < 0 && y >= 0)
		{
			a_int = a2;
			b_int = b2;
		}
		else if(x < 0 && y <= 0)
		{
			a_int = a3;
			b_int = b3;
		}
		else if(x > 0 && y <= 0)
		{
		 	a_int = a4;
		 	b_int = b4;
		}
		else  // happened to be
			opserr << "YieldSurface_BC2D::getDrift(..) - condition not possible, x = " << x << ", y = " << y << endln;

		if(driftDebug)
		{
			opserr << "Equation Internal: a = " << a_int << ", b = " << b_int << "\n";
			opserr << "Equation External: a = " << a_ext << ", b = " << b_ext << "\n";
		}
		x0 = -(b_int - b_ext)/(a_int - a_ext);
		y0 = a_int*x0 + b_int;
	}
	else // x = 0
	{
	 	if(y >= 0)
	 	{
	 		x0 = 0;
	 		y0 = yPos - offset;
	 	}
	 	else
	 	{
	 	 	x0 = 0;
	 	 	y0 = yNeg + offset;
	 	}
	}

	R0 = sqrt(x0*x0 + y0*y0);

	if(driftDebug)
	{
		opserr << "R_xy = " << R_xy << " (x= " << x << ", y= " << y << "), R0 = " << R0;
		opserr << " (x0= " << x0 << ", y0= " << y0 << ")\n";
	}

	// If the R_xy < R0, then the point is inside the
	// internal quad, sdrift is correct
	if(R_xy < R0)
	{
		if(driftDebug)
		{
			opserr << " R_xy < R0, returning sdrift " << sdrift << "\n";
			opserr << "\a";
		}
		return sdrift;
	}
	// Now the point can be between the inner quad and the ys
	// or outside the surface

	if(R0 == 0)
		opserr << "ERROR: YieldSurface_BC2D::getDrift(..) - R0 = 0 (yPos="<<yPos<<", yNeg="<<yNeg<<"\n";

double delx = (x0/R0)*increment; // increment already defined
double dely = (y0/R0)*increment;

double xOld, yOld, xNew, yNew;
double xi, yi, Ri; //incremental

	xOld = x0;
	yOld = y0;

int count = 0;
	while(1)
	{
		Ri = sqrt(xOld*xOld + yOld*yOld);

		// check if point is between the inner quad and
		// the yield surface
		if(R_xy < Ri)
		{
			if(driftDebug)
			{
				opserr << " R_xy < Ri, returning sdrift " << sdrift << "\n";
				opserr << "\a";
			}
			return sdrift;
		}
		xNew = xOld + delx;
		yNew = yOld + dely;

		if(getSurfaceDrift(xNew, yNew) > 0) // we just crossed the surface
		{
			double t = interpolateClose(xOld, yOld, xNew, yNew);
			xi =  xOld + t*delx;
			yi =  yOld + t*dely;

			Ri = sqrt(xi*xi + yi*yi);

			if(driftDebug)
			{
				opserr << " Set to surface at xi = " << xi << ", yi = " << yi << ", Ri = " << Ri << "\n";
				opserr << " Returning drift " << R_xy - Ri << "\n";
				opserr << "\a";
			}
			return (R_xy - Ri);
		}

		xOld = xNew;
		yOld = yNew;
		count++;
		if(count > 100)
		{
			opserr << "ERROR: YieldSurface_BC2D::getDrift(..) - not converging\n";
			opserr << "\a";
		}

	}

	opserr << "YieldSurface_BC2D::getDrift(..) - should not reach here\n";
	opserr << "\a";

	return sdrift;
}




//////////////////////////////////////////////////////////////////////
// Other
//////////////////////////////////////////////////////////////////////

void YieldSurface_BC2D::customizeInterpolate(double &xi, double &yi, double &xj, double &yj)
{
double yCheck = yNeg;

	if(yj > 0)
		yCheck = yPos;

	if(fabs(yj) > fabs(yCheck)) // switch to radial - whatever the case
	{
		xi = 0;
		yi = 0;
	}

}

double	YieldSurface_BC2D::interpolate(double xi, double yi, double xj, double yj)
{
	this->customizeInterpolate(xi, yi, xj, yj);

//check
double di = getDrift(xi, yi);
double dj = getDrift(xj, yj);
	
	if(di > 0 && fabs(di < 1e-7))
		return 0;
	if(dj < 0 && fabs(dj < 1e-7))
		return 1;

	if( di > 0)
	{
	 	opserr << "ERROR - YieldSurface_BC2D::interpolate(xi, yi, xj, yj)\n";
		opserr << "point 1 is outside\n";
		opserr << xi << "," << yi << "  " << xj << "," << yj << " : "<< di<<"\n";
		opserr << "\a";
		return 0;
	}
	else if(   dj <0 )
	{
	 	opserr << "ERROR - YieldSurface_BC2D::interpolate(xi, yi, xj, yj)\n";
		opserr << "point 2 is inside\n";
		opserr << xi << "," << yi << "  " << xj << "," << yj << " : "<< dj<<"\n";
		hModel->Print(opserr);
		opserr << "\a";
		return 0;
	}

double tr, tu, tl, dtu, dtl, dtr=100;
double dy = yj - yi;
double dx = xj - xi;
int count = 0;
	tu = 1; tl =0;

	//double d = getDrift(xi, yi, false);
	//if(d>0) opserr << "WARNING - Orbison2D::interpolate, Drift inside > 0 (" << d << ")\n";

	while(fabs(dtr) > 1.0e-7)
	{
		count++;
		if(count > 1000)
		{
			opserr << "\nYieldSurface_BC2D::Interpolate()-> Error: Unable to converge\n";
			opserr << "xi, yi: " << xi << ","<< yi << "\t xj, yj: " << xj << "," << yj << "\n";
			opserr << "Drift Point j = " << dj << "\n";
			hModel->Print(opserr);
			opserr << "\a";
			return 1;
		}

		dtl = getDrift(xi + tl*dx, yi + tl*dy);
		dtu = getDrift(xi + tu*dx, yi + tu*dy);

		tr	=	tu - (  dtu*(tl - tu)/(dtl - dtu)  );
		dtr = getDrift(xi + tr*dx, yi + tr*dy);

		if(dtr >= 0)
		{
			if(dtu >= 0)
				tu = tr;
			else
				tl = tr;
		}
		else
		{
			if(dtu < 0)
				tu = tr;
			else
				tl = tr;
		}

	}// while
	//opserr << "opserr for iterpolation = " << count << "\n"; 5 ~ 15
	return tr;
}

double	YieldSurface_BC2D::interpolateClose(double xi, double yi, double xj, double yj)
{

//check
double di = getSurfaceDrift(xi, yi);
double dj = getSurfaceDrift(xj, yj);
	if( di > 0)
	{
	 	opserr << "ERROR - YieldSurface_BC2D::interpolateClose(xi, yi, xj, yj)\n";
		opserr << "point 1 is outside\n";
		opserr << xi << "," << yi << "  " << xj << "," << yj << " : "<< di<<"\n";
		opserr << "\a";
		return 0;
	}
	else if(   dj <0 )
	{
	 	opserr << "ERROR - YieldSurface_BC2D::interpolateClose(xi, yi, xj, yj)\n";
		opserr << "point 2 is inside\n";
		opserr << xi << "," << yi << "  " << xj << "," << yj << " : "<< dj<<"\n";
		hModel->Print(opserr);
		opserr << "\a";
		return 0;
	}

double tr, tu, tl, dtu, dtl, dtr=100;
double dy = yj - yi;
double dx = xj - xi;
int count = 0;
	tu = 1; tl =0;

	//double d = getDrift(xi, yi, false);
	//if(d>0) opserr << "WARNING - Orbison2D::interpolate, Drift inside > 0 (" << d << ")\n";

	while(fabs(dtr) > 1.0e-7)
	{
		count++;
		if(count > 1000)
		{
			opserr << "\nYieldSurface_BC2D::InterpolateClose()-> Error: Unable to converge\n";
			opserr << "xi, yi: " << xi << ","<< yi << "\t xj, yj: " << xj << "," << yj << "\n";
			hModel->Print(opserr);
			opserr << "\a";
			return 1;
		}

		dtl = getSurfaceDrift(xi + tl*dx, yi + tl*dy);
		dtu = getSurfaceDrift(xi + tu*dx, yi + tu*dy);

		tr	=	tu - (  dtu*(tl - tu)/(dtl - dtu)  );
		dtr = getSurfaceDrift(xi + tr*dx, yi + tr*dy);

		if(dtr >= 0)
		{
			if(dtu >= 0)
				tu = tr;
			else
				tl = tr;
		}
		else
		{
			if(dtu < 0)
				tu = tr;
			else
				tl = tr;
		}

	}// while
	//opserr << "opserr for iterpolation = " << count << "\n"; 5 ~ 15
	return tr;
}

double YieldSurface_BC2D::setToSurface(Vector &force, int algoType, int color)
{
double x2, y2;
double xi, yi, xj, yj;
double dx, dy, t;
double x, y;

	if(returnDebug)
	{
		opserr << "\nYieldSurface_BC2D::setToSurface(Vector &force, int algoType)\n";
		opserr << "Element system force = " << force << "\n";
	}

	// if force point is already on surface,
	// no need to proceed further
    if(getTrialForceLocation(force) == 0)
    	return 0;


	toLocalSystem(force, x2, y2, true);
	xj = x2;
	yj = y2;

	hModel->toOriginalCoord(xj, yj);

	if(color != 0)
	{
		theView->clearImage();
		this->displaySelf(*theView, 1, 1);
		theView->startImage();
		this->displayForcePoint(false, xj, yj, color);
	}


	if(returnDebug)
	{
		opserr << "Local system force - " << "fx = " << x2 << ",\tfy = " << y2 << "\n";
		opserr << "toOriginalCoord    - " << "fx = " << xj << ",\tfy = " << yj << "\n";
	}

		switch(algoType)
		{
			case 0: //dFReturn:
			{
				xi = fx_hist;
				yi = fy_hist;
				// hModel->toOriginalCoord(xi, yi); - stored as that

				break;
			}

			case 1: //RadialReturn:
			{
				xi = 0;
				yi = 0;

				break;

			}

			case 2: //ConstantXReturn:
			{
				xi = xj;
				yi = 0; 					//if point is outside the ys

				if(getDrift(xj, yj) < 0) 	//point is inside the ys
				{
				  //opserr << "ConstantXReturn - Point is inside: " << getDrift(xj, yj) << "\n";
					if(yj < 0)				//point is below x-axis
						yj = yj - 1;
					else
						yj = yj + 1;
				}
				break;
			}


			case 3: //ConstantYReturn:
			{
				xi = 0; 					//if point is outside the ys
				yi = yj;

				if(getDrift(xj, yj) < 0) 	//point is inside the ys
				{
				    //opserr << " ConstantYReturn - Point Inside\n"; //
					if(xj < 0)				//point is left of y-axis
						xj = xj - 1;
    				else
						xj = xj + 1;
				}
				break;
			}


			default:
			{
				opserr << "YieldSurface_BC2D: Method not implemented yet\n";
				xi = 0; yi = 0; //revert to radial return
				break;
			}
		}//switch

		dx = xj - xi;
		dy = yj - yi;

		t =  interpolate(xi, yi, xj, yj);
		x =  xi + t*dx;
		y =  yi + t*dy;

		if(color != 0)
		{
			this->displayForcePoint(false, x, y, color);
			theView->doneImage();
			opserr << "\a";
		}

		hModel->toDeformedCoord(x, y);

		toElementSystem(force, x, y, true);

        return t;
}

int YieldSurface_BC2D::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
	if(displayMode == this->SurfOnly)
		return 0;

	hModel->displaySelf(theViewer, this->SurfOnly, fact);
	
Vector p1(3), p2(3);
Vector rgb(3);
// Draw the axis
rgb(0) = 0.8; rgb(1) = 0.8; rgb(2) = 0.8;
 double d = 0.04;

	p1(0) =  -10;
	p1(1) =    0;
	p2(0) =   10;
	p2(1) =    0;
	theViewer.drawLine(p1, p2, rgb, rgb);

	p1(0) =    0;
	p1(1) =  -10;
	p2(0) =    0;
	p2(1) =   10;
	theViewer.drawLine(p1, p2, rgb, rgb);
// Mark every 0.5
	for(double i=-10; i <= 10; i = i+ 0.5)
	  {
		p1(0) = -d;
		p1(1) = i;
		p2(0) = d;
		p2(1) = i;
		theViewer.drawLine(p1, p2, rgb, rgb);
	  }
	for(double j=-10; j <= 10; j = j+0.5)
	  {
		p1(0) = j;
		p1(1) = -d;
		p2(0) = j;
		p2(1) = d;
		theViewer.drawLine(p1, p2, rgb, rgb);
	  }

	this->displayCommitForcePoint(theViewer, displayMode, fact);
	this->displayForcePoint(true, 0.0, 0.0, 0);

	// draw the 4 quadrants - to make sure that the
	// line coef are correct, just need 1 point between
	// tested okay for Orbison, Attalla, Hajjar
	// return here till more surfaces have to be tested

	if(!driftDebug)
	{
		return 0;
	}

	// quad 1
	p1(0) = 0;
	p1(1) = yPos - offset;

	p2(0) = 0.5;
	p2(1) = a1*p2(0) + b1;
	theViewer.drawLine(p1, p2, rgb, rgb);

	p1(0) = xPos - offset;
	p1(1) = 0;
	theViewer.drawLine(p1, p2, rgb, rgb);


	// quad 2
	p1(0) = 0;
	p1(1) = yPos - offset;

	p2(0) = -0.5;
	p2(1) = a2*p2(0) + b2;
	theViewer.drawLine(p1, p2, rgb, rgb);

	p1(0) = xNeg + offset;
	p1(1) = 0;
	theViewer.drawLine(p1, p2, rgb, rgb);

	// quad 3
	p1(0) = 0;
	p1(1) = yNeg + offset;

	p2(0) = -0.5;
	p2(1) = a3*p2(0) + b3;
	theViewer.drawLine(p1, p2, rgb, rgb);

	p1(0) = xNeg + offset;
	p1(1) = 0;
	theViewer.drawLine(p1, p2, rgb, rgb);

	// quad 4
	p1(0) = 0;
	p1(1) = yNeg + offset;

	p2(0) = 0.5;
	p2(1) = a4*p2(0) + b4;
	theViewer.drawLine(p1, p2, rgb, rgb);

	p1(0) = xPos - offset;
	p1(1) = 0;
	theViewer.drawLine(p1, p2, rgb, rgb);

	return 0;
}

int YieldSurface_BC2D::displayCommitForcePoint(Renderer &theViewer, int displayMode, float fact)
{
Vector p1(3), p2(3);
Vector rgb(3);
rgb(0) = 1; rgb(1) = 0; rgb(2) = 0;
//  toOriginalCoord(fx_hist, fy_hist);
double isotropicFactor = hModel->getCommitIsotropicFactor(0);
double del = 0.1*isotropicFactor;//( (fx_hist + fy_hist)/2);
 if(del< 0.05) del = 0.05;

double fx = fx_hist;
double fy = fy_hist;

	hModel->toDeformedCoord(fx, fy);

	p1(0) =  fx - del;
	p1(1) =  fy;
	p2(0) =  fx + del;
	p2(1) =  fy;
	theViewer.drawLine(p1, p2, rgb, rgb);

	p1(0) =  fx;
	p1(1) =  fy - del;
	p2(0) =  fx;
	p2(1) =  fy + del;
	theViewer.drawLine(p1, p2, rgb, rgb);

	return 0;
}


int YieldSurface_BC2D::displayForcePoint(Vector &force, int color)
{
 if(!theView)
		return -1;

	double x, y;
	toLocalSystem(force, x, y, true);

	theView->startImage();
		displayForcePoint(false, x, y, color);
	theView->doneImage();

	return 0;
}
// color can be r, g or b
int YieldSurface_BC2D::displayForcePoint(bool toDeformed, double f_x, double f_y, int color)
{
Vector p1(3), p2(3);
Vector rgb(3);

	if(!theView)
		return -1;

	if(color== 1)
	{
		rgb(0) = 1; rgb(1) = 0; rgb(2) = 0;
	}
	else if(color == 2)
	{
       rgb(0) = 0; rgb(1) = 1; rgb(2) = 0;
	}
    else if(color == 3)
	{
       rgb(0) = 0; rgb(1) = 0; rgb(2) = 1;
	}
	else
	{
		rgb(0) = 0; rgb(1) = 0; rgb(2) = 0;
	}
//  toOriginalCoord(fx_hist, fy_hist);
/*
double isotropicFactor = hModel->getTrialIsotropicFactor(0);
double del = 0.1*isotropicFactor;//( (fx_hist + fy_hist)/2);
 if(del< 0.05) del = 0.05;
*/
double fx = f_x;
double fy = f_y;

	if(toDeformed)
		hModel->toDeformedCoord(fx, fy);

	v2(0) = fx;
	v2(1) = fy;

	theView->drawPoint(v2, rgb, 3);
/*
	p1(0) =  fx - del;
	p1(1) =  fy;
	p2(0) =  fx + del;
	p2(1) =  fy;
	theView->drawLine(p1, p2, rgb, rgb);

	p1(0) =  fx;
	p1(1) =  fy - del;
	p2(0) =  fx;
	p2(1) =  fy + del;
	theView->drawLine(p1, p2, rgb, rgb);
*/
	return 0;
}


