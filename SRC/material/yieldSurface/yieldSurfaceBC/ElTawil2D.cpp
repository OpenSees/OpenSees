// ElTawil2D.cpp: implementation of the YieldSurfaceBC class.
//
//////////////////////////////////////////////////////////////////////

#include "ElTawil2D.h"
#include <math.h>

#define ELTAWIL2D_CLASS_TAG -1

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

ElTawil2D::ElTawil2D(int tag, double xbal, double ybal, double ypos, double yneg,
					YS_Evolution &model, double cz_, double ty_)
		 :YieldSurface_BC2D(tag, ELTAWIL2D_CLASS_TAG, 0, 0, model),
		 xBal(xbal), yBal(ybal), yPosCap(ypos), yNegCap(yneg), cz(cz_),
		 yPosCap_orig(ypos), yNegCap_orig(yneg), ty(ty_), qy(0.005) //0.01
{
	//capY = yPosCap;
	capY    = yPosCap - yBal;

	// set to origin
	yPosCap -= yBal;
	yNegCap -= yBal;

	// set translation
	double transY   = yBal/capY;
	Vector t(2);
	t(0) = 0;
	//t(1) = yBal/capY;
	t(1) = transY;
	hModel->setInitTranslation(t);

//	opserr << "Translation Y = " << transY << endln;

	//capX = xBal*(1 - pow( fabs((transY*capY)/yNegCap) , ty));
	capX    = xBal;

	// bad:
	capX_orig = capX;
	capY_orig = capY;
	capXdim = capX;
	capYdim = capY;
	
}

ElTawil2D::~ElTawil2D()
{

}

//////////////////////////////////////////////////////////////////////
// YS specific methods
//////////////////////////////////////////////////////////////////////
void ElTawil2D::setExtent()
{
	// Extent along the axis
//	xPos =  1;
//	xNeg = -1;

	xPos = xBal/capX;
	xNeg = -xPos;

	yPos =  yPosCap/capY - qy; //0.02  0.98;
	yNeg =  yNegCap/capY + qy; //-0.98;

	ytPos = yPos - 0.005; // 0.01 // 0.03 // 0.95
	ytNeg = yNeg + 0.005;

double yValPos = ytPos*capY;
double yValNeg = ytNeg*capY;

	xtPos = xBal*(1 - pow((yValPos/yPosCap) , cz));
	xtNeg = xBal*(1 - pow( fabs(yValNeg/yNegCap) , ty));

	xtPos = xtPos/capX;
	xtNeg = xtNeg/capX;
/*
	opserr << "ytPos = " << ytPos << ", ytNeg = " << ytNeg << ", xtPos = " << xtPos
	                            << ", xtNeg = " << xtNeg << endln;
*/
}


void ElTawil2D::getGradient(double &gx, double &gy, double x, double y)
{
// check if the point is on the surface
    double drift =  getDrift(x, y);
    double loc   =  forceLocation(drift);
    double capx = capXdim;
    double capy = capYdim;

    if(loc != 0)
    {
     	opserr << "ERROR - ElTawil2D::getGradient(double &gx, double &gy, double x, double y)\n";
        opserr << "Force point not on yield surface, drift = " << drift << " loc = " << loc <<"\n";
        // opserr << "\a";
        gx = 1.0;
        gy = 1.0;
    }
    else
    {
		double a = 10.277;//3.043;//4.29293; //by compatibility, check err in gradient
		// double yt = 0.95;

		if(y > ytPos)
		{
		 	gx = 2*a*x/capx;
			gy = 1;
		}
		else if(y < ytNeg)
		{
		 	gx = 2*a*x/capx;
			gy = -1; //-1*y/capy; <-- why did i write this?// -1 ?
		}
		else
		{
			// double xVal = x*capx;
			double yVal = fabs(y*capy);  //!!
			// double yVal = y*capy;           //!!

			gx = 1/xBal;
			if(x < 0)
				gx = -1*gx;

			if(y < 0)
				gy = -1*(1/pow( fabs(yNegCap), ty))*ty*(pow(yVal,ty-1));
			else
				gy = (1/pow(yPosCap, cz))*cz*(pow(yVal,cz-1));
		}
	}

	// opserr << "gx = " << gx << ", gy = " << gy << ", capy = " << capy << "\n";
	// opserr << "\a";
}

double ElTawil2D::getSurfaceDrift(double x, double y)
{
double phi;
double a = 5;//10.277; //3.043;//4.29293; --> effects convergence


	if(y > ytPos && fabs(x) < fabs(y*xtPos/ytPos) )
	{
		phi = a*x*x + y + qy;
	}
	else if(y < ytNeg && fabs(x) < fabs(y*xtNeg/ytNeg))
	{
		phi = a*x*x - y + qy;
	}
	else
	{
		double xVal = x*capX;
		double yVal = y*capY;

		if(y < 0)
			phi = fabs(xVal/xBal) + pow(fabs(yVal/yNegCap), ty);
		else
			phi = fabs(xVal/xBal) + pow(yVal/yPosCap, cz);
	}

	double drift = phi - 1;

//	opserr << "Eltawil2D - x = " << x << ", y = " << y << ", drift = " << drift << endln;

	return drift;
}

// Always call the base class method first
void ElTawil2D::customizeInterpolate(double &xi, double &yi, double &xj, double &yj)
{
	this->YieldSurface_BC2D::customizeInterpolate(xi, yi, xj, yj);

	if(yj >= ytPos && fabs(xj) <= fabs(yj*xtPos/ytPos) )
	{
		xi = 0;
		yi = 0;
	}
	else if(yj <= ytNeg && fabs(xj) <= fabs(yj*xtNeg/ytNeg))
	{
		xi = 0;
		yi = 0;
	}
	else
		; // values are okay


}


YieldSurface_BC *ElTawil2D::getCopy(void)
{
    ElTawil2D *theCopy = new ElTawil2D(this->getTag(), xBal, yBal, yPosCap_orig, yNegCap_orig, *hModel,
                                       cz, ty);
    //later  copy all the state variables
    return theCopy;
}

int ElTawil2D::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
 	this->YieldSurface_BC2D::displaySelf(theViewer, displayMode, fact);
Vector pOld(3), pCurr(3);
Vector rgb(3);
rgb(0) = 0.1; rgb(1) = 0.5; rgb(2) = 0.5;
	if(displayMode == this->SurfOnly)
	{
		rgb(0) = 0.7; rgb(1) = 0.7; rgb(2) = 1.0;
	}


	double incr = ((yPosCap/capY) - (yNegCap/capY)) / 10;
	//double incr =  fabs(0.33333333*yNegCap/capY);
	if(fact < 1) incr = fact;
 	
 	double xOld = 0;
 	double yOld = yNegCap/capY;
 	// hModel->toDeformedCoord(xOld, yOld);
 	
 	double err = 1e-4;
    for(double yc = yNegCap/capY; yc <= yPosCap/capY + err; yc = yc+incr)
	{
		double y    = yc;	
		double yVal = y*capY;
		double xVal;

		if(y < 0)
		{
			xVal = xBal*(1 - pow( fabs(yVal/yNegCap) , ty));
		}
		else
		{
			xVal = xBal*(1 - pow( (yVal/yPosCap) , cz));
		}

		double x = xVal/capX;

		if(displayMode==100)
			opserr << "(undeformed) x = " << x << ", y = " << y;
 		
 		double x1 = x;
		double y1 = y;
		double x2 = -x;
		double y2 = y;
		
		double x1Old = xOld;
		double y1Old = yOld;
		double x2Old = -xOld;
		double y2Old = yOld;
			

        hModel->toDeformedCoord(x1, y1);
        hModel->toDeformedCoord(x1Old, y1Old);
        hModel->toDeformedCoord(x2, y2);
        hModel->toDeformedCoord(x2Old, y2Old);

//		if(displayMode==100)
//			opserr << " (deformed) x = " << x << ", y = " << y << endln;
 		
		pCurr(0) = x1;
		pCurr(1) = y1;
		pOld(0)  = x1Old;
		pOld(1)  = y1Old;
		
		theViewer.drawLine(pOld, pCurr, rgb, rgb);
		
		pCurr(0) = x2;
		pCurr(1) = y2;
		pOld(0)  = x2Old;
		pOld(1)  = y2Old;
		theViewer.drawLine(pOld, pCurr, rgb, rgb);
		
		xOld = x;
		yOld = y;
	}


	// displayForcePoint(theViewer, displayMode, fact);
	return 0;
}

void ElTawil2D::Print(OPS_Stream &s, int flag)
{
    s << "\nYield Surface No: " << this->getTag() << " type: ElTawil2D\n";
	this->YieldSurface_BC::Print(s, flag);
}


/*
// 	return 0;

Vector pOld(3), pCurr(3);
Vector rgb(3);
rgb(0) = 0.1; rgb(1) = 0.5; rgb(2) = 0.5;
double incr = 0.2;// 0.02,  incr = 0.0001;
double x1, y1, xOld, yOld, x2, y2;
double t;

	if(fact < 1) incr = fact;
    xOld = 0;
    //yOld = 1;

	t = interpolate(0, 0, 0, 1);
	yOld = t*1;

	double x, y = 10;
	double err = 1e-5;

    for(double xc = 0; xc <= xPos+err; xc = xc+incr) // xc <= xPos + err
	{
		if(xc > xPos) xc = int(xPos);
	 	double yc = sqrt(xPos - xc*xc); // eqn of a circle 2.1 -> 1

   		t = interpolate(0, 0, xc, yc); //radial return
     	x = t*xc;
		y = t*yc;

		if(fact >=1 && x >= 0.64) incr = 0.05;
//		if(fact >=1 && x >= 0.84) incr = 0.02;
//		if(fact >=1 && x >= 0.94) incr = 0.01;		

		//if(x < 0.06 || x > 0.9)
        {

            //////////////////////// x>0, y>0
            x1 = x;
            y1 = y;
            hModel->toDeformedCoord(x1, y1);
			if(displayMode==100)
			{
				opserr << " x = " << x << ", y = " << y << " ";
				opserr << " x1 = " << x1 << ", y1 = " << y1 << "\n";
			}
            pCurr(0) = x1;
            pCurr(1) = y1;

            x2 = xOld;
            y2 = yOld;
            hModel->toDeformedCoord(x2, y2);
            pOld(0) = x2;
            pOld(1) = y2;
            theViewer.drawLine(pOld, pCurr, rgb, rgb);
			//opserr << "pOld = " << pOld << " pCurr = " << pCurr << "\n";

            ///////////////////////// x<0, y>0
            x1 = -1*x;
            y1 = y;
            hModel->toDeformedCoord(x1, y1);
            pCurr(0) = x1;
            pCurr(1) = y1;

            x2 = -1*xOld;
            y2 = yOld;
            hModel->toDeformedCoord(x2, y2);
            pOld(0) = x2;
            pOld(1) = y2;

            theViewer.drawLine(pOld, pCurr, rgb, rgb);

            //////////////////////// x>0, y<0
            t = interpolate(0, 0, xc, -yc); //radial return
     		x =    t*xc;
			y = -1*t*yc;

            x1 = x;
            y1 = y;
            hModel->toDeformedCoord(x1, y1);
            pCurr(0) = x1;
            pCurr(1) = y1;

            x2 = xOld;
            y2 = -1*yOld;
            hModel->toDeformedCoord(x2, y2);
            pOld(0) = x2;
            pOld(1) = y2;

            theViewer.drawLine(pOld, pCurr, rgb, rgb);

            //////////////////////// x<0, y<0
            x1 = -1*x;
            y1 = -1*y;
            hModel->toDeformedCoord(x1, y1);
            pCurr(0) = x1;
            pCurr(1) = y1;

            x2 = -1*xOld;
            y2 = -1*yOld;
            hModel->toDeformedCoord(x2, y2);
            pOld(0) = x2;
            pOld(1) = y2;

            theViewer.drawLine(pOld, pCurr, rgb, rgb);


            xOld = x;
            yOld = y;
        }//x > 0

	}//while
*/

