// Orbison2D.cpp: implementation of the YieldSurfaceBC class.
//
//////////////////////////////////////////////////////////////////////

#include "Orbison2D.h"
#include <math.h>
#define ORBISON_CLASS_TAG -1
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Orbison2D::Orbison2D(int tag, double xcap, double ycap,
                     YS_Evolution &model)
:YieldSurface_BC2D(tag, ORBISON_CLASS_TAG, xcap, ycap, model)
{

}

Orbison2D::~Orbison2D()
{

}

//////////////////////////////////////////////////////////////////////
// YS specific methods
//////////////////////////////////////////////////////////////////////

void Orbison2D::setExtent()
{
	// Extent along the axis
	xPos =  1;
	xNeg = -1;
	yPos =  1;
	yNeg = -1;
}


void Orbison2D::getGradient(double &gx, double &gy, double x, double y)
{
// check if the point is on the surface
    double drift =  getDrift(x, y);
    double capx = capXdim;
    double capy = capYdim;
    
    if(forceLocation(drift)!=0)
    {
     	opserr << "ERROR - Orbison2D::getGradient(double &gx, double &gy, double x, double y)\n";
        opserr << "Force point not on the yield surface\n";
		opserr << " fx = " << x << ", fy = " << y  << " drift = " << drift << "\n";
        opserr << "\a";
    }
    else
    {
    	gx = 2*x/(capx) + 7.34*pow(y, 2)*(x/(capx));
    	gy = 2.3*y/(capy) - 0.9*pow(y, 5)/(capy) + 7.34*pow(x, 2)*(y/(capy));
//  	p1 = 2.3*p - 0.9*pow(p, 5) + 7.34*pow(m, 2)*(p);
//  	m1 = 2*m + 7.34*pow(p, 2)*(m);

//      gx = 2*x + 7.34*pow(y, 2)*x;
//      gy = 2.3*y - 0.9*pow(y, 5) + 7.34*pow(x, 2)*y;

    }

}

double Orbison2D::getSurfaceDrift(double x, double y)
{
double phi = 1.15*y*y - 0.15*pow(y, 6) + x*x + 3.67*y*y*x*x;
double drift = phi - 1;
	return drift;
}

YieldSurface_BC *Orbison2D::getCopy(void)
{
    Orbison2D *theCopy = new Orbison2D(this->getTag(), capX, capY, *hModel);
    if(theCopy==0)
    {
     	opserr << "Orbison2D::getCopy(void) - unable to make copy\n";
     	opserr << "\a";
    }
    //later  copy all the state variables
    return theCopy;
}

int Orbison2D::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
 	this->YieldSurface_BC2D::displaySelf(theViewer, displayMode, fact);

Vector pOld(3), pCurr(3);
Vector rgb(3);
rgb(0) = 0; rgb(1) = 0; rgb(2) = 0;
double incr = 0.1;
double x1, y1, xOld, yOld, x2, y2;

    xOld = 1; //x = 1
    yOld = 0.0; //y = 0
    //xOld = (1 - 1.15*yOld*yOld + 0.15*pow(yOld, 6))/(1 + 3.67*yOld*yOld);

	double err = 1e-5;

	if(fact < 1) incr = fact;

    for(double y = 0; y <= 1+err; y = y+incr)
    {
		if(y > 1) y = int(1);
        double x = (1 - 1.15*y*y + 0.15*pow(y, 6))/(1 + 3.67*y*y);
        if (x>0) x = sqrt(x);

		//if( x < 0.2) incr = 0.02;

		//if(x < 0.2 || x > 0.85)
		{
			if(displayMode==100)
				opserr << " x = " << x << ", y = " << y << "\n";

            //////////////////////// x>0, y>0
            x1 = x;
            y1 = y;
            hModel->toDeformedCoord(x1, y1);
            pCurr(0) = x1;
            pCurr(1) = y1;

            x2 = xOld;
            y2 = yOld;
            hModel->toDeformedCoord(x2, y2);
            pOld(0) = x2;
            pOld(1) = y2;
            theViewer.drawLine(pOld, pCurr, rgb, rgb);

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
            x1 = x;
            y1 = -1*y;
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
    }

	// displayForcePoint(theViewer, displayMode, fact);

	return 0;
}


void Orbison2D::Print(OPS_Stream &s, int flag)
{
    s << "\nYield Surface No: " << this->getTag() << " type: Attalla2D\n";
}
