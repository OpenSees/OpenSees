// Hajjar2D.cpp: implementation of the YieldSurfaceBC class.
//
//////////////////////////////////////////////////////////////////////

#include "Hajjar2D.h"
#include <math.h>

#define coefDebug 1
#define HAJJAR_CLASS_TAG -1
//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////
Hajjar2D::Hajjar2D(int tag, double xmax, double ymax,
			  YS_Evolution &model,
			  double centroid_y, double c1_, double c2_, double c3_)
:YieldSurface_BC2D(tag, HAJJAR_CLASS_TAG, xmax, ymax, model),
	centroidY(centroid_y),c1(c1_), c2(c2_), c3(c3_)
{
//!!	translateY = translateY_hist = centroidY;   // fix!
}


Hajjar2D::Hajjar2D(	 int tag,
                     YS_Evolution &model,
					 double D, double b, double t,
					 double fc_, double fy_)
:YieldSurface_BC2D(tag, HAJJAR_CLASS_TAG, 0, 0, model),//would blow up if not..
					 depth(D), width(b), thick(t), fc(fc_), fy(fy_)
{
double fr  = 0.623*sqrt(fc);
double Ast = D*b - (D - 2*t)*(b - 2*t);
double Ac  = (D - 2*t)*(b - 2*t);

double x = D/t;
double y = fc/fy;

	c1 = 1.08 - 0.00265*x + 0.000023*x*x - 1.13*1e-7*x*x*x +
		 0.374*y - 1.3*y*y - 0.0419*y*y*y - 0.0691*x*y +
		 0.000234*x*x*y + 0.0754*x*y*y;

	c2 = 0.628 + 0.0259*x - 0.000367*x*x + 1.99*1e-6*x*x*x +
		 4.5*y - 14.9*y*y + 22.4*y*y*y + 0.164*x*y +
		 -0.000756*x*x*y - 0.126*x*y*y;

	c3 = 0.42 + 0.0892*x - 0.00122*x*x + 5.13*1e-6*x*x*x +
		 4.9*y - 16.5*y*y + 16.2*y*y*y - 0.165*x*y +
		 0.000713*x*x*y + 0.12*x*y*y;

	capY = Ast*fy + Ac*fc;

double num   = fc*(b*t - 2*t*t) + 0.5*fr*(D - t)*(b - 2*t) + fy*(2*D*t);
double denom = fc*(b - 2*t) + 0.5*fr*(b - 2*t) + fy*(4*t);

double xn = num/denom;
	capX  = 	  fc*(0.5*(b - 2*t)*(xn - t)*(xn - t))
				+ 0.5*fr*(0.5*(b - 2*t)*(D - xn - t)*(D - xn - t))
				+ fy*((2*t)*(D*D/2 + xn*xn + t*t - D*t - D*xn) + (b*t)*(D - t));

	centroidY = (Ac*fc - Ac*fr)/2;
	centroidY = centroidY/capY;

//!!	translateY = translateY_hist = centroidY;
    Vector tt(2);
	tt(0) = 0;
	tt(1) = centroidY;
	hModel->setInitTranslation(tt);

	if (coefDebug)
	{
		opserr << " c1 = " << c1 << ", c2 = " << c2 << ", c3 = " << c3 << "\n";
		opserr << " centroidY = " << centroidY << "\n";
		opserr << " capX = " << capX << ", capY = " << capY << "\n";
	}
	// bad:
	capX_orig = capX;
	capY_orig = capY;
	capXdim = capX;
	capYdim = capY;

}

Hajjar2D::~Hajjar2D()
{

}

//////////////////////////////////////////////////////////////////////
// YS specific methods
//////////////////////////////////////////////////////////////////////

void Hajjar2D::setExtent()
{
	// Extent along the axis
	xPos = sqrt(1/fabs(c1));
	xNeg = -xPos;
	yPos = sqrt(1/fabs(c2));
	yNeg = -yPos;

}

void Hajjar2D::getGradient(double &gx, double &gy, double xi, double yi)
{
// check if the point is on the surface
    double drift =  getDrift(xi, yi);
    //!! why is capXdim not here??
    
    /*if(drift < -error)
    {
        opserr << "ERROR - Hajjar2D::getGradient(double &gx, double &gy, double x, double y)\n";
        opserr << "Point inside the yield surface\n";
		opserr << " fx = " << xi << ", fy = " << yi  << " drift = " << drift << "\n";
        opserr << "\a";
    }
    else if(drift > error)
    {
        opserr << "ERROR - Hajjar2D::getGradient(double &gx, double &gy, double x, double y)\n";
        opserr << "Point outside the yield surface\n";
		opserr << " fx = " << xi << ", fy = " << yi  << " drift = " << drift << "\n";
        opserr << "\a";
    }*/
    if(forceLocation(drift)!=0)
    {
     	opserr << "ERROR - Hajjar2D::getGradient(double &gx, double &gy, double x, double y)\n";
        opserr << "Force point not on the yield surface\n";
		opserr << " fx = " << xi << ", fy = " << yi  << " drift = " << drift << "\n";
        opserr << "\a";
    }
    else
    {
		double x = xi;
		double y = yi /*- centroidY*/;

    	gx = 2*c1*x + 2*c3*pow(y, 2)*(x);
        gy = 2*c2*y + 2*c3*pow(x, 2)*(y);
    }

}

double Hajjar2D::getSurfaceDrift(double xi, double yi)
{
double x = xi;
double y = yi/* - centroidY*/;

double phi = c1*x*x + c2*y*y + c3*y*y*x*x;
double drift = phi - 1;
	return drift;
}

YieldSurface_BC *Hajjar2D::getCopy(void)
{
    Hajjar2D *theCopy = new Hajjar2D(this->getTag(), *hModel,
									 depth, width, thick, fc, fy);
    //later  copy all the state variables
    return theCopy;
}

int Hajjar2D::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
 	this->YieldSurface_BC2D::displaySelf(theViewer, displayMode, fact);
// 	return 0;

Vector pOld(3), pCurr(3);
Vector rgb(3);
rgb(0) = 0; rgb(1) = 0; rgb(2) = 0;
double incr = 0.1;
double x1, y1, xOld, yOld, x2, y2;
 double phi = 1;
    xOld = 0; //y = 0
	yOld = sqrt( (phi - c1*xOld*xOld)/(c2 + c3*xOld*xOld) ); //x = 1
	/*yOld += centroidY;*/

	double xmax = sqrt(1/c1);

	opserr << " xmax = " << xmax << ", ymax = " << yOld << "( " << sqrt(1/c2) << ")\n";
	if(fact < 1) incr = fact;

	double err = 0.5*incr;

    //for(double y = 0; y <= 1+err; y = y+incr)
	for(double x = 0; x <= xmax+err; x = x+incr)
    {
		if(x > xmax) x = xmax;
        double y =  (phi - c1*x*x)/(c2 + c3*x*x);
		if(y > 0) y  = sqrt(y);
		/*y += centroidY;*/

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
            y1 = -1*y /*+ 2*centroidY*/;
            hModel->toDeformedCoord(x1, y1);
            pCurr(0) = x1;
            pCurr(1) = y1;

            x2 = xOld;
            y2 = -1*yOld/* + 2*centroidY*/;
            hModel->toDeformedCoord(x2, y2);
            pOld(0) = x2;
            pOld(1) = y2;

            theViewer.drawLine(pOld, pCurr, rgb, rgb);

            //////////////////////// x<0, y<0
            x1 = -1*x;
            y1 = -1*y/* + 2*centroidY*/;
            hModel->toDeformedCoord(x1, y1);
            pCurr(0) = x1;
            pCurr(1) = y1;

            x2 = -1*xOld;
            y2 = -1*yOld/* + 2*centroidY*/;
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


void Hajjar2D::Print(OPS_Stream &s, int flag)
{
    s << "\nYield Surface No: " << this->getTag() << " type: Attalla2D\n";
}
