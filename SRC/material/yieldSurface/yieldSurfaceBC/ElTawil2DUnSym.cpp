// ElTawil2DUnSym.cpp: implementation of the YieldSurfaceBC class.
//
//////////////////////////////////////////////////////////////////////

#include "ElTawil2DUnSym.h"
#include <math.h>

#define ELTAWIL2DUNSYM_CLASS_TAG -1

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

ElTawil2DUnSym::ElTawil2DUnSym
        (int tag, double xPosbal, double yPosbal,
		double xNegbal, double yNegbal,
		double ypos, double yneg,
		YS_Evolution &model,
		double cz_pos, double ty_pos,
		double cz_neg, double ty_neg)

	:YieldSurface_BC2D(tag, ELTAWIL2DUNSYM_CLASS_TAG, 0, 0, model),
	xPosBal(xPosbal), yPosBal(yPosbal),
	xNegBal(xNegbal), yNegBal(yNegbal),
	yPosCap(ypos), yNegCap(yneg),
	yPosCap_orig(ypos), yNegCap_orig(yneg),
	czPos(cz_pos), tyPos(ty_pos), czNeg(cz_neg), tyNeg(ty_neg), qy(0.005)
{

// This tight coupling with base-class can be a diaster to debug!
	if(yPosBal < 0 || yNegBal < 0)
		opserr << "WARNING - ElTawil2DUnSym() - yBalance < 0" << endln;
		
	//capY = yPosCap; // why did I change this?
	yBal = yPosBal;

	if(yNegBal < yBal)
		yBal = yNegBal;

//    opserr << "yBal= " << yBal << ", yPosBal= " << yPosBal
//	     << ", yNegBal= " << yNegBal << endl << endln;

// 	capY    = yPosCap - yBal;//!!
	capY = yPosCap;

	// set to origin
	yPosCap -= yBal;
	yNegCap -= yBal;

	// larger of the two will align with the axis
	// the other one will be above
	yPosBal -= yBal;
	yNegBal -= yBal;
/*	
	opserr << "yBal= " << yBal << ", yPosBal= " << yPosBal
	     << ", yNegBal= " << yNegBal << endln;
*/

	// set translation
	double transY   = yBal/capY;
	Vector t(2);
	t(0) = 0;
	//t(1) = yBal/capY;
	t(1) = transY;
	
	hModel->setInitTranslation(t);

	capX    = xPosBal;
	if( fabs(xNegBal) > capX)
		capX = fabs(xNegBal);

	// bad:
	capX_orig = capX;
	capY_orig = capY;
	capXdim = capX;
	capYdim = capY;
		
}

ElTawil2DUnSym::~ElTawil2DUnSym()
{

}

//////////////////////////////////////////////////////////////////////
// YS specific methods
//////////////////////////////////////////////////////////////////////
void ElTawil2DUnSym::setExtent()
{
	// Extent along the axis
//	xPos =  1;
//	xNeg = -1;

	xPos = min_(xPosBal/capX,
	            (xPosBal*(1 - pow( fabs(yPosBal/(yNegCap - yPosBal)) , tyPos)))/capX);
	// xNeg = xNegBal/capX;
    xNeg = max_(xNegBal/capX,
                (xNegBal*(1 - pow( fabs(yNegBal/(yNegCap - yNegBal)) , tyNeg)))/capX);

	yPos =  yPosCap/capY - qy; //0.02  0.98;
	yNeg =  yNegCap/capY + qy; //-0.98;

	ytPos = yPos - 0.005; // 0.01 // 0.03 // 0.95
	ytNeg = yNeg + 0.005;

double yVal1 = ytPos*capY - yPosBal;
double yVal4 = ytNeg*capY - yPosBal;
double yVal2 = ytPos*capY - yNegBal;
double yVal3 = ytNeg*capY - yNegBal;

/*
char c = ' ';

       opserr << "yVals = " << yVal1 << c << yVal2 << c
                          << yVal3 << c << yVal4 << endln;

*/
double	xt1 = xPosBal*(1 - pow((yVal1/(yPosCap - yPosBal)) , czPos));
double	xt4 = xPosBal*(1 - pow( fabs(yVal4/(yNegCap - yPosBal)) , tyPos));
double	xt2 = xNegBal*(1 - pow((yVal2/(yPosCap - yNegBal)) , czNeg));
double	xt3 = xNegBal*(1 - pow( fabs(yVal3/(yNegCap - yNegBal)) , tyNeg));


	xt1 = xt1/capX;
	xt2 = xt2/capX;
	xt3 = xt3/capX;
	xt4 = xt4/capX;

/*
	opserr << "xPos = " << xPos << ", xNeg = " << xNeg << endl
	     << "ytPos = " << ytPos << ", ytNeg = " << ytNeg << endl
	     << "xt1   = " << xt1 << ", xt2 = " << xt2  << endl
		 << "xt3   = " << xt3 << ", xt4 = " << xt4 << endln;
*/
}


void ElTawil2DUnSym::getGradient(double &gx, double &gy, double x, double y)
{
// check if the point is on the surface
    double drift =  getDrift(x, y);
    double loc   =  forceLocation(drift);
    double capx = capXdim;
    double capy = capYdim;
    
	// opserr << "ElTawil2DUnSym::getGradient drift:" << this->YieldSurface_BC2D::getDrift(x, y) << endln;

    if(loc != 0)
    {
     	opserr << "ERROR - ElTawil2D::getGradient(double &gx, double &gy, double x, double y)\n";
        opserr << "Force point not on yield surface, drift = " << drift << " loc = " << loc <<"\n";
        opserr << "\a";
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
			gy = -1; //*y/capy;
		}
		else
		{
			double xVal = x*capx;
			double yVal = y*capy;

			if(xVal >= 0 && yVal >= yPosBal) // quadrant 1
			{
				gx = 1/xPosBal;
				gy = (1/pow( (yPosCap-yPosBal), czPos))*czPos*(pow(yVal-yPosBal,czPos-1));
			}
			else if(xVal >=0 && yVal < yPosBal) // quadrant 1 or 4
			{
				gx = 1/xPosBal;
				gy = -1*(1/pow( fabs(yNegCap-yPosBal), tyPos))*tyPos*(pow(fabs(yVal-yPosBal),tyPos-1));
			}
			else if(xVal< 0 && yVal >= yNegBal) // quadrant 2
			{
				gx = 1/xNegBal;  // xNegBal should be < 0
				gy = (1/pow( (yPosCap-yNegBal), czNeg))*czNeg*(pow(yVal-yNegBal,czNeg-1));
			}
			else if(xVal<0 && yVal < yNegBal) // quadrant 2 or 3
			{
				gx = 1/xNegBal;
				gy = -1*(1/pow( fabs(yNegCap-yNegBal), tyNeg))*tyNeg*(pow(fabs(yVal-yNegBal),tyNeg-1));
			}
			else
			{
				opserr << "Eltawil2DUnsym - condition not possible" << endln;
				opserr << "\a";
			}

			/* gx = 1/xBal;
				if(x < 0)
					gx = -1*gx;

				if(y < 0)
					gy = -1*(1/pow( fabs(yNegCap), ty))*ty*(pow(yVal,ty-1));
				else
					gy = (1/pow(yPosCap, cz))*cz*(pow(yVal,cz-1));
			*/
		}
	}

//	opserr << "gx = " << gx << "gy = " << gy << "\n";
//	opserr << "\a";
}

double ElTawil2DUnSym::getSurfaceDrift(double x, double y)
{
double phi;
double a = 5;//10.277; //3.043;//4.29293; --> effects convergence
// why isnt a = 10.77 here or 5 in getGradient
     double capx = capX;
     double capy = capY;
     
	// ignore the small difference between xt1, xt2
	// and xt4, xt3 for determining whether the
	// quadratic should be used
	if(y > ytPos && fabs(x) < fabs(y*xt1/ytPos) )
	{
		phi = a*x*x + y + qy;
	}
	else if(y < ytNeg && fabs(x) < fabs(y*xt4/ytNeg))
	{
		phi = a*x*x - y + qy;
	}
	else
	{
		double xVal = x*capx;
		double yVal = y*capy;

		if(xVal >=0 && yVal >= yPosBal) 	// quad 1
		{
		 	phi = fabs(xVal/xPosBal) + pow((yVal - yPosBal)/(yPosCap - yPosBal), czPos);
		}
		else if(xVal >=0 && yVal < yPosBal)	// quad 1 or 4
		{
		 	phi = fabs(xVal/xPosBal) + pow(fabs((yVal - yPosBal)/(yNegCap - yPosBal)), tyPos);
		}
		else if(xVal < 0 && yVal >= yNegBal)// quad 2
		{
		 	phi = fabs(xVal/xNegBal) + pow((yVal - yNegBal)/(yPosCap - yNegBal), czNeg);
		}
		else if(xVal < 0 && yVal < yNegBal)	// quad 2 or 3
		{
		
		 	phi = fabs(xVal/xNegBal) + pow(fabs((yVal - yNegBal)/(yNegCap - yNegBal)), tyNeg);
		}
		else
		{
			opserr << "ElTawil2DUnSym::getSurfaceDrift(..) - cond not possible\n";
			opserr << "x=" << x << ", y=" << y << ", capx=" << capx << ", capy=" << capy << endln;
			opserr << "xVal = " << xVal << ", yVal = " << yVal << endln;
			opserr << "\a";
		}
		/*	
		if(y < 0)
			phi = fabs(xVal/xBal) + pow(fabs(yVal/yNegCap), ty);
		else
			phi = fabs(xVal/xBal) + pow(yVal/yPosCap, cz);
		*/
	}

	double drift = phi - 1;
	return drift;
}

// Always call the base class method first
void ElTawil2DUnSym::customizeInterpolate(double &xi, double &yi, double &xj, double &yj)
{
	this->YieldSurface_BC2D::customizeInterpolate(xi, yi, xj, yj);
	
	// again, ignore the small difference between xt1, xt2
	// and xt4, xt3

	if(yj >= ytPos && fabs(xj) <= fabs(yj*xt1/ytPos) )
	{
		xi = 0;
		yi = 0;
	}
	else if(yj <= ytNeg && fabs(xj) <= fabs(yj*xt4/ytNeg))
	{
		xi = 0;
		yi = 0;
	}
	else
		; // values are okay


}


YieldSurface_BC *ElTawil2DUnSym::getCopy(void)
{
    // ElTawil2D *theCopy = new ElTawil2D(this->getTag(), xBal, yBal, yPosCap_orig, yNegCap_orig, *hModel,
    //                                   cz, ty);
    //later  copy all the state variables


    ElTawil2DUnSym *theCopy = new ElTawil2DUnSym(this->getTag(), xPosBal, yPosBal+yBal,
    							  xNegBal, yNegBal+yBal, yPosCap+yBal, yNegCap+yBal,
    							  *hModel, czPos, tyPos, czNeg, tyNeg);

    return theCopy;
}

int ElTawil2DUnSym::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
 	this->YieldSurface_BC2D::displaySelf(theViewer, displayMode, fact);
Vector pOld(3), pCurr(3);
Vector rgb(3);
rgb(0) = 0.1; rgb(1) = 0.5; rgb(2) = 0.5;
// rgb(0) = 0.0; rgb(1) = 0.0; rgb(2) = 0.0;

	if(displayMode == this->SurfOnly)
	{
		rgb(0) = 0.7; rgb(1) = 0.7; rgb(2) = 1.0;
		//opserr << "Displaying self" << endln;
	}


	// double incr = 0.101;
 	double incr =  fabs(0.33333333*yNegCap/capY);
//	double incr =  fabs(0.1*yNegCap/capY);

	if(fact < 1) incr = fact;
 	
 	double xOld = 0;
 	double yOld = yNegCap/capY;
 	hModel->toDeformedCoord(xOld, yOld);
 	
 	double err = 0.01;
 	// quad 1 and 4
	double yc;
    for(yc = yNegCap/capY; yc <= yPosCap/capY + err; yc = yc+incr)
    //int iter = 0;
    //while(1)
	{
		//opserr << "incr = " << incr;
		/*switch (iter)
		{
			case 0:{yc = yNegCap/capY; break;}
			case 1:{yc = 0.333*yNegCap/capY; break;}
			case 2:{yc = 0.666*yNegCap/capY; break;}
			case 2:{yc = 0.0; break;}
			case 3:{yc = yPosBal/capY; break;}
			case 4:{yc = 0.5*yPosCap/capY; break;}
			case 5:{yc = yPosCap/capY; break;}
			default: break;
		}
		if(iter == 6)
			break;
		iter++;
		*/
		
		double y    = yc;
		if(y > yPosCap/capY)
			y = yPosCap/capY;
			
		double yVal = y*capY;
		double xVal;

		if(yVal < yPosBal)
		{
			xVal = xPosBal*(1 - pow( fabs((yVal- yPosBal)/(yNegCap - yPosBal)) , tyPos));
			// xVal = 0;
		}
		else
		{
			xVal = xPosBal*(1 - pow( ((yVal - yPosBal)/(yPosCap - yPosBal)) , czPos));
		}

		double x = xVal/capX;

		if(displayMode==100)
			opserr << "(undeformed) x = " << x << ", y = " << y;

        hModel->toDeformedCoord(x, y);
		if(displayMode==100)
			opserr << " (deformed) x = " << x << ", y = " << y << endln;
 		
		pCurr(0) = x;
		pCurr(1) = y;
		pOld(0)  = xOld;
		pOld(1)  = yOld;

		// x >= 0
		theViewer.drawLine(pOld, pCurr, rgb, rgb);
		xOld = x;
		yOld = y;
	}

 	xOld = 0;
 	yOld = yNegCap/capY;
 	hModel->toDeformedCoord(xOld, yOld);

//	opserr << "Plotting negative side\n";
	// xNegBal itself is < 0
	// quad 2 and 3
	//iter = 0;
    for(yc = yNegCap/capY; yc <= yPosCap/capY + err; yc = yc+incr)
    //while(1)
	{
		/*switch (iter)
		{
			case 0:{yc = yNegCap/capY; break;}
			case 1:{yc = 0.5*yNegCap/capY; break;}
			case 2:{yc = 0.0; break;}
			case 3:{yc = yNegBal/capY; break;}
			case 4:{yc = 0.50*yPosCap/capY; break;}
			case 5:{yc = yPosCap/capY; break;}
			default: break;
		}
		if(iter == 6)
			break;
		iter++;
      */
		double y    = yc;
		if(y > yPosCap/capY)
			y = yPosCap/capY;

		double yVal = y*capY;
		double xVal;

		if(yVal < yNegBal)
		{
			xVal = xNegBal*(1 - pow( fabs((yVal-yNegBal)/(yNegCap-yNegBal)) , tyNeg));
		}
		else
		{
			xVal = xNegBal*(1 - pow( ((yVal-yNegBal)/(yPosCap - yNegBal)) , czNeg));
		}

		double x = xVal/capX;

		if(displayMode==100)
			opserr << "(undeformed) x = " << x << ", y = " << y;

        hModel->toDeformedCoord(x, y);
		if(displayMode==100)
			opserr << " (deformed) x = " << x << ", y = " << y << endln;
 		
		pCurr(0) = x;
		pCurr(1) = y;
		pOld(0)  = xOld;
		pOld(1)  = yOld;
		
		// x >= 0
		theViewer.drawLine(pOld, pCurr, rgb, rgb);
		xOld = x;
		yOld = y;
	}

	//displayForcePoint(theViewer, displayMode, fact);
	return 0;
}

void ElTawil2DUnSym::Print(OPS_Stream &s, int flag)
{
    s << "\nYield Surface No: " << this->getTag() << " type: ElTawil2DUnSym\n";
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

