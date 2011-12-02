// NullYS2D.cpp: implementation of the YieldSurfaceBC class.
//
//////////////////////////////////////////////////////////////////////

#include "NullYS2D.h"
#include <math.h>
#include <NullEvolution.h>

NullEvolution NullYS2D::evolModel(0, 1e+20, 1e+20);
 
#define NULLYS2D_CLASS_TAG -1

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

NullYS2D::NullYS2D(int tag)
:YieldSurface_BC2D(tag, NULLYS2D_CLASS_TAG, 1, 1, evolModel)

{

}

NullYS2D::~NullYS2D()
{

}

//////////////////////////////////////////////////////////////////////
// YS specific methods
//////////////////////////////////////////////////////////////////////
void NullYS2D::setExtent()
{
	// Extent along the axis
	xPos =  1;
	xNeg = -1;
	yPos =  1;
	yNeg = -1;
}


void NullYS2D::getGradient(double &gx, double &gy, double x, double y)
{
    	opserr << "ERROR - NullYS2D::getGradient function should not be called\n";
		gx = 1;
		gy = 1;
}

double NullYS2D::getSurfaceDrift(double x, double y)
{
	return -1;
}

YieldSurface_BC *NullYS2D::getCopy(void)
{
    NullYS2D *theCopy = new NullYS2D(this->getTag());
    return theCopy;
}

int NullYS2D::displaySelf(Renderer &theViewer, int displayMode, float fact)
{
	return 0;
}

void NullYS2D::Print(OPS_Stream &s, int flag)
{
    s << "\nYield Surface No: " << this->getTag() << " type: NullYS2D\n";
}

