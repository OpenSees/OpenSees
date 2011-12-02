// Hajjar2D.h: interface for the Hajjar2D class.
//
//////////////////////////////////////////////////////////////////////

#if !defined HAJJAR2D_H
#define HAJJAR2D_H

#include "YieldSurface_BC2D.h"
#include <UniaxialMaterial.h>

class Hajjar2D : public YieldSurface_BC2D
{
public:
    Hajjar2D( int tag, double xmax, double ymax, YS_Evolution &model,
			  double centroid_y, double c1, double c2, double c3);

    Hajjar2D( int tag, YS_Evolution &model,
			  double D, double b, double t, double fc_, double fy_);
	virtual ~Hajjar2D();

	virtual YieldSurface_BC *getCopy(void);
	virtual int		displaySelf(Renderer &theViewer, int displayMode, float fact);

	virtual void	Print(OPS_Stream &s, int flag =0);

//protected:
//  For the following 2 methods, x, y already non-dimensionalized
    virtual void 	getGradient(double &gx, double &gy, double x, double y);
    virtual double 	getSurfaceDrift(double x, double y);
    virtual void	setExtent();

	double depth, width, thick, fc, fy, centroidY;
	double c1, c2, c3;
};

#endif
