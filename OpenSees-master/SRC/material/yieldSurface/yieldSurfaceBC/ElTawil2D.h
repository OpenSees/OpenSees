// ElTawil2D.h: interface for the ElTawil class.
//
//////////////////////////////////////////////////////////////////////

#if !defined ELTAWIL2D_H
#define ELTAWIL2D_H

#include "YieldSurface_BC2D.h"

class ElTawil2D : public YieldSurface_BC2D
{
public:
    ElTawil2D(int tag, double xbal, double ybal, double ypos, double yneg,
				YS_Evolution &model, double cz=1.6, double ty=1.9);

	virtual ~ElTawil2D();

	virtual YieldSurface_BC *getCopy(void);
	virtual int		displaySelf(Renderer &theViewer, int displayMode, float fact);
	virtual void	Print(OPS_Stream &s, int flag =0);

//protected:
//  For the following 2 methods, x, y already non-dimensionalized
    virtual void 	getGradient(double &gx, double &gy, double x, double y);
    virtual double 	getSurfaceDrift(double x, double y);
    virtual void	setExtent();
	virtual void	customizeInterpolate(double &xi, double &yi, double &xj, double &yj);
protected:

	double xBal, yBal;
	double yPosCap, yNegCap;
	double yPosCap_orig, yNegCap_orig;
	double cz, ty;
	double ytPos, ytNeg;
	double xtPos, xtNeg;
	double qy;

};

#endif
