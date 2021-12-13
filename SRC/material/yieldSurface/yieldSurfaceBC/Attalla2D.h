// Attalla2D.h: interface for the Attalla2D class.
//
//////////////////////////////////////////////////////////////////////

#if !defined ATTALLA2D_H
#define ATTALLA2D_H

#include "YieldSurface_BC2D.h"

class Attalla2D : public YieldSurface_BC2D
{
public:
    Attalla2D(int tag, double xmax, double ymax, YS_Evolution &model,
              double a01=0.19,  double a02=0.54, double a03=-1.4,
              double a04=-1.64, double a05=2.21, double a06=2.10);
	virtual ~Attalla2D();

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
    double a1, a2, a3, a4, a5, a6;
	int    driftAlgo;
};

#endif
