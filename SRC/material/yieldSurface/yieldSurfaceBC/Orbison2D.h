// Orbison2D.h: interface for the Orbison2D class.
//
//////////////////////////////////////////////////////////////////////

#if !defined ORBISON2D_H
#define ORBISON2D_H

#include "YieldSurface_BC2D.h"
#include <UniaxialMaterial.h>

class Orbison2D : public YieldSurface_BC2D
{
public:
    Orbison2D(int tag, double xmax, double ymax, YS_Evolution &model);
			
	virtual ~Orbison2D();

	virtual YieldSurface_BC *getCopy(void);
	virtual int		displaySelf(Renderer &theViewer, int displayMode, float fact);

	virtual void	Print(OPS_Stream &s, int flag =0);

//protected:
//  For the following 2 methods, x, y already non-dimensionalized
    virtual void 	getGradient(double &gx, double &gy, double x, double y);
    virtual double 	getSurfaceDrift(double x, double y);
    virtual void	setExtent();

};

#endif
