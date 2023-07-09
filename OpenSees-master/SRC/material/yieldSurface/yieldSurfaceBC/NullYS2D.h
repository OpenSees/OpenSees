// NullYS2D.h: interface for the NullYS2D class.
//
//////////////////////////////////////////////////////////////////////

#if !defined NULLYS2D_H
#define NULLYS2D_H

#include "YieldSurface_BC2D.h"

class NullEvolution;
class NullYS2D : public YieldSurface_BC2D
{
public:
    NullYS2D(int tag);
	virtual ~NullYS2D();

	virtual YieldSurface_BC *getCopy(void);
	virtual int		displaySelf(Renderer &theViewer, int displayMode, float fact);
	virtual void	Print(OPS_Stream &s, int flag =0);

//protected:
//  For the following 2 methods, x, y already non-dimensionalized
    virtual void 	getGradient(double &gx, double &gy, double x, double y);
    virtual double 	getSurfaceDrift(double x, double y);
    virtual void	setExtent();

private:
 static NullEvolution evolModel;    
};

#endif
