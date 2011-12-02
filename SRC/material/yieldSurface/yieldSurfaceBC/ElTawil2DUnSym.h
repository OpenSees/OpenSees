// ElTawil2DUnSym.h: interface for the ElTawil class.
//
//////////////////////////////////////////////////////////////////////

#if !defined ELTAWIL2DUNSYM_H
#define ELTAWIL2DUNSYM_H

#include "YieldSurface_BC2D.h"

class ElTawil2DUnSym : public YieldSurface_BC2D
{
public:
    ElTawil2DUnSym(int tag, double xPosBal, double yPosBal,
				double xNegBal, double yNegBal,
				double ypos, double yneg,
				YS_Evolution &model,
				double czPos=1.6, double tyPos=1.9,
				double czNeg=1.6, double tyNeg=1.9);

	virtual ~ElTawil2DUnSym();

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
	double xPosBal, yPosBal;
	double xNegBal, yNegBal;
	double yPosCap, yNegCap;
	double yPosCap_orig, yNegCap_orig;
	double czPos, tyPos;
	double czNeg, tyNeg;
	double ytPos, ytNeg;
	double xt1, xt2, xt3, xt4;
	double qy;
};

#endif

/*
             ytPos
         xt2 _  xt1
            /  \
           /    \
          /      )
       <------------>
          \     /
       xt3 \___/ xt4
            ytNeg



*/
