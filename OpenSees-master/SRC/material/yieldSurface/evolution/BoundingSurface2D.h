// BoundingSurface2D.h
//
//////////////////////////////////////////////////////////////////////

#if !defined BOUNDINGSURFACE2D_H
#define BOUNDINGSURFACE2D_H

#include "YS_Evolution2D.h"

class BoundingSurface2D : public YS_Evolution2D
{
public:
    BoundingSurface2D(int tag, int classTag, double min_iso_factor, 
				double iso_ratio, double kin_ratio,
				PlasticHardeningMaterial &kpx,
				PlasticHardeningMaterial &kpy,
				YieldSurface_BC &bound_surface);

	virtual ~BoundingSurface2D();
	virtual int 	commitState();
	virtual int		revertToLastCommit(void);
	virtual const   Vector &getEquiPlasticStiffness(void)=0;

	virtual void	Print(OPS_Stream &s, int flag =0);
	virtual YS_Evolution *getCopy(void)=0;

protected:
	virtual void	setTrialPlasticStrains(double ep, const Vector &f, const Vector &g)=0;
	virtual double getIsoPlasticStiffness(int dir)=0;
	virtual double getKinPlasticStiffness(int dir)=0;
	virtual Vector& getEvolDirection(Vector &f_new)=0;


protected:
//	double minIsoFactor;
	PlasticHardeningMaterial  *kpMatX, *kpMatY;
	YieldSurface_BC   *boundSurface;
};

#endif
