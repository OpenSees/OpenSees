//  YS_Evolution2D.h
//
//////////////////////////////////////////////////////////////////////

#if !defined YS_EVOLUTION2D_H
#define YS_EVOLUTION2D_H

#include "YS_Evolution.h"
#include "PlasticHardeningMaterial.h"

class YS_Evolution2D : public YS_Evolution
{
public:
    YS_Evolution2D(int tag, int classTag,
						double min_iso_factor,
                        double iso_ratio,
						double kin_ratio);

	virtual ~YS_Evolution2D();

//  Methods inherited
	virtual void	Print(OPS_Stream &s, int flag =0) =0;
	virtual YS_Evolution *getCopy(void) = 0;
    virtual Response *setResponse(char **argv, int argc, OPS_Stream &output);
	virtual int 	getResponse(int responseID, Information &info);
	virtual int	 displaySelf(Renderer &theViewer, int displayMode, float fact);
	
	virtual int	update(int flag);

    virtual int 	sendSelf(int commitTag, Channel &theChannel){return -1;}
    virtual int 	recvSelf(int commitTag, Channel &theChannel,
					FEM_ObjectBroker &theBroker){return -1;}

	virtual int 	commitState();
	virtual int		revertToLastCommit(void);
	virtual const   Vector &getEquiPlasticStiffness(void)=0;

	// no checks are performed
	virtual int 	evolveSurface(YieldSurface_BC *ys, double magPlasticDefo, 
								  Vector &G, Vector &F_Surface, int flag=0);

protected:
	virtual void	setTrialPlasticStrains(double ep, const Vector &f, const Vector &g)=0;
	virtual double getIsoPlasticStiffness(int dir)=0;
	virtual double getKinPlasticStiffness(int dir)=0;
	virtual Vector& getEvolDirection(Vector &f_new)=0;

protected:
//	double sumPlasticDeformX, sumPlasticDeformX_hist;
//	double sumPlasticDeformY, sumPlasticDeformY_hist;
	bool   softening;
	static Vector v2;
	double minIsoFactor;
	YieldSurface_BC *tmpYSPtr;
};

#endif

