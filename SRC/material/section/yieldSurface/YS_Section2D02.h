// @ rkaul@stanford.edu
// @ ggd@stanford.edu

#include <YieldSurfaceSection2d.h>

#ifndef YS_Section2D02_h
#define YS_Section2D02_h

class YS_Section2D02: public YieldSurfaceSection2d
{
  public:
	YS_Section2D02 ( int tag,
					 double E, double A, double I,
					 double theta_p_max,
					 YieldSurface_BC *ptrys,
					 bool use_kr=true);
    YS_Section2D02 (void);
	~YS_Section2D02 (void);

	const Matrix &getInitialTangent(void);  
	SectionForceDeformation *getCopy (void);
	void Print (OPS_Stream &s, int flag =0);
    int commitState (void);

  protected:
	void getSectionStiffness(Matrix &Ks);

  private:
    double E, A, I;
    double maxPlstkRot, peakPlstkRot, iFactor;
};

#endif
