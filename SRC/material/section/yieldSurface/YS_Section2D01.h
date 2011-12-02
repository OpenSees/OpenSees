// @ rkaul@stanford.edu
// @ ggd@stanford.edu

#include <YieldSurfaceSection2d.h>

#ifndef YS_Section2D01_h
#define YS_Section2D01_h

class YS_Section2D01: public YieldSurfaceSection2d
{
  public:
	YS_Section2D01 ( int tag, double E, double A, double I,
					 YieldSurface_BC *ptrys,
					 bool use_kr=true);
    YS_Section2D01 (void);
	~YS_Section2D01 (void);
	
	const Matrix &getInitialTangent(void);
	SectionForceDeformation *getCopy (void);
	void Print (OPS_Stream &s, int flag =0);

  protected:
	void getSectionStiffness(Matrix &Ks);

  private:
    double E, A, I;
};

#endif
