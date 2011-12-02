// @ rkaul@stanford.edu
// @ ggd@stanford.edu

#include <math.h>

#include <YieldSurfaceSection2d.h>
#include <Matrix.h>
#include <Vector.h>
#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <MatrixUtil.h>
#include <stdlib.h>

#include <classTags.h>

#define SEC_TAG_YS2d -1

ID      YieldSurfaceSection2d::code(2);
Vector  YieldSurfaceSection2d::dele(2);
Vector  YieldSurfaceSection2d::surfaceForce(2);
Matrix  YieldSurfaceSection2d::G(2,1);
Matrix  YieldSurfaceSection2d::Ktp(2,2);

YieldSurfaceSection2d::YieldSurfaceSection2d(void)
  :SectionForceDeformation(0, SEC_TAG_YS2d),
   use_Kr_orig(true), ys(0), e(2), s(2),
   eCommit(2), sCommit(2), ks(2,2),
   use_Kr(true), split_step(false)
{
  code(0) = SECTION_RESPONSE_P;	// P is the first quantity
  code(1) = SECTION_RESPONSE_MZ;  // Mz is the second
}

YieldSurfaceSection2d::YieldSurfaceSection2d
(int tag, int classtag, YieldSurface_BC *ptrys, bool use_kr)
  :SectionForceDeformation(tag, classtag),
   use_Kr_orig(use_kr), ys(0), e(2), s(2),
   eCommit(2), sCommit(2), ks(2,2),
   use_Kr(use_kr), split_step(false)
{
  code(0) = SECTION_RESPONSE_P;	// P is the first quantity
  code(1) = SECTION_RESPONSE_MZ;	// Mz is the second
  
  if(ptrys==0)
    {
      opserr << "WARNING - InelasticYS2DGNL(): ys1 = 0" << endln;
    }
  else
    {
      ys =  ptrys->getCopy();
      ys->setTransformation(1, 0, 1, -1);   // x-axis is Mz, y-axis is P
      //		ys->setEleInfo(getTag(), 1);
    }
}


YieldSurfaceSection2d::~YieldSurfaceSection2d(void)
{
  if (ys != 0)
    delete ys;
}

int
YieldSurfaceSection2d::commitState(void)
{
  eCommit = e;
  sCommit = s;
  ys->commitState(s);
  split_step = false;
  
  return 0;
}

int
YieldSurfaceSection2d::revertToLastCommit(void)
{
  e = eCommit;
  s = sCommit;
  ys->revertToLastCommit();
  
  return 0;
}

int
YieldSurfaceSection2d::revertToStart(void)
{
  eCommit.Zero();
  sCommit.Zero();
  //	ys->revertToStart();
  
  return 0;
}

int
YieldSurfaceSection2d::setTrialSectionDeformation (const Vector &def)
{
  ys->update(); // important
  use_Kr = use_Kr_orig;
  //	split_step = false;
  //  once determined, leave it till convergence
  
  e = def;
  dele = e - eCommit;

  getSectionStiffness(ks);
  double EA = ks(0,0);
  double EI = ks(1,1);
  
  s(0) = sCommit(0) + EA*dele(0);
  s(1) = sCommit(1) + EI*dele(1);
  
  if(ys->getTrialForceLocation(s) <= 0)
    return 0;
  
  
  // case: if it shoots through
  //         use dF return to surface
  int driftOld = ys->getCommitForceLocation();
  
  if(driftOld < 0) // Inside
    {
      use_Kr = false;
      split_step = true;
      
      surfaceForce = s;
      ys->setToSurface(surfaceForce, ys->dFReturn);  //dFReturn, ConstantYReturn, RadialReturn
      ys->getTrialGradient(G, surfaceForce);
    }
  // Now we know that force point has drifted from the surface
  else if(driftOld == 0) // On surface
    {
      ys->getCommitGradient(G);
      surfaceForce =  sCommit;
    }
  else // driftOld outside - bug or bad constraints or continued from not converged state
    {
      opserr << "WARNING: YieldSurfaceSection2d::setTrialSectionDeformation, driftOld outside [" << getTag()<<"]\n";
    }
  
  double dF_in0 = s(0) - surfaceForce(0);
  double dF_in1 = s(1) - surfaceForce(1);
  
  double g0 = G(0,0);
  double g1 = G(1,0);
  
  Ktp(0,0) = EA;
  Ktp(1,1) = EI;
  ys->addPlasticStiffness(Ktp);
  
  double inv = 1/(Ktp(0,0)*g0*g0 + Ktp(1,1)*g1*g1);
  
  double lamda = inv*(g0*dF_in0 + g1*dF_in1);
  if(fabs(lamda) < 1e-8) lamda = 0.0; // to get rid of -1e-15 etc
  
  if(lamda < 0)
    {
      use_Kr = false;
      lamda = 0.0;
    }
  
  int grow = ys->modifySurface(lamda, surfaceForce, G);
  
  // used to do: (not tested shrinking yet)
  //	if(grow < 0)
  //		algo = ys->ConstantYReturn;
  //	else
  //		algo = algo_orig;
  
  if(use_Kr)
    {
      ks(0,0) = EA - EA*EA*g0*g0*inv;
      ks(0,1) = -1*EA*g0*g1*EI*inv;
      ks(1,0) = ks(0,1);
      ks(1,1) = EI - EI*EI*g1*g1*inv;
    }
  if(split_step)
    {
      s(0) = s(0) - EA*g0*lamda;
      s(1) = s(1) - EI*g1*lamda;
    }
  else
    {
      if(use_Kr)
	s = surfaceForce + ks*dele;
    }
  
  ys->setToSurface(s, ys->ConstantYReturn);
  // used to do centroid return
  // then force-balance using ConstantYReturn
  // comp/tension issue: always use constantP
  
  return 0;
}

const Vector &
YieldSurfaceSection2d::getSectionDeformation (void)
{
  return e;
}

const Vector &
YieldSurfaceSection2d::getStressResultant (void)
{
  return s;
}

const Matrix &
YieldSurfaceSection2d::getSectionTangent(void)
{
  return ks;
}

const Matrix &
YieldSurfaceSection2d::getSectionFlexibility (void)
{
  static Matrix fs(2,2);
  invert2by2Matrix(ks, fs);
  
  return fs;
}

const ID&
YieldSurfaceSection2d::getType ()
{
  return code;
}

int
YieldSurfaceSection2d::getOrder () const
{
  return 2;
}

int
YieldSurfaceSection2d::sendSelf(int commitTag, Channel &theChannel)
{
  return -1;
}

int
YieldSurfaceSection2d::recvSelf(int commitTag, Channel &theChannel,
				FEM_ObjectBroker &theBroker)
{
  return -1;
}

void
YieldSurfaceSection2d::Print(OPS_Stream &s, int flag)
{
  s << "YieldSurfaceSection2d, tag: " << this->getTag() << endln;
  s << "\tYield Surface:" << *ys << endln;
  s << "\tSection Force:" << sCommit;
  s << "\tSection Defom:" << eCommit;
}
