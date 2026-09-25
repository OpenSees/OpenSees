/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
** ****************************************************************** */

// Helpers shared by the elastic beam-column elements.

#ifndef ElasticBeamCommon_h
#define ElasticBeamCommon_h

// Ratio of bending to shear flexibility, phi = 12EI/(G Av L^2).  Zero when the
// section has no shear flexibility, which makes every formula below collapse
// to its Euler-Bernoulli form.
inline double
elasticBeamShearFactor(double E, double I, double G, double Av, double L)
{
  if (G <= 0.0 || Av <= 0.0)
    return 0.0;

  return 12.0*E*I/(G*Av*L*L);
}

// Turns the Euler-Bernoulli fixed-end moments of a transverse load into the
// shear-deformable ones and applies an end release, in place.
//
// The moment map is M_T = k_T * k_EB^-1 * M_EB, which is exact for an
// arbitrary load, not just a uniform one: in the simply supported primary
// system the shear strain integrates to M(L) - M(0) = 0, so shear leaves the
// primary end rotations alone and only the stiffness changes.  Condensing out
// a hinged end then uses the Timoshenko carry-over factor (2 - phi)/(4 + phi).
//
// Both maps are linear and odd, so callers may pass whichever sign convention
// their fixed-end force expressions use, as long as MI is the node I end.
inline void
elasticBeamEndMoments(double &MI, double &MJ, double phi, int release)
{
  double MIe = MI;
  double MJe = MJ;
  MI = (MIe + 0.5*phi*(MIe - MJe))/(1.0 + phi);
  MJ = (MJe + 0.5*phi*(MJe - MIe))/(1.0 + phi);

  double carryOver = (2.0 - phi)/(4.0 + phi);

  switch (release) {
  case 1:
    MJ -= carryOver*MI;
    MI = 0.0;
    break;
  case 2:
    MI -= carryOver*MJ;
    MJ = 0.0;
    break;
  case 3:
    MI = 0.0;
    MJ = 0.0;
    break;
  default:
    break;
  }
}

// Fixed-end forces, in the basic system, of the temperature change of a
// Beam2dTempLoad, which is linear through the depth d and linear along the
// length between its node I and node J values.
//
// Taking the centroid at mid-depth, the free thermal strain is an axial part
// alpha*Tc(x), Tc = (Ttop + Tbot)/2, plus a curvature
// k0(x) = -alpha*(Ttop - Tbot)/d, a hotter top bending the member concave
// down.  With no transverse load the moment is linear in x, and so is k0, so
// M(x) = -EI*k0(x) cancels the curvature everywhere and the clamped member
// stays straight, meeting its end conditions.  The solution being unique,
// this is exact.  With M(x) = MI*(x/L - 1) + MJ*x/L it gives MI = EI*k0(0)
// and MJ = -EI*k0(L), so each end moment depends on the gradient at that end
// only.  The axial force that cancels the mean free strain is
// N = -EA*alpha*(TcI + TcJ)/2.
//
// These moments ignore shear flexibility, which elasticBeamEndMoments then
// adds exactly: the primary system carries no moment under this load, hence
// no shear.  Nor does the load apply any external force, so it adds no
// reactions in the primary system; the shear (MI + MJ)/L of the clamped
// member is derived from the end moments by the coordinate transformation.
//
// Without a depth the gradient cannot be resolved, and the member takes the
// axial force only.
inline void
elasticBeamTemperatureForces(double EA, double EI, double alpha, double d,
                             double TtopI, double TbotI,
                             double TtopJ, double TbotJ,
                             double &N, double &MI, double &MJ)
{
  N = -0.25*EA*alpha*(TtopI + TbotI + TtopJ + TbotJ);

  MI = 0.0;
  MJ = 0.0;
  if (d > 0.0) {
    MI = -EI*alpha*(TtopI - TbotI)/d;
    MJ = EI*alpha*(TtopJ - TbotJ)/d;
  }
}

#endif
