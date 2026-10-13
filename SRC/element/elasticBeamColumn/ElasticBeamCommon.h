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

#endif
