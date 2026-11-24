/* ****************************************************************** **
**    OpenSees - Open System for Earthquake Engineering Simulation    **
**          Pacific Earthquake Engineering Research Center            **
**                                                                    **
**                                                                    **
** (C) Copyright 1999, The Regents of the University of California    **
** All Rights Reserved.                                               **
**                                                                    **
** Commercial use of this program without express permission of the   **
** University of California, Berkeley, is strictly prohibited.  See   **
** file 'COPYRIGHT'  in main directory for information on usage and   **
** redistribution,  and for a DISCLAIMER OF ALL WARRANTIES.           **
**                                                                    **
** Developed by:                                                      **
**   Frank McKenna (fmckenna@ce.berkeley.edu)                         **
**   Gregory lo. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// Purpose: implementation of Corot3DBattini, a 3d corotational coordinate
// transformation with the Battini & Pacoste (2002) element frame, true
// logarithm local rotations and an exact consistent tangent. See
// Corot3DBattini.h for a description.
//
// Notation follows Battini (2002), chapter 4 [B]; the element ends 1 and 2
// of [B] are written I and J here (e.g. q1 -> qI, theta_bar_2 -> thetaBarJ).
// All vectors are in global components unless noted.
//   global dofs d = [uI(3), dwI(3), uJ(3), dwJ(3)], dw = spatial spin
//   xJI = chord vector, ln = |xJI|, lo = undeformed length, r1 = xJI/ln
//   RtI, RtJ = orientation of the nodal triads (R_g1 R_o, R_g2 R_o in [B])
//   qI, qJ = second columns of RtI, RtJ;  q = (qI + qJ)/2
//   r3 = r1 x q/|r1 x q|,  r2 = r3 x r1,  Rr = [r1 r2 r3]
//   frame spin  deltaEpsR = GTg d = r1 x dr1 + deltaEpsR1 r1,
//               deltaEpsR1 = (r3.dq - qe1 (r3.dr1)) / qe2,  qe1 = q.r1, qe2 = q.r2
//   local rotations thetaBarI = log(Rr^T RtI), and
//               dthetaBarI = T_s^-1(thetaBarI) Rr^T (dwI - deltaEpsR)
//   with T_s^-1(t) = I - S(t)/2 + etaTs(|t|) S(t)^2
//
// Internal force  pg = B^T pb, B = d(ub)/d(d)
// Tangent         K  = B^T kb B + Kgeo, Kgeo = d(B^T pb)/d(d) at fixed pb,
//                 assembled one column at a time from directional derivatives.
//
// References:
//   [B]  Battini, J.-M. (2002) "Co-rotational beam elements in instability
//        problems", Doctoral thesis, Dept. of Mechanics, KTH, Stockholm.
//        Chapter 4 (3D co-rotational beam element). Equation numbers in the
//        comments below, e.g. [B 4.27], refer to this thesis.
//   [BP] Battini, J.-M. and Pacoste, C. (2002) "Co-rotational beam elements
//        with warping effects in instability problems", CMAME 191:1755-1789.
//   [PE] Pacoste, C. and Eriksson, A. (1997) "Beam elements in instability
//        problems", CMAME 144:163-197.
//
// Correspondence with [B] (notation there -> here):
//   r1, r2, r3, R_r                     -> r1, r2, r3, Rr          [B 4.24-4.28]
//   l_o, l_n                            -> lo, ln                  [B 4.26, 4.29]
//   R_o                                 -> Ro                      [B 4.27]
//   R_gi R_o (nodal triad t_i)          -> RtI, RtJ                [B Fig. 4.2]
//   q_i = R_gi R_o [0 1 0]^T, q         -> qI, qJ, q               [B 4.27]
//   R_bar_i = R_r^T R_gi R_o            -> RbarI, RbarJ            [B 4.30]
//   theta_bar_i = log(R_bar_i)          -> thetaBarI, thetaBarJ    [B 4.31]
//   T_s^-1(theta_bar), eta, mu          -> TsInvVec, etaTs, dEtaTsOverTheta [B 4.14, 4.37, 4.46]
//   K_h                                 -> dTsInvTVec              [B 4.43-4.45]
//   local q_1, q_2 of R_r^T q; eta      -> qe1, qe2; eta           [B 4.67, 4.70]
//   delta eps_r (frame spin), G^T       -> deltaEpsR, GTg          [B 4.58-4.73]
//   delta w_bar_i, delta theta_bar_i    -> dwBarI, dthetaBarI      [B 4.54, 4.37]
//   n, m_1, m_2 (local forces)          -> n, mI, mJ               [B 4.33]
//   B = B_a B_g                         -> B                       [B 4.39, 4.47, 4.57-4.59]
//   f_g = B_g^T f_a                     -> pg = B^T pb             [B 4.74]
//   K_g = B_g^T K_a B_g + K_m           -> B^T kb B + Kgeo         [B 4.75-4.94]
//   eccentricity v_o, R_g v_o (links)   -> voI/voJ, vI/vJ; A, K_off [B 4.2.4, 4.95-]
// The basic system here is the usual OpenSees one (axial, theta_z i/j,
// theta_y i/j, twist) rather than the 7 local dofs of [B 4.32]; the tangent
// is the exact derivative of the internal force and is equivalent to [B 4.94];
// it is assembled column by column from analytic directional derivatives
// instead of from the closed-form K_m expression.

#include <math.h>
#include <string.h>
#include <Vector.h>
#include <Matrix.h>
#include <Node.h>
#include <Channel.h>
#include <elementAPI.h>
#include <string>
#include <Corot3DBattini.h>

Matrix Corot3DBattini::Kg(12,12);
Vector Corot3DBattini::Pg(12);
Matrix Corot3DBattini::Tlg(12,12);

//
// small 3-vector / 3x3 helpers
//
namespace {

inline double dot3(const double *a, const double *b) {return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}

inline void cross3(const double *a, const double *b, double *c)
{
  double c0 = a[1]*b[2]-a[2]*b[1];
  double c1 = a[2]*b[0]-a[0]*b[2];
  double c2 = a[0]*b[1]-a[1]*b[0];
  c[0] = c0; c[1] = c1; c[2] = c2;
}

inline double norm3(const double *a) {return sqrt(dot3(a,a));}

// y = M^T x
inline void matTvec3(const double M[3][3], const double *x, double *y)
{
  double y0 = M[0][0]*x[0]+M[1][0]*x[1]+M[2][0]*x[2];
  double y1 = M[0][1]*x[0]+M[1][1]*x[1]+M[2][1]*x[2];
  double y2 = M[0][2]*x[0]+M[1][2]*x[1]+M[2][2]*x[2];
  y[0] = y0; y[1] = y1; y[2] = y2;
}

// y = M x
inline void matvec3(const double M[3][3], const double *x, double *y)
{
  double y0 = M[0][0]*x[0]+M[0][1]*x[1]+M[0][2]*x[2];
  double y1 = M[1][0]*x[0]+M[1][1]*x[1]+M[1][2]*x[2];
  double y2 = M[2][0]*x[0]+M[2][1]*x[1]+M[2][2]*x[2];
  y[0] = y0; y[1] = y1; y[2] = y2;
}

// quaternions (w, x, y, z), Hamilton product
void quatFromRotVec(const double *v, double *qt)
{
  double t = norm3(v);
  if (t < 1.0e-14) {
    qt[0] = 1.0; qt[1] = 0.5*v[0]; qt[2] = 0.5*v[1]; qt[3] = 0.5*v[2];
  } else {
    double s = sin(0.5*t)/t;
    qt[0] = cos(0.5*t); qt[1] = s*v[0]; qt[2] = s*v[1]; qt[3] = s*v[2];
  }
  double n = sqrt(qt[0]*qt[0]+qt[1]*qt[1]+qt[2]*qt[2]+qt[3]*qt[3]);
  for (int i=0; i<4; i++) qt[i] /= n;
}

void quatProduct(const double *a, const double *b, double *c)
{
  double w = a[0]*b[0] - a[1]*b[1] - a[2]*b[2] - a[3]*b[3];
  double x = a[0]*b[1] + b[0]*a[1] + a[2]*b[3] - a[3]*b[2];
  double y = a[0]*b[2] + b[0]*a[2] + a[3]*b[1] - a[1]*b[3];
  double z = a[0]*b[3] + b[0]*a[3] + a[1]*b[2] - a[2]*b[1];
  c[0] = w; c[1] = x; c[2] = y; c[3] = z;
}

void quatNormalize(double *qt)
{
  double n = sqrt(qt[0]*qt[0]+qt[1]*qt[1]+qt[2]*qt[2]+qt[3]*qt[3]);
  for (int i=0; i<4; i++) qt[i] /= n;
}

void matFromQuat(const double *qt, double R[3][3])
{
  double w = qt[0], x = qt[1], y = qt[2], z = qt[3];
  R[0][0] = 1 - 2*(y*y + z*z); R[0][1] = 2*(x*y - w*z);     R[0][2] = 2*(x*z + w*y);
  R[1][0] = 2*(x*y + w*z);     R[1][1] = 1 - 2*(x*x + z*z); R[1][2] = 2*(y*z - w*x);
  R[2][0] = 2*(x*z - w*y);     R[2][1] = 2*(y*z + w*x);     R[2][2] = 1 - 2*(x*x + y*y);
}

// Shepperd's method; returns unit quaternion with w >= 0
void quatFromMat(const double R[3][3], double *qt)
{
  double tr = R[0][0] + R[1][1] + R[2][2];
  if (tr >= R[0][0] && tr >= R[1][1] && tr >= R[2][2]) {
    double s = 2.0*sqrt(1.0 + tr);
    qt[0] = 0.25*s;
    qt[1] = (R[2][1] - R[1][2])/s;
    qt[2] = (R[0][2] - R[2][0])/s;
    qt[3] = (R[1][0] - R[0][1])/s;
  } else {
    int i = 0;
    if (R[1][1] > R[i][i]) i = 1;
    if (R[2][2] > R[i][i]) i = 2;
    int j = (i+1)%3, k = (i+2)%3;
    double s = 2.0*sqrt(1.0 + R[i][i] - R[j][j] - R[k][k]);
    double v[3];
    v[i] = 0.25*s;
    v[j] = (R[j][i] + R[i][j])/s;
    v[k] = (R[k][i] + R[i][k])/s;
    qt[0] = (R[k][j] - R[j][k])/s;
    qt[1] = v[0]; qt[2] = v[1]; qt[3] = v[2];
  }
  if (qt[0] < 0.0)
    for (int l=0; l<4; l++) qt[l] = -qt[l];
  quatNormalize(qt);
}

// rotation vector (angle in [0, pi]) of a rotation matrix
void logSO3(const double R[3][3], double *theta)
{
  double qt[4];
  quatFromMat(R, qt);
  double s = sqrt(qt[1]*qt[1]+qt[2]*qt[2]+qt[3]*qt[3]);
  if (s < 1.0e-14) {
    theta[0] = 2.0*qt[1]; theta[1] = 2.0*qt[2]; theta[2] = 2.0*qt[3];
    return;
  }
  double f = 2.0*atan2(s, qt[0])/s;
  theta[0] = f*qt[1]; theta[1] = f*qt[2]; theta[2] = f*qt[3];
}

// T_s^-1(t) = I - S(t)/2 + etaTs S(t)^2 [B 4.14], written with
// S(t)^2 = t t^T - |t|^2 I; etaTs is eta of [B 4.46] (not the eta of [B 4.70]):
//   etaTs = (2 sin t - t (1 + cos t)) / (2 t^2 sin t)
double etaTs(double t)
{
  if (t < 1.0e-4)
    return 1.0/12.0 + t*t/720.0;
  return 1.0/(t*t) - (1.0 + cos(t))/(2.0*t*sin(t));
}

// (d etaTs/dt)/t; plays the role of mu in [B 4.46]
double dEtaTsOverTheta(double t)
{
  if (t < 1.0e-3)
    return 1.0/360.0 + t*t/7560.0;
  double h = 0.5*t, sh = sin(h);
  double db = -2.0/(t*t*t) + 1.0/(4.0*t*sh*sh) + cos(h)/sh/(2.0*t*t);
  return db/t;
}

// y = T_s^-1(theta) x:  dtheta_bar = T_s^-1(theta_bar) dw_bar  [B 4.37, 4.39]
void TsInvVec(const double *theta, const double *x, double *y)
{
  double t = norm3(theta), b = etaTs(t);
  double tx[3], ttx[3];
  cross3(theta, x, tx);
  cross3(theta, tx, ttx);
  for (int i=0; i<3; i++)
    y[i] = x[i] - 0.5*tx[i] + b*ttx[i];
}

// y = T_s^-T(theta) m = m + theta x m / 2 + etaTs theta x (theta x m)
// (moments conjugate to the spins: f_a = B_a^T f_l  [B 4.40])
void TsInvTVec(const double *theta, const double *m, double *y)
{
  double t = norm3(theta), b = etaTs(t);
  double tm[3], ttm[3];
  cross3(theta, m, tm);
  cross3(theta, tm, ttm);
  for (int i=0; i<3; i++)
    y[i] = m[i] + 0.5*tm[i] + b*ttm[i];
}

// directional derivative of T_s^-T(theta) m in direction dtheta (m fixed):
// K_hi dtheta of [B 4.43-4.45] (with v = m_i)
void dTsInvTVec(const double *theta, const double *m, const double *dth, double *y)
{
  double t = norm3(theta), b = etaTs(t), dbt = dEtaTsOverTheta(t);
  double dm[3], tm[3], ttm[3], dtm[3], tdm[3];
  cross3(dth, m, dm);              // dth x m
  cross3(theta, m, tm);            // th x m
  cross3(theta, tm, ttm);          // th x (th x m)
  cross3(dth, tm, dtm);            // dth x (th x m)
  cross3(theta, dm, tdm);          // th x (dth x m)
  double c = dbt*dot3(theta, dth);
  for (int i=0; i<3; i++)
    y[i] = 0.5*dm[i] + c*ttm[i] + b*(dtm[i] + tdm[i]);
}

} // namespace


//
// interpreter command (python / OPS_ interface)
//   geomTransf CorotationalBattini tag vecxzX vecxzY vecxzZ
//
void* OPS_Corot3DBattini()
{
  if (OPS_GetNumRemainingInputArgs() < 4) {
    opserr << "insufficient arguments for Corot3DBattini\n";
    return 0;
  }
  int tag;
  int numData = 1;
  if (OPS_GetIntInput(&numData, &tag) < 0) return 0;

  Vector vec(3);
  numData = 3;
  if (OPS_GetDoubleInput(&numData, &vec(0)) < 0) return 0;

  Vector jntOffsetI(3), jntOffsetJ(3);
  bool symmetric = false;
  while (OPS_GetNumRemainingInputArgs() > 0) {
    std::string type = OPS_GetString();
    if (type == "-jntOffset") {
      if (OPS_GetNumRemainingInputArgs() < 6) {
        opserr << "WARNING Corot3DBattini - -jntOffset needs 6 values\n";
        return 0;
      }
      if (OPS_GetDoubleInput(&numData, &jntOffsetI(0)) < 0) return 0;
      if (OPS_GetDoubleInput(&numData, &jntOffsetJ(0)) < 0) return 0;
    } else if (type == "-symmetric") {
      symmetric = true;
    } else {
      opserr << "WARNING Corot3DBattini - unknown option " << type.c_str() << "\n";
      return 0;
    }
  }
  return new Corot3DBattini(tag, vec, jntOffsetI, jntOffsetJ, symmetric);
}


Corot3DBattini::Corot3DBattini(int tag, const Vector &vecInLocXZPlane,
                               const Vector &rigJntOffsetI,
                               const Vector &rigJntOffsetJ,
                               bool symmetric)
  : CrdTransf(tag, CRDTR_TAG_Corot3DBattini),
    nodeIPtr(0), nodeJPtr(0), lo(0.0), ln(0.0),
    ub(6), ubCommit(6), ubPrev(6),
    nodeIInitialDisp(0), nodeJInitialDisp(0), initialDispChecked(false),
    stateRestored(false), symmetricTangent(symmetric)
{
  if (vecInLocXZPlane.Size() != 3) {
    opserr << "Corot3DBattini::Corot3DBattini - vector that defines local xz plane is invalid; using (0,0,1)\n";
    vecxz[0] = 0.0; vecxz[1] = 0.0; vecxz[2] = 1.0;
  } else {
    for (int i=0; i<3; i++) vecxz[i] = vecInLocXZPlane(i);
  }

  // rigid joint offsets [B 4.2.4], global, from node to element end
  for (int i=0; i<3; i++) {
    voI[i] = (rigJntOffsetI.Size() == 3) ? rigJntOffsetI(i) : 0.0;
    voJ[i] = (rigJntOffsetJ.Size() == 3) ? rigJntOffsetJ(i) : 0.0;
  }
  if ((rigJntOffsetI.Size() != 0 && rigJntOffsetI.Size() != 3) ||
      (rigJntOffsetJ.Size() != 0 && rigJntOffsetJ.Size() != 3))
    opserr << "Corot3DBattini::Corot3DBattini - joint offset vectors must have size 3; using zero values\n";
  hasOffsets = (norm3(voI) != 0.0 || norm3(voJ) != 0.0);

  memset(Ro, 0, sizeof(Ro));
  memset(RtI, 0, sizeof(RtI));  memset(RtJ, 0, sizeof(RtJ));
  memset(Rr,  0, sizeof(Rr));
  memset(GTg, 0, sizeof(GTg));
  memset(B, 0, sizeof(B));  memset(B0, 0, sizeof(B0));
  memset(pEndP0, 0, sizeof(pEndP0));
  for (int i=0; i<3; i++) {
    vI[i] = voI[i];  vJ[i] = voJ[i];
    alphaI[i] = alphaJ[i] = 0.0;
    thetaBarI[i] = thetaBarJ[i] = 0.0;
    qI[i] = qJ[i] = q[i] = 0.0;
  }
  quatRtI[0] = quatRtJ[0] = quatRtICommit[0] = quatRtJCommit[0] = 1.0;
  for (int i=1; i<4; i++)
    quatRtI[i] = quatRtJ[i] = quatRtICommit[i] = quatRtJCommit[i] = 0.0;
}


Corot3DBattini::Corot3DBattini()
  : CrdTransf(0, CRDTR_TAG_Corot3DBattini),
    nodeIPtr(0), nodeJPtr(0), lo(0.0), ln(0.0),
    ub(6), ubCommit(6), ubPrev(6),
    nodeIInitialDisp(0), nodeJInitialDisp(0), initialDispChecked(false),
    stateRestored(false), symmetricTangent(false)
{
  vecxz[0] = 0.0; vecxz[1] = 0.0; vecxz[2] = 1.0;
  for (int i=0; i<3; i++) voI[i] = voJ[i] = 0.0;
  hasOffsets = false;
  memset(Ro, 0, sizeof(Ro));
  memset(RtI, 0, sizeof(RtI));  memset(RtJ, 0, sizeof(RtJ));
  memset(Rr,  0, sizeof(Rr));
  memset(GTg, 0, sizeof(GTg));
  memset(B, 0, sizeof(B));  memset(B0, 0, sizeof(B0));
  memset(pEndP0, 0, sizeof(pEndP0));
  for (int i=0; i<3; i++) {
    vI[i] = voI[i];  vJ[i] = voJ[i];
    alphaI[i] = alphaJ[i] = 0.0;
    thetaBarI[i] = thetaBarJ[i] = 0.0;
    qI[i] = qJ[i] = q[i] = 0.0;
  }
  quatRtI[0] = quatRtJ[0] = quatRtICommit[0] = quatRtJCommit[0] = 1.0;
  for (int i=1; i<4; i++)
    quatRtI[i] = quatRtJ[i] = quatRtICommit[i] = quatRtJCommit[i] = 0.0;
}


Corot3DBattini::~Corot3DBattini()
{
  if (nodeIInitialDisp != 0) delete [] nodeIInitialDisp;
  if (nodeJInitialDisp != 0) delete [] nodeJInitialDisp;
}


//
// geometry
//
int
Corot3DBattini::computeInitialGeometry(void)
{
  const Vector &crdI = nodeIPtr->getCrds();
  const Vector &crdJ = nodeJPtr->getCrds();

  double xJI[3];
  for (int i=0; i<3; i++) {
    xJI[i] = crdJ(i) + voJ[i] - crdI(i) - voI[i];             // between element ends
    if (nodeIInitialDisp != 0) xJI[i] -= nodeIInitialDisp[i];
    if (nodeJInitialDisp != 0) xJI[i] += nodeJInitialDisp[i];
  }

  lo = norm3(xJI);
  if (lo == 0.0) {
    opserr << "\nCorot3DBattini::initialize - transfTag = " << this->getTag()
           << ": element has zero length\n";
    return -2;
  }

  double x[3], y[3], z[3];
  for (int i=0; i<3; i++) x[i] = xJI[i]/lo;

  cross3(vecxz, x, y);
  double ynorm = norm3(y);
  if (ynorm == 0.0) {
    opserr << "\nCorot3DBattini::initialize - transfTag = " << this->getTag()
           << ": vector v that defines plane xz is parallel to x axis\n";
    return -3;
  }
  for (int i=0; i<3; i++) y[i] /= ynorm;
  cross3(x, y, z);

  for (int i=0; i<3; i++) {
    Ro[i][0] = x[i];
    Ro[i][1] = y[i];
    Ro[i][2] = z[i];
  }
  return 0;
}


void
Corot3DBattini::getNodalDisp(Node *node, const double *initDisp, double *d) const
{
  const Vector &disp = node->getTrialDisp();
  for (int i=0; i<6; i++) {
    d[i] = disp(i);
    if (initDisp != 0) d[i] -= initDisp[i];
  }
}


// Given the nodal quaternions and the current nodal displacements, compute
// the corotated frame, the local rotations, ub, GTg and B.
int
Corot3DBattini::computeState(void)
{
  double dI[6], dJ[6];
  this->getNodalDisp(nodeIPtr, nodeIInitialDisp, dI);
  this->getNodalDisp(nodeJPtr, nodeJInitialDisp, dJ);

  matFromQuat(quatRtI, RtI);
  matFromQuat(quatRtJ, RtJ);

  // current rigid offsets v = R_g vo, R_g = Rt Ro^T  [B 4.2.4]
  if (hasOffsets) {
    double t[3];
    matTvec3(Ro, voI, t);  matvec3(RtI, t, vI);
    matTvec3(Ro, voJ, t);  matvec3(RtJ, t, vJ);
  }

  // chord between the element ends: x_end = x_node + u_node + v
  const Vector &crdI = nodeIPtr->getCrds();
  const Vector &crdJ = nodeJPtr->getCrds();
  double xJI[3];
  for (int i=0; i<3; i++) {
    xJI[i] = crdJ(i) - crdI(i) + dJ[i] - dI[i] + vJ[i] - vI[i];
    if (nodeIInitialDisp != 0) xJI[i] -= nodeIInitialDisp[i];
    if (nodeJInitialDisp != 0) xJI[i] += nodeJInitialDisp[i];
  }
  ln = norm3(xJI);                                             // [B 4.26]
  if (ln == 0.0) {
    opserr << "\nCorot3DBattini::update - transfTag = " << this->getTag()
           << ": deformed length is zero\n";
    return -2;
  }

  // Battini frame: r1 = chord [B 4.25], q = (q1+q2)/2 [B 4.27],
  // r3 = r1 x q/|r1 x q|, r2 = r3 x r1 [B 4.28]
  double r1[3], r2[3], r3[3];
  for (int i=0; i<3; i++) {
    r1[i] = xJI[i]/ln;
    qI[i] = RtI[i][1];
    qJ[i] = RtJ[i][1];
    q[i] = 0.5*(qI[i] + qJ[i]);
  }
  cross3(r1, q, r3);
  double n3 = norm3(r3);
  if (n3 < 1.0e-12) {
    opserr << "\nCorot3DBattini::update - transfTag = " << this->getTag()
           << ": mean nodal y-axis is parallel to the element chord; frame undefined\n";
    return -3;
  }
  for (int i=0; i<3; i++) r3[i] /= n3;
  cross3(r3, r1, r2);
  for (int i=0; i<3; i++) {
    Rr[i][0] = r1[i];
    Rr[i][1] = r2[i];
    Rr[i][2] = r3[i];
  }

  // local nodal rotations: R_bar_i = Rr^T Rgi Ro [B 4.30], theta = log(R_bar) [B 4.31]
  // (RtI, RtJ here are the nodal triads Rgi Ro)
  double RbarI[3][3], RbarJ[3][3];
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++) {
      RbarI[i][j] = Rr[0][i]*RtI[0][j] + Rr[1][i]*RtI[1][j] + Rr[2][i]*RtI[2][j];
      RbarJ[i][j] = Rr[0][i]*RtJ[0][j] + Rr[1][i]*RtJ[1][j] + Rr[2][i]*RtJ[2][j];
    }
  logSO3(RbarI, thetaBarI);
  logSO3(RbarJ, thetaBarJ);

  // basic displacements; axial u_bar = ln - lo [B 4.29]
  ub(0) = ln - lo;
  ub(1) = thetaBarI[2];
  ub(2) = thetaBarJ[2];
  ub(3) = thetaBarI[1];
  ub(4) = thetaBarJ[1];
  ub(5) = thetaBarJ[0] - thetaBarI[0];

  // frame spin per global dof: deltaEpsR = r1 x dr1 + deltaEpsR1 r1. This is the matrix G
  // of [B 4.58-4.73] in global components: the bending part from
  // dr1 [B 4.63-4.65] and the twist deltaEpsR1 from dq [B 4.66-4.72], with
  // qe1/qe2 = eta [B 4.70].
  const double qe1 = dot3(q, r1);
  const double qe2 = dot3(q, r2);
  for (int j=0; j<12; j++) {
    double ddx[3] = {0,0,0}, dwI[3] = {0,0,0}, dwJ[3] = {0,0,0};
    if (j < 3)            ddx[j]   = -1.0;
    else if (j < 6)       dwI[j-3] =  1.0;
    else if (j < 9)       ddx[j-6] =  1.0;
    else                  dwJ[j-9] =  1.0;

    double dr1[3], c[3], dq[3], t1[3], t2[3];
    double e1ddx = dot3(r1, ddx);
    for (int i=0; i<3; i++) dr1[i] = (ddx[i] - r1[i]*e1ddx)/ln;
    cross3(dwI, qI, t1);
    cross3(dwJ, qJ, t2);
    for (int i=0; i<3; i++) dq[i] = 0.5*(t1[i] + t2[i]);
    double deltaEpsR1 = (dot3(r3, dq) - qe1*dot3(r3, dr1))/qe2;
    cross3(r1, dr1, c);
    for (int i=0; i<3; i++) GTg[i][j] = c[i] + deltaEpsR1*r1[i];
  }

  // B = d(ub)/d(d) = B_a B_g [B 4.39, 4.59]: dtheta_bar_i = T_s^-1 (dw_i - deltaEpsR)
  // in local components [B 4.37, 4.54, 4.57], axial row r [B 4.47]
  for (int j=0; j<12; j++) {
    double deltaEpsR[3] = {GTg[0][j], GTg[1][j], GTg[2][j]};
    double dwRelI[3], dwRelJ[3], dwBarI[3], dwBarJ[3], dthetaBarI[3], dthetaBarJ[3];
    for (int i=0; i<3; i++) {
      dwRelI[i] = ((j >= 3 && j < 6)  && (j-3 == i) ? 1.0 : 0.0) - deltaEpsR[i];
      dwRelJ[i] = ((j >= 9)           && (j-9 == i) ? 1.0 : 0.0) - deltaEpsR[i];
    }
    matTvec3(Rr, dwRelI, dwBarI);
    matTvec3(Rr, dwRelJ, dwBarJ);
    TsInvVec(thetaBarI, dwBarI, dthetaBarI);
    TsInvVec(thetaBarJ, dwBarJ, dthetaBarJ);

    double axial = 0.0;
    if (j < 3)              axial = -r1[j];
    else if (j >= 6 && j < 9) axial =  r1[j-6];

    B[0][j] = axial;
    B[1][j] = dthetaBarI[2];
    B[2][j] = dthetaBarJ[2];
    B[3][j] = dthetaBarI[1];
    B[4][j] = dthetaBarJ[1];
    B[5][j] = dthetaBarJ[0] - dthetaBarI[0];
  }
  return 0;
}


int
Corot3DBattini::initialize(Node *nodeIPointer, Node *nodeJPointer)
{
  nodeIPtr = nodeIPointer;
  nodeJPtr = nodeJPointer;

  if (nodeIPtr == 0 || nodeJPtr == 0) {
    opserr << "\nCorot3DBattini::initialize - invalid pointers to the element nodes\n";
    return -1;
  }

  // initial displacements are treated as stress free
  if (initialDispChecked == false) {
    const Vector &nodeIDisp = nodeIPtr->getDisp();
    const Vector &nodeJDisp = nodeJPtr->getDisp();
    for (int i=0; i<6; i++)
      if (nodeIDisp(i) != 0.0) {
        nodeIInitialDisp = new double [6];
        for (int j=0; j<6; j++) nodeIInitialDisp[j] = nodeIDisp(j);
        break;
      }
    for (int i=0; i<6; i++)
      if (nodeJDisp(i) != 0.0) {
        nodeJInitialDisp = new double [6];
        for (int j=0; j<6; j++) nodeJInitialDisp[j] = nodeJDisp(j);
        break;
      }
    initialDispChecked = true;
  }

  int error = this->computeInitialGeometry();
  if (error != 0)
    return error;

  // B in the undeformed (reference) configuration, for the initial
  // stiffness: Rr = Ro, nodal triads = Ro, theta = 0, ln = lo
  {
    double r1[3] = {Ro[0][0], Ro[1][0], Ro[2][0]};
    double y0[3] = {Ro[0][1], Ro[1][1], Ro[2][1]};
    double z0[3] = {Ro[0][2], Ro[1][2], Ro[2][2]};
    memset(B0, 0, sizeof(B0));
    for (int j=0; j<12; j++) {
      double ddx[3] = {0,0,0}, dwI[3] = {0,0,0}, dwJ[3] = {0,0,0};
      if (j < 3)            ddx[j]   = -1.0;
      else if (j < 6)       dwI[j-3] =  1.0;
      else if (j < 9)       ddx[j-6] =  1.0;
      else                  dwJ[j-9] =  1.0;
      // frame spin at the reference state: q = y0, qe1 = 0, qe2 = 1
      double dr1[3], c[3], t1[3], t2[3], dq[3], deltaEpsR[3];
      double e1ddx = dot3(r1, ddx);
      for (int i=0; i<3; i++) dr1[i] = (ddx[i] - r1[i]*e1ddx)/lo;
      cross3(dwI, y0, t1);
      cross3(dwJ, y0, t2);
      for (int i=0; i<3; i++) dq[i] = 0.5*(t1[i] + t2[i]);
      double deltaEpsR1 = dot3(z0, dq);
      cross3(r1, dr1, c);
      for (int i=0; i<3; i++) deltaEpsR[i] = c[i] + deltaEpsR1*r1[i];
      double dwRelI[3], dwRelJ[3], dthetaBarI[3], dthetaBarJ[3];
      for (int i=0; i<3; i++) {
        dwRelI[i] = (dwI[i]) - deltaEpsR[i];
        dwRelJ[i] = (dwJ[i]) - deltaEpsR[i];
      }
      matTvec3(Ro, dwRelI, dthetaBarI);    // theta = 0 -> T_s^-1 = I
      matTvec3(Ro, dwRelJ, dthetaBarJ);
      double axial = 0.0;
      if (j < 3)                axial = -r1[j];
      else if (j >= 6 && j < 9) axial =  r1[j-6];
      B0[0][j] = axial;
      B0[1][j] = dthetaBarI[2];
      B0[2][j] = dthetaBarJ[2];
      B0[3][j] = dthetaBarI[1];
      B0[4][j] = dthetaBarJ[1];
      B0[5][j] = dthetaBarJ[0] - dthetaBarI[0];
    }
  }

  if (stateRestored) {
    // state was restored by recvSelf(): keep it
    for (int i=0; i<4; i++) {quatRtI[i] = quatRtICommit[i]; quatRtJ[i] = quatRtJCommit[i];}
    ub = ubCommit;
    ubPrev = ubCommit;
    return 0;
  }

  quatFromMat(Ro, quatRtI);
  quatFromMat(Ro, quatRtJ);
  for (int i=0; i<3; i++) alphaI[i] = alphaJ[i] = 0.0;
  ub.Zero(); ubPrev.Zero();
  this->commitState();
  return this->computeState();
}


int
Corot3DBattini::update(void)
{
  double dI[6], dJ[6];
  this->getNodalDisp(nodeIPtr, nodeIInitialDisp, dI);
  this->getNodalDisp(nodeJPtr, nodeJInitialDisp, dJ);

  // rotational dof increments since the last update are spatial spins
  double dAI[3], dAJ[3];
  for (int i=0; i<3; i++) {
    dAI[i] = dI[i+3] - alphaI[i];
    dAJ[i] = dJ[i+3] - alphaJ[i];
    alphaI[i] = dI[i+3];
    alphaJ[i] = dJ[i+3];
  }
  double dqt[4], qtn[4];
  if (norm3(dAI) != 0.0) {
    quatFromRotVec(dAI, dqt);
    quatProduct(dqt, quatRtI, qtn);
    quatNormalize(qtn);
    for (int i=0; i<4; i++) quatRtI[i] = qtn[i];
  }
  if (norm3(dAJ) != 0.0) {
    quatFromRotVec(dAJ, dqt);
    quatProduct(dqt, quatRtJ, qtn);
    quatNormalize(qtn);
    for (int i=0; i<4; i++) quatRtJ[i] = qtn[i];
  }

  ubPrev = ub;
  return this->computeState();
}


double
Corot3DBattini::getInitialLength(void)
{
  return lo;
}


double
Corot3DBattini::getDeformedLength(void)
{
  return ln;
}


int
Corot3DBattini::commitState(void)
{
  for (int i=0; i<4; i++) {
    quatRtICommit[i] = quatRtI[i];
    quatRtJCommit[i] = quatRtJ[i];
  }
  ubCommit = ub;
  return 0;
}


int
Corot3DBattini::revertToLastCommit(void)
{
  for (int i=0; i<4; i++) {
    quatRtI[i] = quatRtICommit[i];
    quatRtJ[i] = quatRtJCommit[i];
  }

  // the nodes have been reverted: their trial rotations are the reference
  double dI[6], dJ[6];
  this->getNodalDisp(nodeIPtr, nodeIInitialDisp, dI);
  this->getNodalDisp(nodeJPtr, nodeJInitialDisp, dJ);
  for (int i=0; i<3; i++) {
    alphaI[i] = dI[i+3];
    alphaJ[i] = dJ[i+3];
  }

  ub = ubCommit;
  ubPrev = ubCommit;
  return this->computeState();
}


int
Corot3DBattini::revertToStart(void)
{
  quatFromMat(Ro, quatRtI);
  quatFromMat(Ro, quatRtJ);
  for (int i=0; i<3; i++) alphaI[i] = alphaJ[i] = 0.0;
  ub.Zero(); ubPrev.Zero(); ubCommit.Zero();
  for (int i=0; i<4; i++) {
    quatRtICommit[i] = quatRtI[i];
    quatRtJCommit[i] = quatRtJ[i];
  }
  return this->update();
}


const Vector &
Corot3DBattini::getBasicTrialDisp(void)
{
  return ub;
}


const Vector &
Corot3DBattini::getBasicIncrDisp(void)
{
  static Vector Dub(6);
  Dub = ub;
  Dub.addVector(1.0, ubCommit, -1.0);
  return Dub;
}


const Vector &
Corot3DBattini::getBasicIncrDeltaDisp(void)
{
  static Vector dub(6);
  dub = ub;
  dub.addVector(1.0, ubPrev, -1.0);
  return dub;
}


// node -> element-end rates for rigid offsets: u_end = u_node - v x w_node
void
Corot3DBattini::nodeToEnd(const Vector &rI, const Vector &rJ, double *re) const
{
  for (int i=0; i<6; i++) {
    re[i]   = rI(i);
    re[i+6] = rJ(i);
  }
  if (!hasOffsets) return;
  double t[3];
  cross3(vI, &re[3], t);  for (int i=0; i<3; i++) re[i]   -= t[i];
  cross3(vJ, &re[9], t);  for (int i=0; i<3; i++) re[6+i] -= t[i];
}


const Vector &
Corot3DBattini::getBasicTrialVel(void)
{
  // vb = B A vg (nodal angular velocity dofs = spatial angular velocities;
  // A maps node to element-end velocities for rigid offsets)
  static Vector vb(6);
  double ve[12];
  this->nodeToEnd(nodeIPtr->getTrialVel(), nodeJPtr->getTrialVel(), ve);
  for (int k=0; k<6; k++) {
    double s = 0.0;
    for (int j=0; j<12; j++)
      s += B[k][j]*ve[j];
    vb(k) = s;
  }
  return vb;
}


const Vector &
Corot3DBattini::getBasicTrialAccel(void)
{
  // ab ~= B A ag  (the d(BA)/dt vg term is neglected)
  static Vector ab(6);
  double ae[12];
  this->nodeToEnd(nodeIPtr->getTrialAccel(), nodeJPtr->getTrialAccel(), ae);
  for (int k=0; k<6; k++) {
    double s = 0.0;
    for (int j=0; j<12; j++)
      s += B[k][j]*ae[j];
    ab(k) = s;
  }
  return ab;
}


// Rigid offsets [B 4.2.4]: with A = [I -S(v); 0 I] at each end,
//   p_node = A^T p_end:  forces unchanged, moments m_node = m_end + v x f_end
void
Corot3DBattini::endToNodeVector(const double *pEnd, double *pNode) const
{
  for (int i=0; i<12; i++) pNode[i] = pEnd[i];
  if (!hasOffsets) return;
  double t[3];
  cross3(vI, &pEnd[0], t);  for (int i=0; i<3; i++) pNode[3+i] += t[i];
  cross3(vJ, &pEnd[6], t);  for (int i=0; i<3; i++) pNode[9+i] += t[i];
}

//   K_node = A^T K_end A + K_off,  K_off(rot,rot) = v f^T - (v.f) I  per end
void
Corot3DBattini::endToNodeMatrix(double K[12][12], const double *pEnd) const
{
  if (!hasOffsets) return;
  double A[12][12];
  memset(A, 0, sizeof(A));
  for (int i=0; i<12; i++) A[i][i] = 1.0;
  const double *v[2] = {vI, vJ};
  for (int e=0; e<2; e++) {
    const double *w = v[e];
    const int r = 6*e, c = 6*e + 3;
    // -S(v)
    A[r+0][c+1] =  w[2];  A[r+0][c+2] = -w[1];
    A[r+1][c+0] = -w[2];  A[r+1][c+2] =  w[0];
    A[r+2][c+0] =  w[1];  A[r+2][c+1] = -w[0];
  }
  double KA[12][12];
  for (int i=0; i<12; i++)
    for (int j=0; j<12; j++) {
      double s = 0.0;
      for (int k=0; k<12; k++) s += K[i][k]*A[k][j];
      KA[i][j] = s;
    }
  for (int i=0; i<12; i++)
    for (int j=0; j<12; j++) {
      double s = 0.0;
      for (int k=0; k<12; k++) s += A[k][i]*KA[k][j];
      K[i][j] = s;
    }
  for (int e=0; e<2; e++) {
    const double *w = v[e];
    const double *f = &pEnd[6*e];
    const double vf = dot3(w, f);
    const int r = 6*e + 3;
    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
        K[r+i][r+j] += w[i]*f[j] - (i == j ? vf : 0.0);
  }
}


const Vector &
Corot3DBattini::getGlobalResistingForce(const Vector &pb, const Vector &p0)
{
  // end forces: pEnd = B^T pb  (f_g = B_g^T B_a^T f_l  [B 4.40, 4.74])
  double pEnd[12];
  for (int j=0; j<12; j++) {
    double s = 0.0;
    for (int k=0; k<6; k++)
      s += B[k][j]*pb(k);
    pEnd[j] = s;
  }

  // element loads: end reactions assumed in the initial local system
  // (same treatment as CorotCrdTransf3d); kept for the offset term K_off
  memset(pEndP0, 0, sizeof(pEndP0));
  if (p0 != 0.0) {
    double pl0I[3] = {p0(0), p0(1), p0(3)};
    double pl0J[3] = {0.0,   p0(2), p0(4)};
    matvec3(Ro, pl0I, &pEndP0[0]);
    matvec3(Ro, pl0J, &pEndP0[6]);
    for (int i=0; i<12; i++) pEnd[i] += pEndP0[i];
  }

  // element ends -> nodes (rigid offsets)
  double pNode[12];
  this->endToNodeVector(pEnd, pNode);
  for (int i=0; i<12; i++) Pg(i) = pNode[i];
  return Pg;
}


void
Corot3DBattini::computeGeometricStiffness(const Vector &pb, double K[12][12]) const
{
  const double n   = pb(0);
  const double mI[3] = {-pb(5), pb(3), pb(1)};   // local moments conjugate to thetaBarI
  const double mJ[3] = { pb(5), pb(4), pb(2)};   // local moments conjugate to thetaBarJ

  const double r1[3] = {Rr[0][0], Rr[1][0], Rr[2][0]};
  const double r2[3] = {Rr[0][1], Rr[1][1], Rr[2][1]};
  const double r3[3] = {Rr[0][2], Rr[1][2], Rr[2][2]};

  const double qe1 = dot3(q, r1);
  const double qe2 = dot3(q, r2);
  const double eta = qe1/qe2;

  // global nodal moments mgI = Rr T_s^-T(thetaBarI) mI [B 4.40], and their sum mgSum
  double maI[3], maJ[3], mgI[3], mgJ[3], mgSum[3];
  TsInvTVec(thetaBarI, mI, maI);
  TsInvTVec(thetaBarJ, mJ, maJ);
  matvec3(Rr, maI, mgI);
  matvec3(Rr, maJ, mgJ);
  for (int i=0; i<3; i++) mgSum[i] = mgI[i] + mgJ[i];
  const double s = dot3(mgSum, r1);

  double mgSumXr1[3];
  cross3(mgSum, r1, mgSumXr1);

  // coefficient vectors of GTg^T mgSum that depend on the configuration
  double aI[3], aJ[3], t[3];
  cross3(qI, r3, t);  for (int i=0; i<3; i++) aI[i] = t[i]/(2.0*qe2);
  cross3(qJ, r3, t);  for (int i=0; i<3; i++) aJ[i] = t[i]/(2.0*qe2);
  double ax[3];
  for (int i=0; i<3; i++) ax[i] = -eta*r3[i]/ln;

  for (int j=0; j<12; j++) {
    // unit direction
    double ddx[3] = {0,0,0}, dwI[3] = {0,0,0}, dwJ[3] = {0,0,0};
    if (j < 3)            ddx[j]   = -1.0;
    else if (j < 6)       dwI[j-3] =  1.0;
    else if (j < 9)       ddx[j-6] =  1.0;
    else                  dwJ[j-9] =  1.0;

    const double dLn = dot3(r1, ddx);
    double dr1[3];
    for (int i=0; i<3; i++) dr1[i] = (ddx[i] - r1[i]*dLn)/ln;

    const double deltaEpsR[3] = {GTg[0][j], GTg[1][j], GTg[2][j]};
    double dr2[3], dr3[3];
    cross3(deltaEpsR, r2, dr2);
    cross3(deltaEpsR, r3, dr3);

    double dqI[3], dqJ[3], dq[3];
    cross3(dwI, qI, dqI);
    cross3(dwJ, qJ, dqJ);
    for (int i=0; i<3; i++) dq[i] = 0.5*(dqI[i] + dqJ[i]);

    // increments of the local rotations
    double dwRelI[3], dwRelJ[3], dwBarI[3], dwBarJ[3], dthetaBarI[3], dthetaBarJ[3];
    for (int i=0; i<3; i++) {
      dwRelI[i] = dwI[i] - deltaEpsR[i];
      dwRelJ[i] = dwJ[i] - deltaEpsR[i];
    }
    matTvec3(Rr, dwRelI, dwBarI);
    matTvec3(Rr, dwRelJ, dwBarJ);
    TsInvVec(thetaBarI, dwBarI, dthetaBarI);
    TsInvVec(thetaBarJ, dwBarJ, dthetaBarJ);

    // d(Rr ma) = deltaEpsR x (Rr ma) + Rr d(ma)
    double dmgI[3], dmgJ[3], tmp[3], dmu[3];
    cross3(deltaEpsR, mgI, dmgI);
    dTsInvTVec(thetaBarI, mI, dthetaBarI, dmu);
    matvec3(Rr, dmu, tmp);
    for (int i=0; i<3; i++) dmgI[i] += tmp[i];
    cross3(deltaEpsR, mgJ, dmgJ);
    dTsInvTVec(thetaBarJ, mJ, dthetaBarJ, dmu);
    matvec3(Rr, dmu, tmp);
    for (int i=0; i<3; i++) dmgJ[i] += tmp[i];
    double dmgSum[3];
    for (int i=0; i<3; i++) dmgSum[i] = dmgI[i] + dmgJ[i];

    double col[12];
    for (int i=0; i<12; i++) col[i] = 0.0;

    // axial: d(n r) = n D dp, D3 = (I - r1 r1^T)/ln  [B 4.77]
    for (int i=0; i<3; i++) {
      col[i]   += -n*dr1[i];
      col[i+6] +=  n*dr1[i];
    }

    // nodal moments
    for (int i=0; i<3; i++) {
      col[i+3] += dmgI[i];
      col[i+9] += dmgJ[i];
    }

    // - GTg^T dmgSum  (with the nodal-moment terms above: the K_h part
    // [B 4.43] and the rotation of Rr, cf. -Rr Q G^T Rr^T [B 4.80-4.83])
    for (int k=0; k<12; k++)
      col[k] -= GTg[0][k]*dmgSum[0] + GTg[1][k]*dmgSum[1] + GTg[2][k]*dmgSum[2];

    // - d(GTg^T) mgSum  (mgSum fixed): variation of G, cf. Rr G a r^T [B 4.84-4.93]
    double mgSumXdr1[3], dc1[3];
    cross3(mgSum, dr1, mgSumXdr1);
    for (int i=0; i<3; i++)
      dc1[i] = (mgSumXdr1[i] - mgSumXr1[i]*dLn/ln)/ln;

    const double ds  = dot3(mgSum, dr1);
    const double dQ1 = dot3(dq, r1) + dot3(q, dr1);
    const double dQ2 = dot3(dq, r2) + dot3(q, dr2);
    const double dn  = (dQ1*qe2 - qe1*dQ2)/(qe2*qe2);

    double daI[3], daJ[3], dax[3], c1[3], c2[3];
    cross3(dqI, r3, c1);
    cross3(qI, dr3, c2);
    for (int i=0; i<3; i++)
      daI[i] = (c1[i] + c2[i])/(2.0*qe2) - aI[i]*dQ2/qe2;
    cross3(dqJ, r3, c1);
    cross3(qJ, dr3, c2);
    for (int i=0; i<3; i++)
      daJ[i] = (c1[i] + c2[i])/(2.0*qe2) - aJ[i]*dQ2/qe2;
    for (int i=0; i<3; i++)
      dax[i] = -dn*r3[i]/ln - eta*dr3[i]/ln + eta*r3[i]*dLn/(ln*ln);

    for (int i=0; i<3; i++) {
      const double tx = ds*ax[i] + s*dax[i];
      col[i]   -= -dc1[i] - tx;
      col[i+6] -=  dc1[i] + tx;
      col[i+3] -=  ds*aI[i] + s*daI[i];
      col[i+9] -=  ds*aJ[i] + s*daJ[i];
    }

    for (int i=0; i<12; i++)
      K[i][j] = col[i];
  }
}


const Matrix &
Corot3DBattini::getGlobalStiffMatrix(const Matrix &kb, const Vector &pb)
{
  // material part: B^T kb B
  double kbB[6][12];
  for (int k=0; k<6; k++)
    for (int j=0; j<12; j++) {
      double s = 0.0;
      for (int l=0; l<6; l++)
        s += kb(k,l)*B[l][j];
      kbB[k][j] = s;
    }

  double Kgeo[12][12];
  this->computeGeometricStiffness(pb, Kgeo);

  double K[12][12];
  for (int i=0; i<12; i++)
    for (int j=0; j<12; j++) {
      double s = 0.0;
      for (int k=0; k<6; k++)
        s += B[k][i]*kbB[k][j];
      K[i][j] = s + Kgeo[i][j];     // K_g = B^T K_l B + K_h + K_m  [B 4.43, 4.94]
    }

  // element ends -> nodes (rigid offsets): A^T K A + K_off [B 4.2.4]
  if (hasOffsets) {
    double pEnd[12];
    for (int j=0; j<12; j++) {
      double s = pEndP0[j];
      for (int k=0; k<6; k++)
        s += B[k][j]*pb(k);
      pEnd[j] = s;
    }
    this->endToNodeMatrix(K, pEnd);
  }

  for (int i=0; i<12; i++)
    for (int j=0; j<12; j++)
      Kg(i,j) = K[i][j];

  // optional: symmetric part of the (nonsymmetric) consistent tangent,
  // cf. [B] remark after Eq. 4.94 and Nour-Omid & Rankin (1991)
  if (symmetricTangent) {
    for (int i=0; i<12; i++)
      for (int j=i+1; j<12; j++) {
        double a = 0.5*(Kg(i,j) + Kg(j,i));
        Kg(i,j) = a;
        Kg(j,i) = a;
      }
  }

  return Kg;
}


const Matrix &
Corot3DBattini::getInitialGlobalStiffMatrix(const Matrix &kb)
{
  double kbB[6][12];
  for (int k=0; k<6; k++)
    for (int j=0; j<12; j++) {
      double s = 0.0;
      for (int l=0; l<6; l++)
        s += kb(k,l)*B0[l][j];
      kbB[k][j] = s;
    }
  double K[12][12];
  for (int i=0; i<12; i++)
    for (int j=0; j<12; j++) {
      double s = 0.0;
      for (int k=0; k<6; k++)
        s += B0[k][i]*kbB[k][j];
      K[i][j] = s;
    }

  // rigid offsets in the initial configuration (v = vo, no initial forces)
  if (hasOffsets) {
    double vSaveI[3], vSaveJ[3], zero[12];
    for (int i=0; i<3; i++) {
      vSaveI[i] = vI[i];  vSaveJ[i] = vJ[i];
    }
    for (int i=0; i<12; i++) zero[i] = 0.0;
    for (int i=0; i<3; i++) {vI[i] = voI[i];  vJ[i] = voJ[i];}
    this->endToNodeMatrix(K, zero);
    for (int i=0; i<3; i++) {vI[i] = vSaveI[i];  vJ[i] = vSaveJ[i];}
  }

  for (int i=0; i<12; i++)
    for (int j=0; j<12; j++)
      Kg(i,j) = K[i][j];
  return Kg;
}


CrdTransf *
Corot3DBattini::getCopy3d(void)
{
  static Vector vz(3), zero(3);
  for (int i=0; i<3; i++) vz(i) = vecxz[i];
  zero.Zero();

  Corot3DBattini *theCopy = new Corot3DBattini(this->getTag(), vz, zero, zero,
                                               symmetricTangent);

  theCopy->nodeIPtr = nodeIPtr;
  theCopy->nodeJPtr = nodeJPtr;
  memcpy(theCopy->Ro, Ro, sizeof(Ro));
  theCopy->lo  = lo;
  theCopy->ln = ln;
  memcpy(theCopy->quatRtI, quatRtI, sizeof(quatRtI));
  memcpy(theCopy->quatRtJ, quatRtJ, sizeof(quatRtJ));
  memcpy(theCopy->quatRtICommit, quatRtICommit, sizeof(quatRtICommit));
  memcpy(theCopy->quatRtJCommit, quatRtJCommit, sizeof(quatRtJCommit));
  memcpy(theCopy->alphaI, alphaI, sizeof(alphaI));
  memcpy(theCopy->alphaJ, alphaJ, sizeof(alphaJ));
  memcpy(theCopy->RtI, RtI, sizeof(RtI));
  memcpy(theCopy->RtJ, RtJ, sizeof(RtJ));
  memcpy(theCopy->Rr, Rr, sizeof(Rr));
  memcpy(theCopy->qI, qI, sizeof(qI));
  memcpy(theCopy->qJ, qJ, sizeof(qJ));
  memcpy(theCopy->q, q, sizeof(q));
  memcpy(theCopy->thetaBarI, thetaBarI, sizeof(thetaBarI));
  memcpy(theCopy->thetaBarJ, thetaBarJ, sizeof(thetaBarJ));
  memcpy(theCopy->GTg, GTg, sizeof(GTg));
  memcpy(theCopy->B, B, sizeof(B));
  memcpy(theCopy->B0, B0, sizeof(B0));
  theCopy->ub       = ub;
  theCopy->ubCommit = ubCommit;
  theCopy->ubPrev   = ubPrev;
  if (nodeIInitialDisp != 0) {
    theCopy->nodeIInitialDisp = new double [6];
    for (int i=0; i<6; i++) theCopy->nodeIInitialDisp[i] = nodeIInitialDisp[i];
  }
  if (nodeJInitialDisp != 0) {
    theCopy->nodeJInitialDisp = new double [6];
    for (int i=0; i<6; i++) theCopy->nodeJInitialDisp[i] = nodeJInitialDisp[i];
  }
  theCopy->initialDispChecked = initialDispChecked;
  theCopy->stateRestored = stateRestored;
  for (int i=0; i<3; i++) {
    theCopy->voI[i] = voI[i];  theCopy->voJ[i] = voJ[i];
    theCopy->vI[i]  = vI[i];   theCopy->vJ[i]  = vJ[i];
  }
  theCopy->hasOffsets = hasOffsets;
  memcpy(theCopy->pEndP0, pEndP0, sizeof(pEndP0));

  return theCopy;
}


//  data layout:
//   0- 5 ubCommit       6- 9 quatRtICommit   10-13 quatRtJCommit
//  14-16 alphaI        17-19 alphaJ        20-22 vecxz
//  23-28 node I initial disp               29-34 node J initial disp
//  35    has I init    36    has J init    37    lo      38 tag
//  39    symmetric tangent flag             40-45 offsets voI, voJ
int
Corot3DBattini::sendSelf(int cTag, Channel &theChannel)
{
  static Vector data(46);
  data.Zero();
  for (int i=0; i<6; i++) data(i) = ubCommit(i);
  for (int i=0; i<4; i++) {
    data(6+i)  = quatRtICommit[i];
    data(10+i) = quatRtJCommit[i];
  }
  for (int i=0; i<3; i++) {
    data(14+i) = alphaI[i];
    data(17+i) = alphaJ[i];
    data(20+i) = vecxz[i];
  }
  if (nodeIInitialDisp != 0) {
    data(35) = 1.0;
    for (int i=0; i<6; i++) data(23+i) = nodeIInitialDisp[i];
  }
  if (nodeJInitialDisp != 0) {
    data(36) = 1.0;
    for (int i=0; i<6; i++) data(29+i) = nodeJInitialDisp[i];
  }
  data(37) = lo;
  data(38) = this->getTag();
  data(39) = symmetricTangent ? 1.0 : 0.0;
  for (int i=0; i<3; i++) {
    data(40+i) = voI[i];
    data(43+i) = voJ[i];
  }

  if (theChannel.sendVector(this->getDbTag(), cTag, data) < 0) {
    opserr << "Corot3DBattini::sendSelf() - data could not be sent\n";
    return -1;
  }
  return 0;
}


int
Corot3DBattini::recvSelf(int cTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  static Vector data(46);
  if (theChannel.recvVector(this->getDbTag(), cTag, data) < 0) {
    opserr << "Corot3DBattini::recvSelf() - data could not be received\n";
    return -1;
  }

  for (int i=0; i<6; i++) ubCommit(i) = data(i);
  for (int i=0; i<4; i++) {
    quatRtICommit[i] = data(6+i);
    quatRtJCommit[i] = data(10+i);
    quatRtI[i] = quatRtICommit[i];
    quatRtJ[i] = quatRtJCommit[i];
  }
  for (int i=0; i<3; i++) {
    alphaI[i] = data(14+i);
    alphaJ[i] = data(17+i);
    vecxz[i]  = data(20+i);
  }
  if (data(35) != 0.0) {
    if (nodeIInitialDisp == 0) nodeIInitialDisp = new double [6];
    for (int i=0; i<6; i++) nodeIInitialDisp[i] = data(23+i);
  }
  if (data(36) != 0.0) {
    if (nodeJInitialDisp == 0) nodeJInitialDisp = new double [6];
    for (int i=0; i<6; i++) nodeJInitialDisp[i] = data(29+i);
  }
  lo = data(37);
  this->setTag((int)data(38));
  symmetricTangent = (data(39) != 0.0);
  for (int i=0; i<3; i++) {
    voI[i] = data(40+i);
    voJ[i] = data(43+i);
    vI[i] = voI[i];
    vJ[i] = voJ[i];
  }
  hasOffsets = (norm3(voI) != 0.0 || norm3(voJ) != 0.0);

  ub = ubCommit;
  ubPrev = ubCommit;
  initialDispChecked = true;
  stateRestored = true;
  return 0;
}


const Matrix &
Corot3DBattini::getGlobalMatrixFromLocal(const Matrix &ml)
{
  Tlg.Zero();
  for (int b=0; b<4; b++)
    for (int i=0; i<3; i++)
      for (int j=0; j<3; j++)
        Tlg(3*b+i, 3*b+j) = Ro[j][i];     // local = Ro^T global
  Kg.addMatrixTripleProduct(0.0, Tlg, ml, 1.0);
  return Kg;
}


const Vector &
Corot3DBattini::getPointGlobalCoordFromLocal(const Vector &xl)
{
  // undeformed configuration
  static Vector xg(3);
  const Vector &crdI = nodeIPtr->getCrds();
  for (int i=0; i<3; i++)
    xg(i) = crdI(i) + voI[i] + Ro[i][0]*xl(0) + Ro[i][1]*xl(1) + Ro[i][2]*xl(2);
  return xg;
}


const Vector &
Corot3DBattini::getPointGlobalDisplFromBasic(double xi, const Vector &uxb)
{
  static Vector uxg(3);
  opserr << "Corot3DBattini::getPointGlobalDisplFromBasic - not implemented\n";
  uxg.Zero();
  return uxg;
}


const Vector &
Corot3DBattini::getPointLocalDisplFromBasic(double xi, const Vector &uxb)
{
  static Vector uxl(3);
  opserr << "Corot3DBattini::getPointLocalDisplFromBasic - not implemented\n";
  uxl.Zero();
  return uxl;
}


int
Corot3DBattini::getLocalAxes(Vector &XAxis, Vector &YAxis, Vector &ZAxis)
{
  for (int i=0; i<3; i++) {
    XAxis(i) = Ro[i][0];
    YAxis(i) = Ro[i][1];
    ZAxis(i) = Ro[i][2];
  }
  return 0;
}


int
Corot3DBattini::getRigidOffsets(Vector &offsets)
{
  for (int i=0; i<3; i++) {
    offsets(i)   = voI[i];
    offsets(3+i) = voJ[i];
  }
  return 0;
}


void
Corot3DBattini::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_CURRENTSTATE) {
    s << "\nCrdTransf: " << this->getTag() << " Type: Corot3DBattini";
    s << "\tvAxis: " << vecxz[0] << " " << vecxz[1] << " " << vecxz[2];
    if (symmetricTangent) s << "\t-symmetric";
    if (hasOffsets)
      s << "\tjntOffset: " << voI[0] << " " << voI[1] << " " << voI[2] << "  "
        << voJ[0] << " " << voJ[1] << " " << voJ[2];
    s << "\n";
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{\"name\": \"" << this->getTag() << "\", \"type\": \"Corot3DBattini\"";
    s << ", \"vecInLocXZPlane\": [" << vecxz[0] << ", " << vecxz[1] << ", " << vecxz[2] << "]";
    if (symmetricTangent) s << ", \"symmetric\": true";
    if (hasOffsets) {
      s << ", \"iOffset\": [" << voI[0] << ", " << voI[1] << ", " << voI[2] << "]";
      s << ", \"jOffset\": [" << voJ[0] << ", " << voJ[1] << ", " << voJ[2] << "]";
    }
    s << "}";
  }
}
