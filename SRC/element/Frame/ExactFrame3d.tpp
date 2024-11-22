
//
//===----------------------------------------------------------------------===//
//
//         Please cite the following resource in any derivative works:
//                 https://doi.org/10.5281/zenodo.10456866
//
//===----------------------------------------------------------------------===//
//
// [1] Simo J.C. (1985): A finite strain beam formulation. The three-dimensional
//     dynamic problem. Part I.
//     Computer Methods in Applied Mechanics and Engineering, 49(1):55–70.
//     https://doi.org/10.1016/0045-7825(85)90050-7
//
// [2] Simo J.C., Vu-Quoc L. (1986): A three-dimensional finite-strain rod model
//     Part II: Computational aspects.
//     Computer Methods in Applied Mechanics and Engineering, 58(1):79–116.
//     https://doi.org/10/b8wd4z
//
// [3] Perez C.M., and Filippou F.C. (2024):
//     "On Nonlinear Geometric Transformations of Finite Elements" 
//     Int. J. Numer. Meth. Engrg.
//
//===----------------------------------------------------------------------===//
//
// Claudio M. Perez
//
// Implementation
#include <Flag.h>
#include <Node.h>
#include <Logging.h>
#include <Lagrange1D.cpp>
#include <quadrature/GaussLegendre1D.hpp>

//

namespace OpenSees {

template<int nen> static void
G_matrix(MatrixND<6,6> &G, 
         const VectorND<6>& s, const Vector3D& dx, 
         double shape[2][nen], 
         int i, int j)
{
  auto sn = Hat(&s[0]);
  auto sm = Hat(&s[3]);
  G.assemble(         sn, 0, 3, -shape[1][i]*shape[0][j]);
  G.assemble(         sn, 3, 0,  shape[1][j]*shape[0][i]);

  G.assemble(         sm, 3, 3, -shape[1][i]*shape[0][j]);
  G.assemble( Hat(dx)*sn, 3, 3,  shape[0][i]*shape[0][j]);
}

template<int nen> static void
B_matrix(MatrixND<6,6> &B, double shape[2][nen], const Vector3D& dx, int n)
{
  //
  // NOTE This is the transpose of B in the paper by Perez and Filippou (2024)
  //
  for (int i=0; i<6; i++)
    B(i,i) = shape[1][n];
  
  //
  // B(1:3, 4:end) = shape*Hat(dx);
  //
  B(0,3) =  0;
  B(0,4) = -shape[0][n]*dx[2];
  B(0,5) =  shape[0][n]*dx[1];

  B(1,3) =  shape[0][n]*dx[2];
  B(1,4) =  0;
  B(1,5) = -shape[0][n]*dx[0];

  B(2,3) = -shape[0][n]*dx[1];
  B(2,4) =  shape[0][n]*dx[0];
  B(2,5) =  0;
}


template<int nen, int nip>
ExactFrame3d<nen,nip>::ExactFrame3d(int tag, 
                                    std::array<int,nen>& nodes, 
                                    FrameSection *section[nip],
                                    FrameTransform3d& transf)
    : FiniteElement<nen, 3, 6>(tag, 0, nodes),
      transform(&transf),
      logarithm(Logarithm::None),
      stencil(nullptr)
{
//  double wt[nip];
//  double xi[nip];
//  beamIntegr->getSectionLocations(numSections, L, xi);
//  beamIntegr->getSectionWeights(numSections, L, wt);

    p.zero();
    K.zero();

    for (int i=0; i<nip; i++) {
      points[i].point  = 0.0;
      points[i].weight = 0.0;
      points[i].material=section[i]->getFrameCopy(scheme);
    }

}


template<int nen, int nip>
ExactFrame3d<nen,nip>::~ExactFrame3d()
{
  for (GaussPoint& point : points)
    if (point.material != nullptr)
      delete point.material;

  if (stencil != nullptr)
    delete stencil;
}

template<int nen, int nip>
int
ExactFrame3d<nen,nip>::setNodes()
{
  auto& theNodes = this->FiniteElement<nen,3,6>::theNodes;

  if (transform->initialize(theNodes[0], theNodes[nen-1]) != 0) {
      opserr << " -- Error initializing coordinate transformation\n";
      return -1;
  }
  const Vector& xi = theNodes[    0]->getCrds();
  const Vector& xj = theNodes[nen-1]->getCrds();
  double L = (xi-xj).Norm();


  // Node locations in local (scalar) coordinate
  double x[nen];
  for (int i=0; i < nen; i++)
    x[i] = i*L/(nen-1);

  GaussLegendre<1, nip>    gauss;
  for (int i=0; i < nip; i++) {
    points[i].point  = (gauss.pts[i] + 1.0)*L/2.0;
    points[i].weight =  gauss.wts[i]*L/2.0;
    lagrange<nen>(points[i].point, x, points[i].shape);
  }

  // Zero out the state of the Gauss points
  this->revertToStart();

  return 0;
}

template<int nen, int nip>
int
ExactFrame3d<nen,nip>::revertToStart()
{
  // Revert the transformation to start
  if (transform->revertToStart() != 0)
    return -2;

  Vector E1(3), E2(3), E3(3);
  transform->getLocalAxes(E1, E2, E3);

  Matrix3D R0;
  for (int i=0; i<ndm; i++) {
    R0(i, 0) =  E1[i];
    R0(i, 1) =  E2[i];
    R0(i, 2) =  E3[i];
  }

  // Revert the of the Gauss points to start
  for (GaussPoint& point : points) {
    point.curvature.zero();
    point.rotation = R0;
    if (point.material->revertToStart() != 0)
      return -1;
  }
  past = points;

  // Revert the element state to start
  // NOTE: This assumes that there are zero initial stresses?
  p.zero();
  K.zero();

  return 0;
}


template<int nen, int nip>
int
ExactFrame3d<nen,nip>::revertToLastCommit()
{
  points = past;

  for (GaussPoint& point : points) {
    FrameSection& section = *point.material;

    if (section.revertToLastCommit() != 0)
      return -1;
  }

  return 0;
}

template<int nen, int nip>
int
ExactFrame3d<nen,nip>::update()
{

  const Vector3D D {1, 0, 0};
  auto& theNodes = this->FiniteElement<nen,3,6>::theNodes;

  //
  // Collect nodal parameters
  //
  VectorND<ndf> ddu[nen];
  for (int i=0; i < nen; i++) {
    const Vector& ddui = theNodes[i]->getIncrDeltaDisp();
    for (int j=0; j<ndf; j++)
      ddu[i][j] = ddui[j];
  }

  // Form displaced node locations xyz
  VectorND<ndm> xyz[nen];
  for (int i=0; i < nen; i++) {
    const Vector& xi = theNodes[i]->getCrds();
    const Vector& ui = theNodes[i]->getTrialDisp();
    for (int j=0; j<ndm; j++)
      xyz[i][j] = xi[j] + ui[j];

  }

  //
  // Gauss loop
  //
  p.zero();
  K.zero();
  for (int i=0; i<nip; i++) {
    //
    // Interpolate
    //
    Vector3D dx     {0.0};
    Vector3D theta  {0.0};
    Vector3D dtheta {0.0};


    for (int j=0; j < nen; j++) {
      for (int l=0; l<3; l++)
        dx[l]     += points[i].shape[1][j]*xyz[j][l];
      for (int l=0; l<3; l++)
        theta[l]  += points[i].shape[0][j]*ddu[j][l+3];
      for (int l=0; l<3; l++)
        dtheta[l] += points[i].shape[1][j]*ddu[j][l+3];
    }

    //
    //
    MatrixND<3,3> dR = ExpSO3(theta);
    Matrix3D R = dR*points[i].rotation;

    points[i].rotation = R;

    Vector3D omega = dR*points[i].curvature;
    // TODO: choose 'R/L'
//  points[i].curvature = omega + TanSO3(theta, 'R')*dtheta;
    points[i].curvature = omega + dExpSO3(theta)*dtheta;

    Vector3D gamma = (R^dx) - D;
    Vector3D kappa = R^points[i].curvature;

    VectorND<6> e {
      gamma[0], gamma[1], gamma[2],
      kappa[0], kappa[1], kappa[2],
    };

    FrameSection& section = *points[i].material;
    section.setTrialState<nsr,scheme>(e);
    VectorND<nsr> s = section.getResultant<nsr,scheme>();
    MatrixND<nsr,nsr> Ks = section.getTangent<nsr,scheme>(State::Pres);

    //
    //
    //
    // A = diag(R, R);
    // Note that this is transposed
    MatrixND<6,6> A {{
      {R(0,0), R(1,0), R(2,0), 0, 0, 0},
      {R(0,1), R(1,1), R(2,1), 0, 0, 0},
      {R(0,2), R(1,2), R(2,2), 0, 0, 0},
      {0, 0, 0, R(0,0), R(1,0), R(2,0)},
      {0, 0, 0, R(0,1), R(1,1), R(2,1)},
      {0, 0, 0, R(0,2), R(1,2), R(2,2)},
    }};

    MatrixND<6,6> B[nen], Bj;
    Bj.zero();
    for (int j=0; j<nen; j++) {
      B_matrix(Bj,  points[i].shape, dx, j);
      B[j] = A^Bj;

      // p += B s w
      VectorND<ndf> pj = B[j]^s;
      for (int l=0; l<ndf; l++)
        p[j*ndf+l] += points[i].weight * pj[l];
    }

    // Material Tangent
    MatrixND<ndf,ndf> Kjk;
    for (int j=0; j<nen; j++) {
      for (int k=0; k<nen; k++) {
        Kjk.addMatrixTripleProduct(0.0, B[j], Ks, B[k], points[i].weight);

        for (int ii=0; ii<ndf; ii++) {
          for (int jj=0; jj<ndf; jj++) {
            K(j*ndf+ii, k*ndf+jj) += Kjk(ii,jj);
          }
        }

      }
    }

    // Geometric Tangent
    MatrixND<ndf,ndf> G;
    for (int j=0; j<nen; j++) {
      for (int k=0; k<nen; k++) {
        G.zero();
        G_matrix(G, s, dx, points[i].shape, j, k);
        K.assemble(G, ndf*j, ndf*k, points[i].weight);
      }
    }
  }

//commitState();

  return OpenSees::Flag::Success;
}

template<int nen, int nip>
const Vector &
ExactFrame3d<nen,nip>::getResistingForce()
{
  thread_local Vector wrapper;
  wrapper.setData(p);
  return wrapper;
}

template<int nen, int nip>
const Matrix &
ExactFrame3d<nen,nip>::getTangentStiff()
{
  thread_local Matrix wrapper;
  wrapper.setData(K);
  return wrapper;
}

template<int nen, int nip>
const Matrix &
ExactFrame3d<nen,nip>::getInitialStiff()
{
  static MatrixND<ndf*nen,ndf*nen> Ki{};
  static Matrix wrapper(Ki);
  return wrapper;
}

template<int nen, int nip>
const Matrix &
ExactFrame3d<nen,nip>::getMass()
{
  // TODO
  static MatrixND<ndf*nen,ndf*nen> M{0};
  static Matrix wrapper(M);
  return wrapper;
}


template<int nen, int nip>
int
ExactFrame3d<nen,nip>::sendSelf(int commitTag, Channel& theChannel)
{
  // TODO
  return -1;
}

template<int nen, int nip>
int
ExactFrame3d<nen,nip>::recvSelf(int commitTag, Channel& theChannel, FEM_ObjectBroker& theBroker)
{
  // TODO
  return -1;
}

template<int nen, int nip>
void
ExactFrame3d<nen,nip>::Print(OPS_Stream& stream, int flag)
{
  const ID& node_tags = this->getExternalNodes();

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {

    stream << OPS_PRINT_JSON_ELEM_INDENT << "{";
    stream << "\"name\": " << this->getTag() << ", ";
    stream << "\"type\": \"" << this->getClassType() << "\", ";
    stream << "\"nodes\": [" << node_tags(0) << ", " 
                             << node_tags(1) << "]";
    stream << ", ";


    stream << "\"sections\": [";
    for (decltype(points.size()) i = 0; i < points.size() - 1; i++)
      stream << points[i].material->getTag() << ", ";
    stream << points[points.size() - 1].material->getTag() << "]";
    stream << ", ";

    stream << "\"crdTransformation\": " << transform->getTag()  ;
    stream << "}";
  }
}

} // namespace OpenSees
