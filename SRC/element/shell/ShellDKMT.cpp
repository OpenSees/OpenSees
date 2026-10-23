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
**   Gregory L. Fenves (fenves@ce.berkeley.edu)                       **
**   Filip C. Filippou (filippou@ce.berkeley.edu)                     **
**                                                                    **
** ****************************************************************** */

// Written: Nobuyuki Takimoto
//
// Three node flat shell element: CST membrane + DKMT plate bending.
//
// Plate bending (Katili 1993): the rotations are
//   beta = sum_i N_i beta_i + sum_k P_k t_k dbeta_k
// with the quadratic edge bubbles P_k and the tangential corrections
// dbeta_k. The transverse shear strain gamma = Ds^-1 div(Db kappa) is
// constant in the element, and the edge constraints
//   L_k t_k.gamma - 2/3 L_k dbeta_k = (w_j - w_i) + L_k/2 t_k.(beta_i + beta_j)
// give dbeta in terms of the nodal dofs. When the plate is thin the shear
// strains vanish and the element is the DKT.
//
// The relation between dbeta and gamma uses the initial tangent of the
// section, so with a nonlinear section the kinematics stay those of the
// initial (elastic) state while the section forces and tangent follow the
// section. Membrane, bending and shear generalized strains are evaluated at
// the three mid-side points (weights A/3), which integrate the bending
// energy of the element exactly.
//
// The drilling dof is tied to the in-plane rigid rotation of the CST,
// omega = (v,x - u,y)/2, by the penalty drillRatio * D(2,2) * A/3 at each
// node (D(2,2) is the in-plane shear stiffness of the section).

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include <ID.h>
#include <Vector.h>
#include <Matrix.h>
#include <Element.h>
#include <Node.h>
#include <SectionForceDeformation.h>
#include <Domain.h>
#include <ErrorHandler.h>
#include <ShellDKMT.h>
#include <Renderer.h>
#include <ElementResponse.h>
#include <Information.h>
#include <Parameter.h>

#include <Channel.h>
#include <FEM_ObjectBroker.h>
#include <elementAPI.h>

#include <ElementalLoad.h>

void *
OPS_ShellDKMT(void)
{
  if (OPS_GetNumRemainingInputArgs() < 5) {
    opserr << "Want: element ShellDKMT $tag $iNode $jNode $kNode $secTag <-drill $alpha>\n";
    return 0;
  }

  int iData[5];
  int numData = 5;
  if (OPS_GetInt(&numData, iData) != 0) {
    opserr << "WARNING invalid integer tag: element ShellDKMT\n";
    return 0;
  }

  double drillRatio = 1.0e-4;
  while (OPS_GetNumRemainingInputArgs() > 0) {
    const char *type = OPS_GetString();
    if (strcmp(type, "-drill") == 0) {
      numData = 1;
      if (OPS_GetNumRemainingInputArgs() < 1 || OPS_GetDoubleInput(&numData, &drillRatio) < 0) {
        opserr << "WARNING invalid -drill value: element ShellDKMT " << iData[0] << "\n";
        return 0;
      }
    } else {
      opserr << "WARNING unknown option " << type << ": element ShellDKMT " << iData[0] << "\n";
      return 0;
    }
  }

  SectionForceDeformation *theSection = OPS_getSectionForceDeformation(iData[4]);
  if (theSection == 0) {
    opserr << "ERROR: element ShellDKMT " << iData[0] << " section " << iData[4] << " not found\n";
    return 0;
  }

  return new ShellDKMT(iData[0], iData[1], iData[2], iData[3], *theSection, drillRatio);
}

//static data
Matrix ShellDKMT::stiff(18, 18);
Vector ShellDKMT::resid(18);
Matrix ShellDKMT::mass(18, 18);

// mid-side points in area coordinates (xi, eta), weight A/3 each
static const double gxi[3]  = {0.5, 0.5, 0.0};
static const double geta[3] = {0.0, 0.5, 0.5};

// edges k = (i, j)
static const int edgeI[3] = {0, 1, 2};
static const int edgeJ[3] = {1, 2, 0};

//null constructor
ShellDKMT::ShellDKMT() :
  Element(0, ELE_TAG_ShellDKMT),
  connectedExternalNodes(3), drillRatio(1.0e-4), area(0.0), kdrill(0.0),
  load(0), applyLoad(0)
{
  for (int i = 0; i < 3; i++) {
    nodePointers[i] = 0;
    appliedB[i] = 0.0;
  }
  for (int i = 0; i < numGauss; i++) {
    sections[i] = 0;
    B[i].resize(8, 18);
    B1[i].resize(8, 18);
  }
}

//full constructor
ShellDKMT::ShellDKMT(int tag, int node1, int node2, int node3,
                     SectionForceDeformation &theSection, double ratio) :
  Element(tag, ELE_TAG_ShellDKMT),
  connectedExternalNodes(3), drillRatio(ratio), area(0.0), kdrill(0.0),
  load(0), applyLoad(0)
{
  connectedExternalNodes(0) = node1;
  connectedExternalNodes(1) = node2;
  connectedExternalNodes(2) = node3;

  for (int i = 0; i < 3; i++) {
    nodePointers[i] = 0;
    appliedB[i] = 0.0;
  }

  for (int i = 0; i < numGauss; i++) {
    sections[i] = theSection.getCopy();
    if (sections[i] == 0) {
      opserr << "ShellDKMT::constructor - failed to get a copy of the section\n";
      exit(-1);
    }
    B[i].resize(8, 18);
    B1[i].resize(8, 18);
  }
}

//destructor
ShellDKMT::~ShellDKMT()
{
  for (int i = 0; i < numGauss; i++) {
    if (sections[i] != 0)
      delete sections[i];
    sections[i] = 0;
  }
  if (load != 0)
    delete load;
}

void ShellDKMT::setDomain(Domain *theDomain)
{
  for (int i = 0; i < 3; i++) {
    nodePointers[i] = theDomain->getNode(connectedExternalNodes(i));
    if (nodePointers[i] == 0) {
      opserr << "ShellDKMT::setDomain - no node " << connectedExternalNodes(i)
             << " exists in the model\n";
      return;
    }
    if (nodePointers[i]->getNumberDOF() != 6) {
      opserr << "ShellDKMT::setDomain - node " << connectedExternalNodes(i)
             << " NEEDS 6 dof\n";
      return;
    }
  }

  for (int i = 0; i < numGauss; i++) {
    if (sections[i]->getOrder() != 8) {
      opserr << "ShellDKMT::setDomain - element " << this->getTag()
             << " needs a shell section with 8 generalized strains\n";
      return;
    }
  }

  computeBasis();
  formB();

  this->DomainComponent::setDomain(theDomain);
}

int ShellDKMT::getNumExternalNodes() const
{
  return 3;
}

const ID &ShellDKMT::getExternalNodes()
{
  return connectedExternalNodes;
}

Node **ShellDKMT::getNodePtrs()
{
  return nodePointers;
}

int ShellDKMT::getNumDOF()
{
  return 18;
}

int ShellDKMT::commitState()
{
  int success = 0;
  if ((success = this->Element::commitState()) != 0)
    opserr << "ShellDKMT::commitState () - failed in base class";
  for (int i = 0; i < numGauss; i++)
    success += sections[i]->commitState();
  return success;
}

int ShellDKMT::revertToLastCommit()
{
  int success = 0;
  for (int i = 0; i < numGauss; i++)
    success += sections[i]->revertToLastCommit();
  return success;
}

int ShellDKMT::revertToStart()
{
  int success = 0;
  for (int i = 0; i < numGauss; i++)
    success += sections[i]->revertToStart();
  return success;
}

void ShellDKMT::Print(OPS_Stream &s, int flag)
{
  if (flag == OPS_PRINT_CURRENTSTATE) {
    s << endln;
    s << "DKMT Three Node Shell (CST membrane)\n";
    s << "Element Number: " << this->getTag() << endln;
    s << "Node 1 : " << connectedExternalNodes(0) << endln;
    s << "Node 2 : " << connectedExternalNodes(1) << endln;
    s << "Node 3 : " << connectedExternalNodes(2) << endln;
    s << "Drilling ratio : " << drillRatio << endln;
    s << "Section Information : \n ";
    sections[0]->Print(s, flag);
    s << endln;
  }

  if (flag == OPS_PRINT_PRINTMODEL_JSON) {
    s << "\t\t\t{";
    s << "\"name\": " << this->getTag() << ", ";
    s << "\"type\": \"ShellDKMT\", ";
    s << "\"nodes\": [" << connectedExternalNodes(0) << ", ";
    s << connectedExternalNodes(1) << ", ";
    s << connectedExternalNodes(2) << "], ";
    s << "\"section\": \"" << sections[0]->getTag() << "\", ";
    s << "\"drill\": " << drillRatio << "}";
  }
}

// local frame: e1 along edge 1-2, e2 in the plane, e3 = e1 x e2
void ShellDKMT::computeBasis()
{
  const Vector &c0 = nodePointers[0]->getCrds();
  const Vector &c1 = nodePointers[1]->getCrds();
  const Vector &c2 = nodePointers[2]->getCrds();

  double v1[3], v2[3], v3[3];
  for (int a = 0; a < 3; a++) {
    v1[a] = c1(a) - c0(a);
    v2[a] = c2(a) - c0(a);
  }
  double len = sqrt(v1[0]*v1[0] + v1[1]*v1[1] + v1[2]*v1[2]);
  for (int a = 0; a < 3; a++)
    v1[a] /= len;
  double d = v2[0]*v1[0] + v2[1]*v1[1] + v2[2]*v1[2];
  for (int a = 0; a < 3; a++)
    v2[a] -= d * v1[a];
  len = sqrt(v2[0]*v2[0] + v2[1]*v2[1] + v2[2]*v2[2]);
  for (int a = 0; a < 3; a++)
    v2[a] /= len;
  v3[0] = v1[1]*v2[2] - v1[2]*v2[1];
  v3[1] = v1[2]*v2[0] - v1[0]*v2[2];
  v3[2] = v1[0]*v2[1] - v1[1]*v2[0];

  for (int a = 0; a < 3; a++) {
    g[0][a] = v1[a];
    g[1][a] = v2[a];
    g[2][a] = v3[a];
  }

  for (int i = 0; i < 3; i++) {
    const Vector &ci = nodePointers[i]->getCrds();
    double dx[3];
    for (int a = 0; a < 3; a++)
      dx[a] = ci(a) - c0(a);
    xl[i][0] = dx[0]*v1[0] + dx[1]*v1[1] + dx[2]*v1[2];
    xl[i][1] = dx[0]*v2[0] + dx[1]*v2[1] + dx[2]*v2[2];
  }
}

// strain-displacement matrices, local dofs per node [u v w thx thy thz]
// rotations of the plate theory: beta_x = thy, beta_y = -thx
void ShellDKMT::formB()
{
  // jacobian of the map (xi, eta) -> (x, y), rows d/dxi, d/deta
  double J[2][2];
  J[0][0] = xl[1][0] - xl[0][0];  J[0][1] = xl[1][1] - xl[0][1];
  J[1][0] = xl[2][0] - xl[0][0];  J[1][1] = xl[2][1] - xl[0][1];
  double detJ = J[0][0]*J[1][1] - J[0][1]*J[1][0];
  area = 0.5 * detJ;
  double Ji[2][2];
  Ji[0][0] =  J[1][1] / detJ;  Ji[0][1] = -J[0][1] / detJ;
  Ji[1][0] = -J[1][0] / detJ;  Ji[1][1] =  J[0][0] / detJ;

  // derivatives of the linear shape functions
  static const double dNxi[3][2] = {{-1.0, -1.0}, {1.0, 0.0}, {0.0, 1.0}};
  double Nx[3], Ny[3];
  for (int i = 0; i < 3; i++) {
    Nx[i] = Ji[0][0]*dNxi[i][0] + Ji[0][1]*dNxi[i][1];
    Ny[i] = Ji[1][0]*dNxi[i][0] + Ji[1][1]*dNxi[i][1];
  }

  // edge length and unit tangent
  double L[3], c[3], s[3];
  for (int k = 0; k < 3; k++) {
    double dx = xl[edgeJ[k]][0] - xl[edgeI[k]][0];
    double dy = xl[edgeJ[k]][1] - xl[edgeI[k]][1];
    L[k] = sqrt(dx*dx + dy*dy);
    c[k] = dx / L[k];
    s[k] = dy / L[k];
  }

  // bending and shear rigidities of the initial state (standard signs)
  const Matrix &D0 = sections[0]->getInitialTangent();
  double Db[3][3];
  for (int a = 0; a < 3; a++)
    for (int b = 0; b < 3; b++)
      Db[a][b] = -D0(3 + a, 3 + b);
  double Ds[2][2] = {{D0(6, 6), D0(6, 7)}, {D0(7, 6), D0(7, 7)}};
  double detDs = Ds[0][0]*Ds[1][1] - Ds[0][1]*Ds[1][0];
  double Dsi[2][2] = {{0.0, 0.0}, {0.0, 0.0}};
  if (detDs > 0.0) {
    Dsi[0][0] =  Ds[1][1] / detDs;  Dsi[0][1] = -Ds[0][1] / detDs;
    Dsi[1][0] = -Ds[1][0] / detDs;  Dsi[1][1] =  Ds[0][0] / detDs;
  }  // no shear rigidity: Kirchhoff limit (DKT)

  // second derivatives of the edge bubbles P0 = 4 N0 N1, P1 = 4 N1 N2, P2 = 4 N2 N0
  static const double Hxi[3][2][2] = {{{-8.0, -4.0}, {-4.0, 0.0}},
                                      {{ 0.0,  4.0}, { 4.0, 0.0}},
                                      {{ 0.0, -4.0}, {-4.0, -8.0}}};

  // shear forces from the bubbles, Gd = div(Db kappa) per unit dbeta_k
  double Gd[2][3];
  for (int k = 0; k < 3; k++) {
    double H[2][2];
    for (int a = 0; a < 2; a++)
      for (int b = 0; b < 2; b++) {
        H[a][b] = 0.0;
        for (int p = 0; p < 2; p++)
          for (int q = 0; q < 2; q++)
            H[a][b] += Ji[a][p] * Hxi[k][p][q] * Ji[b][q];
      }
    double dkx[3] = {c[k]*H[0][0], s[k]*H[0][1], c[k]*H[0][1] + s[k]*H[0][0]};
    double dky[3] = {c[k]*H[0][1], s[k]*H[1][1], c[k]*H[1][1] + s[k]*H[0][1]};
    double ma[3], mb[3];
    for (int a = 0; a < 3; a++) {
      ma[a] = Db[a][0]*dkx[0] + Db[a][1]*dkx[1] + Db[a][2]*dkx[2];
      mb[a] = Db[a][0]*dky[0] + Db[a][1]*dky[1] + Db[a][2]*dky[2];
    }
    Gd[0][k] = ma[0] + mb[2];
    Gd[1][k] = ma[2] + mb[1];
  }

  // shear strain per unit dbeta_k
  double Gg[2][3];
  for (int a = 0; a < 2; a++)
    for (int k = 0; k < 3; k++)
      Gg[a][k] = Dsi[a][0]*Gd[0][k] + Dsi[a][1]*Gd[1][k];

  // edge constraints  Am dbeta = Aw U
  Matrix Am(3, 3), Aw(3, 18), Q(3, 18);
  Am.Zero();
  Aw.Zero();
  for (int k = 0; k < 3; k++) {
    for (int m = 0; m < 3; m++)
      Am(k, m) = L[k] * (c[k]*Gg[0][m] + s[k]*Gg[1][m]);
    Am(k, k) -= 2.0 / 3.0 * L[k];
    int i = edgeI[k], j = edgeJ[k];
    Aw(k, 6*i + 2) -= 1.0;
    Aw(k, 6*j + 2) += 1.0;
    Aw(k, 6*i + 4) += 0.5 * L[k] * c[k];
    Aw(k, 6*j + 4) += 0.5 * L[k] * c[k];
    Aw(k, 6*i + 3) -= 0.5 * L[k] * s[k];
    Aw(k, 6*j + 3) -= 0.5 * L[k] * s[k];
  }
  if (Am.Solve(Aw, Q) < 0)
    opserr << "ShellDKMT::formB - element " << this->getTag()
           << " singular edge constraint\n";

  // constant shear strains
  double Bs[2][18];
  for (int a = 0; a < 2; a++)
    for (int d = 0; d < 18; d++)
      Bs[a][d] = Gg[a][0]*Q(0, d) + Gg[a][1]*Q(1, d) + Gg[a][2]*Q(2, d);

  for (int gp = 0; gp < numGauss; gp++) {
    double xi = gxi[gp], eta = geta[gp];
    double dPxi[3][2] = {{4.0*(1.0 - 2.0*xi - eta), -4.0*xi},
                         {4.0*eta, 4.0*xi},
                         {-4.0*eta, 4.0*(1.0 - xi - 2.0*eta)}};
    double Px[3], Py[3];
    for (int k = 0; k < 3; k++) {
      Px[k] = Ji[0][0]*dPxi[k][0] + Ji[0][1]*dPxi[k][1];
      Py[k] = Ji[1][0]*dPxi[k][0] + Ji[1][1]*dPxi[k][1];
    }

    // curvatures in the standard signs
    double kap[3][18];
    for (int a = 0; a < 3; a++)
      for (int d = 0; d < 18; d++)
        kap[a][d] = 0.0;
    for (int i = 0; i < 3; i++) {
      kap[0][6*i + 4] += Nx[i];
      kap[1][6*i + 3] -= Ny[i];
      kap[2][6*i + 4] += Ny[i];
      kap[2][6*i + 3] -= Nx[i];
    }
    for (int k = 0; k < 3; k++)
      for (int d = 0; d < 18; d++) {
        kap[0][d] += c[k]*Px[k] * Q(k, d);
        kap[1][d] += s[k]*Py[k] * Q(k, d);
        kap[2][d] += (c[k]*Py[k] + s[k]*Px[k]) * Q(k, d);
      }

    Matrix &Bg = B[gp];
    Bg.Zero();
    for (int i = 0; i < 3; i++) {
      Bg(0, 6*i)     = Nx[i];
      Bg(1, 6*i + 1) = Ny[i];
      Bg(2, 6*i)     = Ny[i];
      Bg(2, 6*i + 1) = Nx[i];
    }
    for (int d = 0; d < 18; d++) {
      for (int a = 0; a < 3; a++)
        Bg(3 + a, d) = -kap[a][d];
      Bg(6, d) = Bs[0][d];
      Bg(7, d) = Bs[1][d];
    }

    B1[gp] = Bg;
    for (int a = 3; a < 6; a++)
      for (int d = 0; d < 18; d++)
        B1[gp](a, d) = -Bg(a, d);
  }

  // drilling penalty on theta_z_i - omega
  for (int i = 0; i < 3; i++) {
    for (int d = 0; d < 18; d++)
      Bdrill[i][d] = 0.0;
    Bdrill[i][6*i + 5] = 1.0;
    for (int j = 0; j < 3; j++) {
      Bdrill[i][6*j + 1] -= 0.5 * Nx[j];
      Bdrill[i][6*j]     += 0.5 * Ny[j];
    }
  }
  kdrill = drillRatio * D0(2, 2) * area / 3.0;
}

// nodal displacements in the local frame
void ShellDKMT::localDisp(Vector &ul)
{
  for (int i = 0; i < 3; i++) {
    const Vector &u = nodePointers[i]->getTrialDisp();
    for (int a = 0; a < 3; a++) {
      double t = 0.0, r = 0.0;
      for (int b = 0; b < 3; b++) {
        t += g[a][b] * u(b);
        r += g[a][b] * u(3 + b);
      }
      ul(6*i + a) = t;
      ul(6*i + 3 + a) = r;
    }
  }
}

// rotate local stiffness and force to the global frame (blocks of 3)
void ShellDKMT::localToGlobal(const Matrix &Kl, const Vector &Rl, Matrix &K, Vector &R, bool doK)
{
  for (int I = 0; I < 6; I++) {
    for (int p = 0; p < 3; p++) {
      double t = 0.0;
      for (int a = 0; a < 3; a++)
        t += g[a][p] * Rl(3*I + a);
      R(3*I + p) = t;
    }
  }
  if (!doK)
    return;
  static double tmp[3][3];
  for (int I = 0; I < 6; I++)
    for (int Jb = 0; Jb < 6; Jb++) {
      // tmp = Kl_IJ g
      for (int a = 0; a < 3; a++)
        for (int q = 0; q < 3; q++) {
          double t = 0.0;
          for (int b = 0; b < 3; b++)
            t += Kl(3*I + a, 3*Jb + b) * g[b][q];
          tmp[a][q] = t;
        }
      for (int p = 0; p < 3; p++)
        for (int q = 0; q < 3; q++) {
          double t = 0.0;
          for (int a = 0; a < 3; a++)
            t += g[a][p] * tmp[a][q];
          K(3*I + p, 3*Jb + q) = t;
        }
    }
}

// tangFlag 0: resisting force, 1: force and tangent, 2: initial tangent
void ShellDKMT::formLocal(Matrix &K, Vector &R, int tangFlag)
{
  static Vector ul(18);
  static Vector strain(8);
  static Matrix DB(8, 18);

  K.Zero();
  R.Zero();
  if (tangFlag != 2)
    localDisp(ul);

  double dA = area / 3.0;
  for (int gp = 0; gp < numGauss; gp++) {
    if (tangFlag == 2) {
      const Matrix &D = sections[gp]->getInitialTangent();
      K.addMatrixTripleProduct(1.0, B1[gp], D, B[gp], dA);
      continue;
    }
    strain.addMatrixVector(0.0, B[gp], ul, 1.0);
    sections[gp]->setTrialSectionDeformation(strain);
    const Vector &S = sections[gp]->getStressResultant();
    R.addMatrixTransposeVector(1.0, B1[gp], S, dA);
    if (tangFlag == 1) {
      const Matrix &D = sections[gp]->getSectionTangent();
      K.addMatrixTripleProduct(1.0, B1[gp], D, B[gp], dA);
    }
  }

  for (int i = 0; i < 3; i++) {
    if (tangFlag != 2) {
      double e = 0.0;
      for (int d = 0; d < 18; d++)
        e += Bdrill[i][d] * ul(d);
      for (int d = 0; d < 18; d++)
        R(d) += kdrill * e * Bdrill[i][d];
    }
    if (tangFlag != 0)
      for (int a = 0; a < 18; a++)
        for (int d = 0; d < 18; d++)
          K(a, d) += kdrill * Bdrill[i][a] * Bdrill[i][d];
  }
}

const Matrix &ShellDKMT::getTangentStiff()
{
  static Matrix Kl(18, 18);
  static Vector Rl(18);
  formLocal(Kl, Rl, 1);
  localToGlobal(Kl, Rl, stiff, resid, true);
  return stiff;
}

const Matrix &ShellDKMT::getInitialStiff()
{
  static Matrix Kl(18, 18);
  static Vector Rl(18);
  formLocal(Kl, Rl, 2);
  localToGlobal(Kl, Rl, stiff, resid, true);
  return stiff;
}

// consistent translational mass, rotational inertia neglected
void ShellDKMT::formInertiaTerms(int tangFlag)
{
  mass.Zero();
  double rhoH = sections[0]->getRho();
  if (rhoH == 0.0)
    return;
  for (int i = 0; i < 3; i++)
    for (int j = 0; j < 3; j++) {
      double m = rhoH * area / 12.0 * (i == j ? 2.0 : 1.0);
      for (int p = 0; p < 3; p++)
        mass(6*i + p, 6*j + p) = m;
    }
}

const Matrix &ShellDKMT::getMass()
{
  formInertiaTerms(1);
  return mass;
}

void ShellDKMT::zeroLoad(void)
{
  if (load != 0)
    load->Zero();
  appliedB[0] = 0.0;
  appliedB[1] = 0.0;
  appliedB[2] = 0.0;
  applyLoad = 0;
}

int ShellDKMT::addLoad(ElementalLoad *theLoad, double loadFactor)
{
  int type;
  const Vector &data = theLoad->getData(type, loadFactor);

  if (type == LOAD_TAG_SelfWeight) {
    applyLoad = 1;
    appliedB[0] += loadFactor * data(0);
    appliedB[1] += loadFactor * data(1);
    appliedB[2] += loadFactor * data(2);
    return 0;
  }

  opserr << "ShellDKMT::addLoad - load type unknown for ele with tag: " << this->getTag() << endln;
  return -1;
}

int ShellDKMT::addInertiaLoadToUnbalance(const Vector &accel)
{
  if (sections[0]->getRho() == 0.0)
    return 0;

  static Vector ra(18);
  for (int i = 0; i < 3; i++) {
    const Vector &Raccel = nodePointers[i]->getRV(accel);
    for (int j = 0; j < 6; j++)
      ra(6*i + j) = Raccel(j);
  }

  formInertiaTerms(1);
  if (load == 0)
    load = new Vector(18);
  load->addMatrixVector(1.0, mass, ra, -1.0);

  return 0;
}

const Vector &ShellDKMT::getResistingForce()
{
  static Matrix Kl(18, 18);
  static Vector Rl(18);
  formLocal(Kl, Rl, 0);
  localToGlobal(Kl, Rl, stiff, resid, false);

  // self weight, lumped to the nodes (exact for a constant body force)
  if (applyLoad == 1) {
    double w = sections[0]->getRho() * area / 3.0;
    for (int i = 0; i < 3; i++)
      for (int p = 0; p < 3; p++)
        resid(6*i + p) -= w * appliedB[p];
  }

  if (load != 0)
    resid -= *load;

  return resid;
}

const Vector &ShellDKMT::getResistingForceIncInertia()
{
  static Vector res(18);
  static Vector acc(18);

  res = this->getResistingForce();

  formInertiaTerms(0);
  for (int i = 0; i < 3; i++) {
    const Vector &a = nodePointers[i]->getTrialAccel();
    for (int j = 0; j < 6; j++)
      acc(6*i + j) = a(j);
  }
  res.addMatrixVector(1.0, mass, acc, 1.0);

  if (alphaM != 0.0 || betaK != 0.0 || betaK0 != 0.0 || betaKc != 0.0)
    res += this->getRayleighDampingForces();

  return res;
}

int ShellDKMT::sendSelf(int commitTag, Channel &theChannel)
{
  int res = 0;
  int dataTag = this->getDbTag();

  static ID idData(11);
  for (int i = 0; i < numGauss; i++) {
    idData(i) = sections[i]->getClassTag();
    int secDbTag = sections[i]->getDbTag();
    if (secDbTag == 0) {
      secDbTag = theChannel.getDbTag();
      if (secDbTag != 0)
        sections[i]->setDbTag(secDbTag);
    }
    idData(i + 3) = secDbTag;
  }
  idData(6) = this->getTag();
  idData(7) = connectedExternalNodes(0);
  idData(8) = connectedExternalNodes(1);
  idData(9) = connectedExternalNodes(2);
  idData(10) = applyLoad;

  res += theChannel.sendID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING ShellDKMT::sendSelf() - " << this->getTag() << " failed to send ID\n";
    return res;
  }

  static Vector vectData(8);
  vectData(0) = alphaM;
  vectData(1) = betaK;
  vectData(2) = betaK0;
  vectData(3) = betaKc;
  vectData(4) = drillRatio;
  vectData(5) = appliedB[0];
  vectData(6) = appliedB[1];
  vectData(7) = appliedB[2];

  res += theChannel.sendVector(dataTag, commitTag, vectData);
  if (res < 0) {
    opserr << "WARNING ShellDKMT::sendSelf() - " << this->getTag() << " failed to send Vector\n";
    return res;
  }

  for (int i = 0; i < numGauss; i++) {
    res += sections[i]->sendSelf(commitTag, theChannel);
    if (res < 0) {
      opserr << "WARNING ShellDKMT::sendSelf() - " << this->getTag() << " failed to send its section\n";
      return res;
    }
  }

  return res;
}

int ShellDKMT::recvSelf(int commitTag, Channel &theChannel, FEM_ObjectBroker &theBroker)
{
  int res = 0;
  int dataTag = this->getDbTag();

  static ID idData(11);
  res += theChannel.recvID(dataTag, commitTag, idData);
  if (res < 0) {
    opserr << "WARNING ShellDKMT::recvSelf() - " << this->getTag() << " failed to receive ID\n";
    return res;
  }

  this->setTag(idData(6));
  connectedExternalNodes(0) = idData(7);
  connectedExternalNodes(1) = idData(8);
  connectedExternalNodes(2) = idData(9);
  applyLoad = idData(10);

  static Vector vectData(8);
  res += theChannel.recvVector(dataTag, commitTag, vectData);
  if (res < 0) {
    opserr << "WARNING ShellDKMT::recvSelf() - " << this->getTag() << " failed to receive Vector\n";
    return res;
  }
  alphaM = vectData(0);
  betaK = vectData(1);
  betaK0 = vectData(2);
  betaKc = vectData(3);
  drillRatio = vectData(4);
  appliedB[0] = vectData(5);
  appliedB[1] = vectData(6);
  appliedB[2] = vectData(7);

  for (int i = 0; i < numGauss; i++) {
    int secClassTag = idData(i);
    int secDbTag = idData(i + 3);
    if (sections[i] == 0 || sections[i]->getClassTag() != secClassTag) {
      if (sections[i] != 0)
        delete sections[i];
      sections[i] = theBroker.getNewSection(secClassTag);
      if (sections[i] == 0) {
        opserr << "ShellDKMT::recvSelf() - Broker could not create Section of class type " << secClassTag << endln;
        return -1;
      }
    }
    sections[i]->setDbTag(secDbTag);
    res += sections[i]->recvSelf(commitTag, theChannel, theBroker);
    if (res < 0) {
      opserr << "ShellDKMT::recvSelf() - section " << i << " failed to recv itself\n";
      return res;
    }
  }

  return res;
}

Response *ShellDKMT::setResponse(const char **argv, int argc, OPS_Stream &output)
{
  Response *theResponse = 0;

  output.tag("ElementOutput");
  output.attr("eleType", "ShellDKMT");
  output.attr("eleTag", this->getTag());
  static char nodeData[32];
  for (int i = 0; i < 3; i++) {
    sprintf(nodeData, "node%d", i + 1);
    output.attr(nodeData, connectedExternalNodes(i));
  }

  if (strcmp(argv[0], "force") == 0 || strcmp(argv[0], "forces") == 0 ||
      strcmp(argv[0], "globalForce") == 0 || strcmp(argv[0], "globalForces") == 0) {
    for (int i = 0; i < 18; i++) {
      sprintf(nodeData, "P%d", i + 1);
      output.tag("ResponseType", nodeData);
    }
    theResponse = new ElementResponse(this, 1, resid);
  }
  else if (strcmp(argv[0], "material") == 0 || strcmp(argv[0], "Material") == 0 ||
           strcmp(argv[0], "section") == 0) {
    if (argc < 2) {
      opserr << "ShellDKMT::setResponse() - need to specify more data\n";
      output.endTag();
      return 0;
    }
    int pointNum = atoi(argv[1]);
    if (pointNum > 0 && pointNum <= numGauss) {
      output.tag("GaussPoint");
      output.attr("number", pointNum);
      output.attr("eta", gxi[pointNum - 1]);
      output.attr("neta", geta[pointNum - 1]);
      theResponse = sections[pointNum - 1]->setResponse(&argv[2], argc - 2, output);
      output.endTag();
    }
  }
  else if (strcmp(argv[0], "stresses") == 0 || strcmp(argv[0], "strains") == 0) {
    bool isStress = strcmp(argv[0], "stresses") == 0;
    static const char *sNames[8] = {"p11", "p22", "p12", "m11", "m22", "m12", "q1", "q2"};
    static const char *eNames[8] = {"eps11", "eps22", "gamma12", "theta11", "theta22", "theta33", "gamma13", "gamma23"};
    for (int i = 0; i < numGauss; i++) {
      output.tag("GaussPoint");
      output.attr("number", i + 1);
      output.attr("eta", gxi[i]);
      output.attr("neta", geta[i]);
      output.tag("SectionForceDeformation");
      output.attr("classType", sections[i]->getClassTag());
      output.attr("tag", sections[i]->getTag());
      for (int j = 0; j < 8; j++)
        output.tag("ResponseType", isStress ? sNames[j] : eNames[j]);
      output.endTag(); // SectionForceDeformation
      output.endTag(); // GaussPoint
    }
    theResponse = new ElementResponse(this, isStress ? 2 : 3, Vector(8 * numGauss));
  }

  output.endTag();
  return theResponse;
}

int ShellDKMT::getResponse(int responseID, Information &eleInfo)
{
  static Vector values(8 * numGauss);

  switch (responseID) {
  case 1: // global forces
    return eleInfo.setVector(this->getResistingForce());

  case 2: // section forces
    for (int i = 0; i < numGauss; i++) {
      const Vector &sigma = sections[i]->getStressResultant();
      for (int j = 0; j < 8; j++)
        values(8*i + j) = sigma(j);
    }
    return eleInfo.setVector(values);

  case 3: // section deformations
    for (int i = 0; i < numGauss; i++) {
      const Vector &eps = sections[i]->getSectionDeformation();
      for (int j = 0; j < 8; j++)
        values(8*i + j) = eps(j);
    }
    return eleInfo.setVector(values);

  default:
    return -1;
  }
}

int ShellDKMT::setParameter(const char **argv, int argc, Parameter &param)
{
  int res = -1;
  for (int i = 0; i < numGauss; i++) {
    int secRes = sections[i]->setParameter(argv, argc, param);
    if (secRes != -1)
      res = secRes;
  }
  return res;
}

int ShellDKMT::displaySelf(Renderer &theViewer, int displayMode, float fact,
                           const char **displayModes, int numModes)
{
  static Vector v1(3);
  static Vector v2(3);
  static Vector v3(3);
  nodePointers[0]->getDisplayCrds(v1, fact, displayMode);
  nodePointers[1]->getDisplayCrds(v2, fact, displayMode);
  nodePointers[2]->getDisplayCrds(v3, fact, displayMode);

  static Matrix coords(3, 3);
  for (int i = 0; i < 3; i++) {
    coords(0, i) = v1(i);
    coords(1, i) = v2(i);
    coords(2, i) = v3(i);
  }

  // display mode 1-8 -> stress resultant at the mid-side points
  static Vector values(3);
  for (int i = 0; i < 3; i++)
    values(i) = 0.0;
  if (displayMode > 0 && displayMode <= 8)
    for (int i = 0; i < 3; i++)
      values(i) = sections[i]->getStressResultant()(displayMode - 1);

  return theViewer.drawPolygon(coords, values, this->getTag());
}
