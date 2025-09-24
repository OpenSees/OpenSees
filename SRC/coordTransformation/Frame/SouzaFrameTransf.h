//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//
//===----------------------------------------------------------------------===//
//
//        OpenSees - Open System for Earthquake Engineering Simulation
//
//===----------------------------------------------------------------------===//
//
// Copyright (c) 2025, Claudio M. Perez
// All rights reserved.  No warranty, explicit or implicit, is provided.
//
// This source code is licensed under the BSD 2-Clause License.
// See LICENSE file or https://opensource.org/licenses/BSD-2-Clause
//
// Please cite the following resource in any derivative works:
//
// Perez, C.M., and Filippou F.C.. "On Nonlinear Geometric Transformations
//     of Finite Elements" Int. J. Numer. Meth. Engrg. 2024; 
//     https://doi.org/10.1002/nme.7506
//
//===----------------------------------------------------------------------===//

//
// Copyright @ 1999-2020 The Regents of the University of California (The Regents). 
// All Rights Reserved.
//
// The Regents grants permission, without fee and without a written license
// agreement, for (a) use, reproduction, modification, and distribution of this
// software and its documentation by educational, research, and non-profit
// entities for noncommercial purposes only; and (b) use, reproduction and
// modification of this software by other entities for internal purposes only.
// The above copyright notice, this paragraph and the following three
// paragraphs must appear in all copies and modifications of the software
// and/or documentation.
//
// Permission to incorporate this software into products for commercial 
// distribution may be obtained by contacting the University of California 
// Office of Technology Licensing 
// 2150 Shattuck Avenue #510, 
// Berkeley, CA 94720-1620, 
// (510) 643-7201.
//
// This software program and documentation are copyrighted by The Regents of
// the University of California. The Regents does not warrant that the
// operation of the program will be uninterrupted or error-free. The end-user
// understands that the program was developed for research purposes and is
// advised not to rely exclusively on the program for any reason.
//
// IN NO EVENT SHALL REGENTS BE LIABLE TO ANY PARTY FOR DIRECT, INDIRECT,
// SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
// ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF
// REGENTS HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.  REGENTS GRANTS
// NO EXPRESS OR IMPLIED LICENSE IN ANY PATENT RIGHTS OF REGENTS BUT HAS
// IMPLEMENTED AN INDIVIDUAL CONTRIBUTOR LICENSE AGREEMENT FOR THE OPENSEES
// PROJECT AT THE UNIVERISTY OF CALIFORNIA, BERKELEY TO BENEFIT THE END USER.
//
// REGENTS SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED
// TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR
// PURPOSE. THE SOFTWARE AND ACCOMPANYING DOCUMENTATION, IF ANY, PROVIDED
// HEREUNDER IS PROVIDED "AS IS". REGENTS HAS NO OBLIGATION TO PROVIDE
// MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
//
//===----------------------------------------------------------------------===//

//
// Description: SouzaFrameTransf implements the formulation 
// of Crisfield (1990) with the objective of maintaining the
// original "Corotational" implementation by Remo Magalhaes de Souza, within
// the new framework proposed by Perez and Filippou (2024).
//
// Written: cmp
// Created: March 2024
//
// Adapted from work by: Remo Magalhaes de Souza (rmsouza@ce.berkeley.edu)
//
// [1] Perez, C.M., and Filippou F.C.. "On Nonlinear Geometric Transformations
//     of Finite Elements" Int. J. Numer. Meth. Engrg. 2024; 
//     https://doi.org/10.1002/nme.7506
//
// [2] Crisfield, M.A. (1990) "A consistent co-rotational formulation for
//     non-linear, three-dimensional, beam-elements", Computer Methods in Applied
//     Mechanics and Engineering, 81(2), pp. 131–150. Available at:
//     https://doi.org/10.1016/0045-7825(90)90106-V.
//
#ifndef SouzaFrameTransf_hpp
#define SouzaFrameTransf_hpp

#include <array>
#include "FrameTransform.h"
#include <Vector.h>
#include <Matrix.h>
#include <Versor.h>
#include <Matrix3D.h>
#include <Vector3D.h>
#include "Isometry/CrisfieldIsometry.h"


namespace OpenSees {

template <int nn, int ndf>
class SouzaFrameTransf: public FrameTransform<nn,ndf>
{
public:
  SouzaFrameTransf(int tag, 
                   const Vector3D &vecxz,
                   const std::array<Vector3D, nn> *offset=nullptr,
                   int offset_flags = 0);

  ~SouzaFrameTransf();

  const char *getClassType() const {
    return "SouzaFrameTransf";
  }

  // NOTE: maybe add arg for rotation parameterization
  FrameTransform<nn,ndf> *getCopy() const final;

  int initialize(std::array<Node*, nn>& new_nodes) final;
  int update() final;
  int commit() final;
  int revertToLastCommit() final;        
  int revertToStart() final;
  int getLocalAxes(Vector3D &x, Vector3D &y, Vector3D &z) const final;
  const std::array<Vector3D,nn> *getRigidOffsets() const final { return offsets; }

  double getInitialLength() final;
  double getDeformedLength() final;

  VectorND<nn*ndf> getStateVariation() final;
  Vector3D getNodePosition(int tag) final;
  Vector3D getNodeRotationLogarithm(int tag) final;

  using Operation = typename FrameTransform<nn,ndf>::Operation;
  int push(VectorND<nn*ndf>&pl, Operation) final;
  int push(MatrixND<nn*ndf,nn*ndf>& kl, const VectorND<nn*ndf>& pl, Operation) final;

  // Sensitivity
  double getLengthGrad() final;
  const Vector &getBasicDisplTotalGrad(int grad); //  final;
  const Vector &getBasicDisplFixedGrad();
  const Vector &getGlobalResistingForceShapeSensitivity(const Vector &pb, const Vector &p0, int gradNumber);

  // Tagged Object
  void Print(OPS_Stream &s, int flag) final;

private:
  constexpr static int n = nn*ndf;

  int  addTangent(MatrixND<12,12>& M, const VectorND<12>& pl, const VectorND<12>&ul);

  enum {
    inx= 0, // axial
    iny= 1, // Vy
    inz= 2, // Vz
    imx= 3, // torsion
    imy= 4, // rot y I
    imz= 5, // rot z I

    jnx= 6, // axial
    jny= 7,
    jnz= 8,
    jmx= 9, // torsion
    jmy=10, // rot y J
    jmz=11, // rot z J
  };

  //
  // Member data
  //
  std::array<Node*, nn> nodes;

  Vector3D xAxis;                              // local x axis
  Vector3D vz;                                 // Vector that lies in local plane xz
  Vector3D dX;
  
  std::array<Vector3D, nn> *offsets;

  double *nodeIInitialDisp, *nodeJInitialDisp;
  bool  initialDispChecked;                    

  double L;                       // initial element length
  double Ln;                      // current element length (at trial state)

  Versor Q_past[nn];              // commited rotations
  Versor Q_pres[nn];              // trial rotations

  Vector3D alphaI;                // last trial rotations end i
  Vector3D alphaJ;                // last trial rotatations end j

  VectorND<n> ul;                 // local displacements
  Vector3D    vr[nn];             //
  VectorND<n> ulcommit;           // commited local displacements
  VectorND<n> ulpr;               // previous local displacements

  OpenSees::MatrixND<n,n> T;     // transformation from local to global system

  OpenSees::Matrix3D R0;         // rotation from local to global coordinates
  CrisfieldIsometry<2,false> crs;

};

} // namespace OpenSees

#include "SouzaFrameTransf.tpp"
#endif
