//===----------------------------------------------------------------------===//
//
//                                   xara
//                              https://xara.so
//===----------------------------------------------------------------------===//
//
//                                  FEDEASLab
//       Finite Elements for Design Evaluation and Analysis of Structures
//
//===----------------------------------------------------------------------===//
//
// Description: This file contains the class definition for
// EuclidFrameTransf.h. EuclidFrameTransf provides the
// abstraction of a linear transformation for a spatial frame
// between the global and basic coordinate systems
//
// Written: cmp
// Created: 04/2025
//
#ifndef EuclidFrameTransf_hpp
#define EuclidFrameTransf_hpp

#include <array>
#include <FrameTransform.h>
#include <Vector3D.h>
#include <MatrixND.h>

template <int nn, int ndf, typename BasisT>
class EuclidFrameTransf: public FrameTransform<nn,ndf>
{
public:
    constexpr static int n = nn*ndf;

    EuclidFrameTransf(int tag, 
                      const Vector3D &vecxz,
                      const std::array<Vector3D, nn> *offset=nullptr,
                      int offset_flags = 0);

    ~EuclidFrameTransf();

    const char *getClassType() const {return "EuclidFrameTransf";}
    
    virtual int getLocalAxes(Vector3D &x, Vector3D &y, Vector3D &z) const;
    
    virtual FrameTransform<nn,ndf> *getCopy() const;

    double getInitialLength() final;
    double getDeformedLength() final;
    const std::array<Vector3D,nn> *getRigidOffsets() const final {return offsets;}
    
    int initialize(std::array<Node*, nn>& new_nodes) final;
    int update() final;
    int commit() final;
    int revertToLastCommit() final;
    int revertToStart() final;

    VectorND<nn*ndf> getStateVariation() final;
    Vector3D getNodePosition(int tag) final;
    Versor   getNodeRotation(int tag) /* final */;
    Vector3D getNodeRotationLogarithm(int tag) final;

    VectorND<nn*ndf>        pushResponse(VectorND<nn*ndf>&pl) final;
    MatrixND<nn*ndf,nn*ndf> pushResponse(MatrixND<nn*ndf,nn*ndf>& kl, const VectorND<nn*ndf>& pl) final;

#if 0
    // method used to rotate consistent mass matrix
    const Matrix &getGlobalMatrixFromLocal(const Matrix &local);
#endif

    // Sensitivity
    //
    bool isShapeSensitivity() final;
    double getLengthGrad() final;
    double getd1overLdh();

    // TaggedObject
    void Print(OPS_Stream &s, int flag = 0) override;

private:

    inline MatrixND<nn*ndf,nn*ndf> 
    getProjection() {

      MatrixND<nn*ndf,nn*ndf> A{};
      A.addDiagonal(1.0);

      // double L = basis.getLength();
      constexpr Vector3D axis{1, 0, 0};
      constexpr Matrix3D ix = Hat(axis);
      MatrixND<3,ndf> Gb{};
      for (int a = 0; a<nn; a++) {
        for (int b = 0; b<nn; b++) {
          
          Gb.template insert<0,0>(basis.getRotationGradient(b), 1.0);
          // TODO(nn>2): Interpolate coordinate?
          A.assemble(ix*Gb, a*ndf  , b*ndf,  double(a)/double(nn-1)*L);
          A.assemble(   Gb, a*ndf+3, b*ndf, -1.0);
        }
      }

      return A;
    }

    int computeElemtLengthAndOrient();

    template<const Vector& (Node::*Getter)()>
    const Vector3D
    pullPosition(int node)
    {
        const Vector &u = (nodes[node]->*Getter)();

        Vector3D v;
        for (int i=0; i<3; i++)
          v[i] = u[i];

        // 1) Offsets
        if (offsets) [[unlikely]] {
          if (!(offset_flags&OffsetLocal))  {
            Vector3D w {u[3], u[4], u[5]};
            v -= offsets->at(node).cross(w);
          }
        }

        // 2) Constant Rotation
        Matrix3D R = basis.getRotation();
        return R^v;
    }

    std::array<Node*, nn> nodes;
    std::array<Vector3D, nn> ur; // rotation vector
    std::array<Vector3D, nn> ux; // displacement vector

    std::array<Vector3D, nn> *offsets;
    int offset_flags;
    Matrix3D R0;
    Vector3D xi, xj, vz;
    double L;           // undeformed element length

    BasisT basis;
};

#include "EuclidFrameTransf.tpp"
#endif

