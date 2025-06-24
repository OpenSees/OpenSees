//===----------------------------------------------------------------------===//
//
//                                   xara
//
//===----------------------------------------------------------------------===//
//                              https://xara.so
//===----------------------------------------------------------------------===//
//
// Description: This file contains the class definition for
// PDeltaFrameTransf.h. PDeltaFrameTransf provides the
// abstraction of a linear transformation for a spatial frame
// between the global and basic coordinate systems
//
// Adapted: Remo Magalhaes de Souza (rmsouza@ce.berkeley.edu)
// Created: 04/2000
//
#ifndef PDeltaFrameTransf_h
#define PDeltaFrameTransf_h

#include <array>
#include <FrameTransform.h>
#include <Vector.h>
#include <Matrix.h>

template <int nn, int ndf>
class PDeltaFrameTransf: public FrameTransform<nn,ndf>
{
public:

    PDeltaFrameTransf(int tag, 
                      const Vector3D &vecxz,
                      const std::array<Vector3D, nn> *offset=nullptr,
                      int offset_flags = 0);

    ~PDeltaFrameTransf();
    
    const char *getClassType() const {return "PDeltaFrameTransf";}
    
    double getInitialLength();
    double getDeformedLength();

    int initialize(std::array<Node*, nn>& new_nodes) final;
    int update() final;
    int commit() final;
    int revertToLastCommit() final;
    int revertToStart() final;

    VectorND<nn*ndf> getStateVariation() final;
    Vector3D getNodePosition(int tag) final;
    Vector3D getNodeRotationLogarithm(int tag) final;
    const std::array<Vector3D,nn> *getRigidOffsets() const final { return linear.getRigidOffsets();}

    virtual VectorND<nn*ndf>    pushResponse(VectorND<nn*ndf>&pl) final;
    virtual MatrixND<nn*ndf,nn*ndf> pushResponse(MatrixND<nn*ndf,nn*ndf>& kl, const VectorND<nn*ndf>& pl) final;

    FrameTransform<nn,ndf> *getCopy() const final;

    int getLocalAxes(Vector3D &x, Vector3D &y, Vector3D &z) const final;

    // Sensitivity
    double getLengthGrad() final;

    // Tagged Object
    void Print(OPS_Stream &s, int flag = 0) final;

    private:
      int offset_flags;
      LinearFrameTransf<nn,ndf> linear;

};
#include "PDeltaFrameTransf3d.tpp"
#endif
