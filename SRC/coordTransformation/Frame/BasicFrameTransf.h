//===----------------------------------------------------------------------===//
//
//                                   xara
//
//===----------------------------------------------------------------------===//
//                              https://xara.so
//===----------------------------------------------------------------------===//
// 
// The purpose of this class is to wrap the more general FrameTransform<>
// templates to reproduce the legacy CrdTransf classes that were derived
// for elements in a "basic" coordinate system.
//
// cmp
//
#ifndef BasicFrameTransf3d_h
#define BasicFrameTransf3d_h

#include <array>
#include <FrameTransform.h>
#include <Vector.h>
#include <Matrix.h>

namespace OpenSees {

template<int ndf=6>
class BasicFrameTransf3d: public CrdTransf
{
public:
    BasicFrameTransf3d(FrameTransform<2,ndf> *t);

    ~BasicFrameTransf3d();

    virtual int getLocalAxes(Vector &x, Vector &y, Vector &z);

    virtual CrdTransf *getCopy3d() final;

    virtual double getInitialLength();
    virtual double getDeformedLength();

    virtual int initialize(Node *ni, Node *nj) final;
    virtual int update() final;
    virtual int commitState() final;
    virtual int revertToLastCommit() final;
    virtual int revertToStart() final;

    virtual const Vector &getBasicTrialDisp() final;
    virtual const Vector &getBasicIncrDisp() final;
    virtual const Vector &getBasicIncrDeltaDisp() final;
    virtual const Vector &getBasicTrialVel() final;

    virtual const Vector &getGlobalResistingForce(const Vector &basicForce, const Vector &p0) final;
    virtual const Matrix &getGlobalStiffMatrix(const Matrix &basicStiff, const Vector &basicForce) final;
    virtual const Matrix &getInitialGlobalStiffMatrix(const Matrix &basicStiff) final;

    // rotate consistent mass matrix
    const Matrix &getGlobalMatrixFromLocal(const Matrix &local);
    
    // methods used in post-processing only
    const Vector &getPointGlobalCoordFromLocal(const Vector &localCoords);
    const Vector &getPointGlobalDisplFromBasic(double xi, const Vector &basicDisps);
    const Vector &getPointLocalDisplFromBasic( double xi, const Vector &basicDisps);    

    //
    // Sensitivity
    //
    const Vector & getBasicDisplFixedGrad();
    const Vector & getBasicDisplTotalGrad(int grad);
    const Vector &getGlobalResistingForceShapeSensitivity (const Vector &basicForce, const Vector &p0, int grad);
    bool isShapeSensitivity();
    double getLengthGrad();
    double getd1overLdh();


    // MovableObject
    virtual int sendSelf(int tag, Channel &);
    virtual int recvSelf(int tag, Channel &, FEM_ObjectBroker &);
    const char *getClassType() const {return "BasicFrameTransf3d";}
    
    // TaggedObject
    void Print(OPS_Stream &s, int flag = 0);


    FrameTransform<2,ndf> &t;
protected:
private:

    constexpr static int NBV = 6;
    constexpr static int NDF = ndf;
    enum : int {
        inx = -12, //  0
        iny = -12, //  1
        inz = -12, //  2
        imx = -12, //  3
        imy =   3, //  4
        imz =   1, //  5
        jnx =   0, //  6
        jny = -12, //  7
        jnz = -12, //  8
        jmx =   5, //  9
        jmy =   4, // 10
        jmz =   2, // 11
    };

    constexpr static int iq[] = {
        // jnx, imz, jmz, imy, jmy, imx
        inx, iny, inz, imx, imy, imz,
        jnx, jny, jnz, jmx, jmy, jmz
    };

};
} // namespace OpenSees
#include "BasicFrameTransf.tpp"
#endif

