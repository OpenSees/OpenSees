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

// Jose Abell (UANDES, github.com/jaabell)
// Massimo Petracca - ASDEA Software, Italy (2022)
//

#ifndef EigenAPI_typedefs_h
#define EigenAPI_typedefs_h

#include "EigenAPI.h"
#include "Eigen/Dense"
#include <Eigen/Eigenvalues>

#include <iostream>
#include <typeinfo>

#include <Vector.h>

// define fixed size 6x1 vectors and 6x6 tensors for voigt notation
namespace Eigen {
typedef Eigen::Matrix<double, 1, 1> Vector1d;
typedef Eigen::Matrix<double, 6, 1> Vector6d;
typedef Eigen::Matrix<double, 6, 6> Matrix6d;
}

// namespace opsEigen {

//     // template<typename Scalar> struct scalar_identity_op_T4 {
//     //     EIGEN_EMPTY_STRUCT_CTOR(scalar_identity_op_T4)
//     //         template<typename IndexType>
//     //     EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar operator() (IndexType row, IndexType col) const { return row < 3 && col < 3 ? Scalar(1) : Scalar(0); }
//     // };

//     // typedef Eigen::CwiseNullaryOp<scalar_identity_op_T4<double>, Eigen::Matrix6d::PlainObject> IdentityReturnTypeT4;

//     template<typename Scalar> struct kronecker_delta {
//         EIGEN_EMPTY_STRUCT_CTOR(kronecker_delta)
//             template<typename IndexType>
//         EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar operator() (IndexType row) const { return row < 3 ? Scalar(1) : Scalar(0); }
//     };

//     typedef Eigen::CwiseNullaryOp<kronecker_delta<double>, Eigen::Vector6d::PlainObject> KroneckerDeltaReturnType;

// }


// inherit from matrix
class VoigtScalar : public Eigen::Vector1d
{
public:


    // empty constructor
    VoigtScalar() : Eigen::Vector1d() {}

    // full constructor
    VoigtScalar(double value)
        : Eigen::Vector1d()
    {
        *this << value;
    }

    // This constructor allows you to construct VoigtScalar from Eigen expressions
    template<typename OtherDerived>
    VoigtScalar(const Eigen::MatrixBase<OtherDerived>& other)
        : Eigen::Vector1d(other)
    {
        // std::cout << "c-tor from EXPR: " << typeid(other).name() << "\n";
    }

    // This method allows you to assign Eigen expressions to VoigtScalar
    template<typename OtherDerived>
    VoigtScalar& operator=(const Eigen::MatrixBase<OtherDerived>& other)
    {
        // std::cout << "assignment from EXPR: " << typeid(other).name() << "\n";
        this->Eigen::Vector1d::operator=(other);
        return *this;
    }
};


class VoigtVector;

EIGEN_STRONG_INLINE VoigtVector kronecker_delta(); //Fwd

// inherit from matrix
class VoigtVector : public Eigen::Vector6d
{
public:

    // empty constructor
    VoigtVector() : Eigen::Vector6d() {}

    // full constructor
    VoigtVector(double v0, double v1, double v2, double v3, double v4, double v5)
        : Eigen::Vector6d()
    {
        *this << v0, v1, v2, v3, v4, v5;
    }

    // from strain (opensees vector)
    static VoigtVector fromStrain(const Vector& strain) {
        return VoigtVector(strain(0), strain(1), strain(2), strain(3) , strain(4) , strain(5) );
    }

    // from stress (opensees vector)
    static VoigtVector fromStress(const Vector& stress) {
        return VoigtVector(stress(0), stress(1), stress(2), stress(3), stress(4) , stress(5) );
    }

    void toStrain(Vector& strain) {
        strain(0) = v11();
        strain(1) = v22();
        strain(2) = v33();
        strain(3) = v12();
        strain(4) = v23();
        strain(5) = v13();
    }

    void toStress(Vector& stress) {
        stress(0) = v11();
        stress(1) = v22();
        stress(2) = v33();
        stress(3) = v12();
        stress(4) = v23();
        stress(5) = v13();
    }


    // This constructor allows you to construct VoigtVector from Eigen expressions
    template<typename OtherDerived>
    VoigtVector(const Eigen::MatrixBase<OtherDerived>& other)
        : Eigen::Vector6d(other)
    {
        // std::cout << "c-tor from EXPR: " << typeid(other).name() << "\n";
    }

    // This method allows you to assign Eigen expressions to VoigtVector
    template<typename OtherDerived>
    VoigtVector& operator=(const Eigen::MatrixBase<OtherDerived>& other)
    {
        // std::cout << "assignment from EXPR: " << typeid(other).name() << "\n";
        this->Eigen::Vector6d::operator=(other);
        return *this;
    }

    // friend ostream & operator << (ostream &out, const VoigtVector &c);
    // friend istream & operator >> (istream &in,  VoigtVector &c);

    /***********************************************************************

    ConsStrainVectorVoigttructors and assignment operators needed for inheriting from Eigen

    @todo: override << operator for copying from OpenSees
    @todo: override >> operator for copying to OpenSees

    ************************************************************************/

public:

    /***********************************************************************

    Specific methods for working with this as if it were a 3x3 tensor

    ************************************************************************/

    EIGEN_STRONG_INLINE double v11()const { return this->operator()(0); }
    EIGEN_STRONG_INLINE double v22()const { return this->operator()(1); }
    EIGEN_STRONG_INLINE double v33()const { return this->operator()(2); }
    EIGEN_STRONG_INLINE double v12()const { return this->operator()(3); }
    EIGEN_STRONG_INLINE double v23()const { return this->operator()(4); }
    EIGEN_STRONG_INLINE double v13()const { return this->operator()(5); }

    EIGEN_STRONG_INLINE double squaredNorm() const {
        return v11() * v11() + v22() * v22() + v33() * v33() +
               2.0 * (v12() * v12() + v23() * v23() + v13() * v13());
    }
    EIGEN_STRONG_INLINE double norm() const { return std::sqrt(this->squaredNorm()); }

    EIGEN_STRONG_INLINE double trace() const { return v11() + v22() + v33(); }

    EIGEN_STRONG_INLINE VoigtVector deviator() const
    {
        double tr_over_3 = this->trace() / 3;
        return *this - tr_over_3 * kronecker_delta();
    };



public:
    Eigen::Vector3d principalStresses() const {
        // Convert Voigt vector to a stress tensor.
        Eigen::Matrix3d stress_tensor;
        stress_tensor << v11(), v12(), v13(),
                      v12(), v22(), v23(),
                      v13(), v23(), v33();

        // Compute the eigenvalues (principal stresses).
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen_solver(stress_tensor);
        return eigen_solver.eigenvalues();
    }

    double meanStress() const {
        return this->trace() / 3.0;
    }

    double stressDeviatorQ() const {
        VoigtVector deviator = this->deviator();
        double J2 = 0.5 * (
                        deviator.v11() * deviator.v11() + deviator.v22() * deviator.v22() + deviator.v33() * deviator.v33()
                        - 2 * deviator.v11() * deviator.v22()
                        - 2 * deviator.v22() * deviator.v33()
                        - 2 * deviator.v33() * deviator.v11()
                        + 6 * (deviator.v12() * deviator.v12() + deviator.v23() * deviator.v23() + deviator.v13() * deviator.v13())
                    );

        return std::sqrt(0.5 * J2);
    }

    double lodeAngle() const {
        Eigen::Vector3d principal_stresses = this->principalStresses();
        double s1 = principal_stresses(0);
        double s2 = principal_stresses(1);
        double s3 = principal_stresses(2);

        double num = std::sqrt(3) * (s1 - s3) ;
        double den = 2 * (s2 - s3) - (s1 - s3);
        if (abs(den) < 1e-20 )
        {
            return 0;
        }
        else
        {
            return std::atan( num / den );
        }
    }

};

EIGEN_STRONG_INLINE VoigtVector kronecker_delta()
{
    return VoigtVector(1.0, 1.0, 1.0, 0, 0, 0);
}


// inherit from matrix
class VoigtMatrix : public Eigen::Matrix6d
{
public:

    /***********************************************************************

    ConsStrainVectorVoigttructors and assignment operators needed for inheriting from Eigen

    @todo: override << operator for copying from OpenSees
    @todo: override >> operator for copying to OpenSees

    ************************************************************************/

    // empty constructor
    VoigtMatrix() : Eigen::Matrix6d() {}

    // full constructor
    VoigtMatrix(double values[36])
        : Eigen::Matrix6d()
    {
        std::cout << "VoigtMatrix CTOR\n";
        for (int i = 0; i < 36; ++i)
        {
            *this << values[i];
        }
    }



    // This constructor allows you to construct VoigtVector from Eigen expressions
    template<typename OtherDerived>
    VoigtMatrix(const Eigen::MatrixBase<OtherDerived>& other)
        : Eigen::Matrix6d(other)
    {
        // std::cout << "c-tor from EXPR: " << typeid(other).name() << "\n";
    }

    // This method allows you to assign Eigen expressions to VoigtVector
    template<typename OtherDerived>
    VoigtMatrix& operator=(const Eigen::MatrixBase<OtherDerived>& other)
    {
        // std::cout << "assignment from EXPR: " << typeid(other).name() << "\n";
        this->Eigen::Matrix6d::operator=(other);
        return *this;
    }

public:
};


#endif // EigenAPI_typedefs_h<
