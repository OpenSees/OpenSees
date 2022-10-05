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

#include <iostream>
#include <typeinfo>

// define fixed size 6x1 vectors and 6x6 tensors for voigt notation
namespace Eigen {
    typedef Eigen::Matrix<double, 6, 1> Vector6d;
    typedef Eigen::Matrix<double, 6, 6> Matrix6d;
}

// namespace opsEigen {

//     template<typename Scalar> struct scalar_identity_op_T4 {
//         EIGEN_EMPTY_STRUCT_CTOR(scalar_identity_op_T4)
//             template<typename IndexType>
//         EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const Scalar operator() (IndexType row, IndexType col) const { return row < 3 && col < 3 ? Scalar(1) : Scalar(0); }
//     };

//     typedef Eigen::CwiseNullaryOp<scalar_identity_op_T4<double>, Eigen::Matrix6d::PlainObject> IdentityReturnTypeT4;

// }


// inherit from matrix
class VoigtVector : public Eigen::Vector6d
{
public:

    /***********************************************************************

    ConsStrainVectorVoigttructors and assignment operators needed for inheriting from Eigen

    @todo: override << operator for copying from OpenSees
    @todo: override >> operator for copying to OpenSees

    ************************************************************************/

    // empty constructor
    VoigtVector() : Eigen::Vector6d() {}

    // full constructor
    VoigtVector(double e11, double e22, double e33, double g12, double g23, double g13)
        : Eigen::Vector6d()
    {
        std::cout << "CTOR\n";
        *this << e11, e22, e33, g12, g23, g13;
    }

    // from strain (opensees vector)
    static VoigtVector fromStrain(const Eigen::Vector6d& op) {
        std::cout << "FACT\n";
        return VoigtVector(op(0), op(1), op(2), op(3) / 2.0, op(4) / 2.0, op(5) / 2.0);
    }

    void toStrain(Eigen::Vector6d& op) {
        // copy to
    }

    // This constructor allows you to construct MyVectorType from Eigen expressions
    template<typename OtherDerived>
    VoigtVector(const Eigen::MatrixBase<OtherDerived>& other)
        : Eigen::Vector6d(other)
    {
        std::cout << "c-tor from EXPR: " << typeid(other).name() << "\n";
    }

    // This method allows you to assign Eigen expressions to MyVectorType
    template<typename OtherDerived>
    VoigtVector& operator=(const Eigen::MatrixBase<OtherDerived>& other)
    {
        std::cout << "assignment from EXPR: " << typeid(other).name() << "\n";
        this->Eigen::Vector6d::operator=(other);
        return *this;
    }

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

    // static EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE const opsEigen::IdentityReturnTypeT4
    //     Identity()
    // {
    //         return Eigen::Matrix6d::NullaryExpr(6, 6, opsEigen::scalar_identity_op_T4<double>());
    // }
};

#endif // EigenAPI_typedefs_h