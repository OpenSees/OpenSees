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

    double maxAbs() const {
        double m = 0.0;
        m = std::max(m, std::abs(this->v11()));
        m = std::max(m, std::abs(this->v22()));
        m = std::max(m, std::abs(this->v33()));
        m = std::max(m, std::abs(this->v12()));
        m = std::max(m, std::abs(this->v23()));
        m = std::max(m, std::abs(this->v13()));
        return m;
    }

public:
    // Eigen::Vector3d principalStresses() const {
    //     // Convert Voigt vector to a stress tensor.
    //     Eigen::Matrix3d stress_tensor;
    //     stress_tensor << v11(), v12(), v13(),
    //                      v12(), v22(), v23(),
    //                      v13(), v23(), v33();

    //     // Compute the eigenvalues (principal stresses).
    //     Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen_solver(stress_tensor);
    //     Eigen::Vector3d eigenvalues = eigen_solver.eigenvalues();

    //     // Sort the eigenvalues in ascending order.
    //     std::sort(eigenvalues.data(), eigenvalues.data() + eigenvalues.size());

    //     return eigenvalues;
    // }

    std::tuple<double, double, double> principalStresses() const {
        // Convert Voigt vector to a stress tensor.
        Eigen::Matrix3d stress_tensor;
        stress_tensor << v11(), v12(), v13(),
                         v12(), v22(), v23(),
                         v13(), v23(), v33();

        // Compute the eigenvalues (principal stresses).
        Eigen::SelfAdjointEigenSolver<Eigen::Matrix3d> eigen_solver(stress_tensor);
        Eigen::Vector3d eigenvalues = eigen_solver.eigenvalues();

        // Sort the eigenvalues in ascending order.
        std::sort(eigenvalues.data(), eigenvalues.data() + eigenvalues.size());

        double s1, s2, s3;
        s1 = eigenvalues[0];
        s2 = eigenvalues[1];
        s3 = eigenvalues[2];

        return std::make_tuple(s1,s2,s3);
    }


    double meanStress() const {
        return this->trace() / 3.0;
    }

    double getI1() const
    {
         return this->trace();
    }

    // double J2_invariant() const {
    //     VoigtVector deviator = this->deviator();
    //     double J2 = 0.5 * (
    //                     deviator.v11() * deviator.v11() + deviator.v22() * deviator.v22() + deviator.v33() * deviator.v33()
    //                     - 2 * deviator.v11() * deviator.v22()
    //                     - 2 * deviator.v22() * deviator.v33()
    //                     - 2 * deviator.v33() * deviator.v11()
    //                     + 6 * (deviator.v12() * deviator.v12() + deviator.v23() * deviator.v23() + deviator.v13() * deviator.v13())
    //                 );
    //     return J2;
    // }

    double getJ2() const {
        VoigtVector deviator = this->deviator();
        double J2 = 0.5 * (
            deviator.v11() * deviator.v11() + 
            deviator.v22() * deviator.v22() + 
            deviator.v33() * deviator.v33()
            + 2 * (deviator.v12() * deviator.v12() + 
            	   deviator.v23() * deviator.v23() + 
            	   deviator.v13() * deviator.v13())
                  );
        return J2;
    }


    double getJ3() const {
        VoigtVector deviator = this->deviator();
        double J3 = deviator.v11() * deviator.v22() * deviator.v33()
                  - deviator.v11() * deviator.v23() * deviator.v23()
                  - deviator.v22() * deviator.v13() * deviator.v13()
                  - deviator.v33() * deviator.v12() * deviator.v12()
                  + 2 * deviator.v12() * deviator.v23() * deviator.v13();
        return J3;
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

    // double lodeAngle() const {
    //     Eigen::Vector3d principal_stresses = this->principalStresses();
    //     double s1 = principal_stresses(0);
    //     double s2 = principal_stresses(1);
    //     double s3 = principal_stresses(2);

    //     double num = std::sqrt(3) * (s1 - s3) ;
    //     double den = 2 * (s2 - s3) - (s1 - s3);
    //     if (abs(den) < 1e-20 )
    //     {
    //         return 0;
    //     }
    //     else
    //     {
    //         return std::atan( num / den );
    //     }
    // }


    // double lodeAngle() const {
    //     Eigen::Vector3d principal_stresses = this->principalStresses();
    //     double s1 = principal_stresses(0);
    //     double s2 = principal_stresses(1);
    //     double s3 = principal_stresses(2);

    //     // Sort principal stresses to ensure s1 >= s2 >= s3
    //     // if (s1 < s3) std::swap(s1, s3);
    //     // if (s1 < s2) std::swap(s1, s2);
    //     // if (s2 < s3) std::swap(s2, s3);

    //     // Calculate the numerator and denominator for the Lode angle
    //     double num = std::sqrt(3.0) * (s2 - s3);
    //     double den = 2.0 * (s2 - s3) - (s1 - s3);

    //     // Use atan2 for better numerical stability
    //     double lode_angle = std::atan2(num, den);

    //     return lode_angle;
    // }

    // double lodeAngle() const {
    //     // Eigen::Vector3d principal_stresses = this->principalStresses();
    //     // double s1 = principal_stresses(0);
    //     // double s2 = principal_stresses(1);
    //     // double s3 = principal_stresses(2);
    //     double J2 = this->getJ2();
    //     double J3 = this->getJ3();

    //     double rLodeAngle;

    //     if (J2 > std::numeric_limits<double>::epsilon()) {
    //         double sint3 = (-3.0 * std::sqrt(3.0) * J3) / (2.0 * J2 * std::sqrt(J2));
    //         if (sint3 < -0.95)
    //             sint3 = -1.0;
    //         else if (sint3 > 0.95)
    //             sint3 = 1.0;
    //         rLodeAngle = std::asin(sint3) / 3.0;
    //     } else {
    //         rLodeAngle = 0.0;
    //     }

    //     return rLodeAngle;
    // }

    double lodeAngle() const {
        double J2 = this->getJ2();
        double J3 = this->getJ3();
        if (J2 <= std::numeric_limits<double>::epsilon()) {
            return 0.0;
        }
        double x = (-3.0 * std::sqrt(3.0) * J3) / (2.0 * J2 * std::sqrt(J2));
        // clamp 
        if (x > 1.0) x = 1.0;
        else if (x < -1.0) x = -1.0;
        return std::asin(x) / 3.0;
    }


    #include <cmath>

    // double lodeAngle() const {
    //     double J2 = this->getJ2();
    //     double J3 = this->getJ3();

    //     // Ensure J2 is positive to avoid division by zero
    //     if (J2 < 0) {
    //         return 0.0; // Handle special case
    //     }

    //     // Calculate the Lode angle
    //     double sqrt_J2_3 = std::pow(J2, 1.5);  // J2 raised to 3/2
    //     double lode_angle = std::asin((3 * std::sqrt(3) * J3) / (2 * sqrt_J2_3)) / 3.0;

    //     return lode_angle;
    // }


};


// Stress-like Frobenius norm 
EIGEN_STRONG_INLINE double tensor_dot_stress_like(const VoigtVector& sigma1, const VoigtVector& sigma2)
{
    double the_dot_product = 
        sigma1.v11() * sigma2.v11() + 
        sigma1.v22() * sigma2.v22() + 
        sigma1.v33() * sigma2.v33()
        + 2 * (sigma1.v12() * sigma2.v12() + 
               sigma1.v23() * sigma2.v23() + 
               sigma1.v13() * sigma2.v13());

    return the_dot_product;
}

// Strain-like Frobenius norm 
EIGEN_STRONG_INLINE double tensor_dot_strain_like(const VoigtVector& epsilon1, const VoigtVector& epsilon2)
{
    double the_dot_product = 
        epsilon1.v11() * epsilon2.v11() + 
        epsilon1.v22() * epsilon2.v22() + 
        epsilon1.v33() * epsilon2.v33()
        + 2 * (epsilon1.v12() * epsilon2.v12() + 
               epsilon1.v23() * epsilon2.v23() + 
               epsilon1.v13() * epsilon2.v13());

    return the_dot_product;
}

// Energy-like Frobenius norm (regular dot product)
EIGEN_STRONG_INLINE double tensor_dot_energy_like(const VoigtVector& sigma, const VoigtVector& epsilon)
{
    double the_dot_product = sigma.dot(epsilon);

    return the_dot_product;
}



EIGEN_STRONG_INLINE VoigtVector kronecker_delta()
{
    return VoigtVector(1.0, 1.0, 1.0, 0, 0, 0);
}

EIGEN_STRONG_INLINE VoigtVector calculate_first_vector()
{
    return VoigtVector(1.0, 1.0, 1.0, 0, 0, 0);
}

EIGEN_STRONG_INLINE VoigtVector calculate_second_vector(const VoigtVector& sigma)
{
    const double J2 = sigma.getJ2();
    const double twosqrtJ2 = 2.0 * std::sqrt(J2);

    VoigtVector rSecondVector;

    if (twosqrtJ2 > 1e-20) {
        rSecondVector = sigma.deviator() / (twosqrtJ2);

        for (int i = 3; i < 6; ++i)
            rSecondVector[i] *= 2.0;
    } 

    return rSecondVector;
}


EIGEN_STRONG_INLINE VoigtVector calculate_third_vector(const VoigtVector& sigma)
{
    const double J2 = sigma.getJ2();
    const double twosqrtJ2 = 2.0 * std::sqrt(J2);
    const double J2thirds = J2 / 3.0;

    VoigtVector rThirdVector;
    VoigtVector rDeviator = sigma.deviator();

    rThirdVector[0] = rDeviator[1] * rDeviator[2] - rDeviator[4] * rDeviator[4] + J2thirds;
    rThirdVector[1] = rDeviator[0] * rDeviator[2] - rDeviator[5] * rDeviator[5] + J2thirds;
    rThirdVector[2] = rDeviator[0] * rDeviator[1] - rDeviator[3] * rDeviator[3] + J2thirds;
    rThirdVector[3] = 2.0 * (rDeviator[4] * rDeviator[5] - rDeviator[3] * rDeviator[2]);
    rThirdVector[4] = 2.0 * (rDeviator[3] * rDeviator[4] - rDeviator[1] * rDeviator[5]);
    rThirdVector[5] = 2.0 * (rDeviator[5] * rDeviator[3] - rDeviator[0] * rDeviator[4]);

    return rThirdVector;
}



// class VoigtMatrix;


// namespace Eigen {
// namespace internal {
// template<> struct traits<VoigtMatrix> : public traits<Eigen::Matrix<double,6,6>> {};
// }}



// class VoigtMatrix : public Eigen::Matrix<double,6,6> {
// public:
//     using Base = Eigen::Matrix<double,6,6>;
//     EIGEN_MAKE_ALIGNED_OPERATOR_NEW
//     EIGEN_DENSE_PUBLIC_INTERFACE(VoigtMatrix)  // o EIGEN_GENERIC_PUBLIC_INTERFACE(VoigtMatrix)
//     using Base::Base; // hereda ctors de Eigen

//     template<typename OtherDerived>
//     VoigtMatrix& operator=(const Eigen::MatrixBase<OtherDerived>& other) {
//         Base::operator=(other);
//         return *this;
//     }

//     // Si de verdad quieres este ctor desde arreglo:
//     // OJO: Eigen es column-major por defecto; si tu arreglo está en row-major, usa RowMajor en el Map.
//     explicit VoigtMatrix(const double* values_row_major)
//       : Base(Eigen::Map<const Eigen::Matrix<double,6,6,Eigen::RowMajor>>(values_row_major)) {}
// };


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
