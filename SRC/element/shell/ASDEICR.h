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

// $Revision: 1.10 $
// $Date: 2020/05/18 22:51:21 $

// Original implementation: Massimo Petracca (ASDEA)
//
// Implementation of the EICR (Element Independent CoRotational) formulation
// 
// References:
//
// Felippa, Carlos A. 
// "A systematic approach to the element-independent corotational 
// dynamics of finite elements". 
// Technical Report CU-CAS-00-03, Center for Aerospace Structures, 2000.
//
// Felippa, Carlos A., and Bjorn Haugen. 
// "A unified formulation of small-strain corotational finite elements: I. Theory." 
// Computer Methods in Applied Mechanics and Engineering 194.21-24 (2005): 2285-2335.
//

#ifndef ASDEICR_h
#define ASDEICR_h

#include <ASDMath.h>
#include <vector>

/** \brief EICR Element Independent CoRotational formulation
*
* E.I.C.R. is a utility class containing static methods related to
* the Element Independent Corotational Formulation.
* This class implements methods that do not depend on the element type,
* and so they can be used by any implementation of a corotational coordinate transformation.
*/
class EICR
{

public:

	typedef std::size_t size_t;

	typedef ASDVector3<double> Vector3Type;

	typedef std::vector<Vector3Type> NodeContainerType;

	typedef Vector VectorType;

	typedef Matrix MatrixType;

	typedef ASDQuaternion<double> QuaternionType;

public:

	/**
	* Computes the Spin of the input vector V, and saves the result into the output matrix S.
	* Note: no check is made on the size of the input-output arguments.
	* @param V the input vector (assumed size: >= 3)
	* @param S the output matrix (assumed size: >= 3x3)
	*/
	template< class TVec, class TMat>
	inline static void Spin(const TVec& V, TMat& S)
	{
		S(0, 0) = 0.00;		S(0, 1) = -V(2);	S(0, 2) = V(1);
		S(1, 0) = V(2);		S(1, 1) = 0.00;		S(1, 2) = -V(0);
		S(2, 0) = -V(1);	S(2, 1) = V(0);		S(2, 2) = 0.00;
	}

	/**
	* Computes the Spin of the input vector V, and saves the result into the output matrix S,
	* at the specified row index.
	* Note: no check is made on the size of the input-output arguments.
	* @param V the input vector (assumed size: >= 3)
	* @param S the output matrix (assumed size: >= 3x3)
	* @param row_index the index of the first row in the output matrix where the spin has to be saved
	*/
	template< class TVec, class TMat>
	inline static void Spin_AtRow(const TVec& V, TMat& S, size_t row_index)
	{
		size_t i0 = row_index;
		size_t i1 = 1 + row_index;
		size_t i2 = 2 + row_index;
		double v0 = V(i0);
		double v1 = V(i1);
		double v2 = V(i2);
		S(i0, 0) = 0.00;	S(i0, 1) = -v2;		S(i0, 2) = v1;
		S(i1, 0) = v2;		S(i1, 1) = 0.00;	S(i1, 2) = -v0;
		S(i2, 0) = -v1;		S(i2, 1) = v0;		S(i2, 2) = 0.00;
	}

	/**
	* Computes the Spin of the input vector V, from the specified index, and saves the result into the output matrix S,
	* at the specified row index.
	* Note: no check is made on the size of the input-output arguments.
	* @param V the input vector (assumed size: >= 3)
	* @param S the output matrix (assumed size: >= 3x3)
	* @param vector_index the index of the first component of the input vector to be used to compute the spin
	* @param row_index the index of the first row in the output matrix where the spin has to be saved
	*/
	template< class TVec, class TMat>
	inline static void Spin_AtRow(const TVec& V, TMat& S, size_t vector_index, size_t matrix_row_index)
	{
		size_t i0 = matrix_row_index;
		size_t i1 = 1 + matrix_row_index;
		size_t i2 = 2 + matrix_row_index;
		double v0 = V(vector_index);
		double v1 = V(vector_index + 1);
		double v2 = V(vector_index + 2);
		S(i0, 0) = 0.00;	S(i0, 1) = -v2;		S(i0, 2) = v1;
		S(i1, 0) = v2;		S(i1, 1) = 0.00;	S(i1, 2) = -v0;
		S(i2, 0) = -v1;		S(i2, 1) = v0;		S(i2, 2) = 0.00;
	}

	/**
	* Computes the Spin of the input vector V, and saves the result into the output matrix S.
	* This version uses a multiplier for the output values.
	* Note: no check is made on the size of the input-output arguments.
	* @param V the input vector (assumed size: >= 3)
	* @param S the output matrix (assumed size: >= 3x3)
	* @param mult the multiplier for the output values
	*/
	template< class TVec, class TMat>
	inline static void Spin(const TVec& V, TMat& S, double mult)
	{
		S(0, 0) = 0.00;			S(0, 1) = -mult * V(2);	S(0, 2) = mult * V(1);
		S(1, 0) = mult * V(2);	S(1, 1) = 0.00;			S(1, 2) = -mult * V(0);
		S(2, 0) = -mult * V(1);	S(2, 1) = mult * V(0);	S(2, 2) = 0.00;
	}

	/**
	* Computes the Spin of the input vector V, and saves the result into the output matrix S,
	* at the specified row index.
	* This version uses a multiplier for the output values.
	* Note: no check is made on the size of the input-output arguments.
	* @param V the input vector (assumed size: >= 3)
	* @param S the output matrix (assumed size: >= 3x3)
	* @param mult the multiplier for the output values
	* @param row_index the index of the first row in the output matrix where the spin has to be saved
	*/
	template< class TVec, class TMat>
	inline static void Spin_AtRow(const TVec& V, TMat& S, double mult, size_t row_index)
	{
		size_t i0 = row_index;
		size_t i1 = 1 + row_index;
		size_t i2 = 2 + row_index;
		double v0 = mult * V(i0);
		double v1 = mult * V(i1);
		double v2 = mult * V(i2);
		S(i0, 0) = 0.00;	S(i0, 1) = -v2;		S(i0, 2) = v1;
		S(i1, 0) = v2;		S(i1, 1) = 0.00;	S(i1, 2) = -v0;
		S(i2, 0) = -v1;		S(i2, 1) = v0;		S(i2, 2) = 0.00;
	}

	/**
	* Computes the Spin of the input vector V, from the specified index, and saves the result into the output matrix S,
	* at the specified row index.
	* This version uses a multiplier for the output values.
	* Note: no check is made on the size of the input-output arguments.
	* @param V the input vector (assumed size: >= 3)
	* @param S the output matrix (assumed size: >= 3x3)
	* @param mult the multiplier for the output values
	* @param vector_index the index of the first component of the input vector to be used to compute the spin
	* @param row_index the index of the first row in the output matrix where the spin has to be saved
	*/
	template< class TVec, class TMat>
	inline static void Spin_AtRow(const TVec& V, TMat& S, double mult, size_t vector_index, size_t matrix_row_index)
	{
		size_t i0 = matrix_row_index;
		size_t i1 = 1 + matrix_row_index;
		size_t i2 = 2 + matrix_row_index;
		double v0 = mult * V(vector_index);
		double v1 = mult * V(vector_index + 1);
		double v2 = mult * V(vector_index + 2);
		S(i0, 0) = 0.00;	S(i0, 1) = -v2;		S(i0, 2) = v1;
		S(i1, 0) = v2;		S(i1, 1) = 0.00;	S(i1, 2) = -v0;
		S(i2, 0) = -v1;		S(i2, 1) = v0;		S(i2, 2) = 0.00;
	}

	/**
	* Sets the input matrix to the zero matrix of requested size. Resize is done if necessary
	* @param n the number of rows
	* @param n the number of columns
	* @param I the input/output matrix
	*/
	inline static void SetZero(size_t n, size_t m, MatrixType& I)
	{
		if ((I.noRows() != n) || (I.noCols() != m))
			I.resize(n, m);
		I.Zero();
	}

	/**
	* Sets the input matrix to the identity of requested size. Resize is done if necessary
	* @param n the number of rows and columns
	* @param I the input/output matrix
	* @param value (optional, default = 1) the value on the diagonal terms
	*/
	inline static void SetIdentity(size_t n, MatrixType& I, double value = 1.0)
	{
		SetZero(n, n, I);
		for (size_t j = 0; j < n; ++j)
			I(j, j) = value;
	}

	/**
	* Copies the block [begin:end[ from A to B. B size must be end-begin. 
	* Note: no size check is made!
	* @param A the first vector
	* @param begin the first index
	* @param end the last+1 index
	* @param B the second vector
	*/
	inline static void GetBlock(const VectorType& A, size_t begin, size_t end, VectorType& B)
	{
		size_t n = end - begin;
		for (size_t i = 0; i < n; ++i)
			B(i) = A(i + begin);
	}

	/**
	* Copies the block [begin:end[ from A to B. B size must be end-begin.
	* Note: no size check is made!
	* @param A the first matrix
	* @param begin the first index
	* @param end the last+1 index
	* @param B the second matrix
	*/
	inline static void GetBlock(const MatrixType& A, size_t begin, size_t end, MatrixType& B)
	{
		size_t n = end - begin;
		for (size_t i = 0; i < n; ++i)
			for (size_t j = 0; j < n; ++j)
				B(i, j) = A(i + begin, j + begin);
	}

	/**
	* Copies the block [begin:end[ from B to A. B size must be end-begin.
	* Note: no size check is made!
	* @param A the first vector
	* @param begin the first index
	* @param end the last+1 index
	* @param B the second vector
	*/
	inline static void SetBlock(MatrixType& A, size_t begin, size_t end, const MatrixType& B)
	{
		size_t n = end - begin;
		for (size_t i = 0; i < n; ++i)
			for (size_t j = 0; j < n; ++j)
				A(i + begin, j + begin) = B(i, j);
	}

	/**
	* computes the outer product A x B, in C.
	* Note: no size check is made!
	* @param A the first vector
	* @param B the second vector
	* @param C the output matrix
	*/
	inline static void OuterProd(const VectorType& A, const VectorType& B, MatrixType& C)
	{
		for (size_t i = 0; i < A.Size(); ++i)
			for (size_t j = 0; j < B.Size(); ++j)
				C(i, j) = A(i) * B(j);
	}

public:

	/**
	* Computes the Translational Projector Matrix.
	* The output is a square matrix of size num_nodes*6.
	* Note that 6 Degrees Of Freedom are assumed for each node.
	* @param num_nodes the number of nodes
	* @param P the Translational Projector Matrix
	*/
	inline static void Compute_Pt(size_t num_nodes, MatrixType& P)
	{
		double a = double(num_nodes - 1) / double(num_nodes);
		double b = -1.0 / double(num_nodes);

		size_t num_dofs = num_nodes * 6;

		SetIdentity(num_dofs, P);

		for (size_t i = 0; i < num_nodes; i++)
		{
			size_t j = i * 6;

			// diagonal block
			P(j, j) = a;
			P(j + 1, j + 1) = a;
			P(j + 2, j + 2) = a;

			// out-of-diagonal block
			for (size_t k = i + 1; k < num_nodes; k++)
			{
				size_t w = k * 6;

				P(j, w) = b;
				P(j + 1, w + 1) = b;
				P(j + 2, w + 2) = b;

				P(w, j) = b;
				P(w + 1, j + 1) = b;
				P(w + 2, j + 2) = b;
			}
		}
	}

	/**
	* Computes the Spin Lever Matrix.
	* The output is a rectangular matrix of 3 columns and nodes.size()*6 rows.
	* Note that 6 Degrees Of Freedom are assumed for each node.
	* @param nodes the input nodes
	* @param S the Spin Lever Matrix
	*/
	inline static void Compute_S(const NodeContainerType& nodes, MatrixType& S)
	{
		size_t num_nodes = nodes.size();
		size_t num_dofs = num_nodes * 6;

		SetZero(num_dofs, 3, S);

		for (size_t i = 0; i < num_nodes; i++)
		{
			size_t j = i * 6;

			Spin_AtRow(nodes[i], S, -1.0, 0, j);

			S(j + 3, 0) = 1.0;
			S(j + 4, 1) = 1.0;
			S(j + 5, 2) = 1.0;
		}
	}

	/**
	* Computes the Axial Vector Jacobian.
	* The output is a square matrix of size displacements.size() (which is num_nodes * 6).
	* Note that 6 Degrees Of Freedom are assumed for each node.
	* @param displacements the vector of nodal displacements and rotations in the local corotational coordinate system. (assumed size = num_nodes*6)
	* @return the H matrix
	*/
	inline static void Compute_H(const VectorType& displacements, MatrixType& H)
	{
		size_t num_dofs = displacements.Size();
		size_t num_nodes = num_dofs / 6;

		SetIdentity(num_dofs, H);

		static MatrixType Omega(3, 3);
		static MatrixType Omega2(3, 3);
		static MatrixType Hi(3, 3);
		static VectorType rv(3);

		for (size_t i = 0; i < num_nodes; i++)
		{
			size_t index = i * 6;

			GetBlock(displacements, index + 3, index + 6, rv);

			double angle = rv.Norm();

			if (angle >= 2.0 * M_PI)
				angle = std::fmod(angle, 2.0 * M_PI);
			
                        // compute eta from Felippa et al 2005, Equation (101)
			double eta;
			if (angle < 0.05) {
				double angle2 = angle * angle;
				double angle4 = angle2 * angle2;
				double angle6 = angle4 * angle2;
				eta = 1.0 / 12.0 + 1.0 / 720.0 * angle2 + 1.0 / 30240.0 * angle4 + 1.0 / 1209600.0 * angle6;
			}
			else {
				eta = (1.0 - 0.5 * angle * std::tan(0.5 * M_PI - 0.5 * angle)) / (angle * angle);
			}

			Spin(rv, Omega);
			Omega2.addMatrixProduct(0.0, Omega, Omega, 1.0);

			// Hi = I - 0.5*Omega + eta*Omega*Omega
			SetIdentity(3, Hi);
			Hi.addMatrix(1.0, Omega, -0.5);
			Hi.addMatrix(1.0, Omega2, eta);

			SetBlock(H, index + 3, index + 6, Hi);
		}
	}

	/**
	* Computes the Spin derivative of (Axial Vector Jacobian)^T contracted with the nodal moment vector.
	* The output is a square matrix of size displacements.size() (which is num_nodes * 6).
	* Note that 6 Degrees Of Freedom are assumed for each node.
	* @param displacements the vector of nodal displacements and rotations in the local corotational coordinate system. (assumed size = num_nodes*6)
	* @param forces the vector of nodal forces and moments in the local corotational coordinate system. (assumed size = num_nodes*6)
	* @param H the Axial Vector Jacobian Matrix computed with a previous call to EICR::Compute_H(displacements)
	* @return the L matrix
	*/
	inline static void Compute_L(const VectorType& displacements, const VectorType& forces, const MatrixType& H, MatrixType& L)
	{
		size_t num_dofs = displacements.Size();
		size_t num_nodes = num_dofs / 6;

		SetZero(num_dofs, num_dofs, L);

		static VectorType rotationVector(3);
		static VectorType momentVector(3);
		static MatrixType Omega(3, 3);
		static MatrixType Omega2(3, 3);
		static MatrixType Li(3, 3);
		static MatrixType LiTemp1(3, 3);
		static MatrixType MxR(3, 3);
		static MatrixType RxM(3, 3);
		static MatrixType Hi(3, 3);

		for (size_t i = 0; i < num_nodes; i++)
		{
			size_t index = i * 6;

			GetBlock(displacements, index + 3, index + 6, rotationVector);
			GetBlock(forces, index + 3, index + 6, momentVector);

			double angle = rotationVector.Norm();

			if (angle >= 2.0 * M_PI)
				angle = std::fmod(angle, 2.0 * M_PI);

			double angle2 = angle * angle;
			double angle4 = angle2 * angle2;
			double angle6 = angle4 * angle2;

			double eta;
			double mu;
			if (angle < 0.05) {
				// Equation (101) from Felippa et al 2005.
				eta = 1.0 / 12.0 + angle2 / 720.0 + angle4 / 30240.0 + angle6 / 1209600.0;
				// Equation (103) from Felippa et al 2005.
				mu = 1.0 / 360.0 + angle2 / 7560.0 + angle4 / 201600.0 + angle6 / 5987520.0;
			}
			else {
				eta = (1.0 - 0.5 * angle * std::tan(0.5 * M_PI - 0.5 * angle)) / (angle * angle);
				double sin_h_angle = std::sin(0.5 * angle);
				mu = (angle2 + 4.0 * std::cos(angle) + angle * std::sin(angle) - 4.0) / (4.0 * angle4 * sin_h_angle * sin_h_angle);
			}

			Spin(rotationVector, Omega);
			Omega2.addMatrixProduct(0.0, Omega, Omega, 1.0);

			OuterProd(momentVector, rotationVector, MxR);
			OuterProd(rotationVector, momentVector, RxM);
			
			SetIdentity(3, Li, (rotationVector ^ momentVector));
			Li.addMatrix(1.0, RxM, 1.0);
			Li.addMatrix(1.0, MxR, -1.0);

			LiTemp1.addMatrixProduct(0.0, Omega2, MxR, mu);
			Spin(momentVector, MxR, 0.5);
			LiTemp1.addMatrix(1.0, MxR, -1.0);

			LiTemp1.addMatrix(1.0, Li, eta);

			GetBlock(H, index + 3, index + 6, Hi);
			Li.addMatrixProduct(0.0, LiTemp1, Hi, 1.0);

			SetBlock(L, index + 3, index + 6, Li);
		}
	}


};

#endif // !ASDEICR_h
