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
// Implementation of a local coordinate system for 4-node shells
//

#ifndef ASDShellQ4LocalCoordinateSystem_h
#define ASDShellQ4LocalCoordinateSystem_h

#include <ASDMath.h>
#include <array>
#include <vector>

// if this macro is defined, the local system is alligned according to the
// original paper of C.Felippa, i.e. the local X axis is alligned to side 1-2
// otherwise we use the same convention as in all the other shells in opensees
#define ALLIGN_AS_FELIPPA

/** \brief ASDShellQ4LocalCoordinateSystem
*
* This class represent the local coordinate system of any element whose geometry
* is a 4-node Quadrilateral in 3D space
*/
class ASDShellQ4LocalCoordinateSystem
{

public:

	typedef ASDVector3<double> Vector3Type;

	typedef ASDQuaternion<double> QuaternionType;

	typedef std::vector<Vector3Type> Vector3ContainerType;

	typedef Matrix MatrixType;

public:

	ASDShellQ4LocalCoordinateSystem(
		const Vector3Type& P1global,
		const Vector3Type& P2global,
		const Vector3Type& P3global,
		const Vector3Type& P4global,
		double alpha = 0.0)
		: m_P(4)
		, m_orientation(3, 3)
	{
		// Form the basis vectors alligning the local X direction
		// with the 1-2 side

		// compute the central point
		m_center = P1global;
		m_center += P2global;
		m_center += P3global;
		m_center += P4global;
		m_center *= 0.25;

#ifdef ALLIGN_AS_FELIPPA

		// compute the diagonal vectors
		Vector3Type d13(P3global - P1global);
		Vector3Type d24(P4global - P2global);

		// compute the Normal vector at the element center
		// as the cross product of the 2 diagonals.
		// While normalizing the normal vector save its norm to compute the area.
		// Note: the norm should be divided by 2 because it's computed from the
		// cross product of the diagonals, which gives twice the area!
		Vector3Type e3 = d13.cross(d24);
		m_area = e3.normalize();
		m_area /= 2.0;

		// compute the local X direction parallel to the projection of the side 1-2 onto
		// the local xy plane.
		Vector3Type e1(P2global - P1global);
		double e1_dot_e3 = e1.dot(e3);
		e1 -= e1_dot_e3 * e3;
		if (std::abs(alpha) > DBL_EPSILON)
			QuaternionType::FromAxisAngle(e3(0), e3(1), e3(2), alpha).rotateVector(e1);
		e1.normalize();

		// finally compute the local Y direction to be orthogonal to both X and Z local directions
		Vector3Type e2 = e3.cross(e1);
		e2.normalize();

#else

		// compute the right and left mid nodes
		Vector3Type p23(P3global + P2global);
		p23 *= 0.5;
		Vector3Type p14(P4global + P1global);
		p14 *= 0.5;

		// compute the local x axis
		Vector3Type e1(p23 - p14);

		// compute the bottom and top mid nodes
		Vector3Type p12(P2global + P1global);
		p12 *= 0.5;
		Vector3Type p43(P3global + P4global);
		p43 *= 0.5;

		// compute the local y axis
		Vector3Type e2(p43 - p12);

		// compute the normal vector at the element center and get area
		// from its norm
		Vector3Type e3 = e1.cross(e2);
		m_area = e3.norm();

		// update local y axis to make an orthonormal basis
		e2 = e3.cross(e1);

		// normalize all local vectors
		e1.normalize();
		e2.normalize();
		e3.normalize();

#endif // ALLIGN_AS_FELIPPA

		// form the 3x3 transformation matrix
		for (int i = 0; i < 3; i++)
		{
			m_orientation(0, i) = e1(i);
			m_orientation(1, i) = e2(i);
			m_orientation(2, i) = e3(i);
		}

		// transform global coordinates to the local coordinate system
		for (int i = 0; i < 3; i++)
		{
			m_P[0](i) = m_orientation(i, 0) * (P1global(0) - m_center(0)) + m_orientation(i, 1) * (P1global(1) - m_center(1)) + m_orientation(i, 2) * (P1global(2) - m_center(2));
			m_P[1](i) = m_orientation(i, 0) * (P2global(0) - m_center(0)) + m_orientation(i, 1) * (P2global(1) - m_center(1)) + m_orientation(i, 2) * (P2global(2) - m_center(2));
			m_P[2](i) = m_orientation(i, 0) * (P3global(0) - m_center(0)) + m_orientation(i, 1) * (P3global(1) - m_center(1)) + m_orientation(i, 2) * (P3global(2) - m_center(2));
			m_P[3](i) = m_orientation(i, 0) * (P4global(0) - m_center(0)) + m_orientation(i, 1) * (P4global(1) - m_center(1)) + m_orientation(i, 2) * (P4global(2) - m_center(2));
		}
	}

public:

	inline const Vector3ContainerType& Nodes()const { return m_P; }

	inline const Vector3Type& P1()const { return m_P[0]; }
	inline const Vector3Type& P2()const { return m_P[1]; }
	inline const Vector3Type& P3()const { return m_P[2]; }
	inline const Vector3Type& P4()const { return m_P[3]; }
	inline const Vector3Type& Center()const { return m_center; }

	inline double X1()const { return m_P[0][0]; }
	inline double X2()const { return m_P[1][0]; }
	inline double X3()const { return m_P[2][0]; }
	inline double X4()const { return m_P[3][0]; }

	inline double Y1()const { return m_P[0][1]; }
	inline double Y2()const { return m_P[1][1]; }
	inline double Y3()const { return m_P[2][1]; }
	inline double Y4()const { return m_P[3][1]; }

	inline double Z1()const { return m_P[0][2]; }
	inline double Z2()const { return m_P[1][2]; }
	inline double Z3()const { return m_P[2][2]; }
	inline double Z4()const { return m_P[3][2]; }

	inline double X(size_t i)const { return m_P[i][0]; }
	inline double Y(size_t i)const { return m_P[i][1]; }
	inline double Z(size_t i)const { return m_P[i][2]; }

	inline double Area()const { return m_area; }

	inline const MatrixType& Orientation()const { return m_orientation; }

	inline Vector3Type Vx()const { return Vector3Type(m_orientation(0, 0), m_orientation(0, 1), m_orientation(0, 2)); }
	inline Vector3Type Vy()const { return Vector3Type(m_orientation(1, 0), m_orientation(1, 1), m_orientation(1, 2)); }
	inline Vector3Type Vz()const { return Vector3Type(m_orientation(2, 0), m_orientation(2, 1), m_orientation(2, 2)); }

	inline double WarpageFactor()const { return Z1(); }
	inline bool IsWarped()const { return std::abs(WarpageFactor()) > 0.0; }

	inline void ComputeTotalRotationMatrix(MatrixType& R)const
	{
		size_t mat_size = 24;
		if (R.noRows() != mat_size || R.noCols() != mat_size)
			R.resize(mat_size, mat_size);

		R.Zero();

		for (size_t k = 0; k < 8; k++)
		{
			size_t i = k * 3;
			R(i, i) = m_orientation(0, 0);   R(i, i + 1) = m_orientation(0, 1);   R(i, i + 2) = m_orientation(0, 2);
			R(i + 1, i) = m_orientation(1, 0);   R(i + 1, i + 1) = m_orientation(1, 1);   R(i + 1, i + 2) = m_orientation(1, 2);
			R(i + 2, i) = m_orientation(2, 0);   R(i + 2, i + 1) = m_orientation(2, 1);   R(i + 2, i + 2) = m_orientation(2, 2);
		}
	}

	inline void ComputeTotalWarpageMatrix(MatrixType& W, double wf)const
	{
		size_t mat_size = 24;
		if (W.noRows() != mat_size || W.noCols() != mat_size)
			W.resize(mat_size, mat_size);

		W.Zero();

		W(0, 4) = -wf;
		W(1, 3) = wf;

		W(6, 10) = wf;
		W(7, 9) = -wf;

		W(12, 16) = -wf;
		W(13, 15) = wf;

		W(18, 22) = wf;
		W(19, 21) = -wf;
	}

	inline void ComputeTotalWarpageMatrix(MatrixType& W)const
	{
		ComputeTotalWarpageMatrix(W, WarpageFactor());
	}

private:

	Vector3ContainerType m_P;
	Vector3Type m_center;
	MatrixType m_orientation;
	double m_area;

};

#endif // !ASDShellQ4LocalCoordinateSystem_h
