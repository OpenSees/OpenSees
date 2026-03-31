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
// $Date: 2024/03

// Original implementation: Massimo Petracca (ASDEA)
//
// Implementation of a local coordinate system for 3-node shells
//

#ifndef ASDShellT3LocalCoordinateSystem_h
#define ASDShellT3LocalCoordinateSystem_h

#include <ASDMath.h>
#include <array>
#include <vector>

/** \brief ASDShellT3LocalCoordinateSystem
*
* This class represent the local coordinate system of any element whose geometry
* is a 3-node Triangle in 3D space
*/
class ASDShellT3LocalCoordinateSystem
{

public:

	typedef ASDVector3<double> Vector3Type;

	typedef ASDQuaternion<double> QuaternionType;

	typedef std::vector<Vector3Type> Vector3ContainerType;

	typedef Matrix MatrixType;

public:

	ASDShellT3LocalCoordinateSystem(
		const Vector3Type& P1global,
		const Vector3Type& P2global,
		const Vector3Type& P3global,
		double alpha = 0.0)
		: m_P(3)
		, m_orientation(3, 3)
	{
		// compute the central point
		m_center = P1global;
		m_center += P2global;
		m_center += P3global;
		m_center /= 3.0;
		
		// compute the local X direction
		Vector3Type e1 = P2global - P1global;
		
		// compute the temporary Y direction for the XY plane
		Vector3Type e2 = P3global - P1global;
		
		// compute the local Z direction (normal).
		Vector3Type e3 = e1.cross(e2);
		
		// if user defined local rotation is included rotate the local X direction
		if (std::abs(alpha) > 0.0)
			QuaternionType::FromAxisAngle(e3(0), e3(1), e3(2), alpha).rotateVector(e1);
		
		// finally compute the local Y direction to be orthogonal to both X and Z local directions
		e2 = e3.cross(e1);
		
		// normalize all local direction, save the e3 norm for the area
		e1.normalize();
		e2.normalize();
		m_area = e3.normalize();
		m_area /= 2.0;
		
		// form the 3x3 transformation matrix (the transposed actually...)
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
		}
	}

public:

	inline const Vector3ContainerType& Nodes()const { return m_P; }

	inline const Vector3Type& P1()const { return m_P[0]; }
	inline const Vector3Type& P2()const { return m_P[1]; }
	inline const Vector3Type& P3()const { return m_P[2]; }
	inline const Vector3Type& Center()const { return m_center; }

	inline double X1()const { return m_P[0][0]; }
	inline double X2()const { return m_P[1][0]; }
	inline double X3()const { return m_P[2][0]; }

	inline double Y1()const { return m_P[0][1]; }
	inline double Y2()const { return m_P[1][1]; }
	inline double Y3()const { return m_P[2][1]; }

	inline double Z1()const { return m_P[0][2]; }
	inline double Z2()const { return m_P[1][2]; }
	inline double Z3()const { return m_P[2][2]; }

	inline double X(size_t i)const { return m_P[i][0]; }
	inline double Y(size_t i)const { return m_P[i][1]; }
	inline double Z(size_t i)const { return m_P[i][2]; }

	inline double Area()const { return m_area; }

	inline const MatrixType& Orientation()const { return m_orientation; }

	inline Vector3Type Vx()const { return Vector3Type(m_orientation(0, 0), m_orientation(0, 1), m_orientation(0, 2)); }
	inline Vector3Type Vy()const { return Vector3Type(m_orientation(1, 0), m_orientation(1, 1), m_orientation(1, 2)); }
	inline Vector3Type Vz()const { return Vector3Type(m_orientation(2, 0), m_orientation(2, 1), m_orientation(2, 2)); }

	inline void ComputeTotalRotationMatrix(MatrixType& R)const
	{
		constexpr size_t mat_size = 18;
		if (R.noRows() != mat_size || R.noCols() != mat_size)
			R.resize(mat_size, mat_size);

		R.Zero();

		for (size_t k = 0; k < 6; k++)
		{
			size_t i = k * 3;
			R(i    , i) = m_orientation(0, 0);   R(i    , i + 1) = m_orientation(0, 1);   R(i    , i + 2) = m_orientation(0, 2);
			R(i + 1, i) = m_orientation(1, 0);   R(i + 1, i + 1) = m_orientation(1, 1);   R(i + 1, i + 2) = m_orientation(1, 2);
			R(i + 2, i) = m_orientation(2, 0);   R(i + 2, i + 1) = m_orientation(2, 1);   R(i + 2, i + 2) = m_orientation(2, 2);
		}
	}

private:

	Vector3ContainerType m_P;
	Vector3Type m_center;
	MatrixType m_orientation;
	double m_area;

};

template<class TStream>
inline static TStream& operator << (TStream& s, const ASDShellT3LocalCoordinateSystem& lc)
{
	s << "Points:\n";
	s << "[";
	for (size_t i = 0; i < 3; i++) {
		s << " (" << lc.X(i) << ", " << lc.Y(i) << ") ";
	}
	s << "]";
	s << "Normal: " << lc.Vz() << "\n";
	return s;
}

#endif // !ASDShellT3LocalCoordinateSystem_h
