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
// Implementation of a Quaternion for compact and robust representation
// of spatial rotations
//

#ifndef ASDMath_h
#define ASDMath_h

#include <stdint.h>
#include <math.h>
#include <cmath>
#include <limits>
#include <Vector.h>
#include <Matrix.h>

#ifndef M_PI
#define M_PI 3.1415926535897932384626433832795
#endif // M_PI

/** \brief ASDQuaternion
* A simple fixed size 3D vector
*/
template<class T>
class ASDVector3
{
public:

	/**
	Creates a Zero ASDVector3.
	*/
	ASDVector3()
	{
		mData[0] = mData[1] = mData[2] = 0.0;
	}

	/**
	Creates a ASDVector3 from its coefficients.
	@param x x coefficient
	@param y y coefficient
	@param z z coefficient
	*/
	ASDVector3(T x, T y, T z)
	{
		mData[0] = x;
		mData[1] = y;
		mData[2] = z;
	}

	/**
	Creates a ASDVector3 from a Vector starting from index.
	Note: no check is made on the input data size!
	@param v the vector
	@param index the index (optional, default = 0)
	*/
	ASDVector3(const Vector& v, size_t index = 0)
	{
		mData[0] = v(0 + index);
		mData[1] = v(1 + index);
		mData[2] = v(2 + index);
	}

	/**
	Creates a ASDVector3 from another ASDVector3.
	@param other the other ASDVector3
	*/
	ASDVector3(const ASDVector3& other)
	{
		mData[0] = other.mData[0];
		mData[1] = other.mData[1];
		mData[2] = other.mData[2];
	}

public:

	/**
	Copies a ASDVector3.
	@param other the other ASDVector3
	*/
	ASDVector3& operator= (const ASDVector3& other)
	{
		if (this != &other) {
			mData[0] = other.mData[0];
			mData[1] = other.mData[1];
			mData[2] = other.mData[2];
		}
		return *this;
	}

public:

	/**
	Returns the X coefficient of this vector.
	@return the X coefficient of this vector.
	*/
	inline const T x()const { return mData[0]; }

	/**
	Returns the Y coefficient of this vector.
	@return the Y coefficient of this vector.
	*/
	inline const T y()const { return mData[1]; }

	/**
	Returns the Z coefficient of this vector.
	@return the Z coefficient of this vector.
	*/
	inline const T z()const { return mData[2]; }

	/**
	Returns the i-th coefficient of this vector.
	@return the i-th coefficient of this vector.
	*/
	inline T operator()(size_t i) const { return mData[i]; }

	/**
	Returns the i-th coefficient of this vector.
	@return the i-th coefficient of this vector.
	*/
	inline T& operator()(size_t i) { return mData[i]; }

	/**
	Returns the i-th coefficient of this vector.
	@return the i-th coefficient of this vector.
	*/
	inline T operator[](size_t i) const { return mData[i]; }

	/**
	Returns the i-th coefficient of this vector.
	@return the i-th coefficient of this vector.
	*/
	inline T& operator[](size_t i) { return mData[i]; }

public:

	/**
	Returns the squared norm this vector.
	@return the squared norm of this vector.
	*/
	inline T squaredNorm() const
	{
		return 
			mData[0] * mData[0] +
			mData[1] * mData[1] +
			mData[2] * mData[2];
	}

	/**
	Returns the norm this vector.
	@return the norm of this vector.
	*/
	inline T norm() const
	{
		return std::sqrt(squaredNorm());
	}

	/**
	makes this vector a unit vector.
	@return the norm of this vector.
	*/
	inline T normalize()
	{
		T n = norm();
		if (n > 0.0) {
			mData[0] /= n;
			mData[1] /= n;
			mData[2] /= n;
		}
		return n;
	}

	/**
	Returns the dot product this vector with another vector.
	@param b the other vector
	@return the dot product.
	*/
	inline T dot(const ASDVector3& b) const
	{
		return
			mData[0] * b.mData[0] +
			mData[1] * b.mData[1] +
			mData[2] * b.mData[2];
	}

	/**
	Returns the cross product this vector with another vector.
	@param b the other vector
	@return the cross product.
	*/
	inline ASDVector3 cross(const ASDVector3& b) const
	{
		ASDVector3 c;
		c.mData[0] = mData[1] * b.mData[2] - mData[2] * b.mData[1];
		c.mData[1] = mData[2] * b.mData[0] - mData[0] * b.mData[2];
		c.mData[2] = mData[0] * b.mData[1] - mData[1] * b.mData[0];
		return c;
	}

public:

	inline void operator += (const ASDVector3& b)
	{
		mData[0] += b.mData[0];
		mData[1] += b.mData[1];
		mData[2] += b.mData[2];
	}

	inline ASDVector3 operator + (const ASDVector3& b) const
	{
		ASDVector3 a(*this);
		a += b;
		return a;
	}

	inline void operator -= (const ASDVector3& b)
	{
		mData[0] -= b.mData[0];
		mData[1] -= b.mData[1];
		mData[2] -= b.mData[2];
	}

	inline ASDVector3 operator - (const ASDVector3& b) const
	{
		ASDVector3 a(*this);
		a -= b;
		return a;
	}

	inline void operator *= (T b)
	{
		mData[0] *= b;
		mData[1] *= b;
		mData[2] *= b;
	}

	inline ASDVector3 operator * (T b) const
	{
		ASDVector3 a(*this);
		a *= b;
		return a;
	}

	inline void operator /= (T b)
	{
		mData[0] /= b;
		mData[1] /= b;
		mData[2] /= b;
	}

	inline ASDVector3 operator / (T b) const
	{
		ASDVector3 a(*this);
		a /= b;
		return a;
	}

private:
	T mData[3];
};

template<class T>
inline ASDVector3<T> operator * (T a, const ASDVector3<T>& b)
{
	return b * a;
}

/**
Prints this vector to a input stream
@param s the output stream
@param v the vector
@return the stream
*/
template<class TStream, class T>
inline TStream& operator << (TStream& s, const ASDVector3<T>& v)
{
	return (s << "(" << v.x() << ", " << v.y() << ", " << v.z() << ")");
}

/** \brief ASDQuaternion
* A simple class that implements the main features of quaternion algebra
*/
template<class T>
class ASDQuaternion
{

public:

	/**
	Creates a Zero ASDQuaternion.
	*/
	ASDQuaternion()
		: mX(0.0)
		, mY(0.0)
		, mZ(0.0)
		, mW(0.0)
	{
	}

	/**
	Creates a ASDQuaternion from its coefficients.
	@param w w coefficient
	@param x x coefficient
	@param y y coefficient
	@param z z coefficient
	*/
	ASDQuaternion(T w, T x, T y, T z)
		: mX(x)
		, mY(y)
		, mZ(z)
		, mW(w)
	{
	}

	/**
	Creates a ASDQuaternion from another ASDQuaternion.
	@param other the other ASDQuaternion
	*/
	ASDQuaternion(const ASDQuaternion& other)
		: mX(other.mX)
		, mY(other.mY)
		, mZ(other.mZ)
		, mW(other.mW)
	{
	}

public:

	/**
	Copies a ASDQuaternion.
	@param other the other ASDQuaternion
	*/
	ASDQuaternion& operator= (const ASDQuaternion& other)
	{
		if (this != &other) {
			mX = other.mX;
			mY = other.mY;
			mZ = other.mZ;
			mW = other.mW;
		}
		return *this;
	}

public:

	/**
	Returns the X coefficient of this quaternion.
	@return the X coefficient of this quaternion.
	*/
	inline const T x()const { return mX; }

	/**
	Returns the Y coefficient of this quaternion.
	@return the Y coefficient of this quaternion.
	*/
	inline const T y()const { return mY; }

	/**
	Returns the Z coefficient of this quaternion.
	@return the Z coefficient of this quaternion.
	*/
	inline const T z()const { return mZ; }

	/**
	Returns the W coefficient of this quaternion.
	@return the W coefficient of this quaternion.
	*/
	inline const T w()const { return mW; }

public:

	/**
	Returns the squared norm of this quaternion.
	x*x + y*y + z*z + w*w
	@return the squared norm of this quaternion.
	*/
	inline const T squaredNorm()const
	{
		return mX * mX + mY * mY + mZ * mZ + mW * mW;
	}

	/**
	Returns the norm of this quaternion.
	sqrt(x*x + y*y + z*z + w*w)
	@return the norm of this quaternion.
	*/
	inline const T norm()const
	{
		return std::sqrt(squaredNorm());
	}

	/**
	Makes this ASDQuaternion a Unit ASDQuaternion.
	If this ASDQuaternion is already normalized this is a no-op
	*/
	inline void normalize()
	{
		T n = squaredNorm();
		if (n > 0.0 && n != 1.0) {
			n = std::sqrt(n);
			mX /= n;
			mY /= n;
			mZ /= n;
			mW /= n;
		}
	}

	/**
	Returns the Conjugate of this ASDQuaternion, which represents the opposite rotation
	@return the Conjugate of this ASDQuaternion
	*/
	inline ASDQuaternion conjugate()const
	{
		return ASDQuaternion(mW, -mX, -mY, -mZ);
	}

	/**
	Constructs a Rotation Matrix from this ASDQuaternion.
	The rotation matrix type is the template argument, no check is made on the type of
	this matrix.
	These assumptions are made:
	The matrix should provide an indexed access like m(i, j) where i and j are indices
	from 0 to 2.
	This means that the input matrix is a C-Style 3x3 Matrix.
	All the 9 coefficients are properly set so there's no need to set the matrix to Zero
	before calling this function.
	@param R the output rotation matrix
	*/
	template<class TMatrix3x3>
	inline void toRotationMatrix(TMatrix3x3& R)const
	{
		R(0, 0) = 2.0 * (mW * mW + mX * mX - 0.5);
		R(0, 1) = 2.0 * (mX * mY - mW * mZ);
		R(0, 2) = 2.0 * (mX * mZ + mW * mY);

		R(1, 0) = 2.0 * (mY * mX + mW * mZ);
		R(1, 1) = 2.0 * (mW * mW + mY * mY - 0.5);
		R(1, 2) = 2.0 * (mY * mZ - mX * mW);

		R(2, 0) = 2.0 * (mZ * mX - mW * mY);
		R(2, 1) = 2.0 * (mZ * mY + mW * mX);
		R(2, 2) = 2.0 * (mW * mW + mZ * mZ - 0.5);
	}

	/**
	Extracts the Rotation Vector from this ASDQuaternion
	@param rx the output x component if the rotation vector
	@param ry the output y component if the rotation vector
	@param rz the output z component if the rotation vector
	*/
	inline void toRotationVector(T& rx, T& ry, T& rz)const
	{
		T xx, yy, zz, ww;

		if (mW < 0.0) {
			xx = -mX;
			yy = -mY;
			zz = -mZ;
			ww = -mW;
		}
		else {
			xx = mX;
			yy = mY;
			zz = mZ;
			ww = mW;
		}

		T vNorm = xx * xx + yy * yy + zz * zz;
		if (vNorm == 0.0) {
			rx = 0.0;
			ry = 0.0;
			rz = 0.0;
			return;
		}

		if (vNorm != 1.0)
			vNorm = std::sqrt(vNorm);

		T mult = (vNorm < ww) ? (2.0 / vNorm * std::asin(vNorm)) : (2.0 / vNorm * std::acos(ww));

		rx = xx * mult;
		ry = yy * mult;
		rz = zz * mult;
	}

	/**
	Extracts the Rotation Vector from this ASDQuaternion
	The vector type is the template parameter. No check is made on this type.
	The following assumptions are made:
	The vector type should provide indexing like vector(i) where i goes from 0 to 2.
	(i.e. a C-Style vector of size 3)
	@param v the output rotation vector
	*/
	template<class TVector3>
	inline void toRotationVector(TVector3& v)const
	{
		toRotationVector(v(0), v(1), v(2));
	}

	/**
	Rotates a vector using this quaternion.
	Note: this is faster than constructing the rotation matrix and perform the matrix
	multiplication for a single vector.
	The vector type is the template parameter. No check is made on this type.
	The following assumptions are made:
	The vector type should provide indexing like vector(i) where i goes from 0 to 2.
	(i.e. a C-Style vector of size 3)
	@param a the input source vector
	@param b the output rotated vector
	*/
	template<class TVector3_A, class TVector3_B>
	inline void rotateVector(const TVector3_A& a, TVector3_B& b)const
	{
		// b = 2.0 * cross( this->VectorialPart, a )
		b(0) = 2.0 * (mY * a(2) - mZ * a(1));
		b(1) = 2.0 * (mZ * a(0) - mX * a(2));
		b(2) = 2.0 * (mX * a(1) - mY * a(0));

		// c = cross( this->VectorialPart, b )
		T c0 = mY * b(2) - mZ * b(1);
		T c1 = mZ * b(0) - mX * b(2);
		T c2 = mX * b(1) - mY * b(0);

		// set results
		b(0) = a(0) + b(0) * mW + c0;
		b(1) = a(1) + b(1) * mW + c1;
		b(2) = a(2) + b(2) * mW + c2;
	}

	/**
	Rotates a vector using this quaternion.
	Note: this is faster than constructing the rotation matrix and perform the matrix
	multiplication for a single vector.
	The vector type is the template parameter. No check is made on this type.
	The following assumptions are made:
	The vector type should provide indexing like vector(i) where i goes from 0 to 2.
	(i.e. a C-Style vector of size 3)
	@param a the input source vector - rotated on exit
	*/
	template<class TVector3>
	inline void rotateVector(TVector3& a)const
	{
		// b = 2.0 * cross( this->VectorialPart, a )
		T b0 = 2.0 * (mY * a(2) - mZ * a(1));
		T b1 = 2.0 * (mZ * a(0) - mX * a(2));
		T b2 = 2.0 * (mX * a(1) - mY * a(0));

		// c = cross( this->VectorialPart, b )
		T c0 = mY * b2 - mZ * b1;
		T c1 = mZ * b0 - mX * b2;
		T c2 = mX * b1 - mY * b0;

		// set results
		a(0) += b0 * mW + c0;
		a(1) += b1 * mW + c1;
		a(2) += b2 * mW + c2;
	}

public:

	/**
	Returns the Identity ASDQuaternion (i.e. a ASDQuaternion that represents a Zero rotation)
	@return the Identity ASDQuaternion
	*/
	static inline ASDQuaternion Identity()
	{
		return ASDQuaternion(1.0, 0.0, 0.0, 0.0);
	}

	/**
	Returns a ASDQuaternion that represents a rotation of an angle 'radians' around the axis (x, y, z)
	@param x the x component of the rotation axis
	@param y the y component of the rotation axis
	@param z the z component of the rotation axis
	@param radians the rotation angle in radians
	@return a ASDQuaternion that represents a rotation of an angle 'radians' around the axis (x, y, z)
	*/
	static inline ASDQuaternion FromAxisAngle(T x, T y, T z, T radians)
	{
		T sqnorm = x * x + y * y + z * z;
		if (sqnorm == 0.0)
			return ASDQuaternion::Identity();

		if (sqnorm > 0.0 && sqnorm != 1.0) {
			T norm = std::sqrt(sqnorm);
			x /= norm;
			y /= norm;
			z /= norm;
		}

		T halfAngle = radians * 0.5;

		T s = std::sin(halfAngle);
		T q0 = std::cos(halfAngle);

		ASDQuaternion result(q0, s * x, s * y, s * z);
		result.normalize();

		return result;
	}

	/**
	Returns a ASDQuaternion from a rotation vector
	@param rx the x component of the source rotation vector
	@param ry the y component of the source rotation vector
	@param rz the z component of the source rotation vector
	@return a ASDQuaternion from a rotation vector
	*/
	static inline ASDQuaternion FromRotationVector(T rx, T ry, T rz)
	{
		T rModulus = rx * rx + ry * ry + rz * rz;
		if (rModulus == 0.0)
			return ASDQuaternion::Identity();

		if (rModulus != 1.0) {
			rModulus = std::sqrt(rModulus);
			rx /= rModulus;
			ry /= rModulus;
			rz /= rModulus;
		}

		T halfAngle = rModulus * 0.5;

		T q0 = std::cos(halfAngle);
		T s = std::sin(halfAngle);

		ASDQuaternion result(q0, rx * s, ry * s, rz * s);
		result.normalize();

		return result;
	}

	/**
	Returns a ASDQuaternion from a rotation vector.
	The vector type is the template parameter. No check is made on this type.
	The following assumptions are made:
	The vector type should provide indexing like vector(i) where i goes from 0 to 2.
	(i.e. a C-Style vector of size 3)
	@param v the source rotation vector
	@return a ASDQuaternion from a rotation vector
	*/
	template<class TVector3>
	static inline ASDQuaternion FromRotationVector(const TVector3& v)
	{
		return ASDQuaternion::FromRotationVector(v(0), v(1), v(2));
	}

	/**
	Returns a ASDQuaternion from a Rotation Matrix.
	The rotation matrix type is the template argument, no check is made on the type of
	this matrix.
	These assumptions are made:
	The matrix should provide an indexed access like m(i, j) where i and j are indices
	from 0 to 2.
	This means that the input matrix is a C-Style 3x3 Matrix.
	@param m the source rotation matrix
	@return a ASDQuaternion from a Rotation Matrix
	*/
	template<class TMatrix3x3>
	static inline ASDQuaternion FromRotationMatrix(const TMatrix3x3& m)
	{
		T xx = m(0, 0);
		T yy = m(1, 1);
		T zz = m(2, 2);
		T tr = xx + yy + zz;
		ASDQuaternion Q;
		if ((tr > xx) && (tr > yy) && (tr > zz))
		{
			T S = std::sqrt(tr + 1.0) * 2.0;
			Q = ASDQuaternion(
				0.25 * S,
				(m(2, 1) - m(1, 2)) / S,
				(m(0, 2) - m(2, 0)) / S,
				(m(1, 0) - m(0, 1)) / S
			);
		}
		else if ((xx > yy) && (xx > zz))
		{
			T S = std::sqrt(1.0 + xx - yy - zz) * 2.0;
			Q = ASDQuaternion(
				(m(2, 1) - m(1, 2)) / S,
				0.25 * S,
				(m(0, 1) + m(1, 0)) / S,
				(m(0, 2) + m(2, 0)) / S
			);
		}
		else if (yy > zz)
		{
			T S = std::sqrt(1.0 + yy - xx - zz) * 2.0;
			Q = ASDQuaternion(
				(m(0, 2) - m(2, 0)) / S,
				(m(0, 1) + m(1, 0)) / S,
				0.25 * S,
				(m(1, 2) + m(2, 1)) / S
			);
		}
		else
		{
			T S = std::sqrt(1.0 + zz - xx - yy) * 2.0;
			Q = ASDQuaternion(
				(m(1, 0) - m(0, 1)) / S,
				(m(0, 2) + m(2, 0)) / S,
				(m(1, 2) + m(2, 1)) / S,
				0.25 * S
			);
		}


		Q.normalize();
		return Q;
	}

private:
	T mX;
	T mY;
	T mZ;
	T mW;

};

/**
Performs a ASDQuaternion-ASDQuaternion product and returns the concatenation
of the two input quaternions (i.e. the compound rotation)
@param a the first quaternion
@param b the second quaternion
@return the compound quaternion
*/
template<class T>
inline ASDQuaternion<T> operator* (const ASDQuaternion<T>& a, const ASDQuaternion<T>& b)
{
	return ASDQuaternion<T>(
		a.w() * b.w() - a.x() * b.x() - a.y() * b.y() - a.z() * b.z(),
		a.w() * b.x() + a.x() * b.w() + a.y() * b.z() - a.z() * b.y(),
		a.w() * b.y() + a.y() * b.w() + a.z() * b.x() - a.x() * b.z(),
		a.w() * b.z() + a.z() * b.w() + a.x() * b.y() - a.y() * b.x()
		);
}

/**
Prints this quaternion to a input stream
@param s the output stream
@param q the quaternion
@return the stream
*/
template<class TStream, class T>
inline TStream& operator << (TStream& s, const ASDQuaternion<T>& q)
{
	return (s << "[(" << q.x() << ", " << q.y() << ", " << q.z() << "), " << q.w() << "]");
}

#endif // !ASDMath_h


