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
// Utility header file to interface with the Eigen library
//

#ifndef EigenAPI_converters_h
#define EigenAPI_converters_h

#include "EigenAPI.h"

//Create a new OpenSees Matrix and copy the components of an eigen matrix (or expression) to it.
template <typename Derived>
Matrix copyToNewMatrix(const Eigen::DenseBase<Derived> &a)
{
	Matrix ops_matrix(a.rows(), a.cols());

	for (int i = 0; i < a.rows(); ++i)
	{
		for (int j = 0; j < a.cols(); ++j)
		{
			ops_matrix(i, j) = a(i, j);
		}
	}
	return ops_matrix;
}

//Copy the components of an eigen matrix (or expression) to a reference to an OpenSees matrix.
template <typename Derived>
void copyToMatrixReference(const Eigen::DenseBase<Derived> &a, Matrix& ops_matrix)
{
	// Check if the opensees matrix has an appropriate size, otherwise resize.
	// might produce a new call to malloc.
	if (ops_matrix.noRows() != a.rows() || ops_matrix.noCols() != a.cols())
	{
		ops_matrix.resize(a.rows(), a.cols());
	}

	for (int i = 0; i < a.rows(); ++i)
	{
		for (int j = 0; j < a.cols(); ++j)
		{
			ops_matrix(i, j) = a(i, j);
		}
	}
	return ;
}

//Create a new OpenSees Vector and copy the components of an eigen vector (or expression) to it.
template <typename Derived>
Vector copyToNewVector(const Eigen::DenseBase<Derived> &a)
{
	Vector ops_vector(a.rows());

	for (int i = 0; i < a.rows(); ++i)
	{
		ops_vector(i) = a(i);
	}
	return ops_vector;
}

//Copy the components of an eigen vector (or expression) to a reference to an OpenSees vector.
template <typename Derived>
void copyToVectorReference(const Eigen::DenseBase<Derived> &a, Vector& ops_vector)
{
	// Check if the opensees matrix has an appropriate size, otherwise resize.
	// might produce a new call to malloc.
	if (ops_vector.Size() != a.rows() )
	{
		ops_vector.resize(a.rows());
	}

	for (int i = 0; i < a.rows(); ++i)
	{

		ops_vector(i) = a(i);

	}
	return ;
}

#endif // EigenAPI_converters_h