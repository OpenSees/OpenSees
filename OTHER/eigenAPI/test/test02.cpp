#include <iostream>
#include "../EigenAPI.h"
#include <StandardStream.h>

StandardStream sserr;
OPS_Stream *opserrPtr = &sserr;


using Eigen::MatrixXd;
using Eigen::Matrix2d;
using Eigen::Matrix3d;
using Eigen::Matrix4d;

template <typename Derived>
void test_copyToNewMatrix(const Eigen::DenseBase<Derived> &m)
{
	std::cout << std::endl;

	std::cout << "Eigen matrix :\n";
	std::cout << m << std::endl;

	std::cout << std::endl;

	Matrix ops_mat = copyToNewMatrix(m);

	*opserrPtr << "OPS matrix:" ;
	*opserrPtr << ops_mat ;

	return ;
}


int main()
{
	test_copyToNewMatrix(Matrix2d::Random());
	test_copyToNewMatrix(Matrix3d::Random());
	test_copyToNewMatrix(Matrix4d::Random());
	test_copyToNewMatrix(MatrixXd::Random(3,8));  
	
	typedef Eigen::Matrix<double, 3, 8> myCustomMatrix;
	test_copyToNewMatrix(myCustomMatrix::Random());
}