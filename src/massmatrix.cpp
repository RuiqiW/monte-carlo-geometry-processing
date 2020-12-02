#include "massmatrix.h"
#include <igl/doublearea.h>
#include <iostream>
using namespace std;
using namespace Eigen;

void massmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::DiagonalMatrix<double,Eigen::Dynamic> & M)
{
	// Assume M is sized correctly, easier to do outside since we don't have #V here



	VectorXd v = VectorXd::Zero(M.rows());

	for (int i = 0; i < F.rows(); i++) {
		double a = l(i, 0);
		double b = l(i, 1);
		double c = l(i, 2);
		double s = (a + b + c) / 2;
		double area = sqrt(s * (s - a) * (s - b) * (s - c));
		
		v(F(i, 0)) += area / 3;
		v(F(i, 1)) += area / 3;
		v(F(i, 2)) += area / 3;
	}	



	M.diagonal() = v;



	//M = v.array().matrix().asDiagonal();	



}

