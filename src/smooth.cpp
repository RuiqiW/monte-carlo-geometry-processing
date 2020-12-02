#include "smooth.h"
#include "massmatrix.h"
#include "cotmatrix.h"

#include <igl/massmatrix.h>
#include <igl/cotmatrix.h>
#include <igl/cotmatrix_entries.h>
#include <igl/edge_lengths.h>

#include <iostream>

using namespace Eigen;
using namespace std;

void smooth(
	const Eigen::MatrixXd & V,
	const Eigen::MatrixXi & F,
	const Eigen::MatrixXd & G,
	double lambda,
	Eigen::MatrixXd & U)
{
	// Replace with your code
	//U = G;

	// random notes 
	// u_i = (u_i + u_(i+1)) / 2 is not a good averaging scheme as it depends on what i ranges on (if you have a coarse grid, amount of sampling
	// is greater compared to a sparse grid
	// forward euler -> u_i = u_i + lambda * u'_(i+1)(t) is bad because error propagates
	// did you know that u * partial^2 u / partial x^2 is quadratic in u, as 
	// u is linear in u, and partial^2 u / partial x^2 is linear in u, so the product is quadratic

	//MatrixXd l;
	//igl::edge_lengths(V, F, l);

	//DiagonalMatrix<double, Eigen::Dynamic> M_pre(V.rows());
	//massmatrix(l, F, M_pre);

	//VectorXd v = M_pre.diagonal();
	//
	//// Is there another way to use the DiagonalMatrix class, this seems clumsy. Can you subtract MatrixXd from DiagonalMatrix?
	//SparseMatrix<double> M(V.rows(), V.rows());
	//for (int i = 0; i < V.rows(); i++) {
	//	M.insert(i, i) = v(i);
	//}
	//
	//SparseMatrix<double> L(V.rows(), V.rows());

	//cotmatrix(l, F, L);

	//SimplicialLLT<SparseMatrix<double>> solver;
	//solver.compute(M - lambda * L);

	//U = solver.solve(M * G);
	//// Why doesn't this work as well?
	//// MatrixXd I = MatrixXd::Identity(M.rows(), M.cols())
	//// A = solver.compute(I);
	//// U = A * M * G
}
