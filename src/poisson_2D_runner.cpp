#include <Eigen/Core>
#include <iostream>
#include "poisson_2D_runner.h"
#include <fstream>
#include "poisson_2D.h"

using namespace Eigen;
using namespace std;

double poisson_boundary_2D(Vector2d boundary_point) {
	return cos(2.0 * 3.14 * boundary_point(0)) * sin(2.0 * 3.14 * boundary_point(0));
	//return ((int)floor(6 * boundary_point(0)) + (int)floor(6 * boundary_point(1))) % 2;
}

double source(RowVector3d point) {
	double M_PI = 3.14;
	return 8.f * (M_PI * M_PI) * cos(2.f * M_PI * point(0)) * sin(2.f * M_PI * point(1));
	//return cos(2.0 * 3.14 * point(0)) * sin(2.0 * 3.14 * point(0));
}


int poisson_2D_runner() {
	//MatrixXd V(8, 3);

	MatrixXd V(4, 3);
	V.row(0) = Vector3d(0, 0, 0);
	V.row(1) = Vector3d(0, 1, 0);
	V.row(2) = Vector3d(1, 0, 0);
	V.row(3) = Vector3d(1, 1, 0);

	MatrixXi F(4, 3);
	F.row(0) = Vector3i(0, 2, 2);
	F.row(1) = Vector3i(2, 3, 3);
	F.row(2) = Vector3i(3, 1, 1);
	F.row(3) = Vector3i(1, 0, 0);


	int n = 25;
	MatrixXd P(n * n, 3);
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			P.row(n * i + j) = 1.0 / n * Vector3d(i, j, 0);
		}
	}
	VectorXd total_U = VectorXd::Zero(n * n);

	int NUM_ITERATIONS = 128;

	for (int k = 0; k < NUM_ITERATIONS; k++) {
		VectorXd U;
		poisson_2D(V, F, poisson_boundary_2D, source, P, U);
		total_U += U;
	}

	total_U /= (float)NUM_ITERATIONS;

	std::cout << total_U << std::endl;
	ofstream out("out.csv");

	for (int i = 0; i < total_U.rows(); i++) {
		out << total_U(i) << endl;
	}

	return 0;
}