#include <Eigen/Core>
#include <iostream>
#include "2D_example_runner.h"
#include <fstream>
#include "walk_on_spheres_2D.h"

using namespace Eigen;
using namespace std;

double boundary_2D(Vector2d boundary_point) {
	return ((int)floor(6 * boundary_point(0)) + (int)floor(6 * boundary_point(1))) % 2;
}

int example_for_2D() {
	MatrixXd V(8, 3);

	V.row(0) = Vector3d(0.5, 0.1, 0);
	V.row(1) = Vector3d(0.9, 0.5, 0);
	V.row(2) = Vector3d(0.5, 0.9, 0);
	V.row(3) = Vector3d(0.1, 0.5, 0);
	V.row(4) = Vector3d(0.5, 0.333333, 0);
	V.row(5) = Vector3d(0.5, 0.66666, 0);
	V.row(6) = Vector3d(0.33333, 0.5, 0);
	V.row(7) = Vector3d(0.666666, 0.5, 0);

	MatrixXi F(5, 3);

	F.row(0) = Vector3i(0, 1, 1);
	F.row(1) = Vector3i(2, 3, 3);
	F.row(2) = Vector3i(3, 0, 0);
	F.row(3) = Vector3i(4, 5, 5);
	F.row(4) = Vector3i(7, 6, 6);


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
		walk_on_spheres_2D(V, F, boundary_2D, P, U);
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