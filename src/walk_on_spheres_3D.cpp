#include "../include/walk_on_spheres_3D.h"
#include <igl/AABB.h>
#include <igl/point_mesh_squared_distance.h>
#include <random>
#include <iostream>

using namespace std;
void walk_on_spheres_3D(
	const Eigen::MatrixXd& V,
	const Eigen::MatrixXi& F,
	double (*B)(Eigen::Vector3d),
	const Eigen::MatrixXd& P,
	Eigen::VectorXd& U)
{
	// TODO: take stoping tolerance as parameter
	double eps = 0.01;

	igl::AABB<Eigen::MatrixXd, 3> tree;
	tree.init(V, F);

	Eigen::MatrixXd Q;
	Q.resizeLike(P);

	Eigen::VectorXd sqrD;
	Eigen::VectorXi I;
	Eigen::MatrixXd C;

	int iter = 0;

	for (int i = 0; i < P.rows(); i++) {
		Q.row(i) = P.row(i);
	}
	while (iter < 5) {
		iter++;
		tree.squared_distance(V, F, Q, sqrD, I, C);

		for (int i = 0; i < Q.rows(); i++) {
			// sample a point on B(x)
			double radius = sqrt(sqrD(i));
			Eigen::RowVector3d center = Q.row(i);
			std::random_device rd;
			std::mt19937 generator(rd());
			std::uniform_real_distribution<double> uniform01(0.0, 1.0);
			double theta = 2 * 3.14 * uniform01(generator);
			double phi = 3.14 * uniform01(generator);
			double x = sin(phi) * cos(theta);
			double y = sin(phi) * sin(theta);
			double z = cos(phi);
			Eigen::RowVector3d sample(x, y, z);
			sample = sample * radius + center;
			Q.row(i) = sample;
		} // end looping over points
	} // end while

	// get closest face
	tree.squared_distance(V, F, Q, sqrD, I, C);

	U.resize(P.rows());

	for (int i = 0; i < P.rows(); i++) {
		U(i) = C.row(i).norm();
	}


}