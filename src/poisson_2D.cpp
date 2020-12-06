#include "../include/walk_on_spheres_2D.h"
#include <igl/AABB.h>
#include <igl/point_mesh_squared_distance.h>
#include <random>
#include <iostream>

using namespace std;
using namespace Eigen;


// 2D Harmonic Green's function for poisson pde
double G(RowVector3d x, RowVector3d y, double R) {
	double pi = 3.14;
	double GrR =  1 / (2 * pi) * log(R /(y - x).norm());
	if (isnan(GrR)) return 0;
	return GrR;
}

void poisson_2D(
	const Eigen::MatrixXd& V,
	const Eigen::MatrixXi& F,
	double (*B)(Eigen::Vector2d),
	double (*f)(Eigen::RowVector3d),
	const Eigen::MatrixXd& P,
	Eigen::VectorXd& U)
{
	// TODO: take stoping tolerance as parameter
	double eps = 1e-8;

	igl::AABB<Eigen::MatrixXd, 3> tree;
	tree.init(V, F);

	Eigen::MatrixXd Q;
	Q.resizeLike(P);


	Eigen::VectorXd sqrD;
	Eigen::VectorXi I;
	Eigen::MatrixXd C;

	int iter = 0;

	U.resizeLike(Q);

	// TODO: check for sqrD and early terminate the loop
	for (int i = 0; i < P.rows(); i++) {
		Q.row(i) = P.row(i);
	}
	vector<int> terminated = vector<int>();
	while (iter < 100) {
		iter++;
		// igl::point_mesh_squared_distance(Q, V, F, sqrD, I, C);
		tree.squared_distance(V, F, Q, sqrD, I, C);

		for (int i = 0; i < Q.rows(); i++) {
			// iterator == terminated.end() means i not in terminated 
			if ((find(terminated.begin(), terminated.end(), i) == terminated.end())) {

				// sample a point on B(x)
				double radius = sqrt(sqrD(i));
				Eigen::RowVector3d center = Q.row(i);

				std::random_device rd;
				std::mt19937 generator(rd());

				// sampling for x_k
				std::uniform_real_distribution<double> uniform01(0.0, 1.0);
				double theta1 = 2 * 3.14 * uniform01(generator);
				Eigen::RowVector3d sample1(cos(theta1), sin(theta1), 0);
				sample1 = sample1 * radius + center;
				Q.row(i) = sample1;
				
				// sampling for y_k
				std::uniform_real_distribution<double> uniform02(0.0, 1.0);
				double theta2 = 2 * 3.14 * uniform01(generator);
				//Eigen::RowVector3d y_k(cos(theta2), sin(theta2), 0);
				Eigen::RowVector3d y_k(cos(theta2), sin(theta2), 0);
				y_k = y_k * radius + center;


				// source f contribution is equal to |B(x_k)| * f(y_k) * G(x_k, y_k) 
				U(i) += 3.14 * radius * radius * f(y_k) * G(sample1, y_k, radius);

				if (sqrD(i) < eps) {
					terminated.push_back(i);
					Q.row(i) = C.row(i);
				}
			}
		}
	}

	// get closest face
	tree.squared_distance(V, F, Q, sqrD, I, C);

	U.resize(P.rows());

	for (int i = 0; i < P.rows(); i++) {
		U(i) += B(C.row(i));
	}


}