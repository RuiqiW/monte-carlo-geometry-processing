#include "../include/walk_on_spheres_biharmonic.h"
#include <igl/AABB.h>
#include <igl/point_mesh_squared_distance.h>
#include <random>
#include <iostream>


void walk_on_spheres_biharmonic(
	const Eigen::MatrixXd& V,
	const Eigen::MatrixXi& F,
	double (*B)(Eigen::Vector3d),
	double (*h)(Eigen::Vector3d),
	const Eigen::MatrixXd& P,
	Eigen::VectorXd& U
	)
{
	// TODO: parameters
	double eps = 0.01;
    Eigen::RowVector3d sourcePoint(0.37, 64, 0);
	// Eigen::RowVector3d sourcePoint(0.5, 0.5, 0.5);

	igl::AABB<Eigen::MatrixXd, 3> tree;
	tree.init(V, F);

	Eigen::MatrixXd Q;
	Q.resizeLike(P);
    for (int i = 0; i < P.rows(); i++) {
		Q.row(i) = P.row(i);
	}

	Eigen::VectorXd sqrD;
	Eigen::VectorXi I;
	Eigen::MatrixXd C;

    Eigen::VectorXd sqrD_y;
	Eigen::VectorXi I_y;
	Eigen::MatrixXd C_y;

    U = Eigen::VectorXd::Zero(P.rows());


	int iter = 0;

	// TODO: check for sqrD and early terminate the loop
	// vector<int> terminated = vector<int>();
	while (iter < 5) {
		iter++;
		tree.squared_distance(V, F, Q, sqrD, I, C);
		
		for (int i = 0; i < Q.rows(); i++) {

				// sample a point on bd B(x)
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

                // sample a point in B(x)
                // r/R
                double k = std::cbrt(uniform01(generator));
                double theta_y = 2 * 3.14 * uniform01(generator);
                double phi_y = 3.14 * uniform01(generator);
                double xx = sin(phi_y) * cos(theta_y);
                double yy = sin(phi_y) * cos(theta_y);
                double zz = cos(phi_y); 

                Eigen::RowVector3d sample_y(xx, yy, zz);
                sample_y = sample_y * k * radius + center;
                

				double volume = 1.0/3.0 * radius * radius; //  4 * pi *radius canceled by G(x, y)
                // (R-r)/rR = (1-k)R/kRR = (1-k)/(k * radius)
                U(i) += volume * (1-k) / k;
		} // end looping over points
	} // end while

	// get closest face
	tree.squared_distance(V, F, Q, sqrD, I, C);


	for (int i = 0; i < P.rows(); i++) {
		// since the walk from yi connects to xi+1 -> xk, estimate after iteration
		double Vy = C.row(i)(0);
		U(i) *= Vy;
		// U(i) *= h(C.row(i)); 
        U(i) += (C.row(i)-sourcePoint).squaredNorm();
	}


}