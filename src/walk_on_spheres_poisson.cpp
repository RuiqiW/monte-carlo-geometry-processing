#include "../include/walk_on_spheres_poisson.h"
#include <igl/AABB.h>
#include <igl/point_mesh_squared_distance.h>
#include <random>
#include <iostream>

using namespace std;
void walk_on_spheres_poisson(
	const Eigen::MatrixXd& V,
	const Eigen::MatrixXi& F,
	double (*B)(Eigen::Vector3d),
	double (*f)(Eigen::Vector3d),
	const Eigen::MatrixXd& P,
	Eigen::VectorXd& U,
	const Eigen::RowVector3d sourcePoint,
	const double c,
	const bool use_importance
	)
{
	// TODO: parameters
	double eps = 0.01;
	//bool use_importance = true;
	//double c = 10000.0;
	// Eigen::RowVector3d sourcePoint(-0.02, 0.09, -0.002);
	//Eigen::RowVector3d sourcePoint(0.5, 0.5, 0.5);
	//c = 10000.0;

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

    U = Eigen::VectorXd::Zero(P.rows());


	int iter = 0;

	// TODO: check for sqrD and early terminate the loop
	// vector<int> terminated = vector<int>();
	while (iter < 5) {
		iter++;
		// igl::point_mesh_squared_distance(Q, V, F, sqrD, I, C);
		tree.squared_distance(V, F, Q, sqrD, I, C);
		
		for (int i = 0; i < Q.rows(); i++) {
			// iterator == terminated.end() means i not in terminated 
			// if ((find(terminated.begin(), terminated.end(), i) == terminated.end())) {

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
				if(use_importance){
					double k = 0;
					double fy = c;
					double volume = 1.0/3.0 * radius * radius;

					// source point inside B(x)
					double dist = (sourcePoint - center).norm();
					if(dist < radius){
						k = dist/radius;
						U(i) += volume * fy * (1-k) / k;

					// not inside B(x)
					}else{
						k = std::cbrt(uniform01(generator));

						double theta_y = 2 * 3.14 * uniform01(generator);
                		double phi_y = 3.14 * uniform01(generator);
                		double xx = sin(phi_y) * cos(theta_y);
                		double yy = sin(phi_y) * cos(theta_y);
                		double zz = cos(phi_y); 

                		Eigen::RowVector3d sample_y(xx, yy, zz);
                		sample_y = sample_y * k * radius + center;

						double r2 = (sample_y - sourcePoint).squaredNorm();

						//fy = c * std::pow(exp(1.0), -r2);
      //          		U(i) += volume * fy * (1-k) / k;
						double volume = 1.0 / 3.0 * radius * radius; //  4 * pi *radius canceled by g(x, y)
						u(i) += volume * f(sample_y.transpose()) * (1 - k) / k;
					}
				}else{
                	// r/R
                	double k = std::cbrt(uniform01(generator));
                	double theta_y = 2 * 3.14 * uniform01(generator);
                	double phi_y = 3.14 * uniform01(generator);
                	double xx = sin(phi_y) * cos(theta_y);
                	double yy = sin(phi_y) * cos(theta_y);
                	double zz = cos(phi_y); 

                	Eigen::RowVector3d sample_y(xx, yy, zz);
                	sample_y = sample_y * k * radius + center;
                
					// double fy = 1;
					//double r2 = (sample_y - sourcePoint).squaredNorm();
					//double fy = c * std::pow(exp(1.0), -r2);

					//double volume = 1.0/3.0 * radius * radius; //  4 * pi *radius canceled by G(x, y)
     //           	// (R-r)/rR = (1-k)R/kRR = (1-k)/(k * radius)
     //           	U(i) += volume * fy * (1-k) / k;
					double volume = 1.0 / 3.0 * radius * radius; //  4 * pi *radius canceled by G(x, y)
					U(i) += volume * f(sample_y) * (1 - k) / k;
				}
				// if (sqrD(i) < eps) {
				// 	terminated.push_back(i);
				// 	Q.row(i) = C.row(i);
				// }
			// }
		} // end looping over points
	} // end while

	// get closest face
	tree.squared_distance(V, F, Q, sqrD, I, C);


	for (int i = 0; i < P.rows(); i++) {

		U(i) += B(C.row(i));
		 //U(i) += 1.0 / (C.row(i) - sourcePoint).norm();
		// U(i) += 1.0 / C.row(i).norm();
	}


}