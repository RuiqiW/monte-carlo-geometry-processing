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

	// TODO: check for sqrD and early terminate the loop
	for (int i = 0; i < P.rows(); i++) {
		Q.row(i) = P.row(i);
	}
	// vector<int> terminated = vector<int>();
	while (iter < 5) {
		iter++;
		// igl::point_mesh_squared_distance(Q, V, F, sqrD, I, C);
		tree.squared_distance(V, F, Q, sqrD, I, C);
		
		for (int i = 0; i < Q.rows(); i++) {
			// iterator == terminated.end() means i not in terminated 
			// if ((find(terminated.begin(), terminated.end(), i) == terminated.end())) {

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

				// double x = cos(theta);
				// double y = sin(theta);
				// Eigen::RowVector3d sample(x, y, 0);
				
				
				sample = sample * radius + center;
				Q.row(i) = sample;
				// if (sqrD(i) < eps) {
				// 	terminated.push_back(i);
				// 	Q.row(i) = C.row(i);
				// }
			// }
		} // end looping over points
	} // end while

	// get closest face
	// igl::point_mesh_squared_distance(Q, V, F, sqrD, I, C);
	tree.squared_distance(V, F, Q, sqrD, I, C);

	U.resize(P.rows());

	for (int i = 0; i < P.rows(); i++) {

		// U(i) = ((int)floor(6 * C(i, 0)) + (int) floor(6 * C(i, 1))) % 2;
		U(i) = C.row(i).norm();

		//// choose the closest vertex (TODO: maybe interpolate the values at the vertices of the closest face?)
		//Eigen::RowVector3i face = F.row(I(i));
		//int closest_idx = face(0);
		//double closest_dist = (V.row(face(0)) - P.row(i)).squaredNorm();
		//double dist = (V.row(face(1)) - P.row(i)).squaredNorm();
		//if (dist < closest_dist) {
		//	closest_idx = face(1);
		//	closest_dist = dist;
		//}
		//// dist = (V.row(face(2)) - P.row(i)).squaredNorm();
		////  if( dist < closest_dist){
		////   closest_idx = face(2);
		////   closest_dist = dist;
		//// }
		//// get the boundary value
		//U(i) = B(closest_idx);
		//// std::cout << closest_dist << " " << closest_idx << " " << U(i) << std::endl;
	}


}