#ifndef WALK_ON_SPHERES_H
#define WALK_ON_SPHERES_H
#include <Eigen/Core>
#include <Eigen/Sparse>
// Solve ∆u = 0 over space at given poinst P subject to B on the given boundary
// mesh (V,F)
//
// Inputs:
//   V  #V by 3 list of surface mesh vertex positions
//   F  #F by 3 list of triangles 
//   B  boundary condition function of the form B:R^2 -> R
//   f  source function of the form B:R^2 -> R
//   P  #P by 3 list of query positions
// Outputs:
//   U  #P by 1 list of values at query positions

void poisson_2D(
	const Eigen::MatrixXd& V,
	const Eigen::MatrixXi& F,
	double (*B)(Eigen::Vector2d),
	double (*f)(Eigen::RowVector3d),
	const Eigen::MatrixXd& P,
	Eigen::VectorXd& U);


#endif