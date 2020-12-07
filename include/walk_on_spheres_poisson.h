#ifndef WALK_ON_SPHERES_POISSONH
#define WALK_ON_SPHERES_POISSONH
#include <Eigen/Core>
#include <Eigen/Sparse>
// Solve ∆u = 0 over space at given poinst P subject to B on the given boundary
// mesh (V,F)
//
// Inputs:
//   V  #V by 3 list of surface mesh vertex positions
//   F  #F by 3 list of triangles 
//   B  boundary condition function of the form B:R^3 -> R
//   f  source function of the form B:R^3 -> R
//   P  #P by 3 list of query positions
// Outputs:
//   U  #P by 1 list of values at query positions

void walk_on_spheres_poisson(
	const Eigen::MatrixXd& V,
	const Eigen::MatrixXi& F,
	double (*B)(Eigen::Vector3d),
	double (*f)(Eigen::Vector3d),
	const Eigen::MatrixXd& P,
	Eigen::VectorXd& U,
	const Eigen::RowVector3d sourcePoint = Eigen::RowVector3d(0, 0, 0),
	const double c = 1,
	bool use_importance = true);

#endif