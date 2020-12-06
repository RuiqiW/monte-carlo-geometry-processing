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
//   B  #V by 1 list of Dirichlet boundary conditions
//   P  #P by 3 list of query positions
// Outputs:
//   U  #P by 1 list of values at query positions

void walk_on_spheres_2D(
	const Eigen::MatrixXd& V,
	const Eigen::MatrixXi& F,
	double (*B)(Eigen::Vector2d),
    const Eigen::MatrixXd & P,
  Eigen::VectorXd & U);


#endif