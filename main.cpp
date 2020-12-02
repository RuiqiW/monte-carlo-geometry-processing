#include "smooth.h"
#include <igl/read_triangle_mesh.h>
#include <igl/readDMAT.h>
#include <igl/parula.h>
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Core>
#include <string>
#include <walk_on_spheres.h>
#include <iostream>

//vector<Segment> scene = {
//   {{ Vec2D(0.5, 0.1), Vec2D(0.9, 0.5) }},
//   {{ Vec2D(0.5, 0.9), Vec2D(0.1, 0.5) }},
//   {{ Vec2D(0.1, 0.5), Vec2D(0.5, 0.1) }},
//   {{ Vec2D(0.5, 0.33333333), Vec2D(0.5, 0.6666666) }},
//   {{ Vec2D(0.33333333, 0.5), Vec2D(0.6666666, 0.5) }}
//};
using namespace Eigen;


int main(int argc, char *argv[])
{
//int main(int argc, char** argv) {

    Matrix3d V(8, 3);

    V.row(0) = Vector3d(0.5, 0.1, 0);
    V.row(1) = Vector3d(0.9, 0.5, 0);
    V.row(2) = Vector3d(0.5, 0.9, 0);
    V.row(3) = Vector3d(0.1, 0.5, 0);
    V.row(4) = Vector3d(0.5, 0.333333, 0);
    V.row(5) = Vector3d(0.5, 0.66666, 0);
    V.row(6) = Vector3d(0.33333, 0.5, 0);
    V.row(7) = Vector3d(0.666666, 0.5, 0);

    MatrixXd F(6, 2);

    F.row(0) = Vector2d(0, 1);
    F.row(1) = Vector2d(2, 3);
    F.row(2) = Vector2d(3, 0);
    F.row(3) = Vector2d(4, 5);
    F.row(4) = Vector2d(7, 6);

    MatrixXd B = V;

    int n = 5;
    MatrixXd P(n, 3);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            P.row(5 * i + j) = 1.0 / n * Vector3d(i, j, 0);
        }
    }

    VectorXd U;
    walk_on_spheres(V, F, B, P, U);

    //srand(time(NULL));
    //ofstream out("out.csv");

    //int s = 128; // image size
    //for (int j = 0; j < s; j++)
    //{
    //    cerr << "row " << j << " of " << s << endl;
    //    for (int i = 0; i < s; i++)
    //    {
    //        Vec2D x0((float)i / (float)s, (float)j / (float)s);
    //        float u = solve(x0, scene, checker);
    //        out << u;
    //        if (i < s - 1) out << ",";
    //    }
    //    out << endl;
    //}
    //return 0;
}


//  // Load input meshes
//  Eigen::MatrixXd OV,V,U;
//  Eigen::MatrixXi F;
//  double lambda = 1e-5;
//  igl::read_triangle_mesh(
//    (argc>1?argv[1]:"../data/sphere-noisy.obj"),OV,F);
//
//
//  // Load data into MatrixXd rather than VectorXd for simpler `smooth` API
//  // Just use y-coordinates as data to be smoothed
//
//
//
//
//
//
//
//  Eigen::MatrixXd G = OV.col(1);
//  if(argc>2)
//  {
//    if(argv[2][0] == 'n')
//    {
//      // Corrupt with a bit of noise
//      G += 0.1*(G.maxCoeff()-G.minCoeff())*
//        Eigen::MatrixXd::Random(G.rows(),G.cols());
//    }else
//    {
//      igl::readDMAT(argv[2],G);
//    }
//  }
//
//  for (int i = 0; i < V.rows(); i++) {
//      V(i, 2) = 0;
//  }
//
//  igl::opengl::glfw::Viewer viewer;
//  std::cout<<R"(
//  D,d  smooth data
//  K    decamate(?) lambda
//  k    decimate lambda
//  L    Toggle lighting
//  M,m  smooth mesh geometry
//  R,r  reset mesh geometry and data
//  L    lighting
//)";
//  const auto & update = [&]()
//  {
//    if((V.array() != V.array()).any())
//    {
//      std::cout<<"Too degenerate to keep smoothing. Better reset"<<std::endl;
//    }
//    viewer.data().set_mesh(V,F);
//    viewer.data().compute_normals();
//    Eigen::MatrixXd C;
//    igl::parula(U,G.minCoeff(),G.maxCoeff(),C);
//    viewer.data().set_colors(C);
//  };
//  const auto & reset = [&]()
//  {
//    V = OV;
//    U = G;
//  };
//  viewer.callback_key_pressed = 
//    [&](igl::opengl::glfw::Viewer &, unsigned int key, int)
//  {
//    switch(key)
//    {
//      case 'D':
//      case 'd':
//        //////////////////////////////////////////////////////////////////////
//        // Smooth data
//        //////////////////////////////////////////////////////////////////////
//        // Use copy constructor to fake in-place API (may be overly
//        // conservative depending on your implementation)
//        smooth(V,F,Eigen::MatrixXd(U),lambda,U);
//
//        //walk_on_spheres(V, F, );
//        
//        break;
//      case 'K':
//      case 'k':
//        lambda = (key=='K'?10.0:0.1)*lambda;
//        std::cout<<"lambda: "<<lambda<<std::endl;
//        break;
//      case 'L':
//        // Toggle lighting
//        viewer.core().lighting_factor = 1.0- viewer.core().lighting_factor;
//        break;
//      case 'M':
//      case 'm':
//      {
//        //////////////////////////////////////////////////////////////////////
//        // Smooth mesh geometry. 
//        //////////////////////////////////////////////////////////////////////
//        // "Linearize" simply by conducting smooth step assuming that vertex
//        // data is a signal defined over current surface: copy is needed to
//        // prevent memory "aliasing"
//        Eigen::MatrixXd Vcpy(V);
//        smooth(Vcpy,F,Vcpy,lambda,V);
//        break;
//      }
//      case 'R':
//      case 'r':
//        reset();
//        break;
//      default:
//        return false;
//    }
//    update();
//    return true;
//  };
//
//  reset();
//  update();
//  viewer.launch();
//
//  return EXIT_SUCCESS;
//}
