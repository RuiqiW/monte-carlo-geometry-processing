#include <igl/read_triangle_mesh.h>
#include <igl/readDMAT.h>
#include <igl/parula.h>
#include <igl/opengl/glfw/Viewer.h>
#include <Eigen/Core>
#include <string>
#include <iostream>
#include <walk_on_spheres_2D.h>
#include <walk_on_spheres_3D.h>
#include <fstream>
#include <omp.h>
#include <igl/fast_winding_number.h>
#include <igl/read_triangle_mesh.h>
#include <igl/slice_mask.h>
#include <Eigen/Geometry>
#include <igl/octree.h>
#include <igl/barycenter.h>
#include <igl/knn.h>
#include <igl/random_points_on_mesh.h>
#include <igl/bounding_box_diagonal.h>
#include <igl/per_face_normals.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/get_seconds.h>
#include <cstdlib>
#include "walk_on_spheres_poisson.h"
#include "walk_on_spheres_screened_poisson.h"
#include "walk_on_spheres_biharmonic.h"


using namespace std;
using namespace Eigen;

double c = 10000.0;

Vector3d sourcePoint;

double laplacian_boundary_3D(Vector3d boundary_point) {
    return boundary_point.norm();
}
double poisson_boundary_3D(Vector3d boundary_point) {
    return 1/ ((boundary_point - sourcePoint).norm());
}

double B(Vector3d point) {
    return (point - sourcePoint).squaredNorm();
}

double h(Vector3d point) {
    return point(0);
}



double source(Vector3d point) {

    Eigen::Vector3d sourcePoint(0.5, 0.5, 0.5);
    double r2 = (point - sourcePoint).squaredNorm();
    return c * std::pow(exp(1.0), -r2);
}

int main(int argc, char* argv[])
{
    const auto time = [](std::function<void(void)> func)->double
    {
        const double t_before = igl::get_seconds();
        func();
        const double t_after = igl::get_seconds();
        return t_after - t_before;
    };


    // Read Mesh
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    //igl::read_triangle_mesh(argc > 1 ? argv[1] : "../data/bunny.off", V, F);
    igl::read_triangle_mesh(argc > 1 ? argv[1] : "../data/cactus.obj", V, F);

    // Sample points inside mesh: https://github.com/libigl/libigl/blob/master/tutorial/717_FastWindingNumber

    // Sample mesh for point cloud
    Eigen::MatrixXd P, N;
    {
        Eigen::VectorXi I;
        Eigen::SparseMatrix<double> B;
        igl::random_points_on_mesh(10000, V, F, B, I);
        P = B * V;
        Eigen::MatrixXd FN;
        igl::per_face_normals(V, F, FN);
        N.resize(P.rows(), 3);
        for (int p = 0; p < I.rows(); p++)
        {
            N.row(p) = FN.row(I(p));
        }
    }
    // Build octree
    std::vector<std::vector<int > > O_PI;
    Eigen::MatrixXi O_CH;
    Eigen::MatrixXd O_CN;
    Eigen::VectorXd O_W;
    igl::octree(P, O_PI, O_CH, O_CN, O_W);
    {
        Eigen::MatrixXi I;
        igl::knn(P, 20, O_PI, O_CH, O_CN, O_W, I);
    }

    if (argc <= 1)
    {
        // corrupt mesh
        Eigen::MatrixXd BC;
        igl::barycenter(V, F, BC);
        Eigen::MatrixXd OV = V;
        V.resize(F.rows() * 3, 3);
        for (int f = 0; f < F.rows(); f++)
        {
            for (int c = 0; c < 3; c++)
            {
                int v = f + c * F.rows();
                // random rotation about barycenter
                Eigen::AngleAxisd R(
                    0.5 * static_cast <double> (rand()) / static_cast <double> (RAND_MAX),
                    Eigen::Vector3d::Random(3, 1));
                V.row(v) = (OV.row(F(f, c)) - BC.row(f)) * R.matrix() + BC.row(f);
                F(f, c) = v;
            }
        }
    }

    // Generate a list of random query points in the bounding box
    Eigen::MatrixXd Q = Eigen::MatrixXd::Random(1000, 3);
    const Eigen::RowVector3d Vmin = V.colwise().minCoeff();
    const Eigen::RowVector3d Vmax = V.colwise().maxCoeff();
    const Eigen::RowVector3d Vdiag = Vmax - Vmin;

    sourcePoint = 0.5 * (Vmax + Vmin);


    for (int q = 0; q < Q.rows(); q++)
    {
        Q.row(q) = (Q.row(q).array() * 0.5 + 0.5) * Vdiag.array() + Vmin.array();
    }

    // Positions of points inside of triangle soup (V,F)
    Eigen::MatrixXd QiV;
    {
        igl::FastWindingNumberBVH fwn_bvh;
        printf("triangle soup precomputation (% 8ld triangles): %g secs\n",
            F.rows(),
            time([&]() {igl::fast_winding_number(V.cast<float>().eval(), F, 2, fwn_bvh); }));
        Eigen::VectorXf WiV;
        printf("triangle soup evaluation  (% 8ld queries):   %g secs\n",
            Q.rows(),
            time([&]() {igl::fast_winding_number(fwn_bvh, 2, Q.cast<float>().eval(), WiV); }));
        igl::slice_mask(Q, WiV.array() > 0.5, 1, QiV);
    }


    // Visualization
    igl::opengl::glfw::Viewer viewer;
    // For dislpaying normals as little line segments
    Eigen::MatrixXd PN(2 * P.rows(), 3);
    Eigen::MatrixXi E(P.rows(), 2);
    const double bbd = igl::bounding_box_diagonal(V);
    for (int p = 0; p < P.rows(); p++)
    {
        E(p, 0) = 2 * p;
        E(p, 1) = 2 * p + 1;
        PN.row(E(p, 0)) = P.row(p);
        PN.row(E(p, 1)) = P.row(p) + bbd * 0.01 * N.row(p);
    }

    bool show_P = false;
    int show_Q = 0;

    int query_data = 0;
    viewer.data_list[query_data].set_mesh(V, F);
    viewer.data_list[query_data].clear();
    viewer.data_list[query_data].point_size = 10;
    viewer.append_mesh();
    int object_data = 1;
    viewer.data_list[object_data].set_mesh(V, F);
    viewer.data_list[object_data].point_size = 5;

    int NUM_ITERATIONS = 64;
    const auto update = [&]()
    {
        viewer.data_list[query_data].clear();
        switch (show_Q)
        {
        case 2:
            // show all Q
            viewer.data_list[query_data].set_points(Q, Eigen::RowVector3d(0.996078, 0.760784, 0.760784));
            break;
        case 3:
            // show all Q inside
        {
            VectorXd total_U_laplacian = VectorXd::Zero(QiV.rows());
            for (int k = 0; k < NUM_ITERATIONS; k++) {
                VectorXd U;
                walk_on_spheres_3D(V, F, laplacian_boundary_3D, QiV, U);
                total_U_laplacian += U;
            }
            VectorXd average_U_laplacian = total_U_laplacian / NUM_ITERATIONS;
            Eigen::MatrixXd CM;
            igl::colormap(igl::COLOR_MAP_TYPE_MAGMA, average_U_laplacian, average_U_laplacian.minCoeff(), average_U_laplacian.maxCoeff(), CM);
            viewer.data_list[query_data].set_points(QiV, CM);
            break;
        }
        case 4:
        {
            // Poisson without importance sampling
            VectorXd total_U_poisson_without_importance = VectorXd::Zero(QiV.rows());
            for (int k = 0; k < NUM_ITERATIONS; k++) {
                VectorXd U;
                walk_on_spheres_poisson(V, F, poisson_boundary_3D, source, QiV, U, sourcePoint, c, false);
                total_U_poisson_without_importance += U;
            }
            VectorXd average_U_poisson_without_importance = total_U_poisson_without_importance / NUM_ITERATIONS;
            Eigen::MatrixXd CM;
            igl::colormap(igl::COLOR_MAP_TYPE_MAGMA, average_U_poisson_without_importance, average_U_poisson_without_importance.minCoeff(), average_U_poisson_without_importance.maxCoeff(), CM);
            viewer.data_list[query_data].set_points(QiV, CM);
            break;
        }
        case 5:
        {
            //cout << "test" << endl;
            // Poisson with importance sampling
            VectorXd total_U_poisson_with_importance = VectorXd::Zero(QiV.rows());
            for (int k = 0; k < NUM_ITERATIONS; k++) {
                VectorXd U;
                walk_on_spheres_poisson(V, F, poisson_boundary_3D, source, QiV, U, sourcePoint, c, true);
                total_U_poisson_with_importance += U;
            }
            VectorXd average_U_poisson_with_importance = total_U_poisson_with_importance / NUM_ITERATIONS;
            Eigen::MatrixXd CM;
            igl::colormap(igl::COLOR_MAP_TYPE_MAGMA, average_U_poisson_with_importance, average_U_poisson_with_importance.minCoeff(), average_U_poisson_with_importance.maxCoeff(), CM);
            viewer.data_list[query_data].set_points(QiV, CM);
            break;
        }
        case 6:
        {
            // Screened poisson without importance
            VectorXd total_U_screened_poisson_without_importance = VectorXd::Zero(QiV.rows());
            for (int k = 0; k < NUM_ITERATIONS; k++) {
                VectorXd U;
                walk_on_spheres_screened_poisson(V, F, poisson_boundary_3D, source, QiV, U, sourcePoint, c, false);
                total_U_screened_poisson_without_importance += U;
            }
            VectorXd average_U_screened_poisson_without_importance = total_U_screened_poisson_without_importance / NUM_ITERATIONS;
            Eigen::MatrixXd CM;
            igl::colormap(igl::COLOR_MAP_TYPE_MAGMA, average_U_screened_poisson_without_importance, average_U_screened_poisson_without_importance.minCoeff(), average_U_screened_poisson_without_importance.maxCoeff(), CM);
            viewer.data_list[query_data].set_points(QiV, CM);
            break;

        }
        case 7:
        {
            // Screened poisson with importance
            VectorXd total_U_screened_poisson_with_importance = VectorXd::Zero(QiV.rows());
            for (int k = 0; k < NUM_ITERATIONS; k++) {
                VectorXd U;
                walk_on_spheres_screened_poisson(V, F, poisson_boundary_3D, source, QiV, U, sourcePoint, c, true);
                total_U_screened_poisson_with_importance += U;
            }
            VectorXd average_U_screened_poisson_with_importance = total_U_screened_poisson_with_importance / NUM_ITERATIONS;
            Eigen::MatrixXd CM;
            igl::colormap(igl::COLOR_MAP_TYPE_MAGMA, average_U_screened_poisson_with_importance, average_U_screened_poisson_with_importance.minCoeff(), average_U_screened_poisson_with_importance.maxCoeff(), CM);
            viewer.data_list[query_data].set_points(QiV, CM);
            break;

        }
        case 8:
        {
            // biharmonic 
            VectorXd total_U_biharmonic = VectorXd::Zero(QiV.rows());
            for (int k = 0; k < NUM_ITERATIONS; k++) {
                VectorXd U;
                walk_on_spheres_biharmonic(V, F, poisson_boundary_3D, source, QiV, U);
                total_U_biharmonic += U;
            }
            VectorXd average_U_biharmonic = total_U_biharmonic / NUM_ITERATIONS;
            Eigen::MatrixXd CM;
            igl::colormap(igl::COLOR_MAP_TYPE_MAGMA, average_U_biharmonic, average_U_biharmonic.minCoeff(), average_U_biharmonic.maxCoeff(), CM);
            viewer.data_list[query_data].set_points(QiV, CM);
            break;

        }
        }
        viewer.data_list[object_data].clear();
        if (show_P)
        {
            viewer.data_list[object_data].set_points(P, Eigen::RowVector3d(1, 1, 1));
            viewer.data_list[object_data].set_edges(PN, E, Eigen::RowVector3d(0.8, 0.8, 0.8));
        }
        else
        {
            viewer.data_list[object_data].set_mesh(V, F);
        }
    };

    viewer.callback_key_pressed =
        [&](igl::opengl::glfw::Viewer&, unsigned int key, int mod)
    {
        switch (key)
        {
        default:
            return false;
        case '1':
            show_P = !show_P;
            break;
        case '2':
            //show_Q = (show_Q + 1) % 3
            show_Q = 2;
            break;
        case '3':
            //show_Q = (show_Q + 1) % 3;
            show_Q = 3;
            break;
        case '4':
            //show_Q = (show_Q + 1) % 3;
            show_Q = 4;
            break;
        case '5':
            //show_Q = (show_Q + 1) % 3;
            show_Q = 5;
            break;
        case '6':
            show_Q = 6;
            break;
        case '7':
            show_Q = 7;
            break;
        case '8':
            show_Q = 8;
            break;
        case 'k':
            if (c < 10000) {
                c *= 10;
            }
            break;
        case 'j':
            if (c > 1) {
                c /= 10;
            }
            break;
        }
        update();
        //cout << "c: " << c << endl;
        return true;
    };

    std::cout << R"(
Monte Carlo PDE computations
  1  Shows Sample Points of for Point Cloud of mesh
  2  Shows Query Points
  3  Solves Laplace Equation for Query Points
  4  Solves Poisson Equation for Query Points
  5  Solves Poisson Equation for Query Points with Importance Sampling
  6  Solves Screened Poisson Equation for Query Points
  7  Solves Screened Poisson Equation for Query Points with Importance Sampling
  8  Solves Biharmonic Equation for Query Points
  k  multiplies c by 10 (max 10000)
  j  divides c by 10 (min 1)
)";


    update();
    viewer.launch();
}
