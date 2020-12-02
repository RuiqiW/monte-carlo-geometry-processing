#include "cotmatrix.h"

#include <igl/edge_lengths.h>

using namespace std;
using namespace Eigen;

void cotmatrix(
  const Eigen::MatrixXd & l,
  const Eigen::MatrixXi & F,
  Eigen::SparseMatrix<double> & L)
{

	typedef Eigen::Triplet<double> T;

	vector<Triplet<double>> tripletlist;

	for (int i = 0; i < F.rows(); i++) {
		
		int k0 = F(i, 0);
		int k1 = F(i, 1);
		int k2 = F(i, 2);

		double a = l(i, 0);
		double b = l(i, 1);
		double c = l(i, 2);
		double s = (a + b + c) / 2;
		double area = sqrt(s * (s - a) * (s - b) * (s - c));

		double sin_alpha = 2 * area / (b * c);
		double sin_beta = 2 * area / (a * c);
		double sin_gamma = 2 * area / (a * b);

		double cos_alpha = (pow(a, 2) - (pow(b, 2) + pow(c, 2))) / (- 2 * b * c);
		double cos_beta = (pow(b, 2) - (pow(a, 2) + pow(c, 2))) / (- 2 * a * c);
		double cos_gamma = (pow(c, 2) - (pow(a, 2) + pow(b, 2))) / (- 2 * a * b);

		double cot_alpha = cos_alpha / sin_alpha;
		double cot_beta = cos_beta / sin_beta;
		double cot_gamma = cos_gamma / sin_gamma;
		
		tripletlist.push_back(T(k0, k1, 1 / 2.0 * cot_alpha));
		tripletlist.push_back(T(k0, k0, - 1 / 2.0 * cot_alpha));
		tripletlist.push_back(T(k1, k2, 1 / 2.0 * cot_beta ));
		tripletlist.push_back(T(k1, k1, -1 / 2.0 * cot_beta));
		tripletlist.push_back(T(k0, k2, 1 / 2.0 * cot_gamma));
		tripletlist.push_back(T(k2, k2, -1 / 2.0 * cot_gamma));


	}
	L.setFromTriplets(tripletlist.begin(), tripletlist.end());

}

