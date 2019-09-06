#include "LinearEigen.h"

namespace Triangulation {

cv::Point3d LinearEigen::triangulate(const cv::Point2d& p0, const cv::Point2d& p1) const
{
	cv::Mat A = cv::Mat::zeros(4, 4, CV_64F);
	A.row(0) = p0.x * P0.row(2) - P0.row(0);
	A.row(1) = p0.y * P0.row(2) - P0.row(1);
	A.row(2) = p1.x * P1.row(2) - P1.row(0);
	A.row(3) = p1.y * P1.row(2) - P1.row(1);
	cv::Mat eigen_vectors;
	std::vector<double> eigen_values;
	cv::eigen(A.t() * A, eigen_values, eigen_vectors);
	cv::Vec4d tmp = eigen_vectors.row(eigen_vectors.rows - 1);
	return cv::Point3d(tmp[0]/tmp[3], tmp[1]/tmp[3], tmp[2]/tmp[3]);
}

}
