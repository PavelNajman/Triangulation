#include "LinearLS.h"

namespace Triangulation {

cv::Point3d LinearLS::triangulate(const cv::Point2d& p0, const cv::Point2d& p1) const
{
	cv::Mat A = cv::Mat::zeros(4, 4, CV_64F);
	A.row(0) = p0.x * P0.row(2) - P0.row(0);
	A.row(1) = p0.y * P0.row(2) - P0.row(1);
	A.row(2) = p1.x * P1.row(2) - P1.row(0);
	A.row(3) = p1.y * P1.row(2) - P1.row(1);
	cv::Mat result;
	cv::SVD::solveZ(A, result);
	return cv::Point3d(result.at<double>(0) / result.at<double>(3), result.at<double>(1) / result.at<double>(3), result.at<double>(2) / result.at<double>(3));
}

}
