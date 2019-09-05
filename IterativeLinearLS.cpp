#include "IterativeLinearLS.h"
#include "LinearLS.h"

namespace Triangulation {

cv::Point3d IterativeLinearLS::triangulate(const cv::Point2d& p0, const cv::Point2d& p1) const
{
	double EPS = 1e-5;
	double w0 = 1, w1 = 1;
	cv::Mat x; 
	for (int i = 0; i < 10; i++) { // 10 iterations should be enough
		cv::Mat A = cv::Mat::zeros(4, 4, CV_64F);
		A.row(0) = (p0.x * P0.row(2) - P0.row(0))/w0;
		A.row(1) = (p0.y * P0.row(2) - P0.row(1))/w0;
		A.row(2) = (p1.x * P1.row(2) - P1.row(0))/w1;
		A.row(3) = (p1.y * P1.row(2) - P1.row(1))/w1;
		cv::SVD::solveZ(A, x);

		double new_w0 = cv::Mat(P0.row(2)*x).at<double>(0);
		double new_w1 = cv::Mat(P1.row(2)*x).at<double>(0);

		if(std::abs(w0 - new_w0) <= EPS && std::abs(w1 - new_w1) <= EPS)
		{
			break;
		}

		w0 = new_w0;
		w1 = new_w1;
	}
	return cv::Point3d(x.at<double>(0) / x.at<double>(3), x.at<double>(1) / x.at<double>(3), x.at<double>(2) / x.at<double>(3));
}

}
