#include "TriangulationBase.h"

namespace Triangulation {

TriangulationBase::TriangulationBase(const cv::Mat& P0, const cv::Mat& P1)
	: P0(P0), P1(P1)
{}

std::vector<cv::Point3d> TriangulationBase::triangulate(const std::vector<cv::Point2d>& p0, const std::vector<cv::Point2d>& p1) const
{
	std::vector<cv::Point3d> result;
	for (size_t i = 0; i < p0.size(); ++i)
	{
		result.emplace_back(triangulate(p0[i], p1[i]));
	}
	return result;
}

}
