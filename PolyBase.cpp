#include "PolyBase.h"

namespace Triangulation {

PolyBase::PolyBase(const PolyBase::Fundamental& F)
	: TriangulationBase(cv::Mat::eye(3, 4, CV_64F), CameraProjectionMatrixFromFundamentalMatrix(F)),
	F(F), LS(P0, P1)
{}

PolyBase::PolyBase(const cv::Mat& P0, const cv::Mat& P1)
	: TriangulationBase(P0, P1), F(ComputeFundamentalMatrix(P0, P1)), LS(P0, P1)
{}

PolyBase::PolyBase(const cv::Mat& P0, const cv::Mat& P1, const PolyBase::Fundamental& F)
	: TriangulationBase(P0, P1), F(F), LS(P0, P1)
{}

std::tuple<PolyBase::Intrinsic, PolyBase::Intrinsic, cv::Mat, cv::Mat> PolyBase::SetOriginToCamera(const cv::Mat& P0, const cv::Mat& P1) const
{
	Intrinsic K0(3, 3, CV_64F);
	Intrinsic K1(3, 3, CV_64F);
	cv::Mat R0(3, 3, CV_64F);
	cv::Mat R1(3, 3, CV_64F);
	cv::Mat T0(4, 1, CV_64F);
	cv::Mat T1(4, 1, CV_64F);
	cv::decomposeProjectionMatrix(P0, K0, R0, T0);
	cv::decomposeProjectionMatrix(P1, K1, R1, T1);

	cv::Mat M = cv::Mat::eye(4, 4, CV_64F);
	M(cv::Rect(0,0,3,3)) = R0.inv();
	M.at<double>(0, 3) = (T0.at<double>(0) / T0.at<double>(3));
	M.at<double>(1, 3) = (T0.at<double>(1) / T0.at<double>(3));
	M.at<double>(2, 3) = (T0.at<double>(2) / T0.at<double>(3));

	// K0.inv() * P0 * M  - should be identity
	cv::Mat tmp = cv::Mat(K1).inv() * P1 * M;

	cv::Mat R = tmp(cv::Rect(0,0,3,3));
	cv::Mat T = tmp(cv::Rect(3,0,1,3));

	return std::make_tuple(K0, K1, R,T);
}

PolyBase::Fundamental PolyBase::ComputeFundamentalMatrix(const cv::Mat& P0, const cv::Mat& P1) const
{
	cv::Mat R, T;
	Intrinsic K0, K1;
	std::tie(K0, K1, R,T) = SetOriginToCamera(P0, P1);

	cv::Mat A = cv::Mat(K0) * R.t() * T;

	cv::Mat C = cv::Mat::zeros(3, 3, CV_64F);
	C.at<double>(0, 1) = -A.at<double>(2);
	C.at<double>(0, 2) = A.at<double>(1);
	C.at<double>(1, 0) = A.at<double>(2);
	C.at<double>(1, 2) = -A.at<double>(0);
	C.at<double>(2, 0) = -A.at<double>(1);
	C.at<double>(2, 1) = A.at<double>(0);

	return cv::Mat_<double>(cv::Mat(K1).inv().t() * R * cv::Mat(K0).t() * C);
}

cv::Point3d PolyBase::triangulate(const cv::Point2d& p0, const cv::Point2d& p1) const
{
	cv::Point2d x0, x1;
	std::tie(x0, x1) = ComputeCorrectedCorrespondences(p0, p1);
	//std::cout << x0 << std::endl;
	//std::cout << x1 << std::endl;
	return LS.triangulate(x0, x1);
}

std::pair<cv::Point2d, cv::Point2d> PolyBase::ComputeCorrectedCorrespondences(const cv::Point2d& p0, const cv::Point2d& p1) const
{
	cv::Mat T0 = TranslateToOrigin(p0);
	cv::Mat T1 = TranslateToOrigin(p1);

	cv::Mat f = T1.t() * cv::Mat(F) * T0;

	Epipole e0 = ComputeRightEpipole(f);
	Epipole e1 = ComputeLeftEpipole(f);

	cv::Mat R0 = FormRotationMatrix(e0);
	cv::Mat R1 = FormRotationMatrix(e1);

	f = R1 * f * R0.t();

	PolyParams params = std::make_tuple(f.at<double>(1, 1), f.at<double>(1, 2), f.at<double>(2, 1), f.at<double>(2, 2), e0.z, e1.z);
	Roots roots = Solve(params);

	double t = EvaluateRoots(roots, params);
	Line l0, l1;
	std::tie(l0, l1) = ConstructLines(t, params);

	cv::Point3d x0 = FindPointOnLineClosestToOrigin(l0);
	cv::Point3d x1 = FindPointOnLineClosestToOrigin(l1);

	x0 = TransferPointToOriginalCoordinates(x0, R0, T0);
	x1 = TransferPointToOriginalCoordinates(x1, R1, T1);

	return std::make_pair(cv::Point2d(x0.x / x0.z, x0.y / x0.z), cv::Point2d(x1.x / x1.z, x1.y / x1.z));
}

cv::Mat PolyBase::TranslateToOrigin(const cv::Point2d& p) const
{
	cv::Mat result = cv::Mat::eye(3, 3, CV_64F);
	result.at<double>(0, 2) = p.x;
	result.at<double>(1, 2) = p.y;
	return result;
}

PolyBase::Epipole PolyBase::ComputeLeftEpipole(const cv::Mat& F) const
{
	cv::Mat W, U, VT;
	cv::SVD::compute(F.t(), W, U, VT);
	return Epipole(VT.row(2));
}

PolyBase::Epipole PolyBase::ComputeRightEpipole(const cv::Mat& F) const
{
	cv::Mat W, U, VT;
	cv::SVD::compute(F, W, U, VT);
	return Epipole(VT.row(2));
}

cv::Mat PolyBase::FormRotationMatrix(const Epipole& e) const
{
	cv::Mat result = cv::Mat::eye(3, 3, CV_64F);
	result.at<double>(0,0) = e.x;
	result.at<double>(0,1) = e.y;
	result.at<double>(1,0) = -e.y;
	result.at<double>(1,1) = e.x;
	return result;
}

int PolyBase::FindPolynomialOrder(const std::vector<double>& coeffs) const
{
	for (int n = static_cast<int>(coeffs.size())-1; n >= 0; --n)
	{
		if (coeffs[n] != 0)
		{
			return n;
		}
	}
	return -1;
}

PolyBase::Roots PolyBase::Solve(const PolyParams& params) const
{
	std::vector<double> coeffs = PreparePolyCoeffs(params);
	if (coeffs.size() <= 1)
	{
		return {0};
	}
	std::vector<cv::Vec2d> roots;
	cv::solvePoly(coeffs, roots);

	Roots result(roots.size());
	for (size_t i = 0u; i < roots.size(); ++i)
	{
		result[i] = std::complex<double>(roots[i][0], roots[i][1]);
	}
	return result;
}

double PolyBase::EvaluateRoots(const PolyBase::Roots& roots, const PolyParams& params) const
{
	std::vector<double> costs = EvaluateRootsCosts(roots, params);
	return roots[std::min_element(costs.begin(), costs.end()) - costs.begin()].real();
}

std::pair<PolyBase::Line, PolyBase::Line> PolyBase::ConstructLines(double t, const PolyBase::PolyParams& params) const
{
	double a, b, c, d, e, f;
	std::tie(a, b, c, d, e, f) = params;
	Line l0(t*e, 1, -t);
	Line l1(-f * (c*t+d), a*t+b, c*t+d);
	return std::make_pair(l0, l1);
}

cv::Point3d PolyBase::FindPointOnLineClosestToOrigin(const PolyBase::Line& l) const
{
	return cv::Point3d(-l[0]*l[2], -l[1]*l[2], l[0]*l[0] + l[1]*l[1]);
}

cv::Point3d PolyBase::TransferPointToOriginalCoordinates(const cv::Point3d& p, const cv::Mat& R, const cv::Mat& T) const
{
	cv::Mat x(3, 1, CV_64F);
	x.at<double>(0) = p.x;
	x.at<double>(1) = p.y;
	x.at<double>(2) = p.z;
	x = T * R.t() * x;
	return cv::Point3d(x.at<double>(0), x.at<double>(1), x.at<double>(2));
}

cv::Mat PolyBase::CameraProjectionMatrixFromFundamentalMatrix(const PolyBase::Fundamental& F) const
{
	cv::Mat e2;
	cv::SVD::solveZ(F.t(), e2);
	cv::Matx33d e2x =
	{
		0, -e2.at<double>(2), e2.at<double>(1),
		e2.at<double>(2), 0, -e2.at<double>(0),
		-e2.at<double>(1), e2.at<double>(0), 0
	};
	cv::Mat result = cv::Mat::eye(3, 4, CV_64FC1);
	result(cv::Rect(0, 0, 3, 3)) = cv::Mat(e2x * F);
	result.at<double>(0, 3) = e2.at<double>(0);
	result.at<double>(1, 3) = e2.at<double>(1);
	result.at<double>(2, 3) = e2.at<double>(2);
	return result;
}

}
