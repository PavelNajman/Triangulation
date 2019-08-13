#include "Poly.h"

#include <gsl/gsl_poly.h>

Poly::Poly(const Poly::Fundamental& F)
	: P0(cv::Mat::eye(3, 4, CV_64F)),
	P1(CameraProjectionMatrixFromFundamentalMatrix(F)),
	F(F)
{}
	
Poly::Poly(const cv::Mat& P0, const cv::Mat& P1)
	: P0(P0), P1(P1), F(ComputeFundamentalMatrix(P0, P1))
{}
	
Poly::Poly(const cv::Mat& P0, const cv::Mat& P1, const Poly::Fundamental& F)
	: P0(P0), P1(P1), F(F)
{}

std::tuple<Poly::Intrinsic, Poly::Intrinsic, cv::Mat, cv::Mat> Poly::SetOriginToCamera(const cv::Mat& P0, const cv::Mat& P1)
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
	
Poly::Fundamental Poly::ComputeFundamentalMatrix(const cv::Mat& P0, const cv::Mat& P1)
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

cv::Point3d Poly::triangulate(const cv::Point2d& p0, const cv::Point2d& p1)
{
	cv::Point2d x0, x1;
	std::tie(x0, x1) = ComputeCorrectedCorrespondences(p0, p1);
	return TriangulateDLT(x0, x1);	
}

cv::Point3d Poly::TriangulateDLT(const cv::Point2d& p0, const cv::Point2d& p1)
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

std::pair<cv::Point2d, cv::Point2d> Poly::ComputeCorrectedCorrespondences(const cv::Point2d& p0, const cv::Point2d& p1)
{
	cv::Mat T0 = TranslateToOrigin(p0);
	cv::Mat T1 = TranslateToOrigin(p1);

	cv::Mat f = T1.inv().t() * cv::Mat(F) * T0.inv();

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
	
cv::Mat Poly::TranslateToOrigin(const cv::Point2d& p) const
{
	cv::Mat result = cv::Mat::eye(3, 3, CV_64F);
	result.at<double>(0, 2) = -p.x;
	result.at<double>(1, 2) = -p.y;
	return result;
}
	
Poly::Epipole Poly::ComputeLeftEpipole(const cv::Mat& F)
{
	cv::Mat values, vectors;
	cv::eigen(F.inv(), values, vectors);
	return Epipole(vectors.row(0));
}

Poly::Epipole Poly::ComputeRightEpipole(const cv::Mat& F)
{
	cv::Mat values, vectors;
	cv::eigen(F, values, vectors);
	return Epipole(vectors.row(0));
}

cv::Mat Poly::FormRotationMatrix(const Epipole& e)
{
	cv::Mat result = cv::Mat::eye(3, 3, CV_64F);
	result.at<double>(0,0) = e.x;
	result.at<double>(0,1) = e.y;
	result.at<double>(1,0) = -e.y;
	result.at<double>(1,1) = e.x;
	return result;
}

std::vector<double> Poly::PreparePolyCoeffs(const Poly::PolyParams& params)
{
	double a, b, c, d, e, f;
	std::tie(a, b, c, d, e, f) = params;
	std::vector<double> result = {
		// a*b*d^2 - b^2*c*d
		a * b * d*d - b*b * c * d,
		// - 2*b^2*d^2*f^2*z + a^2*d^2*z - d^4*f^4*z - b^2*c^2*z - b^4*z
		- 2 * b*b * d*d * f*f + a*a * d*d - d*d*d*d * f*f*f*f - b*b * c*c - b*b*b*b,
		// - 4*b^2*c*d*f^2*z^2 - 2*b^2*c*d*e^2*z^2 - 4*a*b*d^2*f^2*z^2 + 2*a*b*d^2*e^2*z^2 + a^2*c*d*z^2 - 4*c*d^3*f^4*z^2 - a*b*c^2*z^2 - 4*a*b^3*z^2
		- 4 * b*b * c * d * f*f - 2 * b*b * c * d * e*e - 4 * a * b * d*d * f*f + 2 * a * b * d*d * e*e + a*a * c * d - 4 * c * d*d*d * f*f*f*f - a * b * c*c - 4 * a * b*b*b,
		// - 8*a*b*c*d*f^2*z^3 - 6*c^2*d^2*f^4*z^3 - 2*b^2*c^2*f^2*z^3 - 2*a^2*d^2*f^2*z^3 - 2*b^2*c^2*e^2*z^3 + 2*a^2*d^2*e^2*z^3 - 6*a^2*b^2*z^3
		- 8 * a * b * c * d * f*f - 6 * c*c * d*d * f*f*f*f - 2 * b*b * c*c * f*f - 2 * a*a * d*d * f*f - 2 * b*b * c*c * e*e + 2 * a*a * d*d * e*e - 6 * a*a * b*b,
		// - 4*a^2*c*d*f^2*z^4 + 2*a^2*c*d*e^2*z^4 - 4*a*b*c^2*f^2*z^4 + a*b*d^2*e^4*z^4 - 2*a*b*c^2*e^2*z^4 - b^2*c*d*e^4*z^4 - 4*c^3*d*f^4*z^4 - 4*a^3*b*z^4
		- 4 * a*a * c*d * f*f + 2 * a*a * c * d * e*e - 4 * a * b * c*c * f*f + a * b * d*d * e*e*e*e - 2 * a * b * c*c * e*e - b*b * c * d * e*e*e*e - 4 * c*c*c * d *  f*f*f*f - 4 * a*a*a * b,
		// a^2*d^2*e^4*z^5 - 2*a^2*c^2*f^2*z^5 - b^2*c^2*e^4*z^5 - c^4*f^4*z^5 - a^4*z^5
		a*a * d*d * e*e*e*e - 2 * a*a * c*c * f*f - b*b * c*c * e*e*e*e - c*c*c*c * f*f*f*f - a*a*a*a,
		// a^2*c*d*e^4*z^6 - a*b*c^2*e^4*z^6
		a*a * c * d * e*e*e*e - a * b * c*c * e*e*e*e
	};
	result.resize(FindPolynomialOrder(result) + 1);
	return result;
}

int Poly::FindPolynomialOrder(const std::vector<double>& coeffs)
{
	for (int n = coeffs.size()-1; n >= 0; --n)
	{
		if (coeffs[n] != 0)
		{
			return n;
		}
	}
	return -1;
}

Poly::Roots Poly::Solve(const PolyParams& params)
{
	std::vector<double> coeffs = PreparePolyCoeffs(params);
	double z[2 * (coeffs.size() - 1)];
	gsl_poly_complex_workspace* w = gsl_poly_complex_workspace_alloc(coeffs.size());
	gsl_poly_complex_solve(&coeffs[0], coeffs.size(), w, z);
	gsl_poly_complex_workspace_free(w);

	Roots result(coeffs.size()-1);
	for (int i = 0; i < coeffs.size() - 1; ++i)
	{
		result[i] = std::complex<double>(z[2*i], z[2*i+1]);
	}
	return result;
}
	
std::vector<double> Poly::EvaluateRootsCosts(const Poly::Roots& roots, const Poly::PolyParams& params)
{
	double a, b, c, d, e, f;
	std::tie(a, b, c, d, e, f) = params;
	auto cost_function = [&](double t)
	{
		return (t*t / (1 + e*e * t*t)) + ((c*t+d) * (c*t+d) / ((a*t+b) * (a*t+b) + f*f * (c*t+d) * (c*t+d)));
	};

	std::vector<double> result(roots.size());
	for (int i = 0; i < roots.size(); ++i)
	{
		result[i] = cost_function(roots[i].real());
	}
	return result;
}

double Poly::EvaluateRoots(const Poly::Roots& roots, const PolyParams& params)
{
	std::vector<double> costs = EvaluateRootsCosts(roots, params);
	return roots[std::min_element(costs.begin(), costs.end()) - costs.begin()].real();
}

std::pair<Poly::Line, Poly::Line> Poly::ConstructLines(double t, const Poly::PolyParams& params)
{
	double a, b, c, d, e, f;
	std::tie(a, b, c, d, e, f) = params;
	Line l0(t*e, 1, -t);
	Line l1(-f * (c*t+d), a*t+b, c*t+d);
	return std::make_pair(l0, l1);
}

cv::Point3d Poly::FindPointOnLineClosestToOrigin(const Poly::Line& l)
{
	return cv::Point3d(-l[0]*l[2], -l[1]*l[2], l[0]*l[0] + l[1]*l[1]);
}
	
cv::Point3d Poly::TransferPointToOriginalCoordinates(const cv::Point3d& p, const cv::Mat& R, const cv::Mat& T)
{
	cv::Mat x(3, 1, CV_64F);
	x.at<double>(0) = p.x;
	x.at<double>(1) = p.y;
	x.at<double>(2) = p.z;
	x = T.inv() * R.t() * x;
	return cv::Point3d(x.at<double>(0), x.at<double>(1), x.at<double>(2));
}
	
cv::Mat Poly::CameraProjectionMatrixFromFundamentalMatrix(const Poly::Fundamental& F)
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
