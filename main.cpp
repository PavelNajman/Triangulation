#include <gtest/gtest.h>
#include "Poly.h"
#include "PolyAbs.h"
#include "LinearLS.h"
#include "IterativeLinearLS.h"

std::pair<cv::Mat, cv::Mat> SetupGeneralCameraConfiguration()
{
	cv::Mat P0 = cv::Mat::eye(3, 4, CV_64F);
	cv::Mat P1 = cv::Mat::eye(3, 4, CV_64F);

	P0.at<double>(0, 0) = 0.999701;
	P0.at<double>(0, 1) = 0.0174497;
	P0.at<double>(0, 2) = -0.017145;
	P0.at<double>(0, 3) = -500;
	P0.at<double>(1, 0) = -0.0171452;
	P0.at<double>(1, 1) = 0.999695;
	P0.at<double>(1, 2) = 0.0177517;
	P0.at<double>(1, 3) = -100;
	P0.at<double>(2, 0) = 0.0174497;
	P0.at<double>(2, 1) = -0.0174524;
	P0.at<double>(2, 2) = 0.999695;
	P0.at<double>(2, 3) = -100;

	P1.at<double>(0, 0) = 0.99969;
	P1.at<double>(0, 1) = -0.0174497;
	P1.at<double>(0, 2) = 0.0177543;
	P1.at<double>(0, 3) = 500;
	P1.at<double>(1, 0) = 0.0177543;
	P1.at<double>(1, 1) = 0.999695;
	P1.at<double>(1, 2) = -0.0171425;
	P1.at<double>(1, 3) = -100;
	P1.at<double>(2, 0) = -0.0174497;
	P1.at<double>(2, 1) = 0.0174524;
	P1.at<double>(2, 2) = 0.999695;
	P1.at<double>(2, 3) = -100;

	cv::Mat K = cv::Mat::eye(3, 3, CV_64F);
	K.at<double>(0, 0) = 7291.67;
	K.at<double>(1, 1) = 7291.67;
	K.at<double>(0, 2) = 639.5;
	K.at<double>(1, 2) = 511.5;

	P0 = K * P0;
	P1 = K * P1;
	return std::make_pair(P0, P1);
}

std::pair<cv::Mat, cv::Mat> SetupSecondCameraRotatedLeftConfiguration()
{
	cv::Mat P0 = cv::Mat::eye(3, 4, CV_64F);
	cv::Mat P1 = cv::Mat::eye(3, 4, CV_64F);

	P0.at<double>(0, 0) = 0.999701;
	P0.at<double>(0, 1) = 0.0174497;
	P0.at<double>(0, 2) = -0.017145;
	P1.at<double>(0, 3) = -1000;
	P0.at<double>(1, 0) = -0.0171452;
	P0.at<double>(1, 1) = 0.999695;
	P0.at<double>(1, 2) = 0.0177517;
	P0.at<double>(1, 3) = 0;
	P0.at<double>(2, 0) = 0.0174497;
	P0.at<double>(2, 1) = -0.0174524;
	P0.at<double>(2, 2) = 0.999695;
	P0.at<double>(2, 3) = 0;

	cv::Mat K = cv::Mat::eye(3, 3, CV_64F);
	K.at<double>(0, 0) = 7291.67;
	K.at<double>(1, 1) = 7291.67;
	K.at<double>(0, 2) = 639.5;
	K.at<double>(1, 2) = 511.5;

	P0 = K * P0;
	P1 = K * P1;
	return std::make_pair(P0, P1);
}

std::pair<cv::Mat, cv::Mat> SetupSecondCameraRotatedRightConfiguration()
{
	cv::Mat P0 = cv::Mat::eye(3, 4, CV_64F);
	cv::Mat P1 = cv::Mat::eye(3, 4, CV_64F);

	P1.at<double>(0, 0) = 0.999701;
	P1.at<double>(0, 1) = 0.0174497;
	P1.at<double>(0, 2) = -0.017145;
	P1.at<double>(0, 3) = -1000;
	P1.at<double>(1, 0) = -0.0171452;
	P1.at<double>(1, 1) = 0.999695;
	P1.at<double>(1, 2) = 0.0177517;
	P1.at<double>(1, 3) = 0;
	P1.at<double>(2, 0) = 0.0174497;
	P1.at<double>(2, 1) = -0.0174524;
	P1.at<double>(2, 2) = 0.999695;
	P1.at<double>(2, 3) = 0;

	cv::Mat K = cv::Mat::eye(3, 3, CV_64F);
	K.at<double>(0, 0) = 7291.67;
	K.at<double>(1, 1) = 7291.67;
	K.at<double>(0, 2) = 639.5;
	K.at<double>(1, 2) = 511.5;

	P0 = K * P0;
	P1 = K * P1;
	return std::make_pair(P0, P1);
}

std::pair<cv::Mat, cv::Mat> SetupHorizontalConfiguration()
{
	cv::Mat P0 = cv::Mat::eye(3, 4, CV_64F);
	cv::Mat P1 = cv::Mat::eye(3, 4, CV_64F);

	P1.at<double>(0, 3) = -1000;

	cv::Mat K = cv::Mat::eye(3, 3, CV_64F);
	K.at<double>(0, 0) = 7291.67;
	K.at<double>(1, 1) = 7291.67;
	K.at<double>(0, 2) = 639.5;
	K.at<double>(1, 2) = 511.5;

	P0 = K * P0;
	P1 = K * P1;
	return std::make_pair(P0, P1);
}

void EvaluateResult(const cv::Point3d& result, const cv::Point3d& expected_result, double max_percentage_error = 0.001)
{
	double tolerance = cv::norm(expected_result) * max_percentage_error;
	double distance = cv::norm(result - expected_result);
	EXPECT_NEAR(distance, 0.0, tolerance);
}

TEST(PolyTest, GeneralSetup)
{
	cv::Mat P0, P1;
	std::tie(P0, P1) = SetupGeneralCameraConfiguration();
	Triangulation::Poly p(P0, P1);
	cv::Point3d result = p.triangulate(cv::Point2d(146, 642.288), cv::Point2d(1137.31, 385.201));
	cv::Point3d expected_result(0.0, 100.0, 10000.0);
	EvaluateResult(result, expected_result);
}

TEST(PolyTest, RotatedLeft)
{
	cv::Mat P0, P1;
	std::tie(P0, P1) = SetupSecondCameraRotatedLeftConfiguration();
	Triangulation::Poly p(P0, P1);
	cv::Point3d result = p.triangulate(cv::Point2d(878.821, 634.619), cv::Point2d(274.917, 511.5));
	cv::Point3d expected_result(500.0, 0.0, 10000.0);
	EvaluateResult(result, expected_result);
}

TEST(PolyTest, RotatedRight)
{
	cv::Mat P0, P1;
	std::tie(P0, P1) = SetupSecondCameraRotatedRightConfiguration();
	Triangulation::Poly p(P0, P1);
	cv::Point3d result = p.triangulate(cv::Point2d(1004.08, 511.5), cv::Point2d(150.068, 634.618));
	cv::Point3d expected_result(500.0, 0.0, 10000.0);
	EvaluateResult(result, expected_result);
}

TEST(PolyTest, HorizontalStereo)
{
	cv::Mat P0, P1;
	std::tie(P0, P1) = SetupHorizontalConfiguration();
	Triangulation::Poly p(P0, P1);
	cv::Point3d result = p.triangulate(cv::Point2d(1004.08, 511.5), cv::Point2d(274.917, 511.5));
	cv::Point3d expected_result(500.0, 0.0, 10000.0);
	EvaluateResult(result, expected_result);
}

TEST(PolyAbsTest, GeneralSetup)
{
	cv::Mat P0, P1;
	std::tie(P0, P1) = SetupGeneralCameraConfiguration();
	Triangulation::PolyAbs p(P0, P1);
	cv::Point3d result = p.triangulate(cv::Point2d(146, 642.288), cv::Point2d(1137.31, 385.201));
	cv::Point3d expected_result(0.0, 100.0, 10000.0);
	EvaluateResult(result, expected_result);
}

TEST(PolyAbsTest, RotatedLeft)
{
	cv::Mat P0, P1;
	std::tie(P0, P1) = SetupSecondCameraRotatedLeftConfiguration();
	Triangulation::PolyAbs p(P0, P1);
	cv::Point3d result = p.triangulate(cv::Point2d(878.821, 634.619), cv::Point2d(274.917, 511.5));
	cv::Point3d expected_result(500.0, 0.0, 10000.0);
	EvaluateResult(result, expected_result);
}

TEST(PolyAbsTest, RotatedRight)
{
	cv::Mat P0, P1;
	std::tie(P0, P1) = SetupSecondCameraRotatedRightConfiguration();
	Triangulation::PolyAbs p(P0, P1);
	cv::Point3d result = p.triangulate(cv::Point2d(1004.08, 511.5), cv::Point2d(150.068, 634.618));
	cv::Point3d expected_result(500.0, 0.0, 10000.0);
	EvaluateResult(result, expected_result);
}

TEST(PolyAbsTest, HorizontalStereo)
{
	cv::Mat P0, P1;
	std::tie(P0, P1) = SetupHorizontalConfiguration();
	Triangulation::PolyAbs p(P0, P1);
	cv::Point3d result = p.triangulate(cv::Point2d(1004.08, 511.5), cv::Point2d(274.917, 511.5));
	cv::Point3d expected_result(500.0, 0.0, 10000.0);
	EvaluateResult(result, expected_result);
}

TEST(LinearLSTest, GeneralSetup)
{
	cv::Mat P0, P1;
	std::tie(P0, P1) = SetupGeneralCameraConfiguration();
	Triangulation::LinearLS p(P0, P1);
	cv::Point3d result = p.triangulate(cv::Point2d(146, 642.288), cv::Point2d(1137.31, 385.201));
	cv::Point3d expected_result(0.0, 100.0, 10000.0);
	EvaluateResult(result, expected_result);
}

TEST(LinearLSTest, RotatedLeft)
{
	cv::Mat P0, P1;
	std::tie(P0, P1) = SetupSecondCameraRotatedLeftConfiguration();
	Triangulation::LinearLS p(P0, P1);
	cv::Point3d result = p.triangulate(cv::Point2d(878.821, 634.619), cv::Point2d(274.917, 511.5));
	cv::Point3d expected_result(500.0, 0.0, 10000.0);
	EvaluateResult(result, expected_result);
}

TEST(LinearLSTest, RotatedRight)
{
	cv::Mat P0, P1;
	std::tie(P0, P1) = SetupSecondCameraRotatedRightConfiguration();
	Triangulation::LinearLS p(P0, P1);
	cv::Point3d result = p.triangulate(cv::Point2d(1004.08, 511.5), cv::Point2d(150.068, 634.618));
	cv::Point3d expected_result(500.0, 0.0, 10000.0);
	EvaluateResult(result, expected_result);
}

TEST(LinearLSTest, HorizontalStereo)
{
	cv::Mat P0, P1;
	std::tie(P0, P1) = SetupHorizontalConfiguration();
	Triangulation::LinearLS p(P0, P1);
	cv::Point3d result = p.triangulate(cv::Point2d(1004.08, 511.5), cv::Point2d(274.917, 511.5));
	cv::Point3d expected_result(500.0, 0.0, 10000.0);
	EvaluateResult(result, expected_result);
}

TEST(IterativeLinearLSTest, GeneralSetup)
{
	cv::Mat P0, P1;
	std::tie(P0, P1) = SetupGeneralCameraConfiguration();
	Triangulation::IterativeLinearLS p(P0, P1);
	cv::Point3d result = p.triangulate(cv::Point2d(146, 642.288), cv::Point2d(1137.31, 385.201));
	cv::Point3d expected_result(0.0, 100.0, 10000.0);
	EvaluateResult(result, expected_result);
}

TEST(IterativeLinearLSTest, RotatedLeft)
{
	cv::Mat P0, P1;
	std::tie(P0, P1) = SetupSecondCameraRotatedLeftConfiguration();
	Triangulation::IterativeLinearLS p(P0, P1);
	cv::Point3d result = p.triangulate(cv::Point2d(878.821, 634.619), cv::Point2d(274.917, 511.5));
	cv::Point3d expected_result(500.0, 0.0, 10000.0);
	EvaluateResult(result, expected_result);
}

TEST(IterativeLinearLSTest, RotatedRight)
{
	cv::Mat P0, P1;
	std::tie(P0, P1) = SetupSecondCameraRotatedRightConfiguration();
	Triangulation::IterativeLinearLS p(P0, P1);
	cv::Point3d result = p.triangulate(cv::Point2d(1004.08, 511.5), cv::Point2d(150.068, 634.618));
	cv::Point3d expected_result(500.0, 0.0, 10000.0);
	EvaluateResult(result, expected_result);
}

TEST(IterativeLinearLSTest, HorizontalStereo)
{
	cv::Mat P0, P1;
	std::tie(P0, P1) = SetupHorizontalConfiguration();
	Triangulation::IterativeLinearLS p(P0, P1);
	cv::Point3d result = p.triangulate(cv::Point2d(1004.08, 511.5), cv::Point2d(274.917, 511.5));
	cv::Point3d expected_result(500.0, 0.0, 10000.0);
	EvaluateResult(result, expected_result);
}

int main(int argc, char ** argv)
{
	::testing::InitGoogleTest(&argc, argv);
	return RUN_ALL_TESTS();
}
