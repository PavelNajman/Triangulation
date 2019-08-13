#ifndef POLY_H_
#define POLY_H_

#include <array>
#include <complex>
#include <opencv2/opencv.hpp>

/**
 *	\brief	Performs non-linear triangulation of two image points by minimizing geometric error
 */
class Poly
{
public:
	typedef cv::Matx33d Intrinsic;
	typedef cv::Matx33d Fundamental;

	/**
	 *	\brief	Constructor
	 *	\param	F	Fundamental matrix.
	 */
	Poly(const Fundamental& F);
	/**
	 *	\brief	Constructor
	 *	\param	P0	Camera matrix of the first camera.
	 *	\param 	P1	Camera matrix of the second camera.
	 */
	Poly(const cv::Mat& P0, const cv::Mat& P1);
	/**
	 *	\brief	Constructor
	 *	\param	P0	Camera matrix of the first camera.
	 *	\param 	P1	Camera matrix of the second camera.
	 *	\param	F	Fundamental matrix.
	 */
	Poly(const cv::Mat& P0, const cv::Mat& P1, const Fundamental& F);
	/**
	 *	\brief	Triangulates image points.
	 *	\param	p0	Point in the image of the first camera.
	 *	\param	p1	Corresponding point in the image of the second camera.
	 *	\return	Triangulated point.
	 */
	cv::Point3d triangulate(const cv::Point2d& p0, const cv::Point2d& p1);
private:
	typedef cv::Vec3d Line;
	typedef cv::Point3d Epipole;
	typedef std::vector<std::complex<double>> Roots;
	typedef std::tuple<double, double, double, double, double, double> PolyParams;
	/**
	 *	\brief	Computes corrected correspondences, that minimize the geometric error.
	 *	\param	p0	Point in the image of the first camera.
	 *	\param	p1	Corresponding point in the image of the second camera.
	 *	\return	Corrected correspondences in homogeneous coordinates.
	 */
	std::pair<cv::Point2d, cv::Point2d> ComputeCorrectedCorrespondences(const cv::Point2d& p0, const cv::Point2d& p1);
	/**
	 *	\brief	Defines translation matrix that translate the given point to origin.
	 *	\param	p	Translated point.
	 *	\return	Translation matrix.
	 */
	cv::Mat TranslateToOrigin(const cv::Point2d& p) const;
	/**
	 *	\brief	Computes epipole e = (e1, e2, e3) such as eF = 0 and e1*e1 + e2*e2 = 1
	 *	\param	F	Fundamental matrix.
	 *	\return	Computed epipole
	 */
	Epipole ComputeLeftEpipole(const cv::Mat& F);
	/**
	 *	\brief	Computes epipole e = (e1, e2, e3) such as Fe = 0 and e1*e1 + e2*e2 = 1
	 *	\param	F	Fundamental matrix.
	 *	\return	Computed epipole
	 */
	Epipole ComputeRightEpipole(const cv::Mat& F);
	/**
	 *	\brief	Defines rotation matrix using given epipole.
	 *	\param	e	Epipole.
	 *	\return	Rotation matrix.
	 */
	cv::Mat FormRotationMatrix(const Epipole& e);
	/**
	 *	\brief	Prepares polynomial coefficients.
	 *	\param	params	Polynomial coefficients params.
	 *	\return Polynomial coefficients.
	 */
	std::vector<double> PreparePolyCoeffs(const Poly::PolyParams& params);
	/**
	 *	\brief	Forms and solves 6th degree polynomial.
	 *	\param	params	Polynomial coefficients params.
	 *	\return Six roots.
	 */
	Roots Solve(const PolyParams& params);
	/**
	 *	\brief	Evaluates cost function for each root.
	 *	\param	roots	Six roots of 6th degree polynomial.
	 *	\param	params	Polynomial coefficients params.
	 *	\return	Array of cost for each root.
	 */
	std::vector<double> EvaluateRootsCosts(const Roots& roots, const PolyParams& params);
	/**
	 *	\brief	Evaluates cost function for roots.
	 *	\param	roots	Six roots of 6th degree polynomial.
	 *	\param	params	Polynomial coefficients params.
	 *	\return	Root that gives smallest cost function value.
	 */
	double EvaluateRoots(const Roots& roots, const PolyParams& params);
	/**
	 *	\brief	Construct two epipolar lines on which the corrected correspondences lie.
	 *	\param	t	Real root for which the cost function gives the smallest value.
	 *	\param	params	Polynomial coefficients params.
	 *	\return	Pair of epipolar lines.
	 */
	std::pair<Line, Line> ConstructLines(double t, const PolyParams& params);
	/**
	 *	\brief	Finds point on given line that is closest to the origin.
	 *	\param	l	Line on which the point lies.
	 *	\return	Point on line closest to the origin.
	 */
	cv::Point3d FindPointOnLineClosestToOrigin(const Line& l);
	/**
	 *	\brief	Transfers point to original coordinates using R and T.
	 *	\param	p	A point to transfer.
	 *	\param	R	Rotational transformation.
	 *	\param	T	Translational transformation.
	 *	\return 	Point in original coordinates.
	 */
	cv::Point3d TransferPointToOriginalCoordinates(const cv::Point3d& p, const cv::Mat& R, const cv::Mat& T);
	/**
	 *	\brief	Sets origin of world coordinate system to first camera.
	 *	\param	P0	First camera projection matrix.
	 *	\param	P1	Second camera projection matrix.
	 *	\return	K0, K1, R and T - camera intrinsics and orientation of second camera in new world coordinates.
	 */
	std::tuple<Intrinsic, Intrinsic, cv::Mat, cv::Mat> SetOriginToCamera(const cv::Mat& P0, const cv::Mat& P1);
	/**
	 *	\brief	Computes fundamental matrix from camera projection matrices.
	 *	\param	P0	First camera projection matrix.
	 *	\param	P1	Second camera projection matrix.
	 *	\return	Computed fundamental matrix.
	 */
	Fundamental ComputeFundamentalMatrix(const cv::Mat& P0, const cv::Mat& P1);
	/**
	 *	\brief	Triangulate point using The Direct Linear Transformation Method.
	 *	\param	p0	Point in the image of the first camera.
	 *	\param	p1	Corresponding point in the image of the second camera.
	 *	\return	Triangulated point.
	 */
	cv::Point3d TriangulateDLT(const cv::Point2d& p0, const cv::Point2d& p1);
	/**
	 *	\brief	Returns the order of the polynomial with given coefficients (highest non-zero coeffs index).
	 *	\param	coeffs	Polynomial coefficients.
	 *	\return	Polynomial order
	 */
	int FindPolynomialOrder(const std::vector<double>& coeffs);
	/**
	 *	\brief	Returns canonic camera projection matrix of second camera computed from given fundamental matrix.
	 *	\param	F	Fundamental matrix.
	 *	\return	Canonic	camera projection matrix.
	 */
	cv::Mat CameraProjectionMatrixFromFundamentalMatrix(const Fundamental& F);

	const cv::Mat P0;
	const cv::Mat P1;
	const Fundamental F;
};

#endif /* POLY_H_ */
