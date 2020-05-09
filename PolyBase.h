#ifndef POLYBASE_H_
#define POLYBASE_H_

#include <complex>
#include <opencv2/opencv.hpp>
#include "TriangulationBase.h"
#include "LinearLS.h"

namespace Triangulation {

class PolyBase : public TriangulationBase
{
public:
	typedef cv::Matx33d Intrinsic;
	typedef cv::Matx33d Fundamental;

	/**
	 *	\brief	Constructor
	 *	\param	F	Fundamental matrix.
	 */
	explicit PolyBase(const Fundamental& F);
	/**
	 *	\brief	Constructor
	 *	\param	P0	Camera matrix of the first camera.
	 *	\param 	P1	Camera matrix of the second camera.
	 */
	PolyBase(const cv::Mat& P0, const cv::Mat& P1);
	/**
	 *	\brief	Constructor
	 *	\param	P0	Camera matrix of the first camera.
	 *	\param 	P1	Camera matrix of the second camera.
	 *	\param	F	Fundamental matrix.
	 */
	PolyBase(const cv::Mat& P0, const cv::Mat& P1, const Fundamental& F);
	/**
	 *	\brief	Triangulates image points.
	 *	\param	p0	Point in the image of the first camera.
	 *	\param	p1	Corresponding point in the image of the second camera.
	 *	\return	Triangulated point.
	 */
	cv::Point3d triangulate(const cv::Point2d& p0, const cv::Point2d& p1) const override;
protected:
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
	std::pair<cv::Point2d, cv::Point2d> ComputeCorrectedCorrespondences(const cv::Point2d& p0, const cv::Point2d& p1) const;
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
	Epipole ComputeLeftEpipole(const cv::Mat& F) const;
	/**
	 *	\brief	Computes epipole e = (e1, e2, e3) such as Fe = 0 and e1*e1 + e2*e2 = 1
	 *	\param	F	Fundamental matrix.
	 *	\return	Computed epipole
	 */
	Epipole ComputeRightEpipole(const cv::Mat& F) const;
	/**
	 *	\brief	Defines rotation matrix using given epipole.
	 *	\param	e	Epipole.
	 *	\return	Rotation matrix.
	 */
	cv::Mat FormRotationMatrix(const Epipole& e) const;
	/**
	 *	\brief	Prepares polynomial coefficients.
	 *	\param	params	Polynomial coefficients params.
	 *	\return Polynomial coefficients.
	 */
	virtual std::vector<double> PreparePolyCoeffs(const PolyParams& params) const = 0;
	/**
	 *	\brief	Forms and solves 6th degree polynomial.
	 *	\param	params	Polynomial coefficients params.
	 *	\return Six roots.
	 */
	Roots Solve(const PolyParams& params) const;
	/**
	 *	\brief	Evaluates cost function for each root.
	 *	\param	roots	Six roots of 6th degree polynomial.
	 *	\param	params	Polynomial coefficients params.
	 *	\return	Array of cost for each root.
	 */
	virtual std::vector<double> EvaluateRootsCosts(const Roots& roots, const PolyParams& params) const = 0;
	/**
	 *	\brief	Evaluates cost function for roots.
	 *	\param	roots	Six roots of 6th degree polynomial.
	 *	\param	params	Polynomial coefficients params.
	 *	\return	Root that gives smallest cost function value.
	 */
	double EvaluateRoots(const Roots& roots, const PolyParams& params) const;
	/**
	 *	\brief	Construct two epipolar lines on which the corrected correspondences lie.
	 *	\param	t	Real root for which the cost function gives the smallest value.
	 *	\param	params	Polynomial coefficients params.
	 *	\return	Pair of epipolar lines.
	 */
	std::pair<Line, Line> ConstructLines(double t, const PolyParams& params) const;
	/**
	 *	\brief	Finds point on given line that is closest to the origin.
	 *	\param	l	Line on which the point lies.
	 *	\return	Point on line closest to the origin.
	 */
	cv::Point3d FindPointOnLineClosestToOrigin(const Line& l) const;
	/**
	 *	\brief	Transfers point to original coordinates using R and T.
	 *	\param	p	A point to transfer.
	 *	\param	R	Rotational transformation.
	 *	\param	T	Translational transformation.
	 *	\return 	Point in original coordinates.
	 */
	cv::Point3d TransferPointToOriginalCoordinates(const cv::Point3d& p, const cv::Mat& R, const cv::Mat& T) const;
	/**
	 *	\brief	Sets origin of world coordinate system to first camera.
	 *	\param	P0	First camera projection matrix.
	 *	\param	P1	Second camera projection matrix.
	 *	\return	K0, K1, R and T - camera intrinsics and orientation of second camera in new world coordinates.
	 */
	std::tuple<Intrinsic, Intrinsic, cv::Mat, cv::Mat> SetOriginToCamera(const cv::Mat& P0, const cv::Mat& P1) const;
	/**
	 *	\brief	Computes fundamental matrix from camera projection matrices.
	 *	\param	P0	First camera projection matrix.
	 *	\param	P1	Second camera projection matrix.
	 *	\return	Computed fundamental matrix.
	 */
	Fundamental ComputeFundamentalMatrix(const cv::Mat& P0, const cv::Mat& P1) const;
	/**
	 *	\brief	Returns the order of the polynomial with given coefficients (highest non-zero coeffs index).
	 *	\param	coeffs	Polynomial coefficients.
	 *	\return	Polynomial order
	 */
	int FindPolynomialOrder(const std::vector<double>& coeffs) const;
	/**
	 *	\brief	Returns canonic camera projection matrix of second camera computed from given fundamental matrix.
	 *	\param	F	Fundamental matrix.
	 *	\return	Canonic	camera projection matrix.
	 */
	cv::Mat CameraProjectionMatrixFromFundamentalMatrix(const Fundamental& F) const;

	const Fundamental F;
	/// LinearLS method used for final triangulation of corrected points.
	const LinearLS LS;
};

}

#endif /* POLYBASE_H_ */
