#ifndef TRIANGULATIONBASE_H_
#define TRIANGULATIONBASE_H_

#include <opencv2/opencv.hpp>

namespace Triangulation {
/**
 * 	\class	TriangulationBase
 *	\brief	
 */
class TriangulationBase
{
public:
	/**
	 *	\brief	Constructor
	 *	\param	P0	Camera matrix of the first camera.
	 *	\param 	P1	Camera matrix of the second camera.
	 */
	TriangulationBase(const cv::Mat& P0, const cv::Mat& P1);
	virtual ~TriangulationBase() = default;
	/**
	 *	\brief	Triangulates image points.
	 *	\param	p0	Point in the image of the first camera.
	 *	\param	p1	Corresponding point in the image of the second camera.
	 *	\return	Triangulated point.
	 */
	virtual cv::Point3d triangulate(const cv::Point2d& p0, const cv::Point2d& p1) const = 0;
	/**
	 *	\brief	Triangulate points using The Direct Linear Transformation Method.
	 *	\param	p0	Points in the image of the first camera.
	 *	\param	p1	Corresponding points in the image of the second camera.
	 *	\return	Triangulated points.
	 */
	std::vector<cv::Point3d> triangulate(const std::vector<cv::Point2d>& p0, const std::vector<cv::Point2d>& p1) const;
protected:
	const cv::Mat P0, P1;
};

}	// Triangulation

#endif /* TRIANGULATIONBASE_H_ */
