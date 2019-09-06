#ifndef ITERATIVELINEARLS_H_
#define ITERATIVELINEARLS_H_

#include "TriangulationBase.h"

namespace Triangulation {

class IterativeLS : public TriangulationBase
{
public:
	using TriangulationBase::TriangulationBase;
	using TriangulationBase::triangulate;
	/**
	 *	\brief	Triangulate point iteratively using The Direct Linear Transformation Method.
	 *	\param	p0	Point in the image of the first camera.
	 *	\param	p1	Corresponding point in the image of the second camera.
	 *	\return	Triangulated point.
	 */
	cv::Point3d triangulate(const cv::Point2d& p0, const cv::Point2d& p1) const override;
};

}

#endif /* ITERATIVELINEARLS_H_ */
