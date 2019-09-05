#ifndef POLY_H_
#define POLY_H_

#include "PolyBase.h"

namespace Triangulation {

/**
 *	\brief	Performs non-linear triangulation of two image points by minimizing geometric error.
 */
class Poly : public PolyBase
{
public:
	using PolyBase::PolyBase;
	using PolyBase::TriangulationBase::triangulate;
private:
	std::vector<double> PreparePolyCoeffs(const PolyParams& params) const override;
	std::vector<double> EvaluateRootsCosts(const Roots& roots, const PolyParams& params) const override;
};

}

#endif /* POLY_H_ */
