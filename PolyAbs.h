#ifndef POLYABS_H_
#define POLYABS_H_

#include "PolyBase.h"

namespace Triangulation {

/**
 *	\brief	Performs non-linear triangulation of two image points by minimizing the sum of absolute values of the distances.
 */
class PolyAbs : public PolyBase
{
public:
	using PolyBase::PolyBase;
	using PolyBase::TriangulationBase::triangulate;
private:
	std::vector<double> PreparePolyCoeffs(const PolyParams& params) const override;
	std::vector<double> EvaluateRootsCosts(const Roots& roots, const PolyParams& params) const override;
};

}

#endif /* POLYABS_H_ */
