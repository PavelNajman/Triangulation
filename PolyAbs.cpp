#include "PolyAbs.h"

namespace Triangulation {

std::vector<double> PolyAbs::PreparePolyCoeffs(const PolyAbs::PolyParams& params) const
{
	double a, b, c, d, e, f;
	std::tie(a, b, c, d, e, f) = params;
	std::vector<double> result = {
		// 2*a*b^3*c*d + 3*b^2*d^4*f^4 + 3*b^4*d^2*f^2 + d^6*f^6 - a^2*b^2*d^2 - b^4*c^2 + b^6
		+ 2*a*b*b*b*c*d + 3*b*b*d*d*d*d*f*f*f*f + 3*b*b*b*b*d*d*f*f + d*d*d*d*d*d*f*f*f*f*f*f - a*a*b*b*d*d - b*b*b*b*c*c + b*b*b*b*b*b,
		// 12*b^2*c*d^3*f^4*z + 12*a*b^3*d^2*f^2*z + 6*b^4*c*d*f^2*z + 6*a*b*d^4*f^4*z + 4*a^2*b^2*c*d*z + 6*c*d^5*f^6*z - 2*a^3*b*d^2*z - 2*a*b^3*c^2*z + 6*a*b^5*z
		+ 12*b*b*c*d*d*d*f*f*f*f + 12*a*b*b*b*d*d*f*f + 6*b*b*b*b*c*d*f*f + 6*a*b*d*d*d*d*f*f*f*f + 4*a*a*b*b*c*d + 6*c*d*d*d*d*d*f*f*f*f*f*f - 2*a*a*a*b*d*d - 2*a*b*b*b*c*c + 6*a*b*b*b*b*b,
		// 24*a*b*c*d^3*f^4*z^2 + 24*a*b^3*c*d*f^2*z^2 + 6*a*b^3*c*d*e^2*z^2 + 18*b^2*c^2*d^2*f^4*z^2 + 18*a^2*b^2*d^2*f^2*z^2 - 3*a^2*b^2*d^2*e^2*z^2 + 2*a^3*b*c*d*z^2 + 15*c^2*d^4*f^6*z^2 + 3*a^2*d^4*f^4*z^2 + 3*b^4*c^2*f^2*z^2 - 3*b^4*c^2*e^2*z^2 + 15*a^2*b^4*z^2 - a^2*b^2*c^2*z^2 - a^4*d^2*z^2
		+ 24*a*b*c*d*d*d*f*f*f*f + 24*a*b*b*b*c*d*f*f + 6*a*b*b*b*c*d*e*e + 18*b*b*c*c*d*d*f*f*f*f + 18*a*a*b*b*d*d*f*f - 3*a*a*b*b*d*d*e*e + 2*a*a*a*b*c*d + 15*c*c*d*d*d*d*f*f*f*f*f*f + 3*a*a*d*d*d*d*f*f*f*f + 3*b*b*b*b*c*c*f*f - 3*b*b*b*b*c*c*e*e + 15*a*a*b*b*b*b - a*a*b*b*c*c - a*a*a*a*d*d,
		// 36*a*b*c^2*d^2*f^4*z^3 + 36*a^2*b^2*c*d*f^2*z^3 + 12*a^2*b^2*c*d*e^2*z^3 + 12*b^2*c^3*d*f^4*z^3 + 12*a^2*c*d^3*f^4*z^3 + 12*a^3*b*d^2*f^2*z^3 + 12*a*b^3*c^2*f^2*z^3 - 6*a^3*b*d^2*e^2*z^3 - 6*a*b^3*c^2*e^2*z^3 + 20*c^3*d^3*f^6*z^3 + 20*a^3*b^3*z^3
		+ 36*a*b*c*c*d*d*f*f*f*f + 36*a*a*b*b*c*d*f*f + 12*a*a*b*b*c*d*e*e + 12*b*b*c*c*c*d*f*f*f*f + 12*a*a*c*d*d*d*f*f*f*f + 12*a*a*a*b*d*d*f*f + 12*a*b*b*b*c*c*f*f - 6*a*a*a*b*d*d*e*e - 6*a*b*b*b*c*c*e*e + 20*c*c*c*d*d*d*f*f*f*f*f*f + 20*a*a*a*b*b*b,
		// 24*a*b*c^3*d*f^4*z^4 + 24*a^3*b*c*d*f^2*z^4 + 6*a*b^3*c*d*e^4*z^4 + 6*a^3*b*c*d*e^2*z^4 + 18*a^2*c^2*d^2*f^4*z^4 + 18*a^2*b^2*c^2*f^2*z^4 - 3*a^2*b^2*d^2*e^4*z^4 - 3*a^2*b^2*c^2*e^2*z^4 + 15*c^4*d^2*f^6*z^4 + 3*b^2*c^4*f^4*z^4 + 3*a^4*d^2*f^2*z^4 - 3*b^4*c^2*e^4*z^4 - 3*a^4*d^2*e^2*z^4 + 15*a^4*b^2*z^4
		+ 24*a*b*c*c*c*d*f*f*f*f + 24*a*a*a*b*c*d*f*f + 6*a*b*b*b*c*d*e*e*e*e + 6*a*a*a*b*c*d*e*e + 18*a*a*c*c*d*d*f*f*f*f + 18*a*a*b*b*c*c*f*f - 3*a*a*b*b*d*d*e*e*e*e - 3*a*a*b*b*c*c*e*e + 15*c*c*c*c*d*d*f*f*f*f*f*f + 3*b*b*c*c*c*c*f*f*f*f + 3*a*a*a*a*d*d*f*f - 3*b*b*b*b*c*c*e*e*e*e - 3*a*a*a*a*d*d*e*e + 15*a*a*a*a*b*b,
		// 12*a^2*b^2*c*d*e^4*z^5 + 12*a^2*c^3*d*f^4*z^5 + 12*a^3*b*c^2*f^2*z^5 - 6*a^3*b*d^2*e^4*z^5 - 6*a*b^3*c^2*e^4*z^5 + 6*a^4*c*d*f^2*z^5 + 6*a*b*c^4*f^4*z^5 + 6*c^5*d*f^6*z^5 + 6*a^5*b*z^5
		+ 12*a*a*b*b*c*d*e*e*e*e + 12*a*a*c*c*c*d*f*f*f*f + 12*a*a*a*b*c*c*f*f - 6*a*a*a*b*d*d*e*e*e*e - 6*a*b*b*b*c*c*e*e*e*e + 6*a*a*a*a*c*d*f*f + 6*a*b*c*c*c*c*f*f*f*f + 6*c*c*c*c*c*d*f*f*f*f*f*f + 6*a*a*a*a*a*b,
		// 6*a^3*b*c*d*e^4*z^6 + 2*a*b^3*c*d*e^6*z^6 - 3*a^2*b^2*c^2*e^4*z^6 - 3*a^4*d^2*e^4*z^6 + 3*a^2*c^4*f^4*z^6 + 3*a^4*c^2*f^2*z^6 - a^2*b^2*d^2*e^6*z^6 + c^6*f^6*z^6 - b^4*c^2*e^6*z^6 + a^6*z^6
		+ 6*a*a*a*b*c*d*e*e*e*e + 2*a*b*b*b*c*d*e*e*e*e*e*e - 3*a*a*b*b*c*c*e*e*e*e - 3*a*a*a*a*d*d*e*e*e*e + 3*a*a*c*c*c*c*f*f*f*f + 3*a*a*a*a*c*c*f*f - a*a*b*b*d*d*e*e*e*e*e*e + c*c*c*c*c*c*f*f*f*f*f*f - b*b*b*b*c*c*e*e*e*e*e*e + a*a*a*a*a*a,
		// 4*a^2*b^2*c*d*e^6*z^7 - 2*a^3*b*d^2*e^6*z^7 - 2*a*b^3*c^2*e^6*z^7
		+ 4*a*a*b*b*c*d*e*e*e*e*e*e - 2*a*a*a*b*d*d*e*e*e*e*e*e - 2*a*b*b*b*c*c*e*e*e*e*e*e,
		// 2*a^3*b*c*d*e^6*z^8 - a^2*b^2*c^2*e^6*z^8 - a^4*d^2*e^6*z^8
		2*a*a*a*b*c*d*e*e*e*e*e*e - a*a*b*b*c*c*e*e*e*e*e*e - a*a*a*a*d*d*e*e*e*e*e*e
	};
	result.resize(FindPolynomialOrder(result) + 1);
	double max_coeff = *std::max_element(result.begin(), result.end());
	if (max_coeff > 0)
	{
		std::transform(result.begin(), result.end(), result.begin(), [&max_coeff](double& a){ return a / max_coeff; });
	}
	else if (max_coeff < 0)
	{
		double min_coeff = *std::min_element(result.begin(), result.end());
		std::transform(result.begin(), result.end(), result.begin(), [&min_coeff](double& a){ return a / min_coeff; });
	}
	return result;
}

std::vector<double> PolyAbs::EvaluateRootsCosts(const PolyAbs::Roots& roots, const PolyAbs::PolyParams& params) const
{
	double a, b, c, d, e, f;
	std::tie(a, b, c, d, e, f) = params;
	auto cost_function = [&](double t)
	{
		return (std::abs(t) / std::sqrt(1 + e*e * t*t)) + (std::abs(c*t+d) / std::sqrt((a*t+b) * (a*t+b) + f*f * (c*t+d) * (c*t+d)));
	};

	std::vector<double> result(roots.size());
	for (size_t i = 0u; i < roots.size(); ++i)
	{
		result[i] = cost_function(roots[i].real());
	}
	return result;
}

}
