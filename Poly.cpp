#include "Poly.h"

namespace Triangulation {

std::vector<double> Poly::PreparePolyCoeffs(const Poly::PolyParams& params) const
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

std::vector<double> Poly::EvaluateRootsCosts(const Poly::Roots& roots, const Poly::PolyParams& params) const
{
	double a, b, c, d, e, f;
	std::tie(a, b, c, d, e, f) = params;
	auto cost_function = [&](double t)
	{
		return ((t*t) / (1 + e*e * t*t)) + ((c*t+d) * (c*t+d)) / ((a*t+b) * (a*t+b) + f*f * (c*t+d) * (c*t+d));
	};

	std::vector<double> result(roots.size());
	for (size_t i = 0u; i < roots.size(); ++i)
	{
		result[i] = cost_function(roots[i].real());
	}
	return result;
}

}
