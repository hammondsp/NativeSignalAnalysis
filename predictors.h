#pragma once

#include <vector>

namespace Signals
{
	namespace Predictors
	{
		// Predict the next value using 1st order dead reckoning
		double DeadReckoning1(const std::vector<double>& input);

		// Predict the next value using 2nd order dead reckoning
		double DeadReckoning2(const std::vector<double>& input);

		// Predict the next N frames usin linear prediction coding
		std::vector<double> LinearPredictioCoding(std::vector<double> input, int numPoles, int frames);
	}
}