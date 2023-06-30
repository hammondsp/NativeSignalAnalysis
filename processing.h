#pragma once

#include <vector>

namespace Signals
{
	namespace Processing
	{
		void Convolve(const std::vector<double>& left, const std::vector<double>& right, std::vector<double>& result);
	}
}