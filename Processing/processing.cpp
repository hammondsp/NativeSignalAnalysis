#include "processing.h"
#include <cmath>
#include <complex>

using namespace std;

namespace Signals
{
	namespace Processing
	{
		void Convolve(const vector<double>& left, const vector<double>& right, vector<double>& result)
		{
			// (left * right)[n] = sum {left[m]right[n - m]}
			for (int n = 0; n < right.size(); ++n)
			{
				for (int m = 0; m <= n && m < left.size(); ++m)
				{
					result[n] += left[m] * right[n - m];
				}
			}
		}

		void DFT(const vector<double>& signal, vector<complex<double>>& result)
		{
			if (signal.size() != result.size())
			{
				result.resize(signal.size());
			}

			static const double PI = 2.0 * acos(0.0);
			
			for (int idx = 0; idx < signal.size(); ++idx)
			{
				result[idx] = complex<double>(0., 0.);
				const double omegaBase = 2. * PI * (double)idx / (double)signal.size();
				for (int walk = 0; walk < result.size(); ++walk)
				{
					double omega = omegaBase * (double)walk;
					result[idx].real(result[idx].real() + sin(omega) * signal[idx]);
					result[idx].imag(result[idx].imag() + cos(omega) * signal[idx]);
				}
			}
		}

		void IDFT(const vector<complex<double>>& signal, vector<double>& result)
		{
			if (signal.size() != result.size())
			{
				result.resize(signal.size());
			}

			static const double PI = 2.0 * acos(0.0);

			for (int idx = 0; idx < signal.size(); ++idx)
			{
				result[idx] = 0.;
				double omegaBase = 2. * PI * (double)idx / (double)signal.size();
				for (int walk = 0; walk < result.size(); ++walk)
				{
					double omega = omegaBase * (double)walk;
					result[idx] += cos(omega) * signal[walk].real() - sin(omega) * signal[walk].imag();
				}

				result[idx] /= (double)signal.size();
			}
		}
	}
}
