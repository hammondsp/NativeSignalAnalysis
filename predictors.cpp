#include "predictors.h"

using namespace std;

namespace
{
    void CalculatePoles(vector<double>& coeffs, const vector<double>& signal)
    {
        // The size of the input vectors are used to determine how many poles will be calculated
        size_t signalSize = signal.size() - 1;
        size_t coeffCount = coeffs.size();

        // Autocorrelation calculation. This is a O(m * n), m = signal size and n = requested poles.
        vector<double> autoCorr(coeffCount + 1, 0.);
        {
            for (size_t i = 0; i <= coeffCount; i++)
            {
                for (size_t j = 0; j <= signalSize - i; j++)
                {
                    autoCorr[i] += signal[j] * signal[j + i];
                }
            }
        }

        // This is where the coefficients will be calculated
        vector<double> Ak(coeffCount + 1, 0.);
        Ak[0] = 1.;

        // Error tracking
        double Ek = autoCorr[0];

        // Levinson-Durbin recursion
        for (size_t k = 0; k < coeffCount; k++)
        {
            // Calculate the scaling factor for this recursion
            double lambda = 0.;
            for (size_t j = 0; j <= k; j++)
            {
                lambda -= Ak[j] * autoCorr[k + 1 - j];
            }
            lambda /= Ek;

            // Update the coefficients
            for (size_t n = 0; n <= (k + 1) / 2; n++)
            {
                double temp = Ak[k + 1 - n] + lambda * Ak[n];
                Ak[n] = Ak[n] + lambda * Ak[k + 1 - n];
                Ak[k + 1 - n] = temp;
            }

            // Update the error calculation
            Ek *= 1.0 - lambda * lambda;
        }

        // Emit the coefficients ignoring the first
        coeffs.assign(++Ak.begin(), Ak.end());
    }

}

namespace Signals
{
	namespace Predictors
	{
		// Predict the next value using 1st order dead reckoning
		double DeadReckoning1(const vector<double>& input)
		{
			const double tMinus1 = input[input.size() - 1];
			const double tMinus2 = input[input.size() - 2];

			const double delta = tMinus1 - tMinus2;

			return tMinus1 + delta;
		}


		// Predict the next value using 2nd order dead reckoning
		double DeadReckoning2(const vector<double>& input)
		{
			const double tMinus1 = input[input.size() - 1];
			const double tMinus2 = input[input.size() - 2];
			const double tMinus3 = input[input.size() - 3];

			const double vel1 = tMinus1 - tMinus2;
			const double vel2 = tMinus2 - tMinus3;

			const double accel = vel1 - vel2;

			return tMinus1 + vel1 + (accel / 2.0);
		}


		// Predict the next N frames usin linear prediction coding
		std::vector<double> LinearPredictioCoding(vector<double> input, int numPoles, int frames)
        {
            vector<double> poles(numPoles);

            // If there hasn't been movement in the last few frames, predict no movement going forward
            if (abs(input[input.size() - 1] - input[input.size() - 2]) < 0.01 &&
                abs(input[input.size() - 2] - input[input.size() - 3]) < 0.01)
            {
                return vector<double>(frames, input[input.size() - 1]);
            }

            // calculate the coefficients used to predict frames ahead
            vector<double> result(frames, 0.);
            CalculatePoles(poles, input);

            // predict the requested frames
            for (int frame = 0; frame < frames; ++frame)
            {
                int inputIdx = static_cast<int>(input.size() - 1);
                int poleIdx = 0;

                // for each coefficient, apply that to the appropriate previous frame
                while (poleIdx < poles.size() && inputIdx >= 0)
                {
                    result[frame] -= input[inputIdx--] * poles[poleIdx++];
                }

                // store the predicted frame to use for predicting future frames
                input.push_back(result[frame]);
            }

            return result;
        }
	}
}