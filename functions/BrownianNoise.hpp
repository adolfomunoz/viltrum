//----------------------------------------------------------------------------------------
// 	Brownian Noise 
//----------------------------------------------------------------------------------------

# pragma once


namespace Noise
{
	class BrownianNoise
	{
	private:
		unsigned int m_N;
		double m_beta;
		std::uint32_t m_seed;

		// Deterministic random parameter
		// Based on Morgan McGuire
		double random(double x, double y) const
		{
			double val = sin(x*(12.9898+m_seed)+y*(78.233+m_seed))*43758.5453123;
			//printf("Random [@BrownianNoise] (%f,%f): %f\n", x, y, val);
			return val - (int)val;
		}

		double randnoise(double x) const
		{
			double ix = (int)x;
			double fx = x-ix;

			double a = random(ix, 0.);
			double b = random(ix + 1.f, 0.);

			double u = fx*fx*(3.-2.*fx);

			return a*(u) + b*(1-u);
		}

	public:
		explicit BrownianNoise(double beta, unsigned int N, std::uint32_t seed = std::default_random_engine::default_seed):
			m_N(N), m_beta(beta), m_seed(seed)
		{	}

		explicit BrownianNoise(std::uint32_t seed = std::default_random_engine::default_seed):
			m_N(8), m_beta(2.f), m_seed(seed)
		{	}

		double noise(double x) const
		{
			double value = 0.;

			for(unsigned int i=0; i<m_N; ++i)
			{
				double amplitude = powf(2*(i+1.), -m_beta/2);
				value += amplitude * randnoise(x);
				x *= 2;
			}
			return value;
		}
	};
}
