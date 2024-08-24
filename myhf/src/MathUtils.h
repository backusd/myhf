#pragma once


namespace myhf
{
static constexpr unsigned long long int factorial(long int n) noexcept
{
	return n ? n * factorial(n - 1) : 1;
}

static constexpr unsigned long long int double_factorial(long int n) noexcept
{
	return (n > 1) ? n * double_factorial(n - 2) : 1;
}

class MathUtils
{
private:
	static constexpr std::array<unsigned long long int, 21> factorialsTable
	{
		factorial(0),
		factorial(1),
		factorial(2),
		factorial(3),
		factorial(4),
		factorial(5),
		factorial(6),
		factorial(7),
		factorial(8),
		factorial(9),
		factorial(10),
		factorial(11),
		factorial(12),
		factorial(13),
		factorial(14),
		factorial(15),
		factorial(16),
		factorial(17),
		factorial(18),
		factorial(19),
		factorial(20)
	};
    static constexpr std::array<unsigned long long int, 21> doubleFactorialsTable
	{
		double_factorial(0),
		double_factorial(1),
		double_factorial(2),
		double_factorial(3),
		double_factorial(4),
		double_factorial(5),
		double_factorial(6),
		double_factorial(7),
		double_factorial(8),
		double_factorial(9),
		double_factorial(10),
		double_factorial(11),
		double_factorial(12),
		double_factorial(13),
		double_factorial(14),
		double_factorial(15),
		double_factorial(16),
		double_factorial(17),
		double_factorial(18),
		double_factorial(19),
		double_factorial(20)
	};

public:
	static constexpr double Factorial(long int n) noexcept
	{
		// using precalculated (at compile time) values
		if (n < 0) return 1;
		else if (n < 21) return static_cast<double>(factorialsTable[n]);

		double val = static_cast<double>(factorialsTable[20]);
		for (long int i = 21; i <= n; ++i)
			val *= i;

		return val;
	}

	static constexpr double DoubleFactorial(long int n) noexcept
	{
		if (n < 0) return 1;
		else if (n < 21) return static_cast<double>(doubleFactorialsTable[n]);

		double res = 1;
		for (; n > 20; n -= 2)
			res *= n;

		return res * static_cast<double>(doubleFactorialsTable[n]);
	}

	static constexpr unsigned int BinomialCoefficient(unsigned int n, unsigned int k) noexcept
	{
		//return Factorial(n) / (Factorial(k) * Factorial(n-k));

		unsigned int numerator = 1;
		unsigned int denominator = 1;

		for (unsigned int i = 1; i <= k; ++i)
		{
			numerator *= n - i + 1;
			denominator *= i;
		}

		return numerator / denominator;
	}


private:
	static double Abscissa(int n, int i) noexcept
	{
		double val = i * std::numbers::pi / (1. + n);
		double sinVal = std::sin(val);
		const double cosVal = std::cos(val);
	
		return (n + 1. - 2. * i) / (n + 1.) +
			2. / std::numbers::pi * (1. + 2. / 3. * sinVal * sinVal) * cosVal * sinVal;
	}


	static double Omega(int n, int i)
	{
		const double sinVal = std::sin(i * std::numbers::pi / (n + 1.));
		return 16. / (3 * (1. + n)) * sinVal * sinVal * sinVal * sinVal;
	}

public:
	class FunctionFunctor {
	public:
		virtual double operator()(double x) const = 0;
	};


	// this is based on the code from Mathematica Journal article about solving nuclear integrals
	// implemented there in... Wolfram, of course
	// and also on 'A simple, efficient and more reliable scheme for automatic numerical integration'
	// Computer Physics Communications � September 1993, DOI: 10.1016/0010-4655(93)90035-B
	// there is a Fortran implementation in there
	// it is not used by the current code but I leave it here just in case

	static double ChebyshevIntegral(double eps, int M, const FunctionFunctor& F)
	{
		double c0 = std::cos(std::numbers::pi / 6.);
		double s0 = 0.5;
		double c1, s1, q, p, chp, c, s, x;

		double err = 10.;

		int n = 3;
		c1 = s0;
		s1 = c0;

		const double res = Abscissa(2, 1);

		q = (F(res) + F(-res)) * Omega(2, 1);
		p = F(0);

		chp = q + p;

		int j = 0;

		while (err > eps && 2. * n * (1. - j) + j * 4. * n / 3. - 1. <= M)
		{
			j = 1 - j;

			c1 = j * c1 + (1. - j) * c0;
			s1 = j * s1 + (1. - j) * s0;
			c0 = j * c0 + (1. - j) * std::sqrt((1. + c0) / 2.);
			s0 = j * s0 + (1. - j) * s0 / (c0 + c0);

			c = c0;
			s = s0;

			for (int i = 1; i < n; i += 2)
			{
				x = 1. + 2. / (3. * std::numbers::pi) * s * c * (3. + 2. * s * s) - ((double)i) / n;

				const int div = (int)((i + 2. * j) / 3.);
				if (3 * div >= i + j)
					chp += (F(-x) + F(x)) * s * s * s * s;

				x = s;
				s = s * c1 + c * s1;
				c = c * c1 - x * s1;
			}

			n *= 1 + j;
			p += (1. - j) * (chp - q);

			err = 16. * std::abs((1. - j) * (q - 3. * p / 2.) + j * (chp - q - q)) / (3. * n);
			q = (1. - j) * q + j * chp;
		}

		chp = 16. * q / (3. * n);

		return chp;
	}


	// taken from here: http://www.crbond.com/math.htm and adapted.
	// Originally from Zhang and Jin, 'Computation of special functions'
	// faster codes use interpolation on some domain, but I won't bother, at least not yet

	// As an important note, Boys functions are the start for integrals calculations, so if this has numerical errors,
	// they propagate and amplify, so for a serious program more care should be given

	static bool IncompleteGamma(double a, double x, double& gin)
	{
		double xam, r, sum, ga, t0;
		int k;

		gin = 0;

		if (a < 0.0 || x < 0) return false;
		else if (x == 0.0) return true;

		xam = -x + a * std::log(x);

		if (xam > 600 || a > 160.0) return false;
		else if (x <= 1.0 + a) {
			sum = 1.0 / a;
			r = sum;

			for (k = 1; k <= 80; ++k) {
				r *= x / (a + k);
				sum += r;
				if (std::fabs(r / sum) <= std::numeric_limits<double>::epsilon()) break;
			}

			gin = std::exp(xam) * sum;
		}
		else {
			t0 = 0.0;

			for (k = 80; k >= 1; --k)
				t0 = (k - a) / (1.0 + k / (x + t0));

			const double gim = std::exp(xam) / (x + t0);
			ga = std::tgamma(a);
			gin = ga - gim;
		}

		return true;
	}
};
}