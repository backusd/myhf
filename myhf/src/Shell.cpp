#include "pch.h"
#include "Shell.h"


namespace myhf
{
	std::string ContractedGaussianShell::GetShellString() const noexcept
	{
		assert(basisFunctions.size());

		std::string res;

		for (auto& orbital : basisFunctions)
		{
			char c = static_cast<char>(toupper(orbital.AtomicOrbital())); 

			if (std::string::npos == res.find(c))
			{
				res += c;
			}
		}

		return res;
	}

	unsigned int ContractedGaussianShell::AdjustOrbitalsCount(char orbital, unsigned int res) noexcept
	{
		switch (tolower(orbital))
		{
		default:
		case 's':
			break;
		case 'p':
			res /= QuantumNumbers::QuantumNumbers::NumOrbitals(1);
			break;
		case 'd':
			res /= QuantumNumbers::QuantumNumbers::NumOrbitals(2);
			break;
		case 'f':
			res /= QuantumNumbers::QuantumNumbers::NumOrbitals(3);
			break;
		case 'g':
			res /= QuantumNumbers::QuantumNumbers::NumOrbitals(4);
			break;
		case 'h':
			res /= QuantumNumbers::QuantumNumbers::NumOrbitals(5);
			break;
		}

		return res;
	}

	unsigned int ContractedGaussianShell::CountOrbitals(char orbitalChar) const noexcept
	{
		unsigned int res = 0;

		for (const auto& contractedOrbital : basisFunctions)
		{
			if (orbitalChar == contractedOrbital.AtomicOrbital())
				res += static_cast<unsigned int>(contractedOrbital.gaussianOrbitals.size());

//			for (const auto& gaussian : contractedOrbital.gaussianOrbitals)
//				if (orbitalChar == gaussian.AtomicOrbital())
//					++res;
		}

		return AdjustOrbitalsCount(orbitalChar, res);
	}

	unsigned int ContractedGaussianShell::CountContractedOrbitals(char orbitalChar) const noexcept
	{
		unsigned int res = 0;

		for (auto const& contractedOrbital : basisFunctions)
		{
			if (orbitalChar == contractedOrbital.AtomicOrbital())
				++res;
		}

		return AdjustOrbitalsCount(orbitalChar, res);
	}

	unsigned int ContractedGaussianShell::CountNumberOfContractedGaussians() const noexcept
	{
		return static_cast<unsigned int>(basisFunctions.size());
	}

	unsigned int ContractedGaussianShell::CountNumberOfGaussians() const noexcept
	{
		unsigned int res = 0;

		for (auto const& contractedOrbital : basisFunctions)
			res += static_cast<unsigned int>(contractedOrbital.gaussianOrbitals.size());

		return res;
	}

	double ContractedGaussianShell::operator()(const Vec3d& r, const Vec3d& center) const noexcept
	{
		double res = 0;

		for (const auto& contractedOrbital : basisFunctions)
			res += contractedOrbital(r, center);

		return res;
	}


	Vec3d ContractedGaussianShell::GetGradient(const Vec3d& r, const Vec3d& center) const noexcept
	{
		Vec3d res(0, 0, 0);

		for (const auto& orbital : basisFunctions)
			res += orbital.GetGradient(r, center);

		return res;
	}


	double ContractedGaussianShell::GetLaplacian(const Vec3d& r, const Vec3d& center) const noexcept
	{
		double res = 0;

		for (const auto& orbital : basisFunctions)
			res += orbital.GetLaplacian(r, center);

		return res;
	}
}
