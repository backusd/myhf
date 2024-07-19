#pragma once

#include "Vec.h"
#include "Shell.h"

namespace myhf
{
enum class ATOM_TYPE : unsigned int
{
	Hydrogen	= 1,
	Helium		= 2,
	Lithium		= 3,
	Beryllium	= 4,
	Boron		= 5,
	Carbon		= 6,
	Nitrogen	= 7,
	Oxygen		= 8,
	Flourine	= 9,
	Neon		= 10,
	Sodium		= 11, 
	Magnesium	= 12,
	Aluminum	= 13,
	Silicon		= 14, 
	Phosphorous	= 15,
	Sulfur		= 16,
	Chlorine	= 17,
	Argon		= 18,
	Potassium	= 19,
	Calcium		= 20
};

struct Atom
{
	ATOM_TYPE	 type = ATOM_TYPE::Hydrogen;
	unsigned int numberOfElectrons = 1;
	Vec3d		 position = {};
	std::vector<ContractedGaussianShell> shells;

	constexpr unsigned int NumberOfContractedGaussians() const noexcept
	{
		unsigned int res = 0;
		for (auto const& shell : shells)
			res += shell.CountNumberOfContractedGaussians();
		return res;
	}
	constexpr QuantumNumbers GetMaxQN(double alpha) const noexcept
	{
		QuantumNumbers maxQN(0, 0, 0);

		for (const auto& shell : shells)
			for (const auto& orbital : shell.basisFunctions)
				for (const auto& gaussian : orbital.gaussianOrbitals)
					if (alpha == gaussian.alpha)
					{
						maxQN.l = std::max(maxQN.l, orbital.angularMomentum.l);
						maxQN.m = std::max(maxQN.m, orbital.angularMomentum.m);
						maxQN.n = std::max(maxQN.m, orbital.angularMomentum.n);
					}

		return maxQN;
	}
};

}