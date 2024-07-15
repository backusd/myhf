#include "pch.h"
#include "Basis.h"


namespace myhf
{
	void Basis::Print() const noexcept
	{
		std::println("Basis: {}", name);
		for (const auto& atom : atoms)
		{
			std::println("AtomWithShells ({}) - electron # = {}, position = {}", static_cast<unsigned int>(atom.type), atom.numberOfElectrons, atom.position);
			for (const auto& shell : atom.shells)
			{
				std::println("    ContractedGaussianShell - # of basis functions = {}", shell.basisFunctions.size());
				for (const auto& f : shell.basisFunctions)
				{
					std::println("        ContractedGaussianOrbital - l/m/n = ({}, {}, {}), ", f.angularMomentum.l, f.angularMomentum.m, f.angularMomentum.n);
					for (const auto& o : f.gaussianOrbitals)
					{
						std::println("            PrimitiveGaussian - alpha = {}, coeff = {}, normFactor = {}, coeffProdNorm = {}", o.alpha, o.coefficient, o.normalizationFactor, o.coeffProdNorm);
					}
				}
			}
		}
	}


}