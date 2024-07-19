#include "pch.h"

#include <Eigen\eigen>
#include <myhf.h>

using namespace myhf;

void PrintBasis(const Basis& b)
{
	std::println("Basis: {}", b.name);
	for (const auto& atom : b.atoms)
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

int main()
{
	//PrintBasis(STO_3G);
	//PrintBasis(STO_6G);

	std::array<Atom, 2> atoms =
	{
		STO_3G.GetAtom(ATOM_TYPE::Hydrogen),
		STO_3G.GetAtom(ATOM_TYPE::Hydrogen)
	};
	atoms[1].position = { 0.5356598775430879, 0, 0 }; //  0.283459 * 2
	Eigen::MatrixXd overlapMatrix = OverlapMatrix(atoms);


	std::print("{}", overlapMatrix);

	//std::cout << overlapMatrix << std::endl;
}