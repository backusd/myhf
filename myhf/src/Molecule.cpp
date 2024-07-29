#include "pch.h"
#include "Molecule.h"
#include "Integrals.h"

using Eigen::MatrixXd;

namespace myhf
{
double Molecule::GetOverlap(const BasisAtom& atom1, const ContractedGaussianOrbital& orbital1, const Vec3d& position1, const BasisAtom& atom2, const ContractedGaussianOrbital& orbital2, const Vec3d& position2) noexcept
{
	double res = 0;

	for (auto& gaussian1 : orbital1.gaussianOrbitals)
	{
		for (auto& gaussian2 : orbital2.gaussianOrbitals)
		{
			QuantumNumbers maxQN1 = atom1.GetMaxQN(gaussian1.alpha);
			QuantumNumbers maxQN2 = atom2.GetMaxQN(gaussian2.alpha);
			GaussianMoment moment(gaussian1.alpha, gaussian2.alpha, position1, position2, maxQN1, maxQN2);
			double overlap = moment.getOverlap(orbital1.angularMomentum, orbital2.angularMomentum);
			res += gaussian1.normalizationFactor * gaussian2.normalizationFactor * gaussian1.coefficient * gaussian2.coefficient * overlap;
		}
	}
	return res;
}
MatrixXd Molecule::OverlapMatrix() noexcept
{
	Eigen::MatrixXd overlapMatrix;

	assert(atoms.size() > 1);
	unsigned int numberOfContractedGaussians = 0;
	for (const Atom& atom : atoms)
		numberOfContractedGaussians += basis.GetAtom(atom.type).NumberOfContractedGaussians();

	overlapMatrix = Eigen::MatrixXd::Zero(numberOfContractedGaussians, numberOfContractedGaussians);

	int i = 0;
	for (const auto& atom1 : atoms)
	{
		const BasisAtom& ba1 = basis.GetAtom(atom1.type);

		for (const auto& shell1 : ba1.shells)
		{
			for (const auto& orbital1 : shell1.basisFunctions)
			{
				int j = 0;

				for (const auto& atom2 : atoms)
				{
					const BasisAtom& ba2 = basis.GetAtom(atom2.type);

					for (const auto& shell2 : ba2.shells)
					{
						for (const auto& orbital2 : shell2.basisFunctions)
						{
							if (j >= i)
							{
								overlapMatrix(i, j) = GetOverlap(ba1, orbital1, atom1.position, ba2, orbital2, atom2.position);
								if (i != j)
									overlapMatrix(j, i) = overlapMatrix(i, j);
							}

							++j;
						}
					}
				}
				++i;
			}
		}
	}

	return overlapMatrix;
}



}