#include "pch.h"
#include "Molecule.h"
#include "Integrals.h"

using Eigen::MatrixXd;

namespace myhf
{
//MatrixXd Molecule::OverlapMatrix() noexcept
//{
//	Eigen::MatrixXd overlapMatrix;
//
//	assert(atoms.size() > 1);
//	unsigned int numberOfContractedGaussians = 0;
//	for (const Atom& atom : atoms)
//		numberOfContractedGaussians += basis.GetAtom(atom.type).NumberOfContractedGaussians();
//
//	overlapMatrix = Eigen::MatrixXd::Identity(numberOfContractedGaussians, numberOfContractedGaussians);
//
//	int i = 0;
//	for (const auto& atom1 : atoms)
//	{
//		const BasisAtom& ba1 = basis.GetAtom(atom1.type);
//
//		for (const auto& shell1 : ba1.shells)
//		{
//			for (const auto& orbital1 : shell1.basisFunctions)
//			{
//				int j = 0;
//
//				for (const auto& atom2 : atoms)
//				{
//					const BasisAtom& ba2 = basis.GetAtom(atom2.type);
//
//					for (const auto& shell2 : ba2.shells)
//					{
//						for (const auto& orbital2 : shell2.basisFunctions)
//						{
//							if (j > i)
//							{
//								overlapMatrix(i, j) = OverlapOfTwoOrbitals(ba1, orbital1, atom1.position, ba2, orbital2, atom2.position);
//								overlapMatrix(j, i) = overlapMatrix(i, j);								
//							}
//
//							++j;
//						}
//					}
//				}
//				++i;
//			}
//		}
//	}
//
//	return overlapMatrix;
//}
//
//double Molecule::OverlapOfTwoOrbitals(const BasisAtom& atom1, const ContractedGaussianOrbital& orbital1, const Vec3d& position1, const BasisAtom& atom2, const ContractedGaussianOrbital& orbital2, const Vec3d& position2) noexcept
//{
//	double res = 0;
//
//	for (auto& gaussian1 : orbital1.gaussianOrbitals)
//	{
//		for (auto& gaussian2 : orbital2.gaussianOrbitals)
//		{
//			double overlap = OverlapOfTwoPrimitiveGaussians(gaussian1.alpha, gaussian2.alpha, position1, position2, orbital1.angularMomentum, orbital2.angularMomentum);
//			res += gaussian1.normalizationFactor * gaussian2.normalizationFactor * gaussian1.coefficient * gaussian2.coefficient * overlap;
//		}
//	}
//	return res;
//}
//
//double Molecule::OverlapOfTwoPrimitiveGaussians(double alpha1, double alpha2, const Vec3d& position1, const Vec3d& position2, const QuantumNumbers& angularMomentum1, const QuantumNumbers& angularMomentum2) noexcept
//{
//	assert(std::max(angularMomentum1.l, angularMomentum2.l) + 1 < 4);
//	assert(std::max(angularMomentum1.m, angularMomentum2.m) + 1 < 4);
//	assert(std::max(angularMomentum1.n, angularMomentum2.n) + 1 < 4);
//
//	unsigned int maxAngularMomentumX = std::max(angularMomentum1.l, angularMomentum2.l);
//	unsigned int maxAngularMomentumY = std::max(angularMomentum1.m, angularMomentum2.m);
//	unsigned int maxAngularMomentumZ = std::max(angularMomentum1.n, angularMomentum2.n);
//
//	auto GetIndexX = [maxX = maxAngularMomentumX](unsigned int x, unsigned int y) -> unsigned int { return y * maxX + x; };
//	auto GetIndexY = [maxY = maxAngularMomentumY](unsigned int x, unsigned int y) -> unsigned int { return y * maxY + x; };
//	auto GetIndexZ = [maxZ = maxAngularMomentumZ](unsigned int x, unsigned int y) -> unsigned int { return y * maxZ + x; };
//
//	double oneDividedByAlpha1PlusAlpha2 = 1 / (alpha1 + alpha2);
//
//	// X ==============================================================================
//	// NOTE: I use alloca here for 2 reasons:
//	//			1. It allows for allocating the exact number of bytes needed, whereas we could have gone with a fixed length 
//	//             2D array, but that would require us to pick the size ahead of time (which might lead to a poor choice)
//	//			2. There is better cache coherency when all the values are as tightly packed as possible
//	double* s_x = static_cast<double*>(alloca(sizeof(double) * (maxAngularMomentumX + 2) * (maxAngularMomentumX + 1)));
//
//	// X: Initial conditions
//	s_x[GetIndexX(0, 0)] = 1.0;
//	double startingValue = -(position1.x - ((alpha1 * position1.x + alpha2 * position2.x) * oneDividedByAlpha1PlusAlpha2));
//	s_x[GetIndexX(1, 0)] = startingValue;
//
//	// X: Recurrence Index
//	unsigned int maxAngularMomentumPlus1 = maxAngularMomentumX + 1;
//	for (unsigned int iii = 2; iii <= maxAngularMomentumPlus1; ++iii)
//		s_x[GetIndexX(iii, 0)] = startingValue * s_x[GetIndexX(iii - 1, 0)] + ((iii - 1) * 0.5 * oneDividedByAlpha1PlusAlpha2) * s_x[GetIndexX(iii - 2, 0)];
//
//	// X: Transfer Equation
//	for (unsigned int iii = 0; iii < maxAngularMomentumX; ++iii)
//		s_x[GetIndexX(iii, iii + 1)] = s_x[GetIndexX(iii + 1, iii)] + (position1.x - position2.x) * s_x[GetIndexX(iii, iii)];
//
//
//	// Y ==============================================================================
//	double* s_y = static_cast<double*>(alloca(sizeof(double) * (maxAngularMomentumY + 2) * (maxAngularMomentumY + 1)));
//
//	// Y: Initial conditions
//	s_y[GetIndexY(0, 0)] = 1.0;
//	startingValue = -(position1.y - ((alpha1 * position1.y + alpha2 * position2.y) * oneDividedByAlpha1PlusAlpha2));
//	s_y[GetIndexY(1, 0)] = startingValue;
//
//	// Y: Recurrence Index
//	maxAngularMomentumPlus1 = maxAngularMomentumY + 1;
//	for (unsigned int iii = 2; iii <= maxAngularMomentumPlus1; ++iii)
//		s_y[GetIndexY(iii, 0)] = startingValue * s_y[GetIndexY(iii - 1, 0)] + ((iii - 1) * 0.5 * oneDividedByAlpha1PlusAlpha2) * s_y[GetIndexY(iii - 2, 0)];
//
//	// Y: Transfer Equation
//	for (unsigned int iii = 0; iii < maxAngularMomentumY; ++iii)
//		s_y[GetIndexY(iii, iii + 1)] = s_y[GetIndexY(iii + 1, iii)] + (position1.y - position2.y) * s_y[GetIndexY(iii, iii)];
//
//
//	// Z ==============================================================================
//	double* s_z = static_cast<double*>(alloca(sizeof(double) * (maxAngularMomentumZ + 2) * (maxAngularMomentumZ + 1)));
//
//	// Z: Initial conditions
//	s_z[GetIndexZ(0, 0)] = 1.0;
//	startingValue = -(position1.z - ((alpha1 * position1.z + alpha2 * position2.z) * oneDividedByAlpha1PlusAlpha2));
//	s_z[GetIndexZ(1, 0)] = startingValue;
//
//	// Z: Recurrence Index
//	maxAngularMomentumPlus1 = maxAngularMomentumZ + 1;
//	for (unsigned int iii = 2; iii <= maxAngularMomentumPlus1; ++iii)
//		s_z[GetIndexZ(iii, 0)] = startingValue * s_z[GetIndexZ(iii - 1, 0)] + ((iii - 1) * 0.5 * oneDividedByAlpha1PlusAlpha2) * s_z[GetIndexZ(iii - 2, 0)];
//
//	// Z: Transfer Equation
//	for (unsigned int iii = 0; iii < maxAngularMomentumZ; ++iii)
//		s_z[GetIndexZ(iii, iii + 1)] = s_z[GetIndexZ(iii + 1, iii)] + (position1.z - position2.z) * s_z[GetIndexZ(iii, iii)];
//
//	// Overlap =======================================================================
//	Vec3d diff = position2 - position1;
//	return std::exp(-1 * alpha1 * alpha2 * oneDividedByAlpha1PlusAlpha2 * diff.Dot(diff)) *
//		std::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2, 1.5) *
//		s_x[GetIndexX(angularMomentum1.l, angularMomentum2.l)] *
//		s_y[GetIndexY(angularMomentum1.m, angularMomentum2.m)] *
//		s_z[GetIndexZ(angularMomentum1.n, angularMomentum2.n)];
//}

}