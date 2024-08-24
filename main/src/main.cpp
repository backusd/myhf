#include "pch.h"

#include <Eigen\eigen>
#include <myhf.h>

using namespace myhf;

#include <chrono>

//void PrintBasis(const Basis& b)
//{
//	std::println("Basis: {}", b.name);
//	for (const auto& atom : b.atoms)
//	{
//		std::println("AtomWithShells ({}) - electron # = {}, position = {}", static_cast<unsigned int>(atom.type), atom.numberOfElectrons, atom.position);
//		for (const auto& shell : atom.shells)
//		{
//			std::println("    ContractedGaussianShell - # of basis functions = {}", shell.basisFunctions.size());
//			for (const auto& f : shell.basisFunctions)
//			{
//				std::println("        ContractedGaussianOrbital - l/m/n = ({}, {}, {}), ", f.angularMomentum.l, f.angularMomentum.m, f.angularMomentum.n);
//				for (const auto& o : f.gaussianOrbitals)
//				{
//					std::println("            PrimitiveGaussian - alpha = {}, coeff = {}, normFactor = {}, coeffProdNorm = {}", o.alpha, o.coefficient, o.normalizationFactor, o.coeffProdNorm);
//				}
//			}
//		}
//	}
//}

struct scope_time
{
	scope_time() : start(std::chrono::system_clock::now()) {}
	~scope_time()
	{
		std::println("Total Time: {}", std::chrono::duration_cast<std::chrono::nanoseconds>(std::chrono::system_clock::now() - start));
	}

	std::chrono::system_clock::time_point start{};
};

int main()
{
	//PrintBasis(STO_3G);
	//PrintBasis(STO_6G);

//	std::array<Atom, 2> atoms =
//	{
//		STO_3G.GetAtom(ATOM_TYPE::Hydrogen),
//		STO_3G.GetAtom(ATOM_TYPE::Hydrogen)
//	};
//	atoms[1].position = { 0.5356598775430879, 0, 0 }; //  0.283459 * 2
//	Eigen::MatrixXd overlapMatrix = OverlapMatrix(atoms);
//
//
//	std::print("{}", overlapMatrix);

	std::vector<Atom> atoms = 
	{ 
		{ ATOM_TYPE::Hydrogen, 1, { 0,  1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Hydrogen, 1, { 0, -1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Oxygen,   8, { 0,           0,  0.24026010 } }
	};
	Molecule molec(std::move(atoms), STO_3G);

	for(int i = 0; i < 10; ++i)
	{
		scope_time s;
		Eigen::MatrixXd overlapMatrix = molec.OverlapMatrix();
	}

	std::println("=-------------------------------------------------------------");

//	for (int i = 0; i < 10; ++i)
//	{
//		scope_time s;
//		Eigen::MatrixXd overlapMatrix = molec.OverlapMatrix2();
//	}
	
	Eigen::MatrixXd overlapMatrix = molec.OverlapMatrix();
	std::println("{}\n", overlapMatrix);

//	Eigen::MatrixXd overlapMatrix2 = molec.OverlapMatrix2();
//	std::println("{}\n", overlapMatrix2);
}