#include "pch.h"

#include <Eigen\eigen>
#include <myhf.h>

using namespace myhf;

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
	{
		std::vector<Atom> atoms =
		{
			{ ATOM_TYPE::Hydrogen, 1, { 0,  1.43233673, -0.96104039 } },
			{ ATOM_TYPE::Oxygen, 8, { 1,  1.43233673, -0.96104039 } },
			{ ATOM_TYPE::Hydrogen, 1, { 0, -1.43233673, -0.96104039 } },
			{ ATOM_TYPE::Oxygen, 8, { -1,  1.43233673, -0.96104039 } },
		};

	//	std::vector<Atom> atoms =
	//	{
	//		{ ATOM_TYPE::Hydrogen, 1, { 0.0, 0.0, 0.0 } },
	//		{ ATOM_TYPE::Magnesium, 12, { 1.0, 0.5, 0.75 } }
	//	};

	//	Eigen::MatrixXd overlapMatrix = OverlapMatrix(atoms, STO_3G);
	//	std::cout << "Overlap 1:\n" << std::setprecision(5) << overlapMatrix << '\n';
	//
	//	Eigen::MatrixXd kineticMatrix = KineticEnergyMatrix(atoms, STO_3G);
	//	std::cout << "\nKinetic 1:\n" << kineticMatrix << '\n';

		Eigen::MatrixXd overlap2, kinetic2;
		OverlapAndKineticEnergyMatrix2(atoms, STO_3G, overlap2, kinetic2);
		std::cout << "\nOverlap 2:\n" << std::setprecision(8) << overlap2 << '\n';
		std::cout << "\nKinetic 2:\n" << std::setprecision(8) << kinetic2 << '\n';

	}

	std::println("-------------------------------------------------------------");



//	for (int i = 0; i < 10; ++i)
//	{
//		scope_time s;
//		Eigen::MatrixXd overlapMatrix = OverlapMatrix2(atoms, STO_3G);
//	}
//
//	
//	Eigen::MatrixXd overlapMatrix = OverlapMatrix(atoms, STO_3G);
//	std::println("{}\n", overlapMatrix);
//
//
//
//	Eigen::MatrixXd overlapMatrix2 = OverlapMatrix2(atoms, STO_3G);
//	std::println("\n\n{}\n", overlapMatrix2);
}