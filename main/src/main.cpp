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
	//	std::vector<Atom> atoms =
	//	{
	//		{ ATOM_TYPE::Hydrogen, 1, { 0,  1.43233673, -0.96104039 } },
	//		{ ATOM_TYPE::Hydrogen, 1, { 0, -1.43233673, -0.96104039 } },
	//		{ ATOM_TYPE::Oxygen,   8, { 0,           0,  0.24026010 } }
	//	};

		std::vector<Atom> atoms =
		{
			{ ATOM_TYPE::Hydrogen, 1, { 0.0, 0.0, 0.0 } },
			{ ATOM_TYPE::Magnesium, 12, { 1.0, 0.5, 0.75 } }
		};

		Eigen::MatrixXd overlapMatrix = OverlapMatrix(atoms, STO_3G);
		std::cout << "Overlap:\n" << std::setprecision(8) << overlapMatrix << '\n';

		// My original: 0.017479982
		// Expected:    0.017480023


	//	Eigen::MatrixXd kineticMatrix = KineticEnergyMatrix(atoms, STO_3G);
	//	std::cout << "Kinetic:\n" << kineticMatrix << '\n';
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