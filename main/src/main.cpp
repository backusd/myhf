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
//		double alpha1 = 3.425250914;
//		double alpha2 = 5.033151319;
//		Vec3d position1{ 0., 1.43233673, -0.96104039 };
//		Vec3d position2{ 0., 0., 0.24026010 };
//		QuantumNumbers angularMomentum1{0, 0, 0};
//		QuantumNumbers angularMomentum2{0, 0, 1};
//
//		double result = NuclearElectronAttractionOfTwoPrimitiveGaussians(alpha1, alpha2, position1, position2, position1, angularMomentum1, angularMomentum2); 
//		int iii = 0;


		std::vector<Atom> atoms =
		{
			{ ATOM_TYPE::Hydrogen, 1, { 0.0,  1.43233673, -0.96104039 } },
			{ ATOM_TYPE::Hydrogen, 1, { 0.0, -1.43233673, -0.96104039 } },
			{ ATOM_TYPE::Oxygen, 8, { 0.0, 0.0, 0.24026010 } },
		};

	//	std::vector<Atom> atoms =
	//	{
	//		{ ATOM_TYPE::Hydrogen, 1, { 0.0, 0.0, 0.0 } },
	//		{ ATOM_TYPE::Magnesium, 12, { 1.0, 0.5, 0.75 } }
	//	};

	//	Eigen::MatrixXd overlapMatrix = OverlapMatrix(atoms, STO_3G);
	//	std::cout << "Overlap 1:\n" << std::setprecision(5) << overlapMatrix << '\n';

	//	Eigen::MatrixXd kineticMatrix = KineticEnergyMatrix(atoms, STO_3G);
	//	std::cout << "\nKinetic 1:\n" << kineticMatrix << '\n';

		std::println("Started...");
		Eigen::MatrixXd nuclearMatrix = NuclearElectronAttractionEnergyMatrix(atoms, STO_3G);
		std::cout << "\nNuclear:\n" << nuclearMatrix << '\n';



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