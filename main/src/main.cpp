#include "pch.h"

#include <Eigen\eigen>
#include <myhf.h>

using namespace myhf;

int main()
{
//	double alpha1 = 3.425250914;
//	double alpha2 = 5.033151319;
//	Vec3d position1{ 0., 1.43233673, -0.96104039 };
//	Vec3d position2{ 0., 0., 0.24026010 };
//	QuantumNumbers angularMomentum1{0, 0, 0};
//	QuantumNumbers angularMomentum2{0, 0, 1};
//
//	double result = NuclearElectronAttractionOfTwoPrimitiveGaussians(alpha1, alpha2, position1, position2, position1, angularMomentum1, angularMomentum2); 
//	int iii = 0;


	// Water
	std::vector<Atom> atoms =
	{
		{ ATOM_TYPE::Hydrogen, 1, { 0.0,  1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Hydrogen, 1, { 0.0, -1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Oxygen, 8, { 0.0, 0.0, 0.24026010 } },
	};

//	atoms =
//	{
//		{ ATOM_TYPE::Hydrogen, 1, { 0.5,  1.43233673, -0.96104039 } },
//		{ ATOM_TYPE::Aluminum, 13, { 0.2, -0.1234, 0.24026010 } },
//	};

	//std::vector<Atom> atoms =
	//{
	//	{ ATOM_TYPE::Hydrogen, 1, { 0.0, 0.0, 0.0 } },
	//	{ ATOM_TYPE::Lithium, 3, { 1.0, 0.5, 0.75 } },
	//};

	//std::vector<Atom> atoms =
	//{
	//	{ ATOM_TYPE::Hydrogen, 1, { 0.0,  1.43233673, -0.96104039 } },
	//	{ ATOM_TYPE::Helium, 2, { 0.0, -1.43233673, -0.96104039 } },
	//	{ ATOM_TYPE::Hydrogen, 1, { 0.0, 0.0, 0.24026010 } },
	//	{ ATOM_TYPE::Helium, 2, { 1.0, -1.43233673, -0.96104039 } },
	//};

//	Eigen::MatrixXd overlapMatrix = OverlapMatrix(atoms, STO_3G);
//	std::cout << "Overlap 1:\n" << std::setprecision(5) << overlapMatrix << '\n';

//	std::vector<Atom> atoms =
//	{
//		{ ATOM_TYPE::Hydrogen, 1, { 0.0, 0.0, 0.0 } },
//		{ ATOM_TYPE::Magnesium, 12, { 1.0, 0.5, 0.75 } }
//	};

//	Eigen::MatrixXd overlapMatrix = OverlapMatrix(atoms, STO_3G);
//	std::cout << "Overlap 1:\n" << std::setprecision(5) << overlapMatrix << '\n';
//
//	Eigen::MatrixXd overlapMatrix2 = OverlapMatrix_2(atoms, STO_3G);
//	std::cout << "Overlap 2:\n" << std::setprecision(5) << overlapMatrix2 << '\n';

//	Eigen::MatrixXd kineticMatrix = KineticEnergyMatrix(atoms, STO_3G);
//	std::cout << "\nKinetic 1:\n" << kineticMatrix << '\n';

//	PROFILE_BEGIN_SESSION("Session 1", "results.json");
//	std::println("Started...");
//	Eigen::MatrixXd nuclearMatrix = NuclearElectronAttractionEnergyMatrix(atoms, STO_3G);
//	std::cout << "\nNuclear:\n" << std::setprecision(8) << nuclearMatrix << '\n';
//	PROFILE_END_SESSION();

	std::println("Started 2...");
	Eigen::MatrixXd nuclearMatrix2 = NuclearElectronAttractionEnergyMatrix_2(atoms, STO_3G);
	std::cout << "\nNuclear:\n" << std::setprecision(5) << nuclearMatrix2 << '\n';

	std::println("\nStarted 3...");
	Eigen::MatrixXd nuclearMatrix3 = NuclearElectronAttractionEnergyMatrix_3(atoms, STO_3G);
	std::cout << "\nNuclear:\n" << std::setprecision(5) << nuclearMatrix3 << '\n';

//	std::println("\nStarted 4...");
//	Eigen::MatrixXd nuclearMatrix4 = NuclearElectronAttractionEnergyMatrix_4(atoms, STO_3G);
//	std::cout << "\nNuclear:\n" << std::setprecision(5) << nuclearMatrix4 << '\n';

	std::println("\nStarted 5...");
	Eigen::MatrixXd nuclearMatrix5 = NuclearElectronAttractionEnergyMatrix_5(atoms, STO_3G);
	std::cout << "\nNuclear:\n" << std::setprecision(5) << nuclearMatrix5 << '\n';

//
//	Eigen::MatrixXd nuclearMatrix2 = NuclearElectronAttractionEnergyMatrix_Par(atoms, STO_3G);
//	std::cout << "\nNuclear Par:\n" << nuclearMatrix2 << '\n';

	

}