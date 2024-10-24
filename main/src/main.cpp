#include "pch.h"

#include <Eigen\eigen>
#include <myhf.h>

using namespace myhf;

int main()
{

	test::test2::STO_3G::Hydrogen h;

	std::println("{}, {}, {}, {}\n{}, {}, {}, {}\n{}, {}, {}, {}\n", 
		h.orbital_1s.func1.primitiveGaussians[0].alpha,
		h.orbital_1s.func1.primitiveGaussians[0].coefficient,
		h.orbital_1s.func1.primitiveGaussians[0].normalizationFactor,
		h.orbital_1s.func1.primitiveGaussians[0].coeffProdNorm, 
		
		h.orbital_1s.func1.primitiveGaussians[1].alpha,
		h.orbital_1s.func1.primitiveGaussians[1].coefficient,
		h.orbital_1s.func1.primitiveGaussians[1].normalizationFactor,
		h.orbital_1s.func1.primitiveGaussians[1].coeffProdNorm, 
		
		h.orbital_1s.func1.primitiveGaussians[2].alpha,
		h.orbital_1s.func1.primitiveGaussians[2].coefficient,
		h.orbital_1s.func1.primitiveGaussians[2].normalizationFactor,
		h.orbital_1s.func1.primitiveGaussians[2].coeffProdNorm);

	// Water
//	std::vector<Atom> atoms =
//	{
//		{ ATOM_TYPE::Hydrogen, 1, { 0.0,  1.43233673, -0.96104039 } },
//		{ ATOM_TYPE::Hydrogen, 1, { 0.0, -1.43233673, -0.96104039 } },
//		{ ATOM_TYPE::Oxygen, 8, { 0.0, 0.0, 0.24026010 } },
//	};

//	atoms =
//	{
//		{ ATOM_TYPE::Hydrogen, 1, { 0.5,  1.43233673, -0.96104039 } },
//		{ ATOM_TYPE::Aluminum, 13, { 0.2, -0.1234, 0.24026010 } },
//	};






//	PROFILE_BEGIN_SESSION("Session 1", "results.json");
//	std::println("Started...");
//	Eigen::MatrixXd nuclearMatrix = NuclearElectronAttractionEnergyMatrix(atoms, STO_3G);
//	std::cout << "\nNuclear:\n" << std::setprecision(5) << nuclearMatrix << '\n';
//	PROFILE_END_SESSION();





	

}