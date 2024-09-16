#include "pch.h"
#include "MatrixTests.h"

using namespace myhf;

int main()
{
	OverlapMatrixTest();

//	{
//		std::vector<Atom> atoms{
//			{ ATOM_TYPE::Hydrogen, 1, { 0.0, 0.0, 0.0 } },
//			{ ATOM_TYPE::Lithium, 3, { 1.0, 0.0, 0.0 } }
//		};
//
//		Eigen::MatrixXd actualOverlap = OverlapMatrix(atoms, STO_3G);
//		std::cout << "Overlap:\n" << actualOverlap << '\n';
//	}

//	{
//		std::vector<Atom> atoms{
//			{ ATOM_TYPE::Hydrogen, 1, { 0.0, 0.0, 0.0 } },
//			{ ATOM_TYPE::Magnesium, 12, { 1.0, 0.0, 0.0 } }
//		};
//
//		Eigen::MatrixXd actualOverlap = OverlapMatrix(atoms, STO_6G);
//		std::cout << "\nOverlap:\n" << actualOverlap << '\n';
//	}
}