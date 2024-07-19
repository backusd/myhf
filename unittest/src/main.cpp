#include "pch.h"
#include <myhf.h>
using namespace myhf;

int main()
{
	std::cout << "Unit test\n";

	std::array<Atom, 2> atoms =
	{
		STO_3G.GetAtom(ATOM_TYPE::Hydrogen),
		STO_3G.GetAtom(ATOM_TYPE::Hydrogen)
	};
	atoms[1].position = { 0.5356598775430879, 0, 0 }; //  0.283459 * 2
	Eigen::MatrixXd overlapMatrix = OverlapMatrix(atoms);


	std::print("{}", overlapMatrix);
}