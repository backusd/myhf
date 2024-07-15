#include "pch.h"

#include <Eigen\eigen>

#include "Orbital.h"
#include "Vec.h"
#include "MathUtils.h"
#include "Shell.h"
#include "Basis.h"
#include "Integrals.h"

using namespace myhf;

int main()
{
	std::array<Atom, 2> atoms = 
	{ 
		STO3G.GetAtom(ATOM_TYPE::Hydrogen), 
		STO3G.GetAtom(ATOM_TYPE::Hydrogen) 
	};
	atoms[1].position = { 0.283459 * 2, 0, 0 };
	Eigen::MatrixXd overlapMatrix = OverlapMatrix(atoms);


	std::cout << overlapMatrix << std::endl;
}