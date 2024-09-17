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
	std::vector<Atom> atoms = 
	{ 
		{ ATOM_TYPE::Hydrogen, 1, { 0,  1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Hydrogen, 1, { 0, -1.43233673, -0.96104039 } },
		{ ATOM_TYPE::Oxygen,   8, { 0,           0,  0.24026010 } }
	};

	for (int i = 0; i < 10; ++i)
	{
		scope_time s;
		Eigen::MatrixXd overlapMatrix = OverlapMatrix(atoms, STO_3G);
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