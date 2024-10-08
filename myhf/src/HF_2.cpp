#include "pch.h"
#include "HF_2.h"
#include "DataDrivenTest.h"

#include "gcem.hpp"

using Eigen::MatrixXd;

namespace myhf
{
namespace test
{
static unsigned int NumberOfContractedGaussians(ATOM_TYPE atomType) noexcept
{
	static constexpr std::array<unsigned int, 20> counts = {
		Hydrogen::NumberOfContractedGaussians,
		Helium::NumberOfContractedGaussians,
		0,
		0,
		0,
		0,
		Nitrogen::NumberOfContractedGaussians,
		Oxygen::NumberOfContractedGaussians,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0,
		0
	};

	return counts[static_cast<unsigned int>(atomType) - 1];
}

static void OverlapAtomAtom(const Atom& atom_1, const Atom& atom_2, unsigned int row, unsigned int col, MatrixXd& overlapMatrix) noexcept
{
	switch (atom_1.type)
	{
	case ATOM_TYPE::Hydrogen:
	{
		switch (atom_2.type)
		{
		case ATOM_TYPE::Hydrogen: OverlapAtomAtom_1S1S<Hydrogen, Hydrogen>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Helium:   OverlapAtomAtom_1S1S<Hydrogen, Helium>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Nitrogen: OverlapAtomAtom_1S2SP<Hydrogen, Nitrogen>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Oxygen:   OverlapAtomAtom_1S2SP<Hydrogen, Oxygen>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		}
		return;
	}
	case ATOM_TYPE::Helium:
	{
		switch (atom_2.type)
		{
		case ATOM_TYPE::Hydrogen: OverlapAtomAtom_1S1S<Helium, Hydrogen>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Helium:   OverlapAtomAtom_1S1S<Helium, Helium>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Nitrogen: OverlapAtomAtom_1S2SP<Helium, Nitrogen>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Oxygen:   OverlapAtomAtom_1S2SP<Helium, Oxygen>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		}
		return;
	}
	case ATOM_TYPE::Nitrogen:
	{
		switch (atom_2.type)
		{
		case ATOM_TYPE::Hydrogen: OverlapAtomAtom_1S2SP<Nitrogen, Hydrogen>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Helium:   OverlapAtomAtom_1S2SP<Nitrogen, Helium>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Nitrogen: OverlapAtomAtom_2SP2SP<Nitrogen, Nitrogen>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Oxygen:   OverlapAtomAtom_2SP2SP<Nitrogen, Oxygen>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		}
	}
	case ATOM_TYPE::Oxygen:
	{
		switch (atom_2.type)
		{
		case ATOM_TYPE::Hydrogen: OverlapAtomAtom_1S2SP<Oxygen, Hydrogen>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Helium:   OverlapAtomAtom_1S2SP<Oxygen, Helium>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Nitrogen: OverlapAtomAtom_2SP2SP<Oxygen, Nitrogen>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Oxygen:   OverlapAtomAtom_2SP2SP<Oxygen, Oxygen>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		}
	}
	}
}



static void OverlapWithSelf(ATOM_TYPE type, unsigned int topLeftIndex, MatrixXd& overlapMatrix)
{
	switch (type)
	{
	case ATOM_TYPE::Hydrogen: OverlapWithSelf_1S<Hydrogen>(topLeftIndex, overlapMatrix); return;
	case ATOM_TYPE::Helium:   OverlapWithSelf_1S<Helium>(topLeftIndex, overlapMatrix); return;
	case ATOM_TYPE::Nitrogen: OverlapWithSelf_2SP<Nitrogen>(topLeftIndex, overlapMatrix); return;
	case ATOM_TYPE::Oxygen:   OverlapWithSelf_2SP<Oxygen>(topLeftIndex, overlapMatrix); return;
	}
}



MatrixXd OverlapMatrix(std::span<Atom> atoms) noexcept
{
	MatrixXd overlapMatrix;

	assert(atoms.size() > 1);
	unsigned int numberOfContractedGaussians = 0;
	for (const Atom& atom : atoms)
		numberOfContractedGaussians += NumberOfContractedGaussians(atom.type);

	overlapMatrix = MatrixXd::Identity(numberOfContractedGaussians, numberOfContractedGaussians);

	// Optimization for dealing with overlap of atoms with themselves:
	//     The overlap of an atom's orbitals with itself will always be the same. For example, 
	//	   in the case of an oxygen atom, there are 5 orbitals, so each oxygen atom in the 
	//     matrix will include a 5x5 block, but each block will be identical because the overlap
	//     of the orbitals with itself is always the same
	unsigned int topLeftIndex = 0;
	for (const Atom& atom : atoms)
	{
		OverlapWithSelf(atom.type, topLeftIndex, overlapMatrix);
		topLeftIndex += NumberOfContractedGaussians(atom.type);
	}

	unsigned int row = 0;
	for (size_t atom1Index = 0; atom1Index + 1 < atoms.size(); ++atom1Index)
	{
		const Atom& atom_1 = atoms[atom1Index];
		unsigned int col = row + NumberOfContractedGaussians(atom_1.type);

		for (size_t atom2Index = atom1Index + 1; atom2Index < atoms.size(); ++atom2Index)
		{
			const Atom& atom_2 = atoms[atom2Index];
			OverlapAtomAtom(atom_1, atom_2, row, col, overlapMatrix);
			col += NumberOfContractedGaussians(atom_2.type);
		}

		row += NumberOfContractedGaussians(atom_1.type);
	}

	return overlapMatrix;
}

}
}