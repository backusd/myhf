#include "pch.h"
#include "HF_Fast.h"
#include "Basis_Fast.h"

#include "gcem.hpp"

using Eigen::MatrixXd;

namespace myhf
{
namespace test
{
class NumberOfContractedGaussians_util
{
public:
	static constexpr unsigned int GetCount(ATOM_TYPE atomType) noexcept
	{
		return m_counts[static_cast<unsigned int>(atomType) - 1];
	}
private:
	static constexpr std::array<unsigned int, 20> m_counts = {
		Hydrogen::NumberOfContractedGaussians,
		Helium::NumberOfContractedGaussians,
		Lithium::NumberOfContractedGaussians,
		Beryllium::NumberOfContractedGaussians,
		Boron::NumberOfContractedGaussians,
		Carbon::NumberOfContractedGaussians,
		Nitrogen::NumberOfContractedGaussians,
		Oxygen::NumberOfContractedGaussians,
		Flourine::NumberOfContractedGaussians,
		Neon::NumberOfContractedGaussians,
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
};

static constexpr unsigned int NumberOfContractedGaussians(ATOM_TYPE atomType) noexcept
{
	// This is a somewhat silly workaround of the fact that you cannot have a static constexpr
	// variable in a constexpr function. (See here: https://stackoverflow.com/questions/62458079/static-constexpr-variables-in-a-constexpr-function)
	// However, if you encapsulate the data into a class like we do here, it works as expected.
	return NumberOfContractedGaussians_util::GetCount(atomType);
}

static void OverlapAtomAtom(const Atom& atom_1, const Atom& atom_2, unsigned int row, unsigned int col, MatrixXd& overlapMatrix) noexcept
{
	switch (atom_1.type)
	{
	case ATOM_TYPE::Hydrogen:
	{
		switch (atom_2.type)
		{
		case ATOM_TYPE::Hydrogen:	OverlapAtomAtom_1S1S<Hydrogen, Hydrogen>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Helium:		OverlapAtomAtom_1S1S<Hydrogen, Helium>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Lithium:	OverlapAtomAtom_1S2SP<Hydrogen, Lithium>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Beryllium:  OverlapAtomAtom_1S2SP<Hydrogen, Beryllium>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Boron:		OverlapAtomAtom_1S2SP<Hydrogen, Boron>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Carbon:		OverlapAtomAtom_1S2SP<Hydrogen, Carbon>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Nitrogen:	OverlapAtomAtom_1S2SP<Hydrogen, Nitrogen>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Oxygen:		OverlapAtomAtom_1S2SP<Hydrogen, Oxygen>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Flourine:   OverlapAtomAtom_1S2SP<Hydrogen, Flourine>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Neon:		OverlapAtomAtom_1S2SP<Hydrogen, Neon>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		}
		return;
	}
	case ATOM_TYPE::Helium:
	{
		switch (atom_2.type)
		{
		case ATOM_TYPE::Hydrogen:	OverlapAtomAtom_1S1S<Helium, Hydrogen>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Helium:		OverlapAtomAtom_1S1S<Helium, Helium>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Lithium:	OverlapAtomAtom_1S2SP<Helium, Lithium>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Beryllium:  OverlapAtomAtom_1S2SP<Helium, Beryllium>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Boron:		OverlapAtomAtom_1S2SP<Helium, Boron>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Carbon:		OverlapAtomAtom_1S2SP<Helium, Carbon>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Nitrogen:	OverlapAtomAtom_1S2SP<Helium, Nitrogen>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Oxygen:		OverlapAtomAtom_1S2SP<Helium, Oxygen>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Flourine:   OverlapAtomAtom_1S2SP<Helium, Flourine>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Neon:		OverlapAtomAtom_1S2SP<Helium, Neon>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		}
		return;
	}
	case ATOM_TYPE::Lithium:
	{
		switch (atom_2.type)
		{
		case ATOM_TYPE::Hydrogen:	OverlapAtomAtom_1S2SP<Lithium, Hydrogen>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Helium:		OverlapAtomAtom_1S2SP<Lithium, Helium>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Lithium:	OverlapAtomAtom_2SP2SP<Lithium, Lithium>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Beryllium:  OverlapAtomAtom_2SP2SP<Lithium, Beryllium>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Boron:		OverlapAtomAtom_2SP2SP<Lithium, Boron>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Carbon:		OverlapAtomAtom_2SP2SP<Lithium, Carbon>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Nitrogen:	OverlapAtomAtom_2SP2SP<Lithium, Nitrogen>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Oxygen:		OverlapAtomAtom_2SP2SP<Lithium, Oxygen>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Flourine:   OverlapAtomAtom_2SP2SP<Lithium, Flourine>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Neon:		OverlapAtomAtom_2SP2SP<Lithium, Neon>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		}
	}
	case ATOM_TYPE::Beryllium:
	{
		switch (atom_2.type)
		{
		case ATOM_TYPE::Hydrogen:	OverlapAtomAtom_1S2SP<Beryllium, Hydrogen>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Helium:		OverlapAtomAtom_1S2SP<Beryllium, Helium>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Lithium:	OverlapAtomAtom_2SP2SP<Beryllium, Lithium>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Beryllium:  OverlapAtomAtom_2SP2SP<Beryllium, Beryllium>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Boron:		OverlapAtomAtom_2SP2SP<Beryllium, Boron>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Carbon:		OverlapAtomAtom_2SP2SP<Beryllium, Carbon>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Nitrogen:	OverlapAtomAtom_2SP2SP<Beryllium, Nitrogen>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Oxygen:		OverlapAtomAtom_2SP2SP<Beryllium, Oxygen>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Flourine:   OverlapAtomAtom_2SP2SP<Beryllium, Flourine>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Neon:		OverlapAtomAtom_2SP2SP<Beryllium, Neon>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		}
	}
	case ATOM_TYPE::Boron:
	{
		switch (atom_2.type)
		{
		case ATOM_TYPE::Hydrogen:	OverlapAtomAtom_1S2SP<Boron, Hydrogen>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Helium:		OverlapAtomAtom_1S2SP<Boron, Helium>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Lithium:	OverlapAtomAtom_2SP2SP<Boron, Lithium>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Beryllium:  OverlapAtomAtom_2SP2SP<Boron, Beryllium>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Boron:		OverlapAtomAtom_2SP2SP<Boron, Boron>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Carbon:		OverlapAtomAtom_2SP2SP<Boron, Carbon>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Nitrogen:	OverlapAtomAtom_2SP2SP<Boron, Nitrogen>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Oxygen:		OverlapAtomAtom_2SP2SP<Boron, Oxygen>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Flourine:   OverlapAtomAtom_2SP2SP<Boron, Flourine>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Neon:		OverlapAtomAtom_2SP2SP<Boron, Neon>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		}
	}
	case ATOM_TYPE::Carbon:
	{
		switch (atom_2.type)
		{
		case ATOM_TYPE::Hydrogen:	OverlapAtomAtom_1S2SP<Carbon, Hydrogen>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Helium:		OverlapAtomAtom_1S2SP<Carbon, Helium>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Lithium:	OverlapAtomAtom_2SP2SP<Carbon, Lithium>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Beryllium:  OverlapAtomAtom_2SP2SP<Carbon, Beryllium>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Boron:		OverlapAtomAtom_2SP2SP<Carbon, Boron>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Carbon:		OverlapAtomAtom_2SP2SP<Carbon, Carbon>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Nitrogen:	OverlapAtomAtom_2SP2SP<Carbon, Nitrogen>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Oxygen:		OverlapAtomAtom_2SP2SP<Carbon, Oxygen>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Flourine:   OverlapAtomAtom_2SP2SP<Carbon, Flourine>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Neon:		OverlapAtomAtom_2SP2SP<Carbon, Neon>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		}
	}
	case ATOM_TYPE::Nitrogen:
	{
		switch (atom_2.type)
		{
		case ATOM_TYPE::Hydrogen:	OverlapAtomAtom_1S2SP<Nitrogen, Hydrogen>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Helium:		OverlapAtomAtom_1S2SP<Nitrogen, Helium>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Lithium:	OverlapAtomAtom_2SP2SP<Nitrogen, Lithium>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Beryllium:  OverlapAtomAtom_2SP2SP<Nitrogen, Beryllium>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Boron:		OverlapAtomAtom_2SP2SP<Nitrogen, Boron>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Carbon:		OverlapAtomAtom_2SP2SP<Nitrogen, Carbon>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Nitrogen:	OverlapAtomAtom_2SP2SP<Nitrogen, Nitrogen>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Oxygen:		OverlapAtomAtom_2SP2SP<Nitrogen, Oxygen>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Flourine:   OverlapAtomAtom_2SP2SP<Nitrogen, Flourine>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Neon:		OverlapAtomAtom_2SP2SP<Nitrogen, Neon>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		}
	}
	case ATOM_TYPE::Oxygen:
	{
		switch (atom_2.type)
		{
		case ATOM_TYPE::Hydrogen:	OverlapAtomAtom_1S2SP<Oxygen, Hydrogen>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Helium:		OverlapAtomAtom_1S2SP<Oxygen, Helium>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Lithium:	OverlapAtomAtom_2SP2SP<Oxygen, Lithium>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Beryllium:  OverlapAtomAtom_2SP2SP<Oxygen, Beryllium>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Boron:		OverlapAtomAtom_2SP2SP<Oxygen, Boron>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Carbon:		OverlapAtomAtom_2SP2SP<Oxygen, Carbon>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Nitrogen:	OverlapAtomAtom_2SP2SP<Oxygen, Nitrogen>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Oxygen:		OverlapAtomAtom_2SP2SP<Oxygen, Oxygen>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Flourine:   OverlapAtomAtom_2SP2SP<Oxygen, Flourine>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Neon:		OverlapAtomAtom_2SP2SP<Oxygen, Neon>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		}
	}
	case ATOM_TYPE::Flourine:
	{
		switch (atom_2.type)
		{
		case ATOM_TYPE::Hydrogen:	OverlapAtomAtom_1S2SP<Flourine, Hydrogen>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Helium:		OverlapAtomAtom_1S2SP<Flourine, Helium>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Lithium:	OverlapAtomAtom_2SP2SP<Flourine, Lithium>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Beryllium:  OverlapAtomAtom_2SP2SP<Flourine, Beryllium>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Boron:		OverlapAtomAtom_2SP2SP<Flourine, Boron>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Carbon:		OverlapAtomAtom_2SP2SP<Flourine, Carbon>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Nitrogen:	OverlapAtomAtom_2SP2SP<Flourine, Nitrogen>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Oxygen:		OverlapAtomAtom_2SP2SP<Flourine, Oxygen>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Flourine:   OverlapAtomAtom_2SP2SP<Flourine, Flourine>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Neon:		OverlapAtomAtom_2SP2SP<Flourine, Neon>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		}
	}
	case ATOM_TYPE::Neon:
	{
		switch (atom_2.type)
		{
		case ATOM_TYPE::Hydrogen:	OverlapAtomAtom_1S2SP<Neon, Hydrogen>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Helium:		OverlapAtomAtom_1S2SP<Neon, Helium>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Lithium:	OverlapAtomAtom_2SP2SP<Neon, Lithium>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Beryllium:  OverlapAtomAtom_2SP2SP<Neon, Beryllium>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Boron:		OverlapAtomAtom_2SP2SP<Neon, Boron>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Carbon:		OverlapAtomAtom_2SP2SP<Neon, Carbon>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Nitrogen:	OverlapAtomAtom_2SP2SP<Neon, Nitrogen>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Oxygen:		OverlapAtomAtom_2SP2SP<Neon, Oxygen>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Flourine:   OverlapAtomAtom_2SP2SP<Neon, Flourine>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Neon:		OverlapAtomAtom_2SP2SP<Neon, Neon>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		}
	}
	}
}



static void OverlapWithSelf(ATOM_TYPE type, unsigned int topLeftIndex, MatrixXd& overlapMatrix)
{
	switch (type)
	{
	case ATOM_TYPE::Hydrogen:	OverlapWithSelf_1S<Hydrogen>(topLeftIndex, overlapMatrix); return;
	case ATOM_TYPE::Helium:		OverlapWithSelf_1S<Helium>(topLeftIndex, overlapMatrix); return;
	case ATOM_TYPE::Lithium:	OverlapWithSelf_2SP<Lithium>(topLeftIndex, overlapMatrix); return;
	case ATOM_TYPE::Beryllium:  OverlapWithSelf_2SP<Beryllium>(topLeftIndex, overlapMatrix); return;
	case ATOM_TYPE::Boron:		OverlapWithSelf_2SP<Boron>(topLeftIndex, overlapMatrix); return;
	case ATOM_TYPE::Carbon:		OverlapWithSelf_2SP<Carbon>(topLeftIndex, overlapMatrix); return;
	case ATOM_TYPE::Nitrogen:	OverlapWithSelf_2SP<Nitrogen>(topLeftIndex, overlapMatrix); return;
	case ATOM_TYPE::Oxygen:		OverlapWithSelf_2SP<Oxygen>(topLeftIndex, overlapMatrix); return;
	case ATOM_TYPE::Flourine:   OverlapWithSelf_2SP<Flourine>(topLeftIndex, overlapMatrix); return;
	case ATOM_TYPE::Neon:		OverlapWithSelf_2SP<Neon>(topLeftIndex, overlapMatrix); return;
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