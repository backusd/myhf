#pragma once

#include "Atom.h"
#include "Basis.h"


namespace myhf
{
class Molecule
{
public:
	std::vector<Atom> atoms;
	const Basis& basis;

	Molecule(const std::vector<Atom>& _atoms, const Basis& _basis = STO_3G) noexcept :
		atoms(_atoms), basis(_basis)
	{}
	Molecule(std::vector<Atom>&& _atoms, const Basis& _basis = STO_3G) noexcept :
		atoms(std::move(_atoms)), basis(_basis)
	{}
	Molecule(const Molecule&) noexcept = default;
	Molecule(Molecule&&) noexcept = default;
	Molecule& operator=(const Molecule&) noexcept = default;
	Molecule& operator=(Molecule&&) noexcept = default;


	Eigen::MatrixXd OverlapMatrix() noexcept;

private:
	double GetOverlap(const BasisAtom& atom1, const ContractedGaussianOrbital& orbital1, const Vec3d& position1, const BasisAtom& atom2, const ContractedGaussianOrbital& orbital2, const Vec3d& position2) noexcept;
};
}