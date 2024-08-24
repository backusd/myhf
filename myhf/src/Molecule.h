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
	static double OverlapOfTwoOrbitals(const BasisAtom& atom1, const ContractedGaussianOrbital& orbital1, const Vec3d& position1, const BasisAtom& atom2, const ContractedGaussianOrbital& orbital2, const Vec3d& position2) noexcept;
	static double OverlapOfTwoPrimitiveGaussians(double alpha1, double alpha2, const Vec3d& position1, const Vec3d& position2, const QuantumNumbers& angularMomentum1, const QuantumNumbers& angularMomentum2) noexcept;
};
}

template <>
struct std::formatter<myhf::Molecule>
{
	constexpr auto parse(std::format_parse_context& ctx) noexcept { return ctx.begin(); }
	auto format(const myhf::Molecule& molec, std::format_context& ctx) const
	{
		std::stringstream oss;
		oss << "\tBasis: " << molec.basis.name << '\n';
		for (const auto& atom : molec.atoms)
			oss << '\t' << ToStringShortName(atom.type) << ' ' << atom.position << '\n';
		return std::format_to(ctx.out(), "{}", oss.str());
	}
};