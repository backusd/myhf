#pragma once

#include "Orbital.h"


namespace myhf
{
// a shell has its center on the atom nucleus, it contains the basis functions with this center
struct Shell
{
	Shell() noexcept = default;
	Shell(const Shell&) noexcept = default;
	Shell(Shell&&) noexcept = default;
	Shell& operator=(const Shell&) noexcept = default;
	Shell& operator=(Shell&&) noexcept = default;
	virtual ~Shell() noexcept {}

	virtual double operator()(const Vec3d& r, const Vec3d& center) const noexcept = 0;
	virtual Vec3d GetGradient(const Vec3d& r, const Vec3d& center) const noexcept = 0;
	virtual double GetLaplacian(const Vec3d& r, const Vec3d& center) const noexcept = 0;
};

struct ContractedGaussianShell : public Shell
{
	std::vector<ContractedGaussianOrbital> basisFunctions;

	ContractedGaussianShell(std::vector<ContractedGaussianOrbital>&& orbitals) noexcept :
		basisFunctions(std::move(orbitals))
	{}
	ContractedGaussianShell(const ContractedGaussianShell&) noexcept = default;
	ContractedGaussianShell(ContractedGaussianShell&&) noexcept = default;
	ContractedGaussianShell& operator=(const ContractedGaussianShell&) noexcept = default;
	ContractedGaussianShell& operator=(ContractedGaussianShell&&) noexcept = default;

	virtual double operator()(const Vec3d& r, const Vec3d& center) const noexcept override;
	virtual Vec3d GetGradient(const Vec3d& r, const Vec3d& center) const noexcept override;
	virtual double GetLaplacian(const Vec3d& r, const Vec3d& center) const noexcept override;

	std::string GetShellString() const noexcept;

	unsigned int CountOrbitals(char orbitalChar) const noexcept;
	unsigned int CountContractedOrbitals(char orbitalChar) const noexcept;

	unsigned int CountNumberOfContractedGaussians() const noexcept;
	unsigned int CountNumberOfGaussians() const noexcept;

protected:
	static unsigned int AdjustOrbitalsCount(char orbital, unsigned int res) noexcept;
};

}
