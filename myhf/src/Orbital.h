#pragma once

#include "MathUtils.h"
#include "Vec.h"
#include "QuantumNumbers.h"

namespace myhf
{
struct Orbital
{
	Orbital() noexcept = default;
	Orbital(const Orbital&) noexcept = default;
	Orbital(Orbital&&) noexcept = default;
	Orbital& operator=(const Orbital&) noexcept = default;
	Orbital& operator=(Orbital&&) noexcept = default;
	virtual ~Orbital() noexcept {}

	virtual char AtomicOrbital() const noexcept = 0;
	virtual double operator()(const Vec3d& r, const Vec3d& center) const noexcept = 0;
	virtual Vec3d GetGradient(const Vec3d& r, const Vec3d& center) const noexcept = 0;
	virtual double GetLaplacian(const Vec3d& r, const Vec3d& center) const noexcept = 0;
};

struct PrimitiveGaussian
{
	double alpha;
	double coefficient;
	double normalizationFactor;
	double coeffProdNorm;

	constexpr Vec3d ProductCenter(const Vec3d& center1, const Vec3d& center2, const PrimitiveGaussian& other) const noexcept
	{
		return (alpha * center1 + other.alpha * center2) / (alpha + other.alpha);
	}

	double operator()(const Vec3d& r, const Vec3d& center, const QuantumNumbers& angularMomentum) const noexcept;
	Vec3d GetGradient(const Vec3d& r, const Vec3d& center, const QuantumNumbers& angularMomentum) const noexcept;
	double GetLaplacian(const Vec3d& r, const Vec3d& center, const QuantumNumbers& angularMomentum) const noexcept;

	inline void Normalize(const QuantumNumbers& angularMomentum)
	{
		normalizationFactor = GetNormalizationFactor(angularMomentum);
		coeffProdNorm = coefficient * normalizationFactor;
	}

protected:
	inline double GetNormalizationFactor(const QuantumNumbers& angularMomentum) const noexcept
	{
		return std::pow(2. * alpha / std::numbers::pi, 3. / 4.) *
			std::pow(4. * alpha, angularMomentum.AngularMomentum() / 2.) /
			std::sqrt( 
				MathUtils::DoubleFactorial(2 * angularMomentum.l - 1) *  
				MathUtils::DoubleFactorial(2 * angularMomentum.m - 1) *  
				MathUtils::DoubleFactorial(2 * angularMomentum.n - 1) 
			);
	}
};

// the contained gaussian orbitals all have the same center and quantum numbers, whence the derivation from Orbital
struct ContractedGaussianOrbital : public Orbital
{
	std::vector<PrimitiveGaussian> gaussianOrbitals;
	QuantumNumbers angularMomentum{};

	ContractedGaussianOrbital(std::vector<PrimitiveGaussian>&& orbitals, const QuantumNumbers& _angularMomentum) noexcept :
		gaussianOrbitals(std::move(orbitals)), 
		angularMomentum(_angularMomentum)
	{
		assert(gaussianOrbitals.size() > 0);
		Normalize();
	}
	ContractedGaussianOrbital(const ContractedGaussianOrbital&) noexcept = default; // No need to normalize, because that has already been done in the rhs object
	ContractedGaussianOrbital(ContractedGaussianOrbital&&) noexcept = default;
	ContractedGaussianOrbital& operator=(const ContractedGaussianOrbital&) noexcept = default;
	ContractedGaussianOrbital& operator=(ContractedGaussianOrbital&&) noexcept = default;
	virtual ~ContractedGaussianOrbital() noexcept override {}

	char AtomicOrbital() const noexcept override { return angularMomentum.AtomicOrbital(); }
	double operator()(const Vec3d& r, const Vec3d& center) const noexcept override;
	Vec3d GetGradient(const Vec3d& r, const Vec3d& center) const noexcept override;
	double GetLaplacian(const Vec3d& r, const Vec3d& center) const noexcept override;

protected:
	void Normalize() noexcept
	{
		for (auto& orbital : gaussianOrbitals)
			orbital.Normalize(angularMomentum);
	}
};


}