#pragma once

#include "MathUtils.h"
#include "QuantumNumbers.h"
#include "Vec.h"


#include "gcem.hpp"

namespace myhf
{
namespace test
{
	struct PrimitiveGaussian
	{
		constexpr PrimitiveGaussian(double _alpha, double _coefficient, const QuantumNumbers& angularMomentum) noexcept :
			alpha(_alpha), 
			coefficient(_coefficient),
			normalizationFactor(ComputeNormalizationFactor(angularMomentum))
		{
			coeffProdNorm = coefficient * normalizationFactor;
		}
		double alpha;
		double coefficient;
		double normalizationFactor;
		double coeffProdNorm; // coefficient * normalizationFactor

//		constexpr Vec3d ProductCenter(const Vec3d& center1, const Vec3d& center2, const PrimitiveGaussian& other) const noexcept
//		{
//			return (alpha * center1 + other.alpha * center2) / (alpha + other.alpha);
//		}
//
//		double operator()(const Vec3d& r, const Vec3d& center, const QuantumNumbers& angularMomentum) const noexcept;
//		Vec3d GetGradient(const Vec3d& r, const Vec3d& center, const QuantumNumbers& angularMomentum) const noexcept;
//		double GetLaplacian(const Vec3d& r, const Vec3d& center, const QuantumNumbers& angularMomentum) const noexcept;
//
//		inline void Normalize(const QuantumNumbers& angularMomentum)
//		{
//			normalizationFactor = GetNormalizationFactor(angularMomentum);
//			coeffProdNorm = coefficient * normalizationFactor;
//		}

	protected:
		constexpr double ComputeNormalizationFactor(const QuantumNumbers& angularMomentum) const noexcept
		{
			return gcem::pow(2. * alpha / std::numbers::pi, 3. / 4.) *
				   gcem::pow(4. * alpha, angularMomentum.AngularMomentum() / 2.0) /
				   gcem::sqrt(
				       MathUtils::DoubleFactorial(2 * angularMomentum.l - 1) *
				       MathUtils::DoubleFactorial(2 * angularMomentum.m - 1) *
				       MathUtils::DoubleFactorial(2 * angularMomentum.n - 1)
				   );
		}
	};

	template<size_t NumPrimitives>
	struct ContractedGaussian
	{
		std::array<PrimitiveGaussian, NumPrimitives> primitiveGaussians;
		QuantumNumbers angularMomentum{};
	};

	template<size_t NumGaussians, size_t NumPrimitives>
	struct ContractedGaussianShell
	{
		std::array<ContractedGaussian<NumPrimitives>, NumGaussians> contractedGaussians;
	};

	namespace STO3G
	{
		using S_Shell = ContractedGaussianShell<1, 3>;
		using SP_Shell = ContractedGaussianShell<4, 3>;
		
		struct Atom
		{
			Vec3d position{};
		};

		struct Hydrogen : public Atom
		{
			static constexpr S_Shell shell1{
				std::array<ContractedGaussian<3>, 1>{
					ContractedGaussian<3>{
						std::array<PrimitiveGaussian, 3>{
							PrimitiveGaussian(0.3425250914E+01, 0.1543289673E+00, QuantumNumbers(0, 0, 0)),
							PrimitiveGaussian(0.6239137298E+00, 0.5353281423E+00, QuantumNumbers(0, 0, 0)),
							PrimitiveGaussian(0.1688554040E+00, 0.4446345422E+00, QuantumNumbers(0, 0, 0))
						},
						QuantumNumbers(0, 0, 0)
					}
				}
			};
		};

		struct Oxygen : public Atom
		{
			static constexpr S_Shell shell1{
				std::array<ContractedGaussian<3>, 1>{
					ContractedGaussian<3>{
						std::array<PrimitiveGaussian, 3>{
							PrimitiveGaussian(0.1307093214E+03, 0.1543289673E+00, QuantumNumbers(0, 0, 0)),
							PrimitiveGaussian(0.2380886605E+02, 0.5353281423E+00, QuantumNumbers(0, 0, 0)),
							PrimitiveGaussian(0.6443608313E+01, 0.4446345422E+00, QuantumNumbers(0, 0, 0))
						},
						QuantumNumbers(0, 0, 0)
					}
				}
			};
			static constexpr SP_Shell shell2{
				std::array<ContractedGaussian<3>, 4>{
					ContractedGaussian<3>{
						std::array<PrimitiveGaussian, 3>{
							PrimitiveGaussian(0.5033151319E+01, -0.9996722919E-01, QuantumNumbers(0, 0, 0)),
							PrimitiveGaussian(0.1169596125E+01,  0.3995128261E+00, QuantumNumbers(0, 0, 0)),
							PrimitiveGaussian(0.3803889600E+00,  0.7001154689E+00, QuantumNumbers(0, 0, 0))
						},
						QuantumNumbers(0, 0, 0)
					},
					ContractedGaussian<3>{
						std::array<PrimitiveGaussian, 3>{
							PrimitiveGaussian(0.5033151319E+01, 0.1559162750E+00, QuantumNumbers(1, 0, 0)),
							PrimitiveGaussian(0.1169596125E+01, 0.6076837186E+00, QuantumNumbers(1, 0, 0)),
							PrimitiveGaussian(0.3803889600E+00, 0.3919573931E+00, QuantumNumbers(1, 0, 0))
						},
						QuantumNumbers(1, 0, 0)
					},
					ContractedGaussian<3>{
						std::array<PrimitiveGaussian, 3>{
							PrimitiveGaussian(0.5033151319E+01, 0.1559162750E+00, QuantumNumbers(0, 1, 0)),
							PrimitiveGaussian(0.1169596125E+01, 0.6076837186E+00, QuantumNumbers(0, 1, 0)),
							PrimitiveGaussian(0.3803889600E+00, 0.3919573931E+00, QuantumNumbers(0, 1, 0))
						},
						QuantumNumbers(0, 1, 0)
					},
					ContractedGaussian<3>{
						std::array<PrimitiveGaussian, 3>{
							PrimitiveGaussian(0.5033151319E+01, 0.1559162750E+00, QuantumNumbers(0, 0, 1)),
							PrimitiveGaussian(0.1169596125E+01, 0.6076837186E+00, QuantumNumbers(0, 0, 1)),
							PrimitiveGaussian(0.3803889600E+00, 0.3919573931E+00, QuantumNumbers(0, 0, 1))
						},
						QuantumNumbers(0, 0, 1)
					}
				}
			};
		};
	}







}
}