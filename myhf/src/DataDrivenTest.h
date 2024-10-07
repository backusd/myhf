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
public:
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

struct ContractedGaussian
{
	constexpr ContractedGaussian(double a1, double c1, double a2, double c2, double a3, double c3, const QuantumNumbers& angularMomentum) noexcept :
		primitiveGaussians{ 
			PrimitiveGaussian(a1, c1, angularMomentum), 
			PrimitiveGaussian(a2, c2, angularMomentum), 
			PrimitiveGaussian(a3, c3, angularMomentum) 
		},
		angularMomentum(angularMomentum)
	{}
	std::array<PrimitiveGaussian, 3> primitiveGaussians;
	QuantumNumbers angularMomentum{};
};

struct Hydrogen
{
	static constexpr unsigned int NumberOfContractedGaussians = 1;
	static constexpr ContractedGaussian orbital_1s = ContractedGaussian(
		0.3425250914E+01, 0.1543289673E+00, 
		0.6239137298E+00, 0.5353281423E+00,
		0.1688554040E+00, 0.4446345422E+00,
		{0, 0, 0});
};

struct Oxygen
{
	static constexpr unsigned int NumberOfContractedGaussians = 5;

	static constexpr ContractedGaussian orbital_1s = ContractedGaussian(
		0.1307093214E+03, 0.1543289673E+00,
		0.2380886605E+02, 0.5353281423E+00,
		0.6443608313E+01, 0.4446345422E+00,
		{ 0, 0, 0 });

	static constexpr ContractedGaussian orbital_2s = ContractedGaussian(
		0.5033151319E+01, -0.9996722919E-01,
		0.1169596125E+01,  0.3995128261E+00,
		0.3803889600E+00,  0.7001154689E+00,
		{ 0, 0, 0 });

	static constexpr ContractedGaussian orbital_2px = ContractedGaussian(
		0.5033151319E+01, 0.1559162750E+00,
		0.1169596125E+01, 0.6076837186E+00,
		0.3803889600E+00, 0.3919573931E+00,
		{ 1, 0, 0 });

	static constexpr ContractedGaussian orbital_2py = ContractedGaussian(
		0.5033151319E+01, 0.1559162750E+00,
		0.1169596125E+01, 0.6076837186E+00,
		0.3803889600E+00, 0.3919573931E+00,
		{ 0, 1, 0 });

	static constexpr ContractedGaussian orbital_2pz = ContractedGaussian(
		0.5033151319E+01, 0.1559162750E+00,
		0.1169596125E+01, 0.6076837186E+00,
		0.3803889600E+00, 0.3919573931E+00,
		{ 0, 0, 1 });
};

} // namespace test
} // namespace myhf

//#include "gcem.hpp"
//
//namespace myhf
//{
//namespace test
//{
//	struct PrimitiveGaussian
//	{
//		constexpr PrimitiveGaussian(double _alpha, double _coefficient, const QuantumNumbers& angularMomentum) noexcept :
//			alpha(_alpha), 
//			coefficient(_coefficient),
//			normalizationFactor(ComputeNormalizationFactor(angularMomentum))
//		{
//			coeffProdNorm = coefficient * normalizationFactor;
//		}
//		double alpha;
//		double coefficient;
//		double normalizationFactor;
//		double coeffProdNorm; // coefficient * normalizationFactor
//
////		constexpr Vec3d ProductCenter(const Vec3d& center1, const Vec3d& center2, const PrimitiveGaussian& other) const noexcept
////		{
////			return (alpha * center1 + other.alpha * center2) / (alpha + other.alpha);
////		}
////
////		double operator()(const Vec3d& r, const Vec3d& center, const QuantumNumbers& angularMomentum) const noexcept;
////		Vec3d GetGradient(const Vec3d& r, const Vec3d& center, const QuantumNumbers& angularMomentum) const noexcept;
////		double GetLaplacian(const Vec3d& r, const Vec3d& center, const QuantumNumbers& angularMomentum) const noexcept;
////
////		inline void Normalize(const QuantumNumbers& angularMomentum)
////		{
////			normalizationFactor = GetNormalizationFactor(angularMomentum);
////			coeffProdNorm = coefficient * normalizationFactor;
////		}
//
//	protected:
//		constexpr double ComputeNormalizationFactor(const QuantumNumbers& angularMomentum) const noexcept
//		{
//			return gcem::pow(2. * alpha / std::numbers::pi, 3. / 4.) *
//				   gcem::pow(4. * alpha, angularMomentum.AngularMomentum() / 2.0) /
//				   gcem::sqrt(
//				       MathUtils::DoubleFactorial(2 * angularMomentum.l - 1) *
//				       MathUtils::DoubleFactorial(2 * angularMomentum.m - 1) *
//				       MathUtils::DoubleFactorial(2 * angularMomentum.n - 1)
//				   );
//		}
//	};
//
//	template<size_t NumPrimitives>
//	struct ContractedGaussian
//	{
//		std::array<PrimitiveGaussian, NumPrimitives> primitiveGaussians;
//		QuantumNumbers angularMomentum{};
//	};
//
//	template<size_t NumOrbitals, size_t NumPrimitives>
//	struct BasisAtom
//	{
//		std::array<ContractedGaussian<NumPrimitives>, NumOrbitals> orbitals;
//	};
//
//	namespace STO_3G
//	{
//		static constexpr BasisAtom<1, 3> Hydrogen = {
//			std::array<ContractedGaussian<3>, 1>{
//				ContractedGaussian<3>{
//					std::array<PrimitiveGaussian, 3>{
//						PrimitiveGaussian(0.3425250914E+01, 0.1543289673E+00, QN_1s),
//						PrimitiveGaussian(0.6239137298E+00, 0.5353281423E+00, QN_1s),
//						PrimitiveGaussian(0.1688554040E+00, 0.4446345422E+00, QN_1s)
//					},
//					QN_1s
//				}
//			}
//		};
//
//		template<ATOM_TYPE Type>
//		constexpr const auto& GetAtom()
//		{
//			if constexpr (Type == ATOM_TYPE::Hydrogen) return STO_3G::Hydrogen;
//		}
//
//		constexpr const BasisAtom& GetAtom(ATOM_TYPE type) noexcept
//		{
//			assert(static_cast<unsigned int>(type) > 0);
//			assert(static_cast<unsigned int>(type) - 1 < atoms.size());
//			return atoms[static_cast<unsigned int>(type) - 1];
//		}
//	};
//
//	namespace STO_3G
//	{
//		static constexpr size_t PrimitiveCount = 3;
//
//		static constexpr QuantumNumbers QN_1s(0, 0, 0);
//		static constexpr QuantumNumbers QN_2s(0, 0, 0);
//		static constexpr QuantumNumbers QN_2px(1, 0, 0);
//		static constexpr QuantumNumbers QN_2py(0, 1, 0);
//		static constexpr QuantumNumbers QN_2pz(0, 0, 1);
//		
//	//	static constexpr BasisAtom<1, 3> Hydrogen = {
//	//		std::array<ContractedGaussian<3>, 1>{
//	//			ContractedGaussian<3>{
//	//				std::array<PrimitiveGaussian, 3>{
//	//					PrimitiveGaussian(0.3425250914E+01, 0.1543289673E+00, QN_1s),
//	//					PrimitiveGaussian(0.6239137298E+00, 0.5353281423E+00, QN_1s),
//	//					PrimitiveGaussian(0.1688554040E+00, 0.4446345422E+00, QN_1s)
//	//				},
//	//				QN_1s
//	//			}
//	//		}
//	//	};
//	}
//
//
//	template<size_t NumGaussians, size_t NumPrimitives>
//	struct ContractedGaussianShell
//	{
//		std::array<ContractedGaussian<NumPrimitives>, NumGaussians> contractedGaussians;
//	};
//
//	namespace STO3G
//	{
//		using S_Shell = ContractedGaussianShell<1, 3>;
//		using SP_Shell = ContractedGaussianShell<4, 3>;
//		
//		struct Atom
//		{
//			constexpr Atom() noexcept = default;
//			constexpr Atom(double posX, double posY, double posZ) noexcept : position{posX, posY, posZ} {}
//			constexpr Atom(const Vec3d& _position) noexcept : position(_position) {}
//
//			Vec3d position{};
//		};
//
//		struct Hydrogen : public Atom
//		{
//			constexpr Hydrogen() noexcept = default;
//			constexpr Hydrogen(double posX, double posY, double posZ) noexcept : Atom(posX, posY, posZ) {}
//			constexpr Hydrogen(const Vec3d& _position) noexcept : Atom(_position) {}
//
//			static constexpr S_Shell shell1{
//				std::array<ContractedGaussian<3>, 1>{
//					ContractedGaussian<3>{
//						std::array<PrimitiveGaussian, 3>{
//							PrimitiveGaussian(0.3425250914E+01, 0.1543289673E+00, QuantumNumbers(0, 0, 0)),
//							PrimitiveGaussian(0.6239137298E+00, 0.5353281423E+00, QuantumNumbers(0, 0, 0)),
//							PrimitiveGaussian(0.1688554040E+00, 0.4446345422E+00, QuantumNumbers(0, 0, 0))
//						},
//						QuantumNumbers(0, 0, 0)
//					}
//				}
//			};
//		};
//
//		struct Oxygen : public Atom
//		{
//			constexpr Oxygen() noexcept = default;
//			constexpr Oxygen(double posX, double posY, double posZ) noexcept : Atom(posX, posY, posZ) {}
//			constexpr Oxygen(const Vec3d& _position) noexcept : Atom(_position) {}
//
//			static constexpr S_Shell shell1{
//				std::array<ContractedGaussian<3>, 1>{
//					ContractedGaussian<3>{
//						std::array<PrimitiveGaussian, 3>{
//							PrimitiveGaussian(0.1307093214E+03, 0.1543289673E+00, QuantumNumbers(0, 0, 0)),
//							PrimitiveGaussian(0.2380886605E+02, 0.5353281423E+00, QuantumNumbers(0, 0, 0)),
//							PrimitiveGaussian(0.6443608313E+01, 0.4446345422E+00, QuantumNumbers(0, 0, 0))
//						},
//						QuantumNumbers(0, 0, 0)
//					}
//				}
//			};
//			static constexpr SP_Shell shell2{
//				std::array<ContractedGaussian<3>, 4>{
//					ContractedGaussian<3>{
//						std::array<PrimitiveGaussian, 3>{
//							PrimitiveGaussian(0.5033151319E+01, -0.9996722919E-01, QuantumNumbers(0, 0, 0)),
//							PrimitiveGaussian(0.1169596125E+01,  0.3995128261E+00, QuantumNumbers(0, 0, 0)),
//							PrimitiveGaussian(0.3803889600E+00,  0.7001154689E+00, QuantumNumbers(0, 0, 0))
//						},
//						QuantumNumbers(0, 0, 0)
//					},
//					ContractedGaussian<3>{
//						std::array<PrimitiveGaussian, 3>{
//							PrimitiveGaussian(0.5033151319E+01, 0.1559162750E+00, QuantumNumbers(1, 0, 0)),
//							PrimitiveGaussian(0.1169596125E+01, 0.6076837186E+00, QuantumNumbers(1, 0, 0)),
//							PrimitiveGaussian(0.3803889600E+00, 0.3919573931E+00, QuantumNumbers(1, 0, 0))
//						},
//						QuantumNumbers(1, 0, 0)
//					},
//					ContractedGaussian<3>{
//						std::array<PrimitiveGaussian, 3>{
//							PrimitiveGaussian(0.5033151319E+01, 0.1559162750E+00, QuantumNumbers(0, 1, 0)),
//							PrimitiveGaussian(0.1169596125E+01, 0.6076837186E+00, QuantumNumbers(0, 1, 0)),
//							PrimitiveGaussian(0.3803889600E+00, 0.3919573931E+00, QuantumNumbers(0, 1, 0))
//						},
//						QuantumNumbers(0, 1, 0)
//					},
//					ContractedGaussian<3>{
//						std::array<PrimitiveGaussian, 3>{
//							PrimitiveGaussian(0.5033151319E+01, 0.1559162750E+00, QuantumNumbers(0, 0, 1)),
//							PrimitiveGaussian(0.1169596125E+01, 0.6076837186E+00, QuantumNumbers(0, 0, 1)),
//							PrimitiveGaussian(0.3803889600E+00, 0.3919573931E+00, QuantumNumbers(0, 0, 1))
//						},
//						QuantumNumbers(0, 0, 1)
//					}
//				}
//			};
//		};
//	}
//
//
//
//
//
//
//
//}
//}