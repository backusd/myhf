#pragma once

#include "Basis.h"
#include "MathUtils.h"
#include "QuantumNumbers.h"
#include "Vec.h"

#include "gcem.hpp"

namespace myhf
{
namespace test
{

namespace test2
{
static constexpr auto basis_str_STO_3G =
"#----------------------------------------------------------------------\
# Basis Set Exchange													\
# Version 0.10															\
# https://www.basissetexchange.org										\
#----------------------------------------------------------------------	\
#   Basis set: STO-3G													\
# Description: STO-3G Minimal Basis (3 functions/AO)					\
#        Role: orbital													\
#     Version: 1  (Data from Gaussian09)								\
#----------------------------------------------------------------------	\
																		\
																		\
BASIS 'ao basis' SPHERICAL PRINT										\
#BASIS SET: (3s) -> [1s]												\
H    S																	\
      0.3425250914E+01       0.1543289673E+00							\
      0.6239137298E+00       0.5353281423E+00							\
      0.1688554040E+00       0.4446345422E+00							\
#BASIS SET: (3s) -> [1s]												\
He    S																	\
      0.6362421394E+01       0.1543289673E+00							\
      0.1158922999E+01       0.5353281423E+00							\
      0.3136497915E+00       0.4446345422E+00							\
#BASIS SET: (6s,3p) -> [2s,1p]											\
Li    S																	\
      0.1611957475E+02       0.1543289673E+00							\
      0.2936200663E+01       0.5353281423E+00							\
      0.7946504870E+00       0.4446345422E+00							\
Li    SP																\
      0.6362897469E+00      -0.9996722919E-01       0.1559162750E+00	\
      0.1478600533E+00       0.3995128261E+00       0.6076837186E+00	\
      0.4808867840E-01       0.7001154689E+00       0.3919573931E+00	\
#BASIS SET: (6s,3p) -> [2s,1p]											\
Be    S																	\
      0.3016787069E+02       0.1543289673E+00							\
      0.5495115306E+01       0.5353281423E+00							\
      0.1487192653E+01       0.4446345422E+00							\
Be    SP																\
      0.1314833110E+01      -0.9996722919E-01       0.1559162750E+00	\
      0.3055389383E+00       0.3995128261E+00       0.6076837186E+00	\
      0.9937074560E-01       0.7001154689E+00       0.3919573931E+00	\
#BASIS SET: (6s,3p) -> [2s,1p]											\
B    S																	\
      0.4879111318E+02       0.1543289673E+00							\
      0.8887362172E+01       0.5353281423E+00							\
      0.2405267040E+01       0.4446345422E+00							\
B    SP																	\
      0.2236956142E+01      -0.9996722919E-01       0.1559162750E+00	\
      0.5198204999E+00       0.3995128261E+00       0.6076837186E+00	\
      0.1690617600E+00       0.7001154689E+00       0.3919573931E+00	\
#BASIS SET: (6s,3p) -> [2s,1p]											\
C    S																	\
      0.7161683735E+02       0.1543289673E+00							\
      0.1304509632E+02       0.5353281423E+00							\
      0.3530512160E+01       0.4446345422E+00							\
C    SP																	\
      0.2941249355E+01      -0.9996722919E-01       0.1559162750E+00	\
      0.6834830964E+00       0.3995128261E+00       0.6076837186E+00	\
      0.2222899159E+00       0.7001154689E+00       0.3919573931E+00	\
#BASIS SET: (6s,3p) -> [2s,1p]											\
N    S																	\
      0.9910616896E+02       0.1543289673E+00							\
      0.1805231239E+02       0.5353281423E+00							\
      0.4885660238E+01       0.4446345422E+00							\
N    SP																	\
      0.3780455879E+01      -0.9996722919E-01       0.1559162750E+00	\
      0.8784966449E+00       0.3995128261E+00       0.6076837186E+00	\
      0.2857143744E+00       0.7001154689E+00       0.3919573931E+00	\
#BASIS SET: (6s,3p) -> [2s,1p]											\
O    S																	\
      0.1307093214E+03       0.1543289673E+00							\
      0.2380886605E+02       0.5353281423E+00							\
      0.6443608313E+01       0.4446345422E+00							\
O    SP																	\
      0.5033151319E+01      -0.9996722919E-01       0.1559162750E+00	\
      0.1169596125E+01       0.3995128261E+00       0.6076837186E+00	\
      0.3803889600E+00       0.7001154689E+00       0.3919573931E+00	\
#BASIS SET: (6s,3p) -> [2s,1p]											\
F    S																	\
      0.1666791340E+03       0.1543289673E+00							\
      0.3036081233E+02       0.5353281423E+00							\
      0.8216820672E+01       0.4446345422E+00							\
F    SP																	\
      0.6464803249E+01      -0.9996722919E-01       0.1559162750E+00	\
      0.1502281245E+01       0.3995128261E+00       0.6076837186E+00	\
      0.4885884864E+00       0.7001154689E+00       0.3919573931E+00	\
#BASIS SET: (6s,3p) -> [2s,1p]											\
Ne    S																	\
      0.2070156070E+03       0.1543289673E+00							\
      0.3770815124E+02       0.5353281423E+00							\
      0.1020529731E+02       0.4446345422E+00							\
Ne    SP																\
      0.8246315120E+01      -0.9996722919E-01       0.1559162750E+00	\
      0.1916266291E+01       0.3995128261E+00       0.6076837186E+00	\
      0.6232292721E+00       0.7001154689E+00       0.3919573931E+00	\
#BASIS SET: (9s,6p) -> [3s,2p]											\
Na    S																	\
      0.2507724300E+03       0.1543289673E+00							\
      0.4567851117E+02       0.5353281423E+00							\
      0.1236238776E+02       0.4446345422E+00							\
Na    SP																\
      0.1204019274E+02      -0.9996722919E-01       0.1559162750E+00	\
      0.2797881859E+01       0.3995128261E+00       0.6076837186E+00	\
      0.9099580170E+00       0.7001154689E+00       0.3919573931E+00	\
Na    SP																\
      0.1478740622E+01      -0.2196203690E+00       0.1058760429E-01	\
      0.4125648801E+00       0.2255954336E+00       0.5951670053E+00	\
      0.1614750979E+00       0.9003984260E+00       0.4620010120E+00	\
#BASIS SET: (9s,6p) -> [3s,2p]											\
Mg    S																	\
      0.2992374137E+03       0.1543289673E+00							\
      0.5450646845E+02       0.5353281423E+00							\
      0.1475157752E+02       0.4446345422E+00							\
Mg    SP																\
      0.1512182352E+02      -0.9996722919E-01       0.1559162750E+00	\
      0.3513986579E+01       0.3995128261E+00       0.6076837186E+00	\
      0.1142857498E+01       0.7001154689E+00       0.3919573931E+00	\
Mg    SP																\
      0.1395448293E+01      -0.2196203690E+00       0.1058760429E-01	\
      0.3893265318E+00       0.2255954336E+00       0.5951670053E+00	\
      0.1523797659E+00       0.9003984260E+00       0.4620010120E+00	\
#BASIS SET: (9s,6p) -> [3s,2p]											\
Al    S																	\
      0.3514214767E+03       0.1543289673E+00							\
      0.6401186067E+02       0.5353281423E+00							\
      0.1732410761E+02       0.4446345422E+00							\
Al    SP																\
      0.1889939621E+02      -0.9996722919E-01       0.1559162750E+00	\
      0.4391813233E+01       0.3995128261E+00       0.6076837186E+00	\
      0.1428353970E+01       0.7001154689E+00       0.3919573931E+00	\
Al    SP																\
      0.1395448293E+01      -0.2196203690E+00       0.1058760429E-01	\
      0.3893265318E+00       0.2255954336E+00       0.5951670053E+00	\
      0.1523797659E+00       0.9003984260E+00       0.4620010120E+00	\
#BASIS SET: (9s,6p) -> [3s,2p]											\
Si    S																	\
      0.4077975514E+03       0.1543289673E+00							\
      0.7428083305E+02       0.5353281423E+00							\
      0.2010329229E+02       0.4446345422E+00							\
Si    SP																\
      0.2319365606E+02      -0.9996722919E-01       0.1559162750E+00	\
      0.5389706871E+01       0.3995128261E+00       0.6076837186E+00	\
      0.1752899952E+01       0.7001154689E+00       0.3919573931E+00	\
Si    SP																\
      0.1478740622E+01      -0.2196203690E+00       0.1058760429E-01	\
      0.4125648801E+00       0.2255954336E+00       0.5951670053E+00	\
      0.1614750979E+00       0.9003984260E+00       0.4620010120E+00	\
#BASIS SET: (9s,6p) -> [3s,2p]											\
P    S																	\
      0.4683656378E+03       0.1543289673E+00							\
      0.8531338559E+02       0.5353281423E+00							\
      0.2308913156E+02       0.4446345422E+00							\
P    SP																	\
      0.2803263958E+02      -0.9996722919E-01       0.1559162750E+00	\
      0.6514182577E+01       0.3995128261E+00       0.6076837186E+00	\
      0.2118614352E+01       0.7001154689E+00       0.3919573931E+00	\
P    SP																	\
      0.1743103231E+01      -0.2196203690E+00       0.1058760429E-01	\
      0.4863213771E+00       0.2255954336E+00       0.5951670053E+00	\
      0.1903428909E+00       0.9003984260E+00       0.4620010120E+00	\
#BASIS SET: (9s,6p) -> [3s,2p]											\
S    S																	\
      0.5331257359E+03       0.1543289673E+00							\
      0.9710951830E+02       0.5353281423E+00							\
      0.2628162542E+02       0.4446345422E+00							\
S    SP																	\
      0.3332975173E+02      -0.9996722919E-01       0.1559162750E+00	\
      0.7745117521E+01       0.3995128261E+00       0.6076837186E+00	\
      0.2518952599E+01       0.7001154689E+00       0.3919573931E+00	\
S    SP																	\
      0.2029194274E+01      -0.2196203690E+00       0.1058760429E-01	\
      0.5661400518E+00       0.2255954336E+00       0.5951670053E+00	\
      0.2215833792E+00       0.9003984260E+00       0.4620010120E+00	\
#BASIS SET: (9s,6p) -> [3s,2p]											\
Cl    S																	\
      0.6013456136E+03       0.1543289673E+00							\
      0.1095358542E+03       0.5353281423E+00							\
      0.2964467686E+02       0.4446345422E+00							\
Cl    SP																\
      0.3896041889E+02      -0.9996722919E-01       0.1559162750E+00	\
      0.9053563477E+01       0.3995128261E+00       0.6076837186E+00	\
      0.2944499834E+01       0.7001154689E+00       0.3919573931E+00	\
Cl    SP																\
      0.2129386495E+01      -0.2196203690E+00       0.1058760429E-01	\
      0.5940934274E+00       0.2255954336E+00       0.5951670053E+00	\
      0.2325241410E+00       0.9003984260E+00       0.4620010120E+00	\
#BASIS SET: (9s,6p) -> [3s,2p]											\
Ar    S																	\
      0.6744465184E+03       0.1543289673E+00							\
      0.1228512753E+03       0.5353281423E+00							\
      0.3324834945E+02       0.4446345422E+00							\
Ar    SP																\
      0.4516424392E+02      -0.9996722919E-01       0.1559162750E+00	\
      0.1049519900E+02       0.3995128261E+00       0.6076837186E+00	\
      0.3413364448E+01       0.7001154689E+00       0.3919573931E+00	\
Ar    SP																\
      0.2621366518E+01      -0.2196203690E+00       0.1058760429E-01	\
      0.7313546050E+00       0.2255954336E+00       0.5951670053E+00	\
      0.2862472356E+00       0.9003984260E+00       0.4620010120E+00	\
END";
static constexpr auto basis_str_STO_6G = " ";


//struct GenericPrimitiveGaussian
//{
//	double alpha;
//	double coefficient;
//};
struct GenericContractedGaussianOrbital
{
	std::vector<double> primitiveData;
	QuantumNumbers angularMomentum{};
};
struct GenericShell
{
	std::vector<GenericContractedGaussianOrbital> functions;
};
struct GenericAtom
{
	std::vector<GenericShell> shells;
};
struct GenericBasis
{
	constexpr GenericBasis(const char* basisStr) noexcept :
		atoms{}
	{
		//GenericPrimitiveGaussian pg1{ 0.3425250914E+01, 0.1543289673E+00 };
		//GenericPrimitiveGaussian pg2{ 0.6239137298E+00, 0.5353281423E+00 };
		//GenericPrimitiveGaussian pg3{ 0.1688554040E+00, 0.4446345422E+00 };

		std::vector<double> data{ 0.3425250914E+01, 0.1543289673E+00, 0.6239137298E+00, 0.5353281423E+00, 0.1688554040E+00, 0.4446345422E+00 };

		GenericContractedGaussianOrbital cgo{ data, {0, 0, 0} };

		GenericShell shell{ {cgo} };

		GenericAtom atom{ {shell} };

		atoms.push_back(atom);
	}

	constexpr std::span<double> GetContractedGaussianData(size_t atomIndex, size_t shellIndex, size_t functionIndex)
	{
		return atoms[atomIndex].shells[shellIndex].functions[functionIndex].primitiveData;
	}

	std::vector<GenericAtom> atoms;
};








struct PrimitiveGaussian
{
public:
	constexpr PrimitiveGaussian() noexcept {};
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

template<size_t N>
struct ContractedGaussian
{
	constexpr ContractedGaussian(std::span<double> values, const QuantumNumbers& _angularMomentum) noexcept :
		primitiveGaussians{},
		angularMomentum(_angularMomentum)
	{
		for (size_t iii = 0; iii < N; ++iii)
			primitiveGaussians[iii] = PrimitiveGaussian(values[2 * iii], values[2 * iii + 1], angularMomentum);
	}

	std::array<PrimitiveGaussian, N> primitiveGaussians;
	QuantumNumbers angularMomentum{};
};
template<size_t NumFunctions, size_t... NPrimitives>
struct Shell
{
	static_assert(false);
};

template <>
struct Shell<1, 3>
{
	constexpr Shell(std::span<double> values, const QuantumNumbers& angularMomentum) :
		func1(values, angularMomentum) {}
	ContractedGaussian<3> func1;
};

template <>
struct Shell<1, 6>
{
	constexpr Shell(std::span<double> values, const QuantumNumbers& angularMomentum) :
		func1(values, angularMomentum) {}
	ContractedGaussian<6> func1;
};

template<size_t N>
struct STO_NG_Basis
{
	static consteval const char* GetBasisString() noexcept
	{
		if constexpr (N == 3)
			return basis_str_STO_3G;
		else if constexpr (N == 6)
			return basis_str_STO_6G;
		else
			static_assert(false);
	}
	//	static consteval const GenericBasis& GetGenericBasis()
	//	{
	//		static GenericBasis basis(GetBasisString());
	//		return basis;
	//	}

//	static constexpr GenericBasis m_basis = GenericBasis(GetBasisString());

	struct Hydrogen
	{
		static constexpr unsigned int NumberOfContractedGaussians = 1;
		static constexpr Shell<1, N> orbital_1s = Shell<1, N>(GenericBasis(GetBasisString()).GetContractedGaussianData(0, 0, 0), {0, 0, 0});
	};
};

using STO_3G = STO_NG_Basis<3>;
using STO_6G = STO_NG_Basis<6>;









} // namespace test2

















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
struct Helium
{
	static constexpr unsigned int NumberOfContractedGaussians = 1;
	static constexpr ContractedGaussian orbital_1s = ContractedGaussian(
		0.6362421394E+01, 0.1543289673E+00,
		0.1158922999E+01, 0.5353281423E+00,
		0.3136497915E+00, 0.4446345422E+00,
		{ 0, 0, 0 });
};
struct Lithium
{
	static constexpr unsigned int NumberOfContractedGaussians = 5;

	static constexpr ContractedGaussian orbital_1s = ContractedGaussian(
		0.1611957475E+02, 0.1543289673E+00,
		0.2936200663E+01, 0.5353281423E+00,
		0.7946504870E+00, 0.4446345422E+00,
		{ 0, 0, 0 });

	static constexpr ContractedGaussian orbital_2s = ContractedGaussian(
		0.6362897469E+00, -0.9996722919E-01,
		0.1478600533E+00,  0.3995128261E+00,
		0.4808867840E-01,  0.7001154689E+00,
		{ 0, 0, 0 });

	static constexpr ContractedGaussian orbital_2px = ContractedGaussian(
		0.6362897469E+00, 0.1559162750E+00,
		0.1478600533E+00, 0.6076837186E+00,
		0.4808867840E-01, 0.3919573931E+00,
		{ 1, 0, 0 });

	static constexpr ContractedGaussian orbital_2py = ContractedGaussian(
		0.6362897469E+00, 0.1559162750E+00,
		0.1478600533E+00, 0.6076837186E+00,
		0.4808867840E-01, 0.3919573931E+00,
		{ 0, 1, 0 });

	static constexpr ContractedGaussian orbital_2pz = ContractedGaussian(
		0.6362897469E+00, 0.1559162750E+00,
		0.1478600533E+00, 0.6076837186E+00,
		0.4808867840E-01, 0.3919573931E+00,
		{ 0, 0, 1 });
};
struct Beryllium
{
	static constexpr unsigned int NumberOfContractedGaussians = 5;

	static constexpr ContractedGaussian orbital_1s = ContractedGaussian(
		0.3016787069E+02, 0.1543289673E+00,
		0.5495115306E+01, 0.5353281423E+00,
		0.1487192653E+01, 0.4446345422E+00,
		{ 0, 0, 0 });

	static constexpr ContractedGaussian orbital_2s = ContractedGaussian(
		0.1314833110E+01, -0.9996722919E-01,
		0.3055389383E+00, 0.3995128261E+00,
		0.9937074560E-01, 0.7001154689E+00, 
		{ 0, 0, 0 });

	static constexpr ContractedGaussian orbital_2px = ContractedGaussian(
		0.1314833110E+01, 0.1559162750E+00,
		0.3055389383E+00, 0.6076837186E+00,
		0.9937074560E-01, 0.3919573931E+00,
		{ 1, 0, 0 });

	static constexpr ContractedGaussian orbital_2py = ContractedGaussian(
		0.1314833110E+01, 0.1559162750E+00,
		0.3055389383E+00, 0.6076837186E+00,
		0.9937074560E-01, 0.3919573931E+00,
		{ 0, 1, 0 });

	static constexpr ContractedGaussian orbital_2pz = ContractedGaussian(
		0.1314833110E+01, 0.1559162750E+00,
		0.3055389383E+00, 0.6076837186E+00,
		0.9937074560E-01, 0.3919573931E+00,
		{ 0, 0, 1 });
};
struct Boron
{
	static constexpr unsigned int NumberOfContractedGaussians = 5;

	static constexpr ContractedGaussian orbital_1s = ContractedGaussian(
		0.4879111318E+02, 0.1543289673E+00,
		0.8887362172E+01, 0.5353281423E+00,
		0.2405267040E+01, 0.4446345422E+00,
		{ 0, 0, 0 });

	static constexpr ContractedGaussian orbital_2s = ContractedGaussian(
		0.2236956142E+01, -0.9996722919E-01,
		0.5198204999E+00, 0.3995128261E+00,
		0.1690617600E+00, 0.7001154689E+00, 
		{ 0, 0, 0 });

	static constexpr ContractedGaussian orbital_2px = ContractedGaussian(
		0.2236956142E+01, 0.1559162750E+00,
		0.5198204999E+00, 0.6076837186E+00,
		0.1690617600E+00, 0.3919573931E+00,
		{ 1, 0, 0 });

	static constexpr ContractedGaussian orbital_2py = ContractedGaussian(
		0.2236956142E+01, 0.1559162750E+00,
		0.5198204999E+00, 0.6076837186E+00,
		0.1690617600E+00, 0.3919573931E+00,
		{ 0, 1, 0 });

	static constexpr ContractedGaussian orbital_2pz = ContractedGaussian(
		0.2236956142E+01, 0.1559162750E+00,
		0.5198204999E+00, 0.6076837186E+00,
		0.1690617600E+00, 0.3919573931E+00,
		{ 0, 0, 1 });
};
struct Carbon
{
	static constexpr unsigned int NumberOfContractedGaussians = 5;

	static constexpr ContractedGaussian orbital_1s = ContractedGaussian(
		0.7161683735E+02, 0.1543289673E+00,
		0.1304509632E+02, 0.5353281423E+00,
		0.3530512160E+01, 0.4446345422E+00,
		{ 0, 0, 0 });

	static constexpr ContractedGaussian orbital_2s = ContractedGaussian(
		0.2941249355E+01, -0.9996722919E-01,
		0.6834830964E+00, 0.3995128261E+00,
		0.2222899159E+00, 0.7001154689E+00, 
		{ 0, 0, 0 });

	static constexpr ContractedGaussian orbital_2px = ContractedGaussian(
		0.2941249355E+01, 0.1559162750E+00,
		0.6834830964E+00, 0.6076837186E+00,
		0.2222899159E+00, 0.3919573931E+00,
		{ 1, 0, 0 });

	static constexpr ContractedGaussian orbital_2py = ContractedGaussian(
		0.2941249355E+01, 0.1559162750E+00,
		0.6834830964E+00, 0.6076837186E+00,
		0.2222899159E+00, 0.3919573931E+00,
		{ 0, 1, 0 });

	static constexpr ContractedGaussian orbital_2pz = ContractedGaussian(
		0.2941249355E+01, 0.1559162750E+00,
		0.6834830964E+00, 0.6076837186E+00,
		0.2222899159E+00, 0.3919573931E+00,
		{ 0, 0, 1 });
};
struct Nitrogen
{
	static constexpr unsigned int NumberOfContractedGaussians = 5;

	static constexpr ContractedGaussian orbital_1s = ContractedGaussian(
		0.9910616896E+02, 0.1543289673E+00,
		0.1805231239E+02, 0.5353281423E+00,
		0.4885660238E+01, 0.4446345422E+00,
		{ 0, 0, 0 });

	static constexpr ContractedGaussian orbital_2s = ContractedGaussian(
		0.3780455879E+01, -0.9996722919E-01,
		0.8784966449E+00,  0.3995128261E+00,
		0.2857143744E+00,  0.7001154689E+00,
		{ 0, 0, 0 });

	static constexpr ContractedGaussian orbital_2px = ContractedGaussian(
		0.3780455879E+01, 0.1559162750E+00,
		0.8784966449E+00, 0.6076837186E+00,
		0.2857143744E+00, 0.3919573931E+00,
		{ 1, 0, 0 });

	static constexpr ContractedGaussian orbital_2py = ContractedGaussian(
		0.3780455879E+01, 0.1559162750E+00,
		0.8784966449E+00, 0.6076837186E+00,
		0.2857143744E+00, 0.3919573931E+00,
		{ 0, 1, 0 });

	static constexpr ContractedGaussian orbital_2pz = ContractedGaussian(
		0.3780455879E+01, 0.1559162750E+00,
		0.8784966449E+00, 0.6076837186E+00,
		0.2857143744E+00, 0.3919573931E+00,
		{ 0, 0, 1 });
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
struct Flourine
{
	static constexpr unsigned int NumberOfContractedGaussians = 5;

	static constexpr ContractedGaussian orbital_1s = ContractedGaussian(
		0.1666791340E+03, 0.1543289673E+00,
		0.3036081233E+02, 0.5353281423E+00,
		0.8216820672E+01, 0.4446345422E+00,
		{ 0, 0, 0 });

	static constexpr ContractedGaussian orbital_2s = ContractedGaussian(
		0.6464803249E+01, -0.9996722919E-01,
		0.1502281245E+01, 0.3995128261E+00,
		0.4885884864E+00, 0.7001154689E+00, 
		{ 0, 0, 0 });

	static constexpr ContractedGaussian orbital_2px = ContractedGaussian(
		0.6464803249E+01, 0.1559162750E+00,
		0.1502281245E+01, 0.6076837186E+00,
		0.4885884864E+00, 0.3919573931E+00,
		{ 1, 0, 0 });

	static constexpr ContractedGaussian orbital_2py = ContractedGaussian(
		0.6464803249E+01, 0.1559162750E+00,
		0.1502281245E+01, 0.6076837186E+00,
		0.4885884864E+00, 0.3919573931E+00,
		{ 0, 1, 0 });

	static constexpr ContractedGaussian orbital_2pz = ContractedGaussian(
		0.6464803249E+01, 0.1559162750E+00,
		0.1502281245E+01, 0.6076837186E+00,
		0.4885884864E+00, 0.3919573931E+00,
		{ 0, 0, 1 });
};
struct Neon
{
	static constexpr unsigned int NumberOfContractedGaussians = 5;

	static constexpr ContractedGaussian orbital_1s = ContractedGaussian(
		0.2070156070E+03, 0.1543289673E+00,
		0.3770815124E+02, 0.5353281423E+00,
		0.1020529731E+02, 0.4446345422E+00,
		{ 0, 0, 0 });

	static constexpr ContractedGaussian orbital_2s = ContractedGaussian(
		0.8246315120E+01, -0.9996722919E-01,
		0.1916266291E+01, 0.3995128261E+00,
		0.6232292721E+00, 0.7001154689E+00, 
		{ 0, 0, 0 });

	static constexpr ContractedGaussian orbital_2px = ContractedGaussian(
		0.8246315120E+01, 0.1559162750E+00,
		0.1916266291E+01, 0.6076837186E+00,
		0.6232292721E+00, 0.3919573931E+00,
		{ 1, 0, 0 });

	static constexpr ContractedGaussian orbital_2py = ContractedGaussian(
		0.8246315120E+01, 0.1559162750E+00,
		0.1916266291E+01, 0.6076837186E+00,
		0.6232292721E+00, 0.3919573931E+00,
		{ 0, 1, 0 });

	static constexpr ContractedGaussian orbital_2pz = ContractedGaussian(
		0.8246315120E+01, 0.1559162750E+00,
		0.1916266291E+01, 0.6076837186E+00,
		0.6232292721E+00, 0.3919573931E+00,
		{ 0, 0, 1 });
};

} // namespace test
} // namespace myhf