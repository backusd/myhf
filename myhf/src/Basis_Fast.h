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