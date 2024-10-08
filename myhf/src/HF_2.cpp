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
		0,
		0,
		0,
		0,
		0,
		0,
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

// Atom-Atom Overlap Functions
template<>
void OverlapAtomAtomImpl<Hydrogen, Hydrogen>(const Vec3d& position1, const Vec3d& position2, unsigned int row, unsigned int col, MatrixXd& overlapMatrix) noexcept
{
	static constexpr std::array<double, 9> factors = GenerateOverlapFactors_1s1s<Hydrogen, Hydrogen>();
	static constexpr std::array<double, 9> expFactors = GenerateOverlapExponentialFactors_1s1s<Hydrogen, Hydrogen>();

	Vec3d diff = position2 - position1;
	double dotProduct = diff.Dot(diff);

	double result = 0.0;
	result += factors[0] * std::exp(expFactors[0] * dotProduct);
	result += factors[1] * std::exp(expFactors[1] * dotProduct);
	result += factors[2] * std::exp(expFactors[2] * dotProduct);
	result += factors[3] * std::exp(expFactors[3] * dotProduct);
	result += factors[4] * std::exp(expFactors[4] * dotProduct);
	result += factors[5] * std::exp(expFactors[5] * dotProduct);
	result += factors[6] * std::exp(expFactors[6] * dotProduct);
	result += factors[7] * std::exp(expFactors[7] * dotProduct);
	result += factors[8] * std::exp(expFactors[8] * dotProduct);

	overlapMatrix(row, col) = result;
	overlapMatrix(col, row) = result;
}
static Eigen::Array<double, 1, 5> OverlapHydrogenOxygen(const Vec3d& position1, const Vec3d& position2)
{
	Eigen::Array<double, 1, 5> results{};

	static constexpr std::array<double, 9> oneDividedByAlpha1PlusAlpha2_1s2px = GenerateOneDividedByAlpha1PlusAlpha2_1s2px<Hydrogen, Oxygen>();
	static constexpr std::array<double, 9> oneDividedByAlpha1PlusAlpha2_1s2py = GenerateOneDividedByAlpha1PlusAlpha2_1s2py<Hydrogen, Oxygen>();
	static constexpr std::array<double, 9> oneDividedByAlpha1PlusAlpha2_1s2pz = GenerateOneDividedByAlpha1PlusAlpha2_1s2pz<Hydrogen, Oxygen>();
	static constexpr std::array<double, 9> factors_1s1s = GenerateOverlapFactors_1s1s<Hydrogen, Oxygen>();
	static constexpr std::array<double, 9> factors_1s2s = GenerateOverlapFactors_1s2s<Hydrogen, Oxygen>();
	static constexpr std::array<double, 9> factors_1s2px = GenerateOverlapFactors_1s2px<Hydrogen, Oxygen>();
	static constexpr std::array<double, 9> factors_1s2py = GenerateOverlapFactors_1s2py<Hydrogen, Oxygen>();
	static constexpr std::array<double, 9> factors_1s2pz = GenerateOverlapFactors_1s2pz<Hydrogen, Oxygen>();
	static constexpr std::array<double, 9> expFactors_1s1s = GenerateOverlapExponentialFactors_1s1s<Hydrogen, Oxygen>();
	static constexpr std::array<double, 9> expFactors_1s2s = GenerateOverlapExponentialFactors_1s2s<Hydrogen, Oxygen>();
	static constexpr std::array<double, 9> expFactors_1s2px = GenerateOverlapExponentialFactors_1s2px<Hydrogen, Oxygen>();
	static constexpr std::array<double, 9> expFactors_1s2py = GenerateOverlapExponentialFactors_1s2py<Hydrogen, Oxygen>();
	static constexpr std::array<double, 9> expFactors_1s2pz = GenerateOverlapExponentialFactors_1s2pz<Hydrogen, Oxygen>();

	Vec3d diff = position2 - position1;
	double dotProduct = diff.Dot(diff);

	// H_1s - O_1s  ==============================================================
	double result = 0.0;
	result += factors_1s1s[0] * std::exp(expFactors_1s1s[0] * dotProduct);
	result += factors_1s1s[1] * std::exp(expFactors_1s1s[1] * dotProduct);
	result += factors_1s1s[2] * std::exp(expFactors_1s1s[2] * dotProduct);
	result += factors_1s1s[3] * std::exp(expFactors_1s1s[3] * dotProduct);
	result += factors_1s1s[4] * std::exp(expFactors_1s1s[4] * dotProduct);
	result += factors_1s1s[5] * std::exp(expFactors_1s1s[5] * dotProduct);
	result += factors_1s1s[6] * std::exp(expFactors_1s1s[6] * dotProduct);
	result += factors_1s1s[7] * std::exp(expFactors_1s1s[7] * dotProduct);
	result += factors_1s1s[8] * std::exp(expFactors_1s1s[8] * dotProduct);

	results[0] = result;

	// H_1s - O_2s  ==============================================================
	result = 0.0;
	result += factors_1s2s[0] * std::exp(expFactors_1s2s[0] * dotProduct);
	result += factors_1s2s[1] * std::exp(expFactors_1s2s[1] * dotProduct);
	result += factors_1s2s[2] * std::exp(expFactors_1s2s[2] * dotProduct);
	result += factors_1s2s[3] * std::exp(expFactors_1s2s[3] * dotProduct);
	result += factors_1s2s[4] * std::exp(expFactors_1s2s[4] * dotProduct);
	result += factors_1s2s[5] * std::exp(expFactors_1s2s[5] * dotProduct);
	result += factors_1s2s[6] * std::exp(expFactors_1s2s[6] * dotProduct);
	result += factors_1s2s[7] * std::exp(expFactors_1s2s[7] * dotProduct);
	result += factors_1s2s[8] * std::exp(expFactors_1s2s[8] * dotProduct);

	results[1] = result;

	// H_1s - O_2px  ==============================================================
	{
		const auto& atom1Prims = Hydrogen::orbital_1s.primitiveGaussians;
		const auto& atom2Prims = Oxygen::orbital_2px.primitiveGaussians;
		result = 0.0;
		result += factors_1s2px[0] * std::exp(expFactors_1s2px[0] * dotProduct) * (-position2.x + ((position1.x * atom1Prims[0].alpha + position2.x * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_1s2px[0]));
		result += factors_1s2px[1] * std::exp(expFactors_1s2px[1] * dotProduct) * (-position2.x + ((position1.x * atom1Prims[0].alpha + position2.x * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_1s2px[1]));
		result += factors_1s2px[2] * std::exp(expFactors_1s2px[2] * dotProduct) * (-position2.x + ((position1.x * atom1Prims[0].alpha + position2.x * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_1s2px[2]));
		result += factors_1s2px[3] * std::exp(expFactors_1s2px[3] * dotProduct) * (-position2.x + ((position1.x * atom1Prims[1].alpha + position2.x * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_1s2px[3]));
		result += factors_1s2px[4] * std::exp(expFactors_1s2px[4] * dotProduct) * (-position2.x + ((position1.x * atom1Prims[1].alpha + position2.x * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_1s2px[4]));
		result += factors_1s2px[5] * std::exp(expFactors_1s2px[5] * dotProduct) * (-position2.x + ((position1.x * atom1Prims[1].alpha + position2.x * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_1s2px[5]));
		result += factors_1s2px[6] * std::exp(expFactors_1s2px[6] * dotProduct) * (-position2.x + ((position1.x * atom1Prims[2].alpha + position2.x * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_1s2px[6]));
		result += factors_1s2px[7] * std::exp(expFactors_1s2px[7] * dotProduct) * (-position2.x + ((position1.x * atom1Prims[2].alpha + position2.x * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_1s2px[7]));
		result += factors_1s2px[8] * std::exp(expFactors_1s2px[8] * dotProduct) * (-position2.x + ((position1.x * atom1Prims[2].alpha + position2.x * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_1s2px[8]));
	}
	results[2] = result;

	// H_1s - O_2py  ==============================================================
	{
		const auto& atom1Prims = Hydrogen::orbital_1s.primitiveGaussians;
		const auto& atom2Prims = Oxygen::orbital_2py.primitiveGaussians;
		result = 0.0;
		result += factors_1s2py[0] * std::exp(expFactors_1s2py[0] * dotProduct) * (-position2.y + ((position1.y * atom1Prims[0].alpha + position2.y * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_1s2py[0]));
		result += factors_1s2py[1] * std::exp(expFactors_1s2py[1] * dotProduct) * (-position2.y + ((position1.y * atom1Prims[0].alpha + position2.y * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_1s2py[1]));
		result += factors_1s2py[2] * std::exp(expFactors_1s2py[2] * dotProduct) * (-position2.y + ((position1.y * atom1Prims[0].alpha + position2.y * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_1s2py[2]));
		result += factors_1s2py[3] * std::exp(expFactors_1s2py[3] * dotProduct) * (-position2.y + ((position1.y * atom1Prims[1].alpha + position2.y * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_1s2py[3]));
		result += factors_1s2py[4] * std::exp(expFactors_1s2py[4] * dotProduct) * (-position2.y + ((position1.y * atom1Prims[1].alpha + position2.y * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_1s2py[4]));
		result += factors_1s2py[5] * std::exp(expFactors_1s2py[5] * dotProduct) * (-position2.y + ((position1.y * atom1Prims[1].alpha + position2.y * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_1s2py[5]));
		result += factors_1s2py[6] * std::exp(expFactors_1s2py[6] * dotProduct) * (-position2.y + ((position1.y * atom1Prims[2].alpha + position2.y * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_1s2py[6]));
		result += factors_1s2py[7] * std::exp(expFactors_1s2py[7] * dotProduct) * (-position2.y + ((position1.y * atom1Prims[2].alpha + position2.y * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_1s2py[7]));
		result += factors_1s2py[8] * std::exp(expFactors_1s2py[8] * dotProduct) * (-position2.y + ((position1.y * atom1Prims[2].alpha + position2.y * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_1s2py[8]));
	}
	results[3] = result;

	// H_1s - O_2pz  ==============================================================
	{
		const auto& atom1Prims = Hydrogen::orbital_1s.primitiveGaussians;
		const auto& atom2Prims = Oxygen::orbital_2pz.primitiveGaussians;
		result = 0.0;
		result += factors_1s2pz[0] * std::exp(expFactors_1s2pz[0] * dotProduct) * (-position2.z + ((position1.z * atom1Prims[0].alpha + position2.z * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_1s2pz[0]));
		result += factors_1s2pz[1] * std::exp(expFactors_1s2pz[1] * dotProduct) * (-position2.z + ((position1.z * atom1Prims[0].alpha + position2.z * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_1s2pz[1]));
		result += factors_1s2pz[2] * std::exp(expFactors_1s2pz[2] * dotProduct) * (-position2.z + ((position1.z * atom1Prims[0].alpha + position2.z * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_1s2pz[2]));
		result += factors_1s2pz[3] * std::exp(expFactors_1s2pz[3] * dotProduct) * (-position2.z + ((position1.z * atom1Prims[1].alpha + position2.z * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_1s2pz[3]));
		result += factors_1s2pz[4] * std::exp(expFactors_1s2pz[4] * dotProduct) * (-position2.z + ((position1.z * atom1Prims[1].alpha + position2.z * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_1s2pz[4]));
		result += factors_1s2pz[5] * std::exp(expFactors_1s2pz[5] * dotProduct) * (-position2.z + ((position1.z * atom1Prims[1].alpha + position2.z * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_1s2pz[5]));
		result += factors_1s2pz[6] * std::exp(expFactors_1s2pz[6] * dotProduct) * (-position2.z + ((position1.z * atom1Prims[2].alpha + position2.z * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_1s2pz[6]));
		result += factors_1s2pz[7] * std::exp(expFactors_1s2pz[7] * dotProduct) * (-position2.z + ((position1.z * atom1Prims[2].alpha + position2.z * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_1s2pz[7]));
		result += factors_1s2pz[8] * std::exp(expFactors_1s2pz[8] * dotProduct) * (-position2.z + ((position1.z * atom1Prims[2].alpha + position2.z * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_1s2pz[8]));
	}
	results[4] = result;

	return results;
}
template<>
void OverlapAtomAtomImpl<Hydrogen, Oxygen>(const Vec3d& position1, const Vec3d& position2, unsigned int row, unsigned int col, MatrixXd& overlapMatrix) noexcept
{
	// NOTE: With the ordering of <Hydrogen, Oxygen>, the block is a horizontal section of 5 values
	overlapMatrix.block<1, 5>(row, col) = OverlapHydrogenOxygen(position1, position2);
	overlapMatrix.block<5, 1>(col, row) = overlapMatrix.block<1, 5>(row, col); 
}
template<>
void OverlapAtomAtomImpl<Oxygen, Hydrogen>(const Vec3d& position1, const Vec3d& position2, unsigned int row, unsigned int col, MatrixXd& overlapMatrix) noexcept
{
	// NOTE: With the ordering of <Oxygen, Hydrogen>, the block is a vertical section of 5 values
	// NOTE: we flip the order of the positions in the call to OverlapHydrogenOxygen because it expects the hydrogen position to come first
	overlapMatrix.block<5, 1>(row, col) = OverlapHydrogenOxygen(position2, position1);
	overlapMatrix.block<1, 5>(col, row) = overlapMatrix.block<5, 1>(row, col);
}

template<>
void OverlapAtomAtomImpl<Oxygen, Oxygen>(const Vec3d& position1, const Vec3d& position2, unsigned int row, unsigned int col, MatrixXd& overlapMatrix) noexcept
{
	static constexpr std::array<double, 9> oneDividedByAlpha1PlusAlpha2_1s2px  = GenerateOneDividedByAlpha1PlusAlpha2_1s2px<Oxygen, Oxygen>();
	static constexpr std::array<double, 9> oneDividedByAlpha1PlusAlpha2_1s2py  = GenerateOneDividedByAlpha1PlusAlpha2_1s2py<Oxygen, Oxygen>();
	static constexpr std::array<double, 9> oneDividedByAlpha1PlusAlpha2_1s2pz  = GenerateOneDividedByAlpha1PlusAlpha2_1s2pz<Oxygen, Oxygen>();
	static constexpr std::array<double, 9> oneDividedByAlpha1PlusAlpha2_2s2px  = GenerateOneDividedByAlpha1PlusAlpha2_2s2px<Oxygen, Oxygen>();
	static constexpr std::array<double, 9> oneDividedByAlpha1PlusAlpha2_2s2py  = GenerateOneDividedByAlpha1PlusAlpha2_2s2py<Oxygen, Oxygen>();
	static constexpr std::array<double, 9> oneDividedByAlpha1PlusAlpha2_2s2pz  = GenerateOneDividedByAlpha1PlusAlpha2_2s2pz<Oxygen, Oxygen>();
	static constexpr std::array<double, 9> oneDividedByAlpha1PlusAlpha2_2px2px = GenerateOneDividedByAlpha1PlusAlpha2_2px2px<Oxygen, Oxygen>();
	static constexpr std::array<double, 9> oneDividedByAlpha1PlusAlpha2_2py2py = GenerateOneDividedByAlpha1PlusAlpha2_2py2py<Oxygen, Oxygen>();
	static constexpr std::array<double, 9> oneDividedByAlpha1PlusAlpha2_2pz2pz = GenerateOneDividedByAlpha1PlusAlpha2_2pz2pz<Oxygen, Oxygen>();

	static constexpr std::array<double, 9> factors_1s1s     = GenerateOverlapFactors_1s1s<Oxygen, Oxygen>();
	static constexpr std::array<double, 9> factors_1s2s     = GenerateOverlapFactors_1s2s<Oxygen, Oxygen>();
	static constexpr std::array<double, 9> factors_1s2px    = GenerateOverlapFactors_1s2px<Oxygen, Oxygen>();
	static constexpr std::array<double, 9> factors_1s2py    = GenerateOverlapFactors_1s2py<Oxygen, Oxygen>();
	static constexpr std::array<double, 9> factors_1s2pz    = GenerateOverlapFactors_1s2pz<Oxygen, Oxygen>();
	static constexpr std::array<double, 9> factors_2s2s     = GenerateOverlapFactors_2s2s<Oxygen, Oxygen>();
	static constexpr std::array<double, 9> factors_2s2px    = GenerateOverlapFactors_2s2px<Oxygen, Oxygen>();
	static constexpr std::array<double, 9> factors_2s2py    = GenerateOverlapFactors_2s2py<Oxygen, Oxygen>();
	static constexpr std::array<double, 9> factors_2s2pz    = GenerateOverlapFactors_2s2pz<Oxygen, Oxygen>();
	static constexpr std::array<double, 9> factors_2px2px   = GenerateOverlapFactors_2px2px<Oxygen, Oxygen>();
	static constexpr std::array<double, 9> factors_2py2py   = GenerateOverlapFactors_2py2py<Oxygen, Oxygen>();
	static constexpr std::array<double, 9> factors_2pz2pz   = GenerateOverlapFactors_2pz2pz<Oxygen, Oxygen>();

	static constexpr std::array<double, 9> expFactors_1s1s   = GenerateOverlapExponentialFactors_1s1s<Oxygen, Oxygen>();
	static constexpr std::array<double, 9> expFactors_1s2s   = GenerateOverlapExponentialFactors_1s2s<Oxygen, Oxygen>();
	static constexpr std::array<double, 9> expFactors_1s2px  = GenerateOverlapExponentialFactors_1s2px<Oxygen, Oxygen>();
	static constexpr std::array<double, 9> expFactors_1s2py  = GenerateOverlapExponentialFactors_1s2py<Oxygen, Oxygen>();
	static constexpr std::array<double, 9> expFactors_1s2pz  = GenerateOverlapExponentialFactors_1s2pz<Oxygen, Oxygen>();
	static constexpr std::array<double, 9> expFactors_2s2s   = GenerateOverlapExponentialFactors_2s2s<Oxygen, Oxygen>();
	static constexpr std::array<double, 9> expFactors_2s2px  = GenerateOverlapExponentialFactors_2s2px<Oxygen, Oxygen>();
	static constexpr std::array<double, 9> expFactors_2s2py  = GenerateOverlapExponentialFactors_2s2py<Oxygen, Oxygen>();
	static constexpr std::array<double, 9> expFactors_2s2pz  = GenerateOverlapExponentialFactors_2s2pz<Oxygen, Oxygen>();
	static constexpr std::array<double, 9> expFactors_2px2px = GenerateOverlapExponentialFactors_2px2px<Oxygen, Oxygen>();
	static constexpr std::array<double, 9> expFactors_2py2py = GenerateOverlapExponentialFactors_2py2py<Oxygen, Oxygen>();
	static constexpr std::array<double, 9> expFactors_2pz2pz = GenerateOverlapExponentialFactors_2pz2pz<Oxygen, Oxygen>();

	Vec3d diff = position2 - position1;
	double dotProduct = diff.Dot(diff);

	// O_1s - O_1s  ==============================================================
	double result = 0.0;
	result += factors_1s1s[0] * std::exp(expFactors_1s1s[0] * dotProduct);
	result += factors_1s1s[1] * std::exp(expFactors_1s1s[1] * dotProduct);
	result += factors_1s1s[2] * std::exp(expFactors_1s1s[2] * dotProduct);
	result += factors_1s1s[3] * std::exp(expFactors_1s1s[3] * dotProduct);
	result += factors_1s1s[4] * std::exp(expFactors_1s1s[4] * dotProduct);
	result += factors_1s1s[5] * std::exp(expFactors_1s1s[5] * dotProduct);
	result += factors_1s1s[6] * std::exp(expFactors_1s1s[6] * dotProduct);
	result += factors_1s1s[7] * std::exp(expFactors_1s1s[7] * dotProduct);
	result += factors_1s1s[8] * std::exp(expFactors_1s1s[8] * dotProduct);

	overlapMatrix(row, col) = result;

	// O_1s - O_2s  ==============================================================
	result = 0.0;
	result += factors_1s2s[0] * std::exp(expFactors_1s2s[0] * dotProduct);
	result += factors_1s2s[1] * std::exp(expFactors_1s2s[1] * dotProduct);
	result += factors_1s2s[2] * std::exp(expFactors_1s2s[2] * dotProduct);
	result += factors_1s2s[3] * std::exp(expFactors_1s2s[3] * dotProduct);
	result += factors_1s2s[4] * std::exp(expFactors_1s2s[4] * dotProduct);
	result += factors_1s2s[5] * std::exp(expFactors_1s2s[5] * dotProduct);
	result += factors_1s2s[6] * std::exp(expFactors_1s2s[6] * dotProduct);
	result += factors_1s2s[7] * std::exp(expFactors_1s2s[7] * dotProduct);
	result += factors_1s2s[8] * std::exp(expFactors_1s2s[8] * dotProduct);

	overlapMatrix(row, col + 1) = result;
	overlapMatrix(row + 1, col) = result; // O_2s - O_1s is the same as O_1s - O_2s

	// O_1s - O_2px  ==============================================================
	{
		const auto& atom1Prims = Oxygen::orbital_1s.primitiveGaussians;
		const auto& atom2Prims = Oxygen::orbital_2px.primitiveGaussians;
		result = 0.0;
		result += factors_1s2px[0] * std::exp(expFactors_1s2px[0] * dotProduct) * (-position2.x + ((position1.x * atom1Prims[0].alpha + position2.x * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_1s2px[0]));
		result += factors_1s2px[1] * std::exp(expFactors_1s2px[1] * dotProduct) * (-position2.x + ((position1.x * atom1Prims[0].alpha + position2.x * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_1s2px[1]));
		result += factors_1s2px[2] * std::exp(expFactors_1s2px[2] * dotProduct) * (-position2.x + ((position1.x * atom1Prims[0].alpha + position2.x * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_1s2px[2]));
		result += factors_1s2px[3] * std::exp(expFactors_1s2px[3] * dotProduct) * (-position2.x + ((position1.x * atom1Prims[1].alpha + position2.x * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_1s2px[3]));
		result += factors_1s2px[4] * std::exp(expFactors_1s2px[4] * dotProduct) * (-position2.x + ((position1.x * atom1Prims[1].alpha + position2.x * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_1s2px[4]));
		result += factors_1s2px[5] * std::exp(expFactors_1s2px[5] * dotProduct) * (-position2.x + ((position1.x * atom1Prims[1].alpha + position2.x * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_1s2px[5]));
		result += factors_1s2px[6] * std::exp(expFactors_1s2px[6] * dotProduct) * (-position2.x + ((position1.x * atom1Prims[2].alpha + position2.x * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_1s2px[6]));
		result += factors_1s2px[7] * std::exp(expFactors_1s2px[7] * dotProduct) * (-position2.x + ((position1.x * atom1Prims[2].alpha + position2.x * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_1s2px[7]));
		result += factors_1s2px[8] * std::exp(expFactors_1s2px[8] * dotProduct) * (-position2.x + ((position1.x * atom1Prims[2].alpha + position2.x * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_1s2px[8]));
	}
	overlapMatrix(row, col + 2) = result;
	overlapMatrix(row + 2, col) = -1 * result; // O_2px - O_1s is just the negative of O_1s - O_2px

	// O_1s - O_2py  ==============================================================
	{
		const auto& atom1Prims = Oxygen::orbital_1s.primitiveGaussians;
		const auto& atom2Prims = Oxygen::orbital_2py.primitiveGaussians;
		result = 0.0;
		result += factors_1s2py[0] * std::exp(expFactors_1s2py[0] * dotProduct) * (-position2.y + ((position1.y * atom1Prims[0].alpha + position2.y * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_1s2py[0]));
		result += factors_1s2py[1] * std::exp(expFactors_1s2py[1] * dotProduct) * (-position2.y + ((position1.y * atom1Prims[0].alpha + position2.y * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_1s2py[1]));
		result += factors_1s2py[2] * std::exp(expFactors_1s2py[2] * dotProduct) * (-position2.y + ((position1.y * atom1Prims[0].alpha + position2.y * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_1s2py[2]));
		result += factors_1s2py[3] * std::exp(expFactors_1s2py[3] * dotProduct) * (-position2.y + ((position1.y * atom1Prims[1].alpha + position2.y * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_1s2py[3]));
		result += factors_1s2py[4] * std::exp(expFactors_1s2py[4] * dotProduct) * (-position2.y + ((position1.y * atom1Prims[1].alpha + position2.y * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_1s2py[4]));
		result += factors_1s2py[5] * std::exp(expFactors_1s2py[5] * dotProduct) * (-position2.y + ((position1.y * atom1Prims[1].alpha + position2.y * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_1s2py[5]));
		result += factors_1s2py[6] * std::exp(expFactors_1s2py[6] * dotProduct) * (-position2.y + ((position1.y * atom1Prims[2].alpha + position2.y * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_1s2py[6]));
		result += factors_1s2py[7] * std::exp(expFactors_1s2py[7] * dotProduct) * (-position2.y + ((position1.y * atom1Prims[2].alpha + position2.y * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_1s2py[7]));
		result += factors_1s2py[8] * std::exp(expFactors_1s2py[8] * dotProduct) * (-position2.y + ((position1.y * atom1Prims[2].alpha + position2.y * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_1s2py[8]));
	}
	overlapMatrix(row, col + 3) = result;
	overlapMatrix(row + 3, col) = -1 * result; // O_2py - O_1s is just the negative of O_1s - O_2py

	// O_1s - O_2pz  ==============================================================
	{
		const auto& atom1Prims = Oxygen::orbital_1s.primitiveGaussians;
		const auto& atom2Prims = Oxygen::orbital_2pz.primitiveGaussians;
		result = 0.0;
		result += factors_1s2pz[0] * std::exp(expFactors_1s2pz[0] * dotProduct) * (-position2.z + ((position1.z * atom1Prims[0].alpha + position2.z * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_1s2pz[0]));
		result += factors_1s2pz[1] * std::exp(expFactors_1s2pz[1] * dotProduct) * (-position2.z + ((position1.z * atom1Prims[0].alpha + position2.z * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_1s2pz[1]));
		result += factors_1s2pz[2] * std::exp(expFactors_1s2pz[2] * dotProduct) * (-position2.z + ((position1.z * atom1Prims[0].alpha + position2.z * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_1s2pz[2]));
		result += factors_1s2pz[3] * std::exp(expFactors_1s2pz[3] * dotProduct) * (-position2.z + ((position1.z * atom1Prims[1].alpha + position2.z * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_1s2pz[3]));
		result += factors_1s2pz[4] * std::exp(expFactors_1s2pz[4] * dotProduct) * (-position2.z + ((position1.z * atom1Prims[1].alpha + position2.z * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_1s2pz[4]));
		result += factors_1s2pz[5] * std::exp(expFactors_1s2pz[5] * dotProduct) * (-position2.z + ((position1.z * atom1Prims[1].alpha + position2.z * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_1s2pz[5]));
		result += factors_1s2pz[6] * std::exp(expFactors_1s2pz[6] * dotProduct) * (-position2.z + ((position1.z * atom1Prims[2].alpha + position2.z * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_1s2pz[6]));
		result += factors_1s2pz[7] * std::exp(expFactors_1s2pz[7] * dotProduct) * (-position2.z + ((position1.z * atom1Prims[2].alpha + position2.z * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_1s2pz[7]));
		result += factors_1s2pz[8] * std::exp(expFactors_1s2pz[8] * dotProduct) * (-position2.z + ((position1.z * atom1Prims[2].alpha + position2.z * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_1s2pz[8]));
	}
	overlapMatrix(row, col + 4) = result;
	overlapMatrix(row + 4, col) = -1 * result; // O_2pz - O_1s is just the negative of O_1s - O_2pz

	// O_2s - O_2s  ==============================================================
	result = 0.0;
	result += factors_2s2s[0] * std::exp(expFactors_2s2s[0] * dotProduct);
	result += factors_2s2s[1] * std::exp(expFactors_2s2s[1] * dotProduct);
	result += factors_2s2s[2] * std::exp(expFactors_2s2s[2] * dotProduct);
	result += factors_2s2s[3] * std::exp(expFactors_2s2s[3] * dotProduct);
	result += factors_2s2s[4] * std::exp(expFactors_2s2s[4] * dotProduct);
	result += factors_2s2s[5] * std::exp(expFactors_2s2s[5] * dotProduct);
	result += factors_2s2s[6] * std::exp(expFactors_2s2s[6] * dotProduct);
	result += factors_2s2s[7] * std::exp(expFactors_2s2s[7] * dotProduct);
	result += factors_2s2s[8] * std::exp(expFactors_2s2s[8] * dotProduct);

	overlapMatrix(row + 1, col + 1) = result;

	// O_2s - O_2px  ==============================================================
	{
		const auto& atom1Prims = Oxygen::orbital_2s.primitiveGaussians;
		const auto& atom2Prims = Oxygen::orbital_2px.primitiveGaussians;
		result = 0.0;
		result += factors_2s2px[0] * std::exp(expFactors_2s2px[0] * dotProduct) * (-position2.x + ((position1.x * atom1Prims[0].alpha + position2.x * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2s2px[0]));
		result += factors_2s2px[1] * std::exp(expFactors_2s2px[1] * dotProduct) * (-position2.x + ((position1.x * atom1Prims[0].alpha + position2.x * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2s2px[1]));
		result += factors_2s2px[2] * std::exp(expFactors_2s2px[2] * dotProduct) * (-position2.x + ((position1.x * atom1Prims[0].alpha + position2.x * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2s2px[2]));
		result += factors_2s2px[3] * std::exp(expFactors_2s2px[3] * dotProduct) * (-position2.x + ((position1.x * atom1Prims[1].alpha + position2.x * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2s2px[3]));
		result += factors_2s2px[4] * std::exp(expFactors_2s2px[4] * dotProduct) * (-position2.x + ((position1.x * atom1Prims[1].alpha + position2.x * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2s2px[4]));
		result += factors_2s2px[5] * std::exp(expFactors_2s2px[5] * dotProduct) * (-position2.x + ((position1.x * atom1Prims[1].alpha + position2.x * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2s2px[5]));
		result += factors_2s2px[6] * std::exp(expFactors_2s2px[6] * dotProduct) * (-position2.x + ((position1.x * atom1Prims[2].alpha + position2.x * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2s2px[6]));
		result += factors_2s2px[7] * std::exp(expFactors_2s2px[7] * dotProduct) * (-position2.x + ((position1.x * atom1Prims[2].alpha + position2.x * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2s2px[7]));
		result += factors_2s2px[8] * std::exp(expFactors_2s2px[8] * dotProduct) * (-position2.x + ((position1.x * atom1Prims[2].alpha + position2.x * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2s2px[8]));
	}
	overlapMatrix(row + 1, col + 2) = result;
	overlapMatrix(row + 2, col + 1) = -1 * result; // O_2px - O_2s is just the negative of O_2s - O_2px

	// O_2s - O_2py  ==============================================================
	{
		const auto& atom1Prims = Oxygen::orbital_2s.primitiveGaussians;
		const auto& atom2Prims = Oxygen::orbital_2py.primitiveGaussians;
		result = 0.0;
		result += factors_2s2py[0] * std::exp(expFactors_2s2py[0] * dotProduct) * (-position2.y + ((position1.y * atom1Prims[0].alpha + position2.y * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2s2py[0]));
		result += factors_2s2py[1] * std::exp(expFactors_2s2py[1] * dotProduct) * (-position2.y + ((position1.y * atom1Prims[0].alpha + position2.y * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2s2py[1]));
		result += factors_2s2py[2] * std::exp(expFactors_2s2py[2] * dotProduct) * (-position2.y + ((position1.y * atom1Prims[0].alpha + position2.y * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2s2py[2]));
		result += factors_2s2py[3] * std::exp(expFactors_2s2py[3] * dotProduct) * (-position2.y + ((position1.y * atom1Prims[1].alpha + position2.y * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2s2py[3]));
		result += factors_2s2py[4] * std::exp(expFactors_2s2py[4] * dotProduct) * (-position2.y + ((position1.y * atom1Prims[1].alpha + position2.y * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2s2py[4]));
		result += factors_2s2py[5] * std::exp(expFactors_2s2py[5] * dotProduct) * (-position2.y + ((position1.y * atom1Prims[1].alpha + position2.y * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2s2py[5]));
		result += factors_2s2py[6] * std::exp(expFactors_2s2py[6] * dotProduct) * (-position2.y + ((position1.y * atom1Prims[2].alpha + position2.y * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2s2py[6]));
		result += factors_2s2py[7] * std::exp(expFactors_2s2py[7] * dotProduct) * (-position2.y + ((position1.y * atom1Prims[2].alpha + position2.y * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2s2py[7]));
		result += factors_2s2py[8] * std::exp(expFactors_2s2py[8] * dotProduct) * (-position2.y + ((position1.y * atom1Prims[2].alpha + position2.y * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2s2py[8]));
	}
	overlapMatrix(row + 1, col + 3) = result;
	overlapMatrix(row + 3, col + 1) = -1 * result; // O_2py - O_2s is just the negative of O_2s - O_2py

	// O_2s - O_2pz  ==============================================================
	{
		const auto& atom1Prims = Oxygen::orbital_2s.primitiveGaussians;
		const auto& atom2Prims = Oxygen::orbital_2pz.primitiveGaussians;
		result = 0.0;
		result += factors_2s2pz[0] * std::exp(expFactors_2s2pz[0] * dotProduct) * (-position2.z + ((position1.z * atom1Prims[0].alpha + position2.z * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2s2pz[0]));
		result += factors_2s2pz[1] * std::exp(expFactors_2s2pz[1] * dotProduct) * (-position2.z + ((position1.z * atom1Prims[0].alpha + position2.z * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2s2pz[1]));
		result += factors_2s2pz[2] * std::exp(expFactors_2s2pz[2] * dotProduct) * (-position2.z + ((position1.z * atom1Prims[0].alpha + position2.z * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2s2pz[2]));
		result += factors_2s2pz[3] * std::exp(expFactors_2s2pz[3] * dotProduct) * (-position2.z + ((position1.z * atom1Prims[1].alpha + position2.z * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2s2pz[3]));
		result += factors_2s2pz[4] * std::exp(expFactors_2s2pz[4] * dotProduct) * (-position2.z + ((position1.z * atom1Prims[1].alpha + position2.z * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2s2pz[4]));
		result += factors_2s2pz[5] * std::exp(expFactors_2s2pz[5] * dotProduct) * (-position2.z + ((position1.z * atom1Prims[1].alpha + position2.z * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2s2pz[5]));
		result += factors_2s2pz[6] * std::exp(expFactors_2s2pz[6] * dotProduct) * (-position2.z + ((position1.z * atom1Prims[2].alpha + position2.z * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2s2pz[6]));
		result += factors_2s2pz[7] * std::exp(expFactors_2s2pz[7] * dotProduct) * (-position2.z + ((position1.z * atom1Prims[2].alpha + position2.z * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2s2pz[7]));
		result += factors_2s2pz[8] * std::exp(expFactors_2s2pz[8] * dotProduct) * (-position2.z + ((position1.z * atom1Prims[2].alpha + position2.z * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2s2pz[8]));
	}
	overlapMatrix(row + 1, col + 4) = result;
	overlapMatrix(row + 4, col + 1) = -1 * result; // O_2pz - O_2s is just the negative of O_2s - O_2pz

	// O_2px - O_2px  ==============================================================
	{
		const auto& atom1Prims = Oxygen::orbital_2px.primitiveGaussians;
		const auto& atom2Prims = Oxygen::orbital_2px.primitiveGaussians;
		result = 0.0;

		double positionDelta = position1.x - position2.x;

		double s_x_1_0 = -1 * (position1.x - ((position1.x * atom1Prims[0].alpha + position2.x * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2px2px[0]));
		result += factors_2px2px[0] * std::exp(expFactors_2px2px[0] * dotProduct) * (s_x_1_0 * s_x_1_0 + 0.5 * oneDividedByAlpha1PlusAlpha2_2px2px[0] + positionDelta * s_x_1_0);
		
		s_x_1_0 = -1 * (position1.x - ((position1.x * atom1Prims[0].alpha + position2.x * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2px2px[1]));
		result += factors_2px2px[1] * std::exp(expFactors_2px2px[1] * dotProduct) * (s_x_1_0 * s_x_1_0 + 0.5 * oneDividedByAlpha1PlusAlpha2_2px2px[1] + positionDelta * s_x_1_0);

		s_x_1_0 = -1 * (position1.x - ((position1.x * atom1Prims[0].alpha + position2.x * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2px2px[2]));
		result += factors_2px2px[2] * std::exp(expFactors_2px2px[2] * dotProduct) * (s_x_1_0 * s_x_1_0 + 0.5 * oneDividedByAlpha1PlusAlpha2_2px2px[2] + positionDelta * s_x_1_0);

		s_x_1_0 = -1 * (position1.x - ((position1.x * atom1Prims[1].alpha + position2.x * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2px2px[3]));
		result += factors_2px2px[3] * std::exp(expFactors_2px2px[3] * dotProduct) * (s_x_1_0 * s_x_1_0 + 0.5 * oneDividedByAlpha1PlusAlpha2_2px2px[3] + positionDelta * s_x_1_0);

		s_x_1_0 = -1 * (position1.x - ((position1.x * atom1Prims[1].alpha + position2.x * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2px2px[4]));
		result += factors_2px2px[4] * std::exp(expFactors_2px2px[4] * dotProduct) * (s_x_1_0 * s_x_1_0 + 0.5 * oneDividedByAlpha1PlusAlpha2_2px2px[4] + positionDelta * s_x_1_0);

		s_x_1_0 = -1 * (position1.x - ((position1.x * atom1Prims[1].alpha + position2.x * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2px2px[5]));
		result += factors_2px2px[5] * std::exp(expFactors_2px2px[5] * dotProduct) * (s_x_1_0 * s_x_1_0 + 0.5 * oneDividedByAlpha1PlusAlpha2_2px2px[5] + positionDelta * s_x_1_0);

		s_x_1_0 = -1 * (position1.x - ((position1.x * atom1Prims[2].alpha + position2.x * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2px2px[6]));
		result += factors_2px2px[6] * std::exp(expFactors_2px2px[6] * dotProduct) * (s_x_1_0 * s_x_1_0 + 0.5 * oneDividedByAlpha1PlusAlpha2_2px2px[6] + positionDelta * s_x_1_0);

		s_x_1_0 = -1 * (position1.x - ((position1.x * atom1Prims[2].alpha + position2.x * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2px2px[7]));
		result += factors_2px2px[7] * std::exp(expFactors_2px2px[7] * dotProduct) * (s_x_1_0 * s_x_1_0 + 0.5 * oneDividedByAlpha1PlusAlpha2_2px2px[7] + positionDelta * s_x_1_0);

		s_x_1_0 = -1 * (position1.x - ((position1.x * atom1Prims[2].alpha + position2.x * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2px2px[8]));
		result += factors_2px2px[8] * std::exp(expFactors_2px2px[8] * dotProduct) * (s_x_1_0 * s_x_1_0 + 0.5 * oneDividedByAlpha1PlusAlpha2_2px2px[8] + positionDelta * s_x_1_0);
	}
	overlapMatrix(row + 2, col + 2) = result;

	// O_2py - O_2py  ==============================================================
	{
		const auto& atom1Prims = Oxygen::orbital_2py.primitiveGaussians;
		const auto& atom2Prims = Oxygen::orbital_2py.primitiveGaussians;
		result = 0.0;

		double positionDelta = position1.y - position2.y;

		double s_y_1_0 = -1 * (position1.y - ((position1.y * atom1Prims[0].alpha + position2.y * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2py2py[0]));
		result += factors_2py2py[0] * std::exp(expFactors_2py2py[0] * dotProduct) * (s_y_1_0 * s_y_1_0 + 0.5 * oneDividedByAlpha1PlusAlpha2_2py2py[0] + positionDelta * s_y_1_0);

		s_y_1_0 = -1 * (position1.y - ((position1.y * atom1Prims[0].alpha + position2.y * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2py2py[1]));
		result += factors_2py2py[1] * std::exp(expFactors_2py2py[1] * dotProduct) * (s_y_1_0 * s_y_1_0 + 0.5 * oneDividedByAlpha1PlusAlpha2_2py2py[1] + positionDelta * s_y_1_0);

		s_y_1_0 = -1 * (position1.y - ((position1.y * atom1Prims[0].alpha + position2.y * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2py2py[2]));
		result += factors_2py2py[2] * std::exp(expFactors_2py2py[2] * dotProduct) * (s_y_1_0 * s_y_1_0 + 0.5 * oneDividedByAlpha1PlusAlpha2_2py2py[2] + positionDelta * s_y_1_0);

		s_y_1_0 = -1 * (position1.y - ((position1.y * atom1Prims[1].alpha + position2.y * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2py2py[3]));
		result += factors_2py2py[3] * std::exp(expFactors_2py2py[3] * dotProduct) * (s_y_1_0 * s_y_1_0 + 0.5 * oneDividedByAlpha1PlusAlpha2_2py2py[3] + positionDelta * s_y_1_0);

		s_y_1_0 = -1 * (position1.y - ((position1.y * atom1Prims[1].alpha + position2.y * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2py2py[4]));
		result += factors_2py2py[4] * std::exp(expFactors_2py2py[4] * dotProduct) * (s_y_1_0 * s_y_1_0 + 0.5 * oneDividedByAlpha1PlusAlpha2_2py2py[4] + positionDelta * s_y_1_0);

		s_y_1_0 = -1 * (position1.y - ((position1.y * atom1Prims[1].alpha + position2.y * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2py2py[5]));
		result += factors_2py2py[5] * std::exp(expFactors_2py2py[5] * dotProduct) * (s_y_1_0 * s_y_1_0 + 0.5 * oneDividedByAlpha1PlusAlpha2_2py2py[5] + positionDelta * s_y_1_0);

		s_y_1_0 = -1 * (position1.y - ((position1.y * atom1Prims[2].alpha + position2.y * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2py2py[6]));
		result += factors_2py2py[6] * std::exp(expFactors_2py2py[6] * dotProduct) * (s_y_1_0 * s_y_1_0 + 0.5 * oneDividedByAlpha1PlusAlpha2_2py2py[6] + positionDelta * s_y_1_0);

		s_y_1_0 = -1 * (position1.y - ((position1.y * atom1Prims[2].alpha + position2.y * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2py2py[7]));
		result += factors_2py2py[7] * std::exp(expFactors_2py2py[7] * dotProduct) * (s_y_1_0 * s_y_1_0 + 0.5 * oneDividedByAlpha1PlusAlpha2_2py2py[7] + positionDelta * s_y_1_0);

		s_y_1_0 = -1 * (position1.y - ((position1.y * atom1Prims[2].alpha + position2.y * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2py2py[8]));
		result += factors_2py2py[8] * std::exp(expFactors_2py2py[8] * dotProduct) * (s_y_1_0 * s_y_1_0 + 0.5 * oneDividedByAlpha1PlusAlpha2_2py2py[8] + positionDelta * s_y_1_0);
	}
	overlapMatrix(row + 3, col + 3) = result;

	// O_2pz - O_2pz  ==============================================================
	{
		const auto& atom1Prims = Oxygen::orbital_2pz.primitiveGaussians;
		const auto& atom2Prims = Oxygen::orbital_2pz.primitiveGaussians;
		result = 0.0;

		double positionDelta = position1.z - position2.z;

		double s_z_1_0 = -1 * (position1.z - ((position1.z * atom1Prims[0].alpha + position2.z * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2pz[0]));
		result += factors_2pz2pz[0] * std::exp(expFactors_2pz2pz[0] * dotProduct) * (s_z_1_0 * s_z_1_0 + 0.5 * oneDividedByAlpha1PlusAlpha2_2pz2pz[0] + positionDelta * s_z_1_0);

		s_z_1_0 = -1 * (position1.z - ((position1.z * atom1Prims[0].alpha + position2.z * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2pz[1]));
		result += factors_2pz2pz[1] * std::exp(expFactors_2pz2pz[1] * dotProduct) * (s_z_1_0 * s_z_1_0 + 0.5 * oneDividedByAlpha1PlusAlpha2_2pz2pz[1] + positionDelta * s_z_1_0);

		s_z_1_0 = -1 * (position1.z - ((position1.z * atom1Prims[0].alpha + position2.z * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2pz[2]));
		result += factors_2pz2pz[2] * std::exp(expFactors_2pz2pz[2] * dotProduct) * (s_z_1_0 * s_z_1_0 + 0.5 * oneDividedByAlpha1PlusAlpha2_2pz2pz[2] + positionDelta * s_z_1_0);

		s_z_1_0 = -1 * (position1.z - ((position1.z * atom1Prims[1].alpha + position2.z * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2pz[3]));
		result += factors_2pz2pz[3] * std::exp(expFactors_2pz2pz[3] * dotProduct) * (s_z_1_0 * s_z_1_0 + 0.5 * oneDividedByAlpha1PlusAlpha2_2pz2pz[3] + positionDelta * s_z_1_0);

		s_z_1_0 = -1 * (position1.z - ((position1.z * atom1Prims[1].alpha + position2.z * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2pz[4]));
		result += factors_2pz2pz[4] * std::exp(expFactors_2pz2pz[4] * dotProduct) * (s_z_1_0 * s_z_1_0 + 0.5 * oneDividedByAlpha1PlusAlpha2_2pz2pz[4] + positionDelta * s_z_1_0);

		s_z_1_0 = -1 * (position1.z - ((position1.z * atom1Prims[1].alpha + position2.z * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2pz[5]));
		result += factors_2pz2pz[5] * std::exp(expFactors_2pz2pz[5] * dotProduct) * (s_z_1_0 * s_z_1_0 + 0.5 * oneDividedByAlpha1PlusAlpha2_2pz2pz[5] + positionDelta * s_z_1_0);

		s_z_1_0 = -1 * (position1.z - ((position1.z * atom1Prims[2].alpha + position2.z * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2pz[6]));
		result += factors_2pz2pz[6] * std::exp(expFactors_2pz2pz[6] * dotProduct) * (s_z_1_0 * s_z_1_0 + 0.5 * oneDividedByAlpha1PlusAlpha2_2pz2pz[6] + positionDelta * s_z_1_0);

		s_z_1_0 = -1 * (position1.z - ((position1.z * atom1Prims[2].alpha + position2.z * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2pz[7]));
		result += factors_2pz2pz[7] * std::exp(expFactors_2pz2pz[7] * dotProduct) * (s_z_1_0 * s_z_1_0 + 0.5 * oneDividedByAlpha1PlusAlpha2_2pz2pz[7] + positionDelta * s_z_1_0);

		s_z_1_0 = -1 * (position1.z - ((position1.z * atom1Prims[2].alpha + position2.z * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2pz[8]));
		result += factors_2pz2pz[8] * std::exp(expFactors_2pz2pz[8] * dotProduct) * (s_z_1_0 * s_z_1_0 + 0.5 * oneDividedByAlpha1PlusAlpha2_2pz2pz[8] + positionDelta * s_z_1_0);
	}
	overlapMatrix(row + 4, col + 4) = result;

	// Copy the values so the matrix is symmetric
	overlapMatrix.block<5, 5>(col, row) = overlapMatrix.block<5, 5>(row, col).transpose();
}

static void OverlapAtomAtom(const Atom& atom_1, const Atom& atom_2, unsigned int row, unsigned int col, MatrixXd& overlapMatrix) noexcept
{
	switch (atom_1.type)
	{
	case ATOM_TYPE::Hydrogen:
	{
		switch (atom_2.type)
		{
		case ATOM_TYPE::Hydrogen: OverlapAtomAtomImpl<Hydrogen, Hydrogen>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Oxygen:   OverlapAtomAtomImpl<Hydrogen, Oxygen>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		}
		return;
	}
	case ATOM_TYPE::Oxygen:
	{
		switch (atom_2.type)
		{
		case ATOM_TYPE::Hydrogen: OverlapAtomAtomImpl<Oxygen, Hydrogen>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		case ATOM_TYPE::Oxygen:   OverlapAtomAtomImpl<Oxygen, Oxygen>(atom_1.position, atom_2.position, row, col, overlapMatrix); return;
		}
	}
	}
}


// Self Overlap Functions
template<>
void OverlapWithSelfImpl<Hydrogen>(unsigned int topLeftIndex, MatrixXd& overlapMatrix) noexcept
{
	overlapMatrix(topLeftIndex, topLeftIndex) = 1.0;
}
template<>
void OverlapWithSelfImpl<Oxygen>(unsigned int topLeftIndex, MatrixXd& overlapMatrix) noexcept
{
	// With oxygen, all diagonal values will be 1, and the only off-diagonal value that is non-zero
	// is the 1s-2s overlap

	// Diagonal values
	overlapMatrix(topLeftIndex, topLeftIndex) = 1.0;
	overlapMatrix(topLeftIndex + 1, topLeftIndex + 1) = 1.0;
	overlapMatrix(topLeftIndex + 2, topLeftIndex + 2) = 1.0;
	overlapMatrix(topLeftIndex + 3, topLeftIndex + 3) = 1.0;
	overlapMatrix(topLeftIndex + 4, topLeftIndex + 4) = 1.0;

	// 1s-2s overlap
	static constexpr std::array<double, 9> factors_1s2s = GenerateOverlapFactors_1s2s<Oxygen, Oxygen>();
	static constexpr double overlap_1s2s = std::accumulate(factors_1s2s.begin(), factors_1s2s.end(), 0.0);
	
	overlapMatrix(topLeftIndex, topLeftIndex + 1) = overlap_1s2s;
	overlapMatrix(topLeftIndex + 1, topLeftIndex) = overlap_1s2s;
}
static void OverlapWithSelf(ATOM_TYPE type, unsigned int topLeftIndex, MatrixXd& overlapMatrix)
{
	switch (type)
	{
	case ATOM_TYPE::Hydrogen: OverlapWithSelfImpl<Hydrogen>(topLeftIndex, overlapMatrix); return;
	case ATOM_TYPE::Oxygen:   OverlapWithSelfImpl<Oxygen>(topLeftIndex, overlapMatrix); return;
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