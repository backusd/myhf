#pragma once

#include "Atom.h"
#include "Basis_Fast.h"


namespace myhf
{
namespace test
{

	[[nodiscard]] Eigen::MatrixXd OverlapMatrix(std::span<Atom> atoms) noexcept;




	// Overlap With Self
	template<typename Atom1>
	void OverlapWithSelf_1S(unsigned int topLeftIndex, Eigen::MatrixXd& overlapMatrix) noexcept
	{
		static_assert(Atom1::NumberOfContractedGaussians == 1);
		overlapMatrix(topLeftIndex, topLeftIndex) = 1.0;
	}
	template<typename Atom1>
	void OverlapWithSelf_2SP(unsigned int topLeftIndex, Eigen::MatrixXd& overlapMatrix) noexcept
	{
		static_assert(Atom1::NumberOfContractedGaussians == 5);

		// With atoms having 1s2s2p orbitals, all diagonal values will be 1, and the only off-diagonal 
		// value that is non-zero is the 1s-2s overlap

		// Diagonal values
		overlapMatrix(topLeftIndex, topLeftIndex) = 1.0;
		overlapMatrix(topLeftIndex + 1, topLeftIndex + 1) = 1.0;
		overlapMatrix(topLeftIndex + 2, topLeftIndex + 2) = 1.0;
		overlapMatrix(topLeftIndex + 3, topLeftIndex + 3) = 1.0;
		overlapMatrix(topLeftIndex + 4, topLeftIndex + 4) = 1.0;

		// 1s-2s overlap
		static constexpr std::array<double, 9> factors_1s2s = GenerateOverlapFactors_1s2s<Atom1, Atom1>();
		static constexpr double overlap_1s2s = std::accumulate(factors_1s2s.begin(), factors_1s2s.end(), 0.0);

		overlapMatrix(topLeftIndex, topLeftIndex + 1) = overlap_1s2s;
		overlapMatrix(topLeftIndex + 1, topLeftIndex) = overlap_1s2s;
	}


	// Overlap Between Atoms
	template<typename Atom1, typename Atom2>
	void OverlapAtomAtom_1S1S(const Vec3d& position1, const Vec3d& position2, unsigned int row, unsigned int col, Eigen::MatrixXd& overlapMatrix)
	{
		static_assert(Atom1::NumberOfContractedGaussians == 1);
		static_assert(Atom2::NumberOfContractedGaussians == 1);

		static constexpr std::array<double, 9> factors = GenerateOverlapFactors_1s1s<Atom1, Atom2>();
		static constexpr std::array<double, 9> expFactors = GenerateOverlapExponentialFactors_1s1s<Atom1, Atom2>();

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
	template<typename Atom1, typename Atom2>
	Eigen::Array<double, 1, 5> OverlapAtomAtom_1S2SP_Impl(const Vec3d& position1, const Vec3d& position2)
	{
		static_assert(Atom1::NumberOfContractedGaussians == 1);
		static_assert(Atom2::NumberOfContractedGaussians == 5);

		Eigen::Array<double, 1, 5> results{};

		static constexpr std::array<double, 9> oneDividedByAlpha1PlusAlpha2_1s2px = GenerateOneDividedByAlpha1PlusAlpha2_1s2px<Atom1, Atom2>();
		static constexpr std::array<double, 9> oneDividedByAlpha1PlusAlpha2_1s2py = GenerateOneDividedByAlpha1PlusAlpha2_1s2py<Atom1, Atom2>();
		static constexpr std::array<double, 9> oneDividedByAlpha1PlusAlpha2_1s2pz = GenerateOneDividedByAlpha1PlusAlpha2_1s2pz<Atom1, Atom2>();
		static constexpr std::array<double, 9> factors_1s1s = GenerateOverlapFactors_1s1s<Atom1, Atom2>();
		static constexpr std::array<double, 9> factors_1s2s = GenerateOverlapFactors_1s2s<Atom1, Atom2>();
		static constexpr std::array<double, 9> factors_1s2px = GenerateOverlapFactors_1s2px<Atom1, Atom2>();
		static constexpr std::array<double, 9> factors_1s2py = GenerateOverlapFactors_1s2py<Atom1, Atom2>();
		static constexpr std::array<double, 9> factors_1s2pz = GenerateOverlapFactors_1s2pz<Atom1, Atom2>();
		static constexpr std::array<double, 9> expFactors_1s1s = GenerateOverlapExponentialFactors_1s1s<Atom1, Atom2>();
		static constexpr std::array<double, 9> expFactors_1s2s = GenerateOverlapExponentialFactors_1s2s<Atom1, Atom2>();
		static constexpr std::array<double, 9> expFactors_1s2px = GenerateOverlapExponentialFactors_1s2px<Atom1, Atom2>();
		static constexpr std::array<double, 9> expFactors_1s2py = GenerateOverlapExponentialFactors_1s2py<Atom1, Atom2>();
		static constexpr std::array<double, 9> expFactors_1s2pz = GenerateOverlapExponentialFactors_1s2pz<Atom1, Atom2>();

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
			const auto& atom1Prims = Atom1::orbital_1s.primitiveGaussians;
			const auto& atom2Prims = Atom2::orbital_2px.primitiveGaussians;
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
			const auto& atom1Prims = Atom1::orbital_1s.primitiveGaussians;
			const auto& atom2Prims = Atom2::orbital_2py.primitiveGaussians;
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
			const auto& atom1Prims = Atom1::orbital_1s.primitiveGaussians;
			const auto& atom2Prims = Atom2::orbital_2pz.primitiveGaussians;
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
	template<typename Atom1, typename Atom2>
	void OverlapAtomAtom_1S2SP(const Vec3d& position1, const Vec3d& position2, unsigned int row, unsigned int col, Eigen::MatrixXd& overlapMatrix)
	{
		if constexpr (Atom1::NumberOfContractedGaussians == 1 && Atom2::NumberOfContractedGaussians == 5)
		{
			// NOTE: With the ordering of <1S, 2SP> (for example: <Hydrogen, Oxygen>), the block is a horizontal section of 5 values
			auto result = OverlapAtomAtom_1S2SP_Impl<Atom1, Atom2>(position1, position2);
			overlapMatrix.block<1, 5>(row, col) = result;
			overlapMatrix.block<5, 1>(col, row) = result;
		}
		else if constexpr (Atom1::NumberOfContractedGaussians == 5 && Atom2::NumberOfContractedGaussians == 1)
		{
			// NOTE: With the ordering of <2SP, 1S> (for example: <Oxygen, Hydrogen>), the block is a vertical section of 5 values
			// NOTE: we flip the order of the positions in the call to OverlapAtomAtom_1S2SP_Impl because it expects the 1S position to come first
			auto result = OverlapAtomAtom_1S2SP_Impl<Atom2, Atom1>(position2, position1);
			overlapMatrix.block<5, 1>(row, col) = result;
			overlapMatrix.block<1, 5>(col, row) = result;
		}
		else
		{
			static_assert(false);
		}
	}
	template<typename Atom1, typename Atom2>
	void OverlapAtomAtom_2SP2SP(const Vec3d& position1, const Vec3d& position2, unsigned int row, unsigned int col, Eigen::MatrixXd& overlapMatrix)
	{
		static_assert(Atom1::NumberOfContractedGaussians == 5);
		static_assert(Atom2::NumberOfContractedGaussians == 5);

		// 1 / (alpha_1 + alpha_2)
		static constexpr std::array<double, 9> oneDividedByAlpha1PlusAlpha2_1s2px = GenerateOneDividedByAlpha1PlusAlpha2_1s2px<Atom1, Atom2>();
		static constexpr std::array<double, 9> oneDividedByAlpha1PlusAlpha2_1s2py = GenerateOneDividedByAlpha1PlusAlpha2_1s2py<Atom1, Atom2>();
		static constexpr std::array<double, 9> oneDividedByAlpha1PlusAlpha2_1s2pz = GenerateOneDividedByAlpha1PlusAlpha2_1s2pz<Atom1, Atom2>();
		static constexpr std::array<double, 9> oneDividedByAlpha1PlusAlpha2_2s2px = GenerateOneDividedByAlpha1PlusAlpha2_2s2px<Atom1, Atom2>();
		static constexpr std::array<double, 9> oneDividedByAlpha1PlusAlpha2_2s2py = GenerateOneDividedByAlpha1PlusAlpha2_2s2py<Atom1, Atom2>();
		static constexpr std::array<double, 9> oneDividedByAlpha1PlusAlpha2_2s2pz = GenerateOneDividedByAlpha1PlusAlpha2_2s2pz<Atom1, Atom2>();

		static constexpr std::array<double, 9> oneDividedByAlpha1PlusAlpha2_2px1s  = GenerateOneDividedByAlpha1PlusAlpha2_1s2px<Atom2, Atom1>(); // Note the swap because we want the 2px on Atom1 and 1s on Atom2
		static constexpr std::array<double, 9> oneDividedByAlpha1PlusAlpha2_2px2s  = GenerateOneDividedByAlpha1PlusAlpha2_2s2px<Atom2, Atom1>(); // Note the swap because we want the 2px on Atom1 and 2s on Atom2
		static constexpr std::array<double, 9> oneDividedByAlpha1PlusAlpha2_2px2px = GenerateOneDividedByAlpha1PlusAlpha2_2px2px<Atom1, Atom2>();
		static constexpr std::array<double, 9> oneDividedByAlpha1PlusAlpha2_2px2py = GenerateOneDividedByAlpha1PlusAlpha2_2px2py<Atom1, Atom2>();
		static constexpr std::array<double, 9> oneDividedByAlpha1PlusAlpha2_2px2pz = GenerateOneDividedByAlpha1PlusAlpha2_2px2pz<Atom1, Atom2>();

		static constexpr std::array<double, 9> oneDividedByAlpha1PlusAlpha2_2py1s  = GenerateOneDividedByAlpha1PlusAlpha2_1s2py<Atom2, Atom1>(); // Note the swap because we want the 2py on Atom1 and 1s on Atom2
		static constexpr std::array<double, 9> oneDividedByAlpha1PlusAlpha2_2py2s  = GenerateOneDividedByAlpha1PlusAlpha2_2s2py<Atom2, Atom1>(); // Note the swap because we want the 2py on Atom1 and 2s on Atom2
		static constexpr std::array<double, 9> oneDividedByAlpha1PlusAlpha2_2py2px = GenerateOneDividedByAlpha1PlusAlpha2_2px2py<Atom2, Atom1>(); // Note the swap because we want the 2py on Atom1 and 2px on Atom2
		static constexpr std::array<double, 9> oneDividedByAlpha1PlusAlpha2_2py2py = GenerateOneDividedByAlpha1PlusAlpha2_2py2py<Atom1, Atom2>();
		static constexpr std::array<double, 9> oneDividedByAlpha1PlusAlpha2_2py2pz = GenerateOneDividedByAlpha1PlusAlpha2_2py2pz<Atom1, Atom2>();

		static constexpr std::array<double, 9> oneDividedByAlpha1PlusAlpha2_2pz1s  = GenerateOneDividedByAlpha1PlusAlpha2_1s2pz<Atom2, Atom1>(); // Note the swap because we want the 2pz on Atom1 and 1s on Atom2
		static constexpr std::array<double, 9> oneDividedByAlpha1PlusAlpha2_2pz2s  = GenerateOneDividedByAlpha1PlusAlpha2_2s2pz<Atom2, Atom1>(); // Note the swap because we want the 2pz on Atom1 and 2s on Atom2
		static constexpr std::array<double, 9> oneDividedByAlpha1PlusAlpha2_2pz2px = GenerateOneDividedByAlpha1PlusAlpha2_2px2pz<Atom2, Atom1>(); // Note the swap because we want the 2pz on Atom1 and 2px on Atom2
		static constexpr std::array<double, 9> oneDividedByAlpha1PlusAlpha2_2pz2py = GenerateOneDividedByAlpha1PlusAlpha2_2py2pz<Atom2, Atom1>(); // Note the swap because we want the 2pz on Atom1 and 2py on Atom2
		static constexpr std::array<double, 9> oneDividedByAlpha1PlusAlpha2_2pz2pz = GenerateOneDividedByAlpha1PlusAlpha2_2pz2pz<Atom1, Atom2>();

		// Factors
		static constexpr std::array<double, 9> factors_1s1s  = GenerateOverlapFactors_1s1s<Atom1, Atom2>();
		static constexpr std::array<double, 9> factors_1s2s  = GenerateOverlapFactors_1s2s<Atom1, Atom2>();
		static constexpr std::array<double, 9> factors_1s2px = GenerateOverlapFactors_1s2px<Atom1, Atom2>();
		static constexpr std::array<double, 9> factors_1s2py = GenerateOverlapFactors_1s2py<Atom1, Atom2>();
		static constexpr std::array<double, 9> factors_1s2pz = GenerateOverlapFactors_1s2pz<Atom1, Atom2>();

		static constexpr std::array<double, 9> factors_2s1s  = GenerateOverlapFactors_1s2s<Atom2, Atom1>(); // Note the swap because we want the 2s on Atom1 and 1s on Atom2
		static constexpr std::array<double, 9> factors_2s2s  = GenerateOverlapFactors_2s2s<Atom1, Atom2>();
		static constexpr std::array<double, 9> factors_2s2px = GenerateOverlapFactors_2s2px<Atom1, Atom2>();
		static constexpr std::array<double, 9> factors_2s2py = GenerateOverlapFactors_2s2py<Atom1, Atom2>();
		static constexpr std::array<double, 9> factors_2s2pz = GenerateOverlapFactors_2s2pz<Atom1, Atom2>();

		static constexpr std::array<double, 9> factors_2px1s  = GenerateOverlapFactors_1s2px<Atom2, Atom1>(); // Note the swap because we want the 2px on Atom1 and 1s on Atom2
		static constexpr std::array<double, 9> factors_2px2s  = GenerateOverlapFactors_2s2px<Atom2, Atom1>(); // Note the swap because we want the 2px on Atom1 and 2s on Atom2
		static constexpr std::array<double, 9> factors_2px2px = GenerateOverlapFactors_2px2px<Atom1, Atom2>();
		static constexpr std::array<double, 9> factors_2px2py = GenerateOverlapFactors_2px2py<Atom1, Atom2>();
		static constexpr std::array<double, 9> factors_2px2pz = GenerateOverlapFactors_2px2pz<Atom1, Atom2>();

		static constexpr std::array<double, 9> factors_2py1s  = GenerateOverlapFactors_1s2py<Atom2, Atom1>(); // Note the swap because we want the 2py on Atom1 and 1s on Atom2
		static constexpr std::array<double, 9> factors_2py2s  = GenerateOverlapFactors_2s2py<Atom2, Atom1>(); // Note the swap because we want the 2py on Atom1 and 2s on Atom2
		static constexpr std::array<double, 9> factors_2py2px = GenerateOverlapFactors_2px2py<Atom2, Atom1>(); // Note the swap because we want the 2py on Atom1 and 2px on Atom2
		static constexpr std::array<double, 9> factors_2py2py = GenerateOverlapFactors_2py2py<Atom1, Atom2>();
		static constexpr std::array<double, 9> factors_2py2pz = GenerateOverlapFactors_2py2pz<Atom1, Atom2>();

		static constexpr std::array<double, 9> factors_2pz1s  = GenerateOverlapFactors_1s2pz<Atom2, Atom1>(); // Note the swap because we want the 2pz on Atom1 and 1s on Atom2
		static constexpr std::array<double, 9> factors_2pz2s  = GenerateOverlapFactors_2s2pz<Atom2, Atom1>(); // Note the swap because we want the 2pz on Atom1 and 2s on Atom2
		static constexpr std::array<double, 9> factors_2pz2px = GenerateOverlapFactors_2px2pz<Atom2, Atom1>(); // Note the swap because we want the 2pz on Atom1 and 2px on Atom2
		static constexpr std::array<double, 9> factors_2pz2py = GenerateOverlapFactors_2py2pz<Atom2, Atom1>(); // Note the swap because we want the 2pz on Atom1 and 2py on Atom2
		static constexpr std::array<double, 9> factors_2pz2pz = GenerateOverlapFactors_2pz2pz<Atom1, Atom2>();

		// Exponential Factors
		static constexpr std::array<double, 9> expFactors_1s1s  = GenerateOverlapExponentialFactors_1s1s<Atom1, Atom2>();
		static constexpr std::array<double, 9> expFactors_1s2s  = GenerateOverlapExponentialFactors_1s2s<Atom1, Atom2>();
		static constexpr std::array<double, 9> expFactors_1s2px = GenerateOverlapExponentialFactors_1s2px<Atom1, Atom2>();
		static constexpr std::array<double, 9> expFactors_1s2py = GenerateOverlapExponentialFactors_1s2py<Atom1, Atom2>();
		static constexpr std::array<double, 9> expFactors_1s2pz = GenerateOverlapExponentialFactors_1s2pz<Atom1, Atom2>();

		static constexpr std::array<double, 9> expFactors_2s1s  = GenerateOverlapExponentialFactors_1s2s<Atom2, Atom1>(); // Note the swap because we want the 2s on Atom1 and 1s on Atom2
		static constexpr std::array<double, 9> expFactors_2s2s  = GenerateOverlapExponentialFactors_2s2s<Atom1, Atom2>();
		static constexpr std::array<double, 9> expFactors_2s2px = GenerateOverlapExponentialFactors_2s2px<Atom1, Atom2>();
		static constexpr std::array<double, 9> expFactors_2s2py = GenerateOverlapExponentialFactors_2s2py<Atom1, Atom2>();
		static constexpr std::array<double, 9> expFactors_2s2pz = GenerateOverlapExponentialFactors_2s2pz<Atom1, Atom2>();

		static constexpr std::array<double, 9> expFactors_2px1s  = GenerateOverlapExponentialFactors_1s2px<Atom2, Atom1>(); // Note the swap because we want the 2px on Atom1 and 1s on Atom2
		static constexpr std::array<double, 9> expFactors_2px2s  = GenerateOverlapExponentialFactors_2s2px<Atom2, Atom1>(); // Note the swap because we want the 2px on Atom1 and 2s on Atom2
		static constexpr std::array<double, 9> expFactors_2px2px = GenerateOverlapExponentialFactors_2px2px<Atom1, Atom2>();
		static constexpr std::array<double, 9> expFactors_2px2py = GenerateOverlapExponentialFactors_2px2py<Atom1, Atom2>();
		static constexpr std::array<double, 9> expFactors_2px2pz = GenerateOverlapExponentialFactors_2px2pz<Atom1, Atom2>();

		static constexpr std::array<double, 9> expFactors_2py1s  = GenerateOverlapExponentialFactors_1s2py<Atom2, Atom1>(); // Note the swap because we want the 2py on Atom1 and 1s on Atom2
		static constexpr std::array<double, 9> expFactors_2py2s  = GenerateOverlapExponentialFactors_2s2py<Atom2, Atom1>(); // Note the swap because we want the 2py on Atom1 and 2s on Atom2
		static constexpr std::array<double, 9> expFactors_2py2px = GenerateOverlapExponentialFactors_2px2py<Atom2, Atom1>(); // Note the swap because we want the 2py on Atom1 and 2px on Atom2
		static constexpr std::array<double, 9> expFactors_2py2py = GenerateOverlapExponentialFactors_2py2py<Atom1, Atom2>();
		static constexpr std::array<double, 9> expFactors_2py2pz = GenerateOverlapExponentialFactors_2py2pz<Atom1, Atom2>();

		static constexpr std::array<double, 9> expFactors_2pz1s  = GenerateOverlapExponentialFactors_1s2pz<Atom2, Atom1>(); // Note the swap because we want the 2pz on Atom1 and 1s on Atom2
		static constexpr std::array<double, 9> expFactors_2pz2s  = GenerateOverlapExponentialFactors_2s2pz<Atom2, Atom1>(); // Note the swap because we want the 2pz on Atom1 and 2s on Atom2
		static constexpr std::array<double, 9> expFactors_2pz2px = GenerateOverlapExponentialFactors_2px2pz<Atom2, Atom1>(); // Note the swap because we want the 2pz on Atom1 and 2px on Atom2
		static constexpr std::array<double, 9> expFactors_2pz2py = GenerateOverlapExponentialFactors_2py2pz<Atom2, Atom1>(); // Note the swap because we want the 2pz on Atom1 and 2py on Atom2
		static constexpr std::array<double, 9> expFactors_2pz2pz = GenerateOverlapExponentialFactors_2pz2pz<Atom1, Atom2>();

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

		// O_1s - O_2px  ==============================================================
		{
			const auto& atom1Prims = Atom1::orbital_1s.primitiveGaussians;
			const auto& atom2Prims = Atom2::orbital_2px.primitiveGaussians;
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

		// O_1s - O_2py  ==============================================================
		{
			const auto& atom1Prims = Atom1::orbital_1s.primitiveGaussians;
			const auto& atom2Prims = Atom2::orbital_2py.primitiveGaussians;
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

		// O_1s - O_2pz  ==============================================================
		{
			const auto& atom1Prims = Atom1::orbital_1s.primitiveGaussians;
			const auto& atom2Prims = Atom2::orbital_2pz.primitiveGaussians;
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

		// O_2s - O_1s  ==============================================================
		result = 0.0;
		result += factors_2s1s[0] * std::exp(expFactors_2s1s[0] * dotProduct);
		result += factors_2s1s[1] * std::exp(expFactors_2s1s[1] * dotProduct);
		result += factors_2s1s[2] * std::exp(expFactors_2s1s[2] * dotProduct);
		result += factors_2s1s[3] * std::exp(expFactors_2s1s[3] * dotProduct);
		result += factors_2s1s[4] * std::exp(expFactors_2s1s[4] * dotProduct);
		result += factors_2s1s[5] * std::exp(expFactors_2s1s[5] * dotProduct);
		result += factors_2s1s[6] * std::exp(expFactors_2s1s[6] * dotProduct);
		result += factors_2s1s[7] * std::exp(expFactors_2s1s[7] * dotProduct);
		result += factors_2s1s[8] * std::exp(expFactors_2s1s[8] * dotProduct);

		overlapMatrix(row + 1, col) = result;

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
			const auto& atom1Prims = Atom1::orbital_2s.primitiveGaussians;
			const auto& atom2Prims = Atom2::orbital_2px.primitiveGaussians;
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

		// O_2s - O_2py  ==============================================================
		{
			const auto& atom1Prims = Atom1::orbital_2s.primitiveGaussians;
			const auto& atom2Prims = Atom2::orbital_2py.primitiveGaussians;
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

		// O_2s - O_2pz  ==============================================================
		{
			const auto& atom1Prims = Atom1::orbital_2s.primitiveGaussians;
			const auto& atom2Prims = Atom2::orbital_2pz.primitiveGaussians;
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

		// O_2px - O_1s  ==============================================================
		{
			const auto& atom1Prims = Atom1::orbital_2px.primitiveGaussians;
			const auto& atom2Prims = Atom2::orbital_1s.primitiveGaussians;
			result = 0.0;
			// Note the ordering of the indices on the factors is rearranged. This is because the function that computes them orders them assuming
			// the 2px atom is the second atom, but it is actually the first, which leads to the weird ordering
			result += factors_2px1s[0] * std::exp(expFactors_2px1s[0] * dotProduct) * (-position1.x + ((position1.x * atom1Prims[0].alpha + position2.x * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2px1s[0]));
			result += factors_2px1s[3] * std::exp(expFactors_2px1s[3] * dotProduct) * (-position1.x + ((position1.x * atom1Prims[0].alpha + position2.x * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2px1s[3]));
			result += factors_2px1s[6] * std::exp(expFactors_2px1s[6] * dotProduct) * (-position1.x + ((position1.x * atom1Prims[0].alpha + position2.x * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2px1s[6]));
			result += factors_2px1s[1] * std::exp(expFactors_2px1s[1] * dotProduct) * (-position1.x + ((position1.x * atom1Prims[1].alpha + position2.x * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2px1s[1]));
			result += factors_2px1s[4] * std::exp(expFactors_2px1s[4] * dotProduct) * (-position1.x + ((position1.x * atom1Prims[1].alpha + position2.x * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2px1s[4]));
			result += factors_2px1s[7] * std::exp(expFactors_2px1s[7] * dotProduct) * (-position1.x + ((position1.x * atom1Prims[1].alpha + position2.x * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2px1s[7]));
			result += factors_2px1s[2] * std::exp(expFactors_2px1s[2] * dotProduct) * (-position1.x + ((position1.x * atom1Prims[2].alpha + position2.x * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2px1s[2]));
			result += factors_2px1s[5] * std::exp(expFactors_2px1s[5] * dotProduct) * (-position1.x + ((position1.x * atom1Prims[2].alpha + position2.x * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2px1s[5]));
			result += factors_2px1s[8] * std::exp(expFactors_2px1s[8] * dotProduct) * (-position1.x + ((position1.x * atom1Prims[2].alpha + position2.x * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2px1s[8]));
		}
		overlapMatrix(row + 2, col) = result;

		// O_2px - O_2s  ==============================================================
		{
			const auto& atom1Prims = Atom1::orbital_2px.primitiveGaussians;
			const auto& atom2Prims = Atom2::orbital_2s.primitiveGaussians;
			result = 0.0;
			// Note the ordering of the indices on the factors is rearranged. This is because the function that computes them orders them assuming
			// the 2px atom is the second atom, but it is actually the first, which leads to the weird ordering
			result += factors_2px2s[0] * std::exp(expFactors_2px2s[0] * dotProduct) * (-position1.x + ((position1.x * atom1Prims[0].alpha + position2.x * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2px2s[0]));
			result += factors_2px2s[3] * std::exp(expFactors_2px2s[3] * dotProduct) * (-position1.x + ((position1.x * atom1Prims[0].alpha + position2.x * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2px2s[3]));
			result += factors_2px2s[6] * std::exp(expFactors_2px2s[6] * dotProduct) * (-position1.x + ((position1.x * atom1Prims[0].alpha + position2.x * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2px2s[6]));
			result += factors_2px2s[1] * std::exp(expFactors_2px2s[1] * dotProduct) * (-position1.x + ((position1.x * atom1Prims[1].alpha + position2.x * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2px2s[1]));
			result += factors_2px2s[4] * std::exp(expFactors_2px2s[4] * dotProduct) * (-position1.x + ((position1.x * atom1Prims[1].alpha + position2.x * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2px2s[4]));
			result += factors_2px2s[7] * std::exp(expFactors_2px2s[7] * dotProduct) * (-position1.x + ((position1.x * atom1Prims[1].alpha + position2.x * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2px2s[7]));
			result += factors_2px2s[2] * std::exp(expFactors_2px2s[2] * dotProduct) * (-position1.x + ((position1.x * atom1Prims[2].alpha + position2.x * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2px2s[2]));
			result += factors_2px2s[5] * std::exp(expFactors_2px2s[5] * dotProduct) * (-position1.x + ((position1.x * atom1Prims[2].alpha + position2.x * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2px2s[5]));
			result += factors_2px2s[8] * std::exp(expFactors_2px2s[8] * dotProduct) * (-position1.x + ((position1.x * atom1Prims[2].alpha + position2.x * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2px2s[8]));
		}
		overlapMatrix(row + 2, col + 1) = result;

		// O_2px - O_2px  ==============================================================
		{
			const auto& atom1Prims = Atom1::orbital_2px.primitiveGaussians;
			const auto& atom2Prims = Atom2::orbital_2px.primitiveGaussians;
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

		// O_2px - O_2py  ==============================================================
		{
			const auto& atom1Prims = Atom1::orbital_2px.primitiveGaussians;
			const auto& atom2Prims = Atom2::orbital_2py.primitiveGaussians;
			result = 0.0;
			result += factors_2px2py[0] * std::exp(expFactors_2px2py[0] * dotProduct) * (-position1.x + ((position1.x * atom1Prims[0].alpha + position2.x * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2px2py[0])) * (-position2.y + ((position1.y * atom1Prims[0].alpha + position2.y * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2px2py[0]));
			result += factors_2px2py[1] * std::exp(expFactors_2px2py[1] * dotProduct) * (-position1.x + ((position1.x * atom1Prims[0].alpha + position2.x * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2px2py[1])) * (-position2.y + ((position1.y * atom1Prims[0].alpha + position2.y * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2px2py[1]));
			result += factors_2px2py[2] * std::exp(expFactors_2px2py[2] * dotProduct) * (-position1.x + ((position1.x * atom1Prims[0].alpha + position2.x * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2px2py[2])) * (-position2.y + ((position1.y * atom1Prims[0].alpha + position2.y * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2px2py[2]));
			result += factors_2px2py[3] * std::exp(expFactors_2px2py[3] * dotProduct) * (-position1.x + ((position1.x * atom1Prims[1].alpha + position2.x * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2px2py[3])) * (-position2.y + ((position1.y * atom1Prims[1].alpha + position2.y * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2px2py[3]));
			result += factors_2px2py[4] * std::exp(expFactors_2px2py[4] * dotProduct) * (-position1.x + ((position1.x * atom1Prims[1].alpha + position2.x * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2px2py[4])) * (-position2.y + ((position1.y * atom1Prims[1].alpha + position2.y * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2px2py[4]));
			result += factors_2px2py[5] * std::exp(expFactors_2px2py[5] * dotProduct) * (-position1.x + ((position1.x * atom1Prims[1].alpha + position2.x * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2px2py[5])) * (-position2.y + ((position1.y * atom1Prims[1].alpha + position2.y * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2px2py[5]));
			result += factors_2px2py[6] * std::exp(expFactors_2px2py[6] * dotProduct) * (-position1.x + ((position1.x * atom1Prims[2].alpha + position2.x * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2px2py[6])) * (-position2.y + ((position1.y * atom1Prims[2].alpha + position2.y * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2px2py[6]));
			result += factors_2px2py[7] * std::exp(expFactors_2px2py[7] * dotProduct) * (-position1.x + ((position1.x * atom1Prims[2].alpha + position2.x * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2px2py[7])) * (-position2.y + ((position1.y * atom1Prims[2].alpha + position2.y * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2px2py[7]));
			result += factors_2px2py[8] * std::exp(expFactors_2px2py[8] * dotProduct) * (-position1.x + ((position1.x * atom1Prims[2].alpha + position2.x * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2px2py[8])) * (-position2.y + ((position1.y * atom1Prims[2].alpha + position2.y * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2px2py[8]));
		}
		overlapMatrix(row + 2, col + 3) = result;

		// O_2px - O_2pz  ==============================================================
		{
			const auto& atom1Prims = Atom1::orbital_2px.primitiveGaussians;
			const auto& atom2Prims = Atom2::orbital_2pz.primitiveGaussians;
			result = 0.0;
			result += factors_2px2pz[0] * std::exp(expFactors_2px2pz[0] * dotProduct) * (-position1.x + ((position1.x * atom1Prims[0].alpha + position2.x * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2px2pz[0])) * (-position2.z + ((position1.z * atom1Prims[0].alpha + position2.z * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2px2pz[0]));
			result += factors_2px2pz[1] * std::exp(expFactors_2px2pz[1] * dotProduct) * (-position1.x + ((position1.x * atom1Prims[0].alpha + position2.x * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2px2pz[1])) * (-position2.z + ((position1.z * atom1Prims[0].alpha + position2.z * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2px2pz[1]));
			result += factors_2px2pz[2] * std::exp(expFactors_2px2pz[2] * dotProduct) * (-position1.x + ((position1.x * atom1Prims[0].alpha + position2.x * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2px2pz[2])) * (-position2.z + ((position1.z * atom1Prims[0].alpha + position2.z * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2px2pz[2]));
			result += factors_2px2pz[3] * std::exp(expFactors_2px2pz[3] * dotProduct) * (-position1.x + ((position1.x * atom1Prims[1].alpha + position2.x * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2px2pz[3])) * (-position2.z + ((position1.z * atom1Prims[1].alpha + position2.z * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2px2pz[3]));
			result += factors_2px2pz[4] * std::exp(expFactors_2px2pz[4] * dotProduct) * (-position1.x + ((position1.x * atom1Prims[1].alpha + position2.x * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2px2pz[4])) * (-position2.z + ((position1.z * atom1Prims[1].alpha + position2.z * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2px2pz[4]));
			result += factors_2px2pz[5] * std::exp(expFactors_2px2pz[5] * dotProduct) * (-position1.x + ((position1.x * atom1Prims[1].alpha + position2.x * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2px2pz[5])) * (-position2.z + ((position1.z * atom1Prims[1].alpha + position2.z * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2px2pz[5]));
			result += factors_2px2pz[6] * std::exp(expFactors_2px2pz[6] * dotProduct) * (-position1.x + ((position1.x * atom1Prims[2].alpha + position2.x * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2px2pz[6])) * (-position2.z + ((position1.z * atom1Prims[2].alpha + position2.z * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2px2pz[6]));
			result += factors_2px2pz[7] * std::exp(expFactors_2px2pz[7] * dotProduct) * (-position1.x + ((position1.x * atom1Prims[2].alpha + position2.x * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2px2pz[7])) * (-position2.z + ((position1.z * atom1Prims[2].alpha + position2.z * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2px2pz[7]));
			result += factors_2px2pz[8] * std::exp(expFactors_2px2pz[8] * dotProduct) * (-position1.x + ((position1.x * atom1Prims[2].alpha + position2.x * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2px2pz[8])) * (-position2.z + ((position1.z * atom1Prims[2].alpha + position2.z * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2px2pz[8]));
		}
		overlapMatrix(row + 2, col + 4) = result;

		// O_2py - O_1s  ==============================================================
		{
			const auto& atom1Prims = Atom1::orbital_2py.primitiveGaussians;
			const auto& atom2Prims = Atom2::orbital_1s.primitiveGaussians;
			result = 0.0;
			// Note the ordering of the indices on the factors is rearranged. This is because the function that computes them orders them assuming
			// the 2py atom is the second atom, but it is actually the first, which leads to the weird ordering
			result += factors_2py1s[0] * std::exp(expFactors_2py1s[0] * dotProduct) * (-position1.y + ((position1.y * atom1Prims[0].alpha + position2.y * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2py1s[0]));
			result += factors_2py1s[3] * std::exp(expFactors_2py1s[3] * dotProduct) * (-position1.y + ((position1.y * atom1Prims[0].alpha + position2.y * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2py1s[3]));
			result += factors_2py1s[6] * std::exp(expFactors_2py1s[6] * dotProduct) * (-position1.y + ((position1.y * atom1Prims[0].alpha + position2.y * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2py1s[6]));
			result += factors_2py1s[1] * std::exp(expFactors_2py1s[1] * dotProduct) * (-position1.y + ((position1.y * atom1Prims[1].alpha + position2.y * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2py1s[1]));
			result += factors_2py1s[4] * std::exp(expFactors_2py1s[4] * dotProduct) * (-position1.y + ((position1.y * atom1Prims[1].alpha + position2.y * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2py1s[4]));
			result += factors_2py1s[7] * std::exp(expFactors_2py1s[7] * dotProduct) * (-position1.y + ((position1.y * atom1Prims[1].alpha + position2.y * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2py1s[7]));
			result += factors_2py1s[2] * std::exp(expFactors_2py1s[2] * dotProduct) * (-position1.y + ((position1.y * atom1Prims[2].alpha + position2.y * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2py1s[2]));
			result += factors_2py1s[5] * std::exp(expFactors_2py1s[5] * dotProduct) * (-position1.y + ((position1.y * atom1Prims[2].alpha + position2.y * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2py1s[5]));
			result += factors_2py1s[8] * std::exp(expFactors_2py1s[8] * dotProduct) * (-position1.y + ((position1.y * atom1Prims[2].alpha + position2.y * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2py1s[8]));
		}
		overlapMatrix(row + 3, col) = result;

		// O_2py - O_2s  ==============================================================
		{
			const auto& atom1Prims = Atom1::orbital_2py.primitiveGaussians;
			const auto& atom2Prims = Atom2::orbital_2s.primitiveGaussians;
			result = 0.0;
			// Note the ordering of the indices on the factors is rearranged. This is because the function that computes them orders them assuming
			// the 2py atom is the second atom, but it is actually the first, which leads to the weird ordering
			result += factors_2py2s[0] * std::exp(expFactors_2py2s[0] * dotProduct) * (-position1.y + ((position1.y * atom1Prims[0].alpha + position2.y * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2py2s[0]));
			result += factors_2py2s[3] * std::exp(expFactors_2py2s[3] * dotProduct) * (-position1.y + ((position1.y * atom1Prims[0].alpha + position2.y * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2py2s[3]));
			result += factors_2py2s[6] * std::exp(expFactors_2py2s[6] * dotProduct) * (-position1.y + ((position1.y * atom1Prims[0].alpha + position2.y * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2py2s[6]));
			result += factors_2py2s[1] * std::exp(expFactors_2py2s[1] * dotProduct) * (-position1.y + ((position1.y * atom1Prims[1].alpha + position2.y * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2py2s[1]));
			result += factors_2py2s[4] * std::exp(expFactors_2py2s[4] * dotProduct) * (-position1.y + ((position1.y * atom1Prims[1].alpha + position2.y * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2py2s[4]));
			result += factors_2py2s[7] * std::exp(expFactors_2py2s[7] * dotProduct) * (-position1.y + ((position1.y * atom1Prims[1].alpha + position2.y * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2py2s[7]));
			result += factors_2py2s[2] * std::exp(expFactors_2py2s[2] * dotProduct) * (-position1.y + ((position1.y * atom1Prims[2].alpha + position2.y * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2py2s[2]));
			result += factors_2py2s[5] * std::exp(expFactors_2py2s[5] * dotProduct) * (-position1.y + ((position1.y * atom1Prims[2].alpha + position2.y * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2py2s[5]));
			result += factors_2py2s[8] * std::exp(expFactors_2py2s[8] * dotProduct) * (-position1.y + ((position1.y * atom1Prims[2].alpha + position2.y * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2py2s[8]));
		}
		overlapMatrix(row + 3, col + 1) = result;

		// O_2py - O_2px  ==============================================================
		{
			const auto& atom1Prims = Atom1::orbital_2py.primitiveGaussians;
			const auto& atom2Prims = Atom2::orbital_2px.primitiveGaussians;
			result = 0.0;
			// Note the ordering of the indices on the factors is rearranged. This is because the function that computes them orders them assuming
			// the 2py atom is the second atom, but it is actually the first, which leads to the weird ordering
			result += factors_2py2px[0] * std::exp(expFactors_2py2px[0] * dotProduct) * (-position2.x + ((position1.x * atom1Prims[0].alpha + position2.x * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2py2px[0])) * (-position1.y + ((position1.y * atom1Prims[0].alpha + position2.y * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2py2px[0]));
			result += factors_2py2px[3] * std::exp(expFactors_2py2px[3] * dotProduct) * (-position2.x + ((position1.x * atom1Prims[0].alpha + position2.x * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2py2px[3])) * (-position1.y + ((position1.y * atom1Prims[0].alpha + position2.y * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2py2px[3]));
			result += factors_2py2px[6] * std::exp(expFactors_2py2px[6] * dotProduct) * (-position2.x + ((position1.x * atom1Prims[0].alpha + position2.x * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2py2px[6])) * (-position1.y + ((position1.y * atom1Prims[0].alpha + position2.y * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2py2px[6]));
			result += factors_2py2px[1] * std::exp(expFactors_2py2px[1] * dotProduct) * (-position2.x + ((position1.x * atom1Prims[1].alpha + position2.x * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2py2px[1])) * (-position1.y + ((position1.y * atom1Prims[1].alpha + position2.y * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2py2px[1]));
			result += factors_2py2px[4] * std::exp(expFactors_2py2px[4] * dotProduct) * (-position2.x + ((position1.x * atom1Prims[1].alpha + position2.x * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2py2px[4])) * (-position1.y + ((position1.y * atom1Prims[1].alpha + position2.y * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2py2px[4]));
			result += factors_2py2px[7] * std::exp(expFactors_2py2px[7] * dotProduct) * (-position2.x + ((position1.x * atom1Prims[1].alpha + position2.x * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2py2px[7])) * (-position1.y + ((position1.y * atom1Prims[1].alpha + position2.y * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2py2px[7]));
			result += factors_2py2px[2] * std::exp(expFactors_2py2px[2] * dotProduct) * (-position2.x + ((position1.x * atom1Prims[2].alpha + position2.x * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2py2px[2])) * (-position1.y + ((position1.y * atom1Prims[2].alpha + position2.y * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2py2px[2]));
			result += factors_2py2px[5] * std::exp(expFactors_2py2px[5] * dotProduct) * (-position2.x + ((position1.x * atom1Prims[2].alpha + position2.x * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2py2px[5])) * (-position1.y + ((position1.y * atom1Prims[2].alpha + position2.y * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2py2px[5]));
			result += factors_2py2px[8] * std::exp(expFactors_2py2px[8] * dotProduct) * (-position2.x + ((position1.x * atom1Prims[2].alpha + position2.x * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2py2px[8])) * (-position1.y + ((position1.y * atom1Prims[2].alpha + position2.y * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2py2px[8]));
		}
		overlapMatrix(row + 3, col + 2) = result;

		// O_2py - O_2py  ==============================================================
		{
			const auto& atom1Prims = Atom1::orbital_2py.primitiveGaussians;
			const auto& atom2Prims = Atom2::orbital_2py.primitiveGaussians;
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


		// O_2py - O_2pz  ==============================================================
		{
			const auto& atom1Prims = Atom1::orbital_2py.primitiveGaussians;
			const auto& atom2Prims = Atom2::orbital_2pz.primitiveGaussians;
			result = 0.0;
			result += factors_2py2pz[0] * std::exp(expFactors_2py2pz[0] * dotProduct) * (-position1.y + ((position1.y * atom1Prims[0].alpha + position2.y * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2py2pz[0])) * (-position2.z + ((position1.z * atom1Prims[0].alpha + position2.z * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2py2pz[0]));
			result += factors_2py2pz[1] * std::exp(expFactors_2py2pz[1] * dotProduct) * (-position1.y + ((position1.y * atom1Prims[0].alpha + position2.y * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2py2pz[1])) * (-position2.z + ((position1.z * atom1Prims[0].alpha + position2.z * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2py2pz[1]));
			result += factors_2py2pz[2] * std::exp(expFactors_2py2pz[2] * dotProduct) * (-position1.y + ((position1.y * atom1Prims[0].alpha + position2.y * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2py2pz[2])) * (-position2.z + ((position1.z * atom1Prims[0].alpha + position2.z * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2py2pz[2]));
			result += factors_2py2pz[3] * std::exp(expFactors_2py2pz[3] * dotProduct) * (-position1.y + ((position1.y * atom1Prims[1].alpha + position2.y * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2py2pz[3])) * (-position2.z + ((position1.z * atom1Prims[1].alpha + position2.z * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2py2pz[3]));
			result += factors_2py2pz[4] * std::exp(expFactors_2py2pz[4] * dotProduct) * (-position1.y + ((position1.y * atom1Prims[1].alpha + position2.y * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2py2pz[4])) * (-position2.z + ((position1.z * atom1Prims[1].alpha + position2.z * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2py2pz[4]));
			result += factors_2py2pz[5] * std::exp(expFactors_2py2pz[5] * dotProduct) * (-position1.y + ((position1.y * atom1Prims[1].alpha + position2.y * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2py2pz[5])) * (-position2.z + ((position1.z * atom1Prims[1].alpha + position2.z * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2py2pz[5]));
			result += factors_2py2pz[6] * std::exp(expFactors_2py2pz[6] * dotProduct) * (-position1.y + ((position1.y * atom1Prims[2].alpha + position2.y * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2py2pz[6])) * (-position2.z + ((position1.z * atom1Prims[2].alpha + position2.z * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2py2pz[6]));
			result += factors_2py2pz[7] * std::exp(expFactors_2py2pz[7] * dotProduct) * (-position1.y + ((position1.y * atom1Prims[2].alpha + position2.y * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2py2pz[7])) * (-position2.z + ((position1.z * atom1Prims[2].alpha + position2.z * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2py2pz[7]));
			result += factors_2py2pz[8] * std::exp(expFactors_2py2pz[8] * dotProduct) * (-position1.y + ((position1.y * atom1Prims[2].alpha + position2.y * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2py2pz[8])) * (-position2.z + ((position1.z * atom1Prims[2].alpha + position2.z * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2py2pz[8]));
		}
		overlapMatrix(row + 3, col + 4) = result;

		// O_2pz - O_1s  ==============================================================
		{
			const auto& atom1Prims = Atom1::orbital_2pz.primitiveGaussians;
			const auto& atom2Prims = Atom2::orbital_1s.primitiveGaussians;
			result = 0.0;
			// Note the ordering of the indices on the factors is rearranged. This is because the function that computes them orders them assuming
			// the 2pz atom is the second atom, but it is actually the first, which leads to the weird ordering
			result += factors_2pz1s[0] * std::exp(expFactors_2pz1s[0] * dotProduct) * (-position1.z + ((position1.z * atom1Prims[0].alpha + position2.z * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2pz1s[0]));
			result += factors_2pz1s[3] * std::exp(expFactors_2pz1s[3] * dotProduct) * (-position1.z + ((position1.z * atom1Prims[0].alpha + position2.z * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2pz1s[3]));
			result += factors_2pz1s[6] * std::exp(expFactors_2pz1s[6] * dotProduct) * (-position1.z + ((position1.z * atom1Prims[0].alpha + position2.z * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2pz1s[6]));
			result += factors_2pz1s[1] * std::exp(expFactors_2pz1s[1] * dotProduct) * (-position1.z + ((position1.z * atom1Prims[1].alpha + position2.z * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2pz1s[1]));
			result += factors_2pz1s[4] * std::exp(expFactors_2pz1s[4] * dotProduct) * (-position1.z + ((position1.z * atom1Prims[1].alpha + position2.z * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2pz1s[4]));
			result += factors_2pz1s[7] * std::exp(expFactors_2pz1s[7] * dotProduct) * (-position1.z + ((position1.z * atom1Prims[1].alpha + position2.z * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2pz1s[7]));
			result += factors_2pz1s[2] * std::exp(expFactors_2pz1s[2] * dotProduct) * (-position1.z + ((position1.z * atom1Prims[2].alpha + position2.z * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2pz1s[2]));
			result += factors_2pz1s[5] * std::exp(expFactors_2pz1s[5] * dotProduct) * (-position1.z + ((position1.z * atom1Prims[2].alpha + position2.z * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2pz1s[5]));
			result += factors_2pz1s[8] * std::exp(expFactors_2pz1s[8] * dotProduct) * (-position1.z + ((position1.z * atom1Prims[2].alpha + position2.z * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2pz1s[8]));
		}
		overlapMatrix(row + 4, col) = result;

		// O_2pz - O_2s  ==============================================================
		{
			const auto& atom1Prims = Atom1::orbital_2pz.primitiveGaussians;
			const auto& atom2Prims = Atom2::orbital_2s.primitiveGaussians;
			result = 0.0;
			// Note the ordering of the indices on the factors is rearranged. This is because the function that computes them orders them assuming
			// the 2pz atom is the second atom, but it is actually the first, which leads to the weird ordering
			result += factors_2pz2s[0] * std::exp(expFactors_2pz2s[0] * dotProduct) * (-position1.z + ((position1.z * atom1Prims[0].alpha + position2.z * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2s[0]));
			result += factors_2pz2s[3] * std::exp(expFactors_2pz2s[3] * dotProduct) * (-position1.z + ((position1.z * atom1Prims[0].alpha + position2.z * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2s[3]));
			result += factors_2pz2s[6] * std::exp(expFactors_2pz2s[6] * dotProduct) * (-position1.z + ((position1.z * atom1Prims[0].alpha + position2.z * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2s[6]));
			result += factors_2pz2s[1] * std::exp(expFactors_2pz2s[1] * dotProduct) * (-position1.z + ((position1.z * atom1Prims[1].alpha + position2.z * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2s[1]));
			result += factors_2pz2s[4] * std::exp(expFactors_2pz2s[4] * dotProduct) * (-position1.z + ((position1.z * atom1Prims[1].alpha + position2.z * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2s[4]));
			result += factors_2pz2s[7] * std::exp(expFactors_2pz2s[7] * dotProduct) * (-position1.z + ((position1.z * atom1Prims[1].alpha + position2.z * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2s[7]));
			result += factors_2pz2s[2] * std::exp(expFactors_2pz2s[2] * dotProduct) * (-position1.z + ((position1.z * atom1Prims[2].alpha + position2.z * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2s[2]));
			result += factors_2pz2s[5] * std::exp(expFactors_2pz2s[5] * dotProduct) * (-position1.z + ((position1.z * atom1Prims[2].alpha + position2.z * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2s[5]));
			result += factors_2pz2s[8] * std::exp(expFactors_2pz2s[8] * dotProduct) * (-position1.z + ((position1.z * atom1Prims[2].alpha + position2.z * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2s[8]));
		}
		overlapMatrix(row + 4, col + 1) = result;

		// O_2pz - O_2px  ==============================================================
		{
			const auto& atom1Prims = Atom1::orbital_2pz.primitiveGaussians;
			const auto& atom2Prims = Atom2::orbital_2px.primitiveGaussians;
			result = 0.0;
			// Note the ordering of the indices on the factors is rearranged. This is because the function that computes them orders them assuming
			// the 2pz atom is the second atom, but it is actually the first, which leads to the weird ordering
			result += factors_2pz2px[0] * std::exp(expFactors_2pz2px[0] * dotProduct) * (-position2.x + ((position1.x * atom1Prims[0].alpha + position2.x * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2px[0])) * (-position1.z + ((position1.z * atom1Prims[0].alpha + position2.z * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2px[0]));
			result += factors_2pz2px[3] * std::exp(expFactors_2pz2px[3] * dotProduct) * (-position2.x + ((position1.x * atom1Prims[0].alpha + position2.x * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2px[3])) * (-position1.z + ((position1.z * atom1Prims[0].alpha + position2.z * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2px[3]));
			result += factors_2pz2px[6] * std::exp(expFactors_2pz2px[6] * dotProduct) * (-position2.x + ((position1.x * atom1Prims[0].alpha + position2.x * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2px[6])) * (-position1.z + ((position1.z * atom1Prims[0].alpha + position2.z * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2px[6]));
			result += factors_2pz2px[1] * std::exp(expFactors_2pz2px[1] * dotProduct) * (-position2.x + ((position1.x * atom1Prims[1].alpha + position2.x * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2px[1])) * (-position1.z + ((position1.z * atom1Prims[1].alpha + position2.z * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2px[1]));
			result += factors_2pz2px[4] * std::exp(expFactors_2pz2px[4] * dotProduct) * (-position2.x + ((position1.x * atom1Prims[1].alpha + position2.x * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2px[4])) * (-position1.z + ((position1.z * atom1Prims[1].alpha + position2.z * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2px[4]));
			result += factors_2pz2px[7] * std::exp(expFactors_2pz2px[7] * dotProduct) * (-position2.x + ((position1.x * atom1Prims[1].alpha + position2.x * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2px[7])) * (-position1.z + ((position1.z * atom1Prims[1].alpha + position2.z * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2px[7]));
			result += factors_2pz2px[2] * std::exp(expFactors_2pz2px[2] * dotProduct) * (-position2.x + ((position1.x * atom1Prims[2].alpha + position2.x * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2px[2])) * (-position1.z + ((position1.z * atom1Prims[2].alpha + position2.z * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2px[2]));
			result += factors_2pz2px[5] * std::exp(expFactors_2pz2px[5] * dotProduct) * (-position2.x + ((position1.x * atom1Prims[2].alpha + position2.x * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2px[5])) * (-position1.z + ((position1.z * atom1Prims[2].alpha + position2.z * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2px[5]));
			result += factors_2pz2px[8] * std::exp(expFactors_2pz2px[8] * dotProduct) * (-position2.x + ((position1.x * atom1Prims[2].alpha + position2.x * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2px[8])) * (-position1.z + ((position1.z * atom1Prims[2].alpha + position2.z * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2px[8]));
		}
		overlapMatrix(row + 4, col + 2) = result;

		// O_2pz - O_2py  ==============================================================
		{
			const auto& atom1Prims = Atom1::orbital_2pz.primitiveGaussians;
			const auto& atom2Prims = Atom2::orbital_2py.primitiveGaussians;
			result = 0.0;
			// Note the ordering of the indices on the factors is rearranged. This is because the function that computes them orders them assuming
			// the 2pz atom is the second atom, but it is actually the first, which leads to the weird ordering
			result += factors_2pz2py[0] * std::exp(expFactors_2pz2py[0] * dotProduct) * (-position2.y + ((position1.y * atom1Prims[0].alpha + position2.y * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2py[0])) * (-position1.z + ((position1.z * atom1Prims[0].alpha + position2.z * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2py[0]));
			result += factors_2pz2py[3] * std::exp(expFactors_2pz2py[3] * dotProduct) * (-position2.y + ((position1.y * atom1Prims[0].alpha + position2.y * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2py[3])) * (-position1.z + ((position1.z * atom1Prims[0].alpha + position2.z * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2py[3]));
			result += factors_2pz2py[6] * std::exp(expFactors_2pz2py[6] * dotProduct) * (-position2.y + ((position1.y * atom1Prims[0].alpha + position2.y * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2py[6])) * (-position1.z + ((position1.z * atom1Prims[0].alpha + position2.z * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2py[6]));
			result += factors_2pz2py[1] * std::exp(expFactors_2pz2py[1] * dotProduct) * (-position2.y + ((position1.y * atom1Prims[1].alpha + position2.y * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2py[1])) * (-position1.z + ((position1.z * atom1Prims[1].alpha + position2.z * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2py[1]));
			result += factors_2pz2py[4] * std::exp(expFactors_2pz2py[4] * dotProduct) * (-position2.y + ((position1.y * atom1Prims[1].alpha + position2.y * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2py[4])) * (-position1.z + ((position1.z * atom1Prims[1].alpha + position2.z * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2py[4]));
			result += factors_2pz2py[7] * std::exp(expFactors_2pz2py[7] * dotProduct) * (-position2.y + ((position1.y * atom1Prims[1].alpha + position2.y * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2py[7])) * (-position1.z + ((position1.z * atom1Prims[1].alpha + position2.z * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2py[7]));
			result += factors_2pz2py[2] * std::exp(expFactors_2pz2py[2] * dotProduct) * (-position2.y + ((position1.y * atom1Prims[2].alpha + position2.y * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2py[2])) * (-position1.z + ((position1.z * atom1Prims[2].alpha + position2.z * atom2Prims[0].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2py[2]));
			result += factors_2pz2py[5] * std::exp(expFactors_2pz2py[5] * dotProduct) * (-position2.y + ((position1.y * atom1Prims[2].alpha + position2.y * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2py[5])) * (-position1.z + ((position1.z * atom1Prims[2].alpha + position2.z * atom2Prims[1].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2py[5]));
			result += factors_2pz2py[8] * std::exp(expFactors_2pz2py[8] * dotProduct) * (-position2.y + ((position1.y * atom1Prims[2].alpha + position2.y * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2py[8])) * (-position1.z + ((position1.z * atom1Prims[2].alpha + position2.z * atom2Prims[2].alpha) * oneDividedByAlpha1PlusAlpha2_2pz2py[8]));
		}
		overlapMatrix(row + 4, col + 3) = result;

		// O_2pz - O_2pz  ==============================================================
		{
			const auto& atom1Prims = Atom1::orbital_2pz.primitiveGaussians;
			const auto& atom2Prims = Atom2::orbital_2pz.primitiveGaussians;
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


	// 1s - 1s
	template<typename Atom1, typename Atom2>
	consteval std::array<double, 9> GenerateOverlapFactors_1s1s() noexcept
	{
		constexpr double oneDividedByAlpha1PlusAlpha2_1 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[0].alpha + Atom2::orbital_1s.primitiveGaussians[0].alpha);
		constexpr double d1 =
			Atom1::orbital_1s.primitiveGaussians[0].normalizationFactor *
			Atom1::orbital_1s.primitiveGaussians[0].coefficient *
			Atom2::orbital_1s.primitiveGaussians[0].normalizationFactor *
			Atom2::orbital_1s.primitiveGaussians[0].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_1, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_2 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[0].alpha + Atom2::orbital_1s.primitiveGaussians[1].alpha);
		constexpr double d2 =
			Atom1::orbital_1s.primitiveGaussians[0].normalizationFactor *
			Atom1::orbital_1s.primitiveGaussians[0].coefficient *
			Atom2::orbital_1s.primitiveGaussians[1].normalizationFactor *
			Atom2::orbital_1s.primitiveGaussians[1].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_2, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_3 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[0].alpha + Atom2::orbital_1s.primitiveGaussians[2].alpha);
		constexpr double d3 =
			Atom1::orbital_1s.primitiveGaussians[0].normalizationFactor *
			Atom1::orbital_1s.primitiveGaussians[0].coefficient *
			Atom2::orbital_1s.primitiveGaussians[2].normalizationFactor *
			Atom2::orbital_1s.primitiveGaussians[2].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_3, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_4 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[1].alpha + Atom2::orbital_1s.primitiveGaussians[0].alpha);
		constexpr double d4 =
			Atom1::orbital_1s.primitiveGaussians[1].normalizationFactor *
			Atom1::orbital_1s.primitiveGaussians[1].coefficient *
			Atom2::orbital_1s.primitiveGaussians[0].normalizationFactor *
			Atom2::orbital_1s.primitiveGaussians[0].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_4, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_5 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[1].alpha + Atom2::orbital_1s.primitiveGaussians[1].alpha);
		constexpr double d5 =
			Atom1::orbital_1s.primitiveGaussians[1].normalizationFactor *
			Atom1::orbital_1s.primitiveGaussians[1].coefficient *
			Atom2::orbital_1s.primitiveGaussians[1].normalizationFactor *
			Atom2::orbital_1s.primitiveGaussians[1].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_5, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_6 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[1].alpha + Atom2::orbital_1s.primitiveGaussians[2].alpha);
		constexpr double d6 =
			Atom1::orbital_1s.primitiveGaussians[1].normalizationFactor *
			Atom1::orbital_1s.primitiveGaussians[1].coefficient *
			Atom2::orbital_1s.primitiveGaussians[2].normalizationFactor *
			Atom2::orbital_1s.primitiveGaussians[2].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_6, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_7 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[2].alpha + Atom2::orbital_1s.primitiveGaussians[0].alpha);
		constexpr double d7 =
			Atom1::orbital_1s.primitiveGaussians[2].normalizationFactor *
			Atom1::orbital_1s.primitiveGaussians[2].coefficient *
			Atom2::orbital_1s.primitiveGaussians[0].normalizationFactor *
			Atom2::orbital_1s.primitiveGaussians[0].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_7, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_8 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[2].alpha + Atom2::orbital_1s.primitiveGaussians[1].alpha);
		constexpr double d8 =
			Atom1::orbital_1s.primitiveGaussians[2].normalizationFactor *
			Atom1::orbital_1s.primitiveGaussians[2].coefficient *
			Atom2::orbital_1s.primitiveGaussians[1].normalizationFactor *
			Atom2::orbital_1s.primitiveGaussians[1].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_8, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_9 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[2].alpha + Atom2::orbital_1s.primitiveGaussians[2].alpha);
		constexpr double d9 =
			Atom1::orbital_1s.primitiveGaussians[2].normalizationFactor *
			Atom1::orbital_1s.primitiveGaussians[2].coefficient *
			Atom2::orbital_1s.primitiveGaussians[2].normalizationFactor *
			Atom2::orbital_1s.primitiveGaussians[2].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_9, 1.5);

		return { d1, d2, d3, d4, d5, d6, d7, d8, d9 };
	}
	template<typename Atom1, typename Atom2>
	consteval std::array<double, 9> GenerateOverlapExponentialFactors_1s1s() noexcept
	{
		constexpr double oneDividedByAlpha1PlusAlpha2_1 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[0].alpha + Atom2::orbital_1s.primitiveGaussians[0].alpha);
		constexpr double f1 = -1 * oneDividedByAlpha1PlusAlpha2_1 * Atom1::orbital_1s.primitiveGaussians[0].alpha * Atom2::orbital_1s.primitiveGaussians[0].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_2 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[0].alpha + Atom2::orbital_1s.primitiveGaussians[1].alpha);
		constexpr double f2 = -1 * oneDividedByAlpha1PlusAlpha2_2 * Atom1::orbital_1s.primitiveGaussians[0].alpha * Atom2::orbital_1s.primitiveGaussians[1].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_3 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[0].alpha + Atom2::orbital_1s.primitiveGaussians[2].alpha);
		constexpr double f3 = -1 * oneDividedByAlpha1PlusAlpha2_3 * Atom1::orbital_1s.primitiveGaussians[0].alpha * Atom2::orbital_1s.primitiveGaussians[2].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_4 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[1].alpha + Atom2::orbital_1s.primitiveGaussians[0].alpha);
		constexpr double f4 = -1 * oneDividedByAlpha1PlusAlpha2_4 * Atom1::orbital_1s.primitiveGaussians[1].alpha * Atom2::orbital_1s.primitiveGaussians[0].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_5 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[1].alpha + Atom2::orbital_1s.primitiveGaussians[1].alpha);
		constexpr double f5 = -1 * oneDividedByAlpha1PlusAlpha2_5 * Atom1::orbital_1s.primitiveGaussians[1].alpha * Atom2::orbital_1s.primitiveGaussians[1].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_6 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[1].alpha + Atom2::orbital_1s.primitiveGaussians[2].alpha);
		constexpr double f6 = -1 * oneDividedByAlpha1PlusAlpha2_6 * Atom1::orbital_1s.primitiveGaussians[1].alpha * Atom2::orbital_1s.primitiveGaussians[2].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_7 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[2].alpha + Atom2::orbital_1s.primitiveGaussians[0].alpha);
		constexpr double f7 = -1 * oneDividedByAlpha1PlusAlpha2_7 * Atom1::orbital_1s.primitiveGaussians[2].alpha * Atom2::orbital_1s.primitiveGaussians[0].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_8 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[2].alpha + Atom2::orbital_1s.primitiveGaussians[1].alpha);
		constexpr double f8 = -1 * oneDividedByAlpha1PlusAlpha2_8 * Atom1::orbital_1s.primitiveGaussians[2].alpha * Atom2::orbital_1s.primitiveGaussians[1].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_9 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[2].alpha + Atom2::orbital_1s.primitiveGaussians[2].alpha);
		constexpr double f9 = -1 * oneDividedByAlpha1PlusAlpha2_9 * Atom1::orbital_1s.primitiveGaussians[2].alpha * Atom2::orbital_1s.primitiveGaussians[2].alpha;

		return { f1, f2, f3, f4, f5, f6, f7, f8, f9 };
	}

	// 1s - 2s
	template<typename Atom1, typename Atom2>
	consteval std::array<double, 9> GenerateOverlapFactors_1s2s() noexcept
	{
		constexpr double oneDividedByAlpha1PlusAlpha2_1 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[0].alpha + Atom2::orbital_2s.primitiveGaussians[0].alpha);
		constexpr double d1 =
			Atom1::orbital_1s.primitiveGaussians[0].normalizationFactor *
			Atom1::orbital_1s.primitiveGaussians[0].coefficient *
			Atom2::orbital_2s.primitiveGaussians[0].normalizationFactor *
			Atom2::orbital_2s.primitiveGaussians[0].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_1, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_2 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[0].alpha + Atom2::orbital_2s.primitiveGaussians[1].alpha);
		constexpr double d2 =
			Atom1::orbital_1s.primitiveGaussians[0].normalizationFactor *
			Atom1::orbital_1s.primitiveGaussians[0].coefficient *
			Atom2::orbital_2s.primitiveGaussians[1].normalizationFactor *
			Atom2::orbital_2s.primitiveGaussians[1].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_2, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_3 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[0].alpha + Atom2::orbital_2s.primitiveGaussians[2].alpha);
		constexpr double d3 =
			Atom1::orbital_1s.primitiveGaussians[0].normalizationFactor *
			Atom1::orbital_1s.primitiveGaussians[0].coefficient *
			Atom2::orbital_2s.primitiveGaussians[2].normalizationFactor *
			Atom2::orbital_2s.primitiveGaussians[2].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_3, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_4 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[1].alpha + Atom2::orbital_2s.primitiveGaussians[0].alpha);
		constexpr double d4 =
			Atom1::orbital_1s.primitiveGaussians[1].normalizationFactor *
			Atom1::orbital_1s.primitiveGaussians[1].coefficient *
			Atom2::orbital_2s.primitiveGaussians[0].normalizationFactor *
			Atom2::orbital_2s.primitiveGaussians[0].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_4, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_5 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[1].alpha + Atom2::orbital_2s.primitiveGaussians[1].alpha);
		constexpr double d5 =
			Atom1::orbital_1s.primitiveGaussians[1].normalizationFactor *
			Atom1::orbital_1s.primitiveGaussians[1].coefficient *
			Atom2::orbital_2s.primitiveGaussians[1].normalizationFactor *
			Atom2::orbital_2s.primitiveGaussians[1].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_5, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_6 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[1].alpha + Atom2::orbital_2s.primitiveGaussians[2].alpha);
		constexpr double d6 =
			Atom1::orbital_1s.primitiveGaussians[1].normalizationFactor *
			Atom1::orbital_1s.primitiveGaussians[1].coefficient *
			Atom2::orbital_2s.primitiveGaussians[2].normalizationFactor *
			Atom2::orbital_2s.primitiveGaussians[2].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_6, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_7 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[2].alpha + Atom2::orbital_2s.primitiveGaussians[0].alpha);
		constexpr double d7 =
			Atom1::orbital_1s.primitiveGaussians[2].normalizationFactor *
			Atom1::orbital_1s.primitiveGaussians[2].coefficient *
			Atom2::orbital_2s.primitiveGaussians[0].normalizationFactor *
			Atom2::orbital_2s.primitiveGaussians[0].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_7, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_8 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[2].alpha + Atom2::orbital_2s.primitiveGaussians[1].alpha);
		constexpr double d8 =
			Atom1::orbital_1s.primitiveGaussians[2].normalizationFactor *
			Atom1::orbital_1s.primitiveGaussians[2].coefficient *
			Atom2::orbital_2s.primitiveGaussians[1].normalizationFactor *
			Atom2::orbital_2s.primitiveGaussians[1].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_8, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_9 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[2].alpha + Atom2::orbital_2s.primitiveGaussians[2].alpha);
		constexpr double d9 =
			Atom1::orbital_1s.primitiveGaussians[2].normalizationFactor *
			Atom1::orbital_1s.primitiveGaussians[2].coefficient *
			Atom2::orbital_2s.primitiveGaussians[2].normalizationFactor *
			Atom2::orbital_2s.primitiveGaussians[2].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_9, 1.5);

		return { d1, d2, d3, d4, d5, d6, d7, d8, d9 };
	}
	template<typename Atom1, typename Atom2>
	consteval std::array<double, 9> GenerateOverlapExponentialFactors_1s2s() noexcept
	{
		constexpr double oneDividedByAlpha1PlusAlpha2_1 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[0].alpha + Atom2::orbital_2s.primitiveGaussians[0].alpha);
		constexpr double f1 = -1 * oneDividedByAlpha1PlusAlpha2_1 * Atom1::orbital_1s.primitiveGaussians[0].alpha * Atom2::orbital_2s.primitiveGaussians[0].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_2 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[0].alpha + Atom2::orbital_2s.primitiveGaussians[1].alpha);
		constexpr double f2 = -1 * oneDividedByAlpha1PlusAlpha2_2 * Atom1::orbital_1s.primitiveGaussians[0].alpha * Atom2::orbital_2s.primitiveGaussians[1].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_3 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[0].alpha + Atom2::orbital_2s.primitiveGaussians[2].alpha);
		constexpr double f3 = -1 * oneDividedByAlpha1PlusAlpha2_3 * Atom1::orbital_1s.primitiveGaussians[0].alpha * Atom2::orbital_2s.primitiveGaussians[2].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_4 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[1].alpha + Atom2::orbital_2s.primitiveGaussians[0].alpha);
		constexpr double f4 = -1 * oneDividedByAlpha1PlusAlpha2_4 * Atom1::orbital_1s.primitiveGaussians[1].alpha * Atom2::orbital_2s.primitiveGaussians[0].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_5 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[1].alpha + Atom2::orbital_2s.primitiveGaussians[1].alpha);
		constexpr double f5 = -1 * oneDividedByAlpha1PlusAlpha2_5 * Atom1::orbital_1s.primitiveGaussians[1].alpha * Atom2::orbital_2s.primitiveGaussians[1].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_6 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[1].alpha + Atom2::orbital_2s.primitiveGaussians[2].alpha);
		constexpr double f6 = -1 * oneDividedByAlpha1PlusAlpha2_6 * Atom1::orbital_1s.primitiveGaussians[1].alpha * Atom2::orbital_2s.primitiveGaussians[2].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_7 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[2].alpha + Atom2::orbital_2s.primitiveGaussians[0].alpha);
		constexpr double f7 = -1 * oneDividedByAlpha1PlusAlpha2_7 * Atom1::orbital_1s.primitiveGaussians[2].alpha * Atom2::orbital_2s.primitiveGaussians[0].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_8 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[2].alpha + Atom2::orbital_2s.primitiveGaussians[1].alpha);
		constexpr double f8 = -1 * oneDividedByAlpha1PlusAlpha2_8 * Atom1::orbital_1s.primitiveGaussians[2].alpha * Atom2::orbital_2s.primitiveGaussians[1].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_9 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[2].alpha + Atom2::orbital_2s.primitiveGaussians[2].alpha);
		constexpr double f9 = -1 * oneDividedByAlpha1PlusAlpha2_9 * Atom1::orbital_1s.primitiveGaussians[2].alpha * Atom2::orbital_2s.primitiveGaussians[2].alpha;

		return { f1, f2, f3, f4, f5, f6, f7, f8, f9 };
	}
	
	// 1s - 2px
	template<typename Atom1, typename Atom2>
	consteval std::array<double, 9> GenerateOverlapFactors_1s2px() noexcept
	{
		constexpr double oneDividedByAlpha1PlusAlpha2_1 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[0].alpha + Atom2::orbital_2px.primitiveGaussians[0].alpha);
		constexpr double d1 =
			Atom1::orbital_1s.primitiveGaussians[0].normalizationFactor *
			Atom1::orbital_1s.primitiveGaussians[0].coefficient *
			Atom2::orbital_2px.primitiveGaussians[0].normalizationFactor *
			Atom2::orbital_2px.primitiveGaussians[0].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_1, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_2 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[0].alpha + Atom2::orbital_2px.primitiveGaussians[1].alpha);
		constexpr double d2 =
			Atom1::orbital_1s.primitiveGaussians[0].normalizationFactor *
			Atom1::orbital_1s.primitiveGaussians[0].coefficient *
			Atom2::orbital_2px.primitiveGaussians[1].normalizationFactor *
			Atom2::orbital_2px.primitiveGaussians[1].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_2, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_3 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[0].alpha + Atom2::orbital_2px.primitiveGaussians[2].alpha);
		constexpr double d3 =
			Atom1::orbital_1s.primitiveGaussians[0].normalizationFactor *
			Atom1::orbital_1s.primitiveGaussians[0].coefficient *
			Atom2::orbital_2px.primitiveGaussians[2].normalizationFactor *
			Atom2::orbital_2px.primitiveGaussians[2].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_3, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_4 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[1].alpha + Atom2::orbital_2px.primitiveGaussians[0].alpha);
		constexpr double d4 =
			Atom1::orbital_1s.primitiveGaussians[1].normalizationFactor *
			Atom1::orbital_1s.primitiveGaussians[1].coefficient *
			Atom2::orbital_2px.primitiveGaussians[0].normalizationFactor *
			Atom2::orbital_2px.primitiveGaussians[0].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_4, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_5 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[1].alpha + Atom2::orbital_2px.primitiveGaussians[1].alpha);
		constexpr double d5 =
			Atom1::orbital_1s.primitiveGaussians[1].normalizationFactor *
			Atom1::orbital_1s.primitiveGaussians[1].coefficient *
			Atom2::orbital_2px.primitiveGaussians[1].normalizationFactor *
			Atom2::orbital_2px.primitiveGaussians[1].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_5, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_6 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[1].alpha + Atom2::orbital_2px.primitiveGaussians[2].alpha);
		constexpr double d6 =
			Atom1::orbital_1s.primitiveGaussians[1].normalizationFactor *
			Atom1::orbital_1s.primitiveGaussians[1].coefficient *
			Atom2::orbital_2px.primitiveGaussians[2].normalizationFactor *
			Atom2::orbital_2px.primitiveGaussians[2].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_6, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_7 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[2].alpha + Atom2::orbital_2px.primitiveGaussians[0].alpha);
		constexpr double d7 =
			Atom1::orbital_1s.primitiveGaussians[2].normalizationFactor *
			Atom1::orbital_1s.primitiveGaussians[2].coefficient *
			Atom2::orbital_2px.primitiveGaussians[0].normalizationFactor *
			Atom2::orbital_2px.primitiveGaussians[0].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_7, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_8 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[2].alpha + Atom2::orbital_2px.primitiveGaussians[1].alpha);
		constexpr double d8 =
			Atom1::orbital_1s.primitiveGaussians[2].normalizationFactor *
			Atom1::orbital_1s.primitiveGaussians[2].coefficient *
			Atom2::orbital_2px.primitiveGaussians[1].normalizationFactor *
			Atom2::orbital_2px.primitiveGaussians[1].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_8, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_9 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[2].alpha + Atom2::orbital_2px.primitiveGaussians[2].alpha);
		constexpr double d9 =
			Atom1::orbital_1s.primitiveGaussians[2].normalizationFactor *
			Atom1::orbital_1s.primitiveGaussians[2].coefficient *
			Atom2::orbital_2px.primitiveGaussians[2].normalizationFactor *
			Atom2::orbital_2px.primitiveGaussians[2].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_9, 1.5);

		return { d1, d2, d3, d4, d5, d6, d7, d8, d9 };
	}
	template<typename Atom1, typename Atom2>
	consteval std::array<double, 9> GenerateOverlapExponentialFactors_1s2px() noexcept
	{
		constexpr double oneDividedByAlpha1PlusAlpha2_1 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[0].alpha + Atom2::orbital_2px.primitiveGaussians[0].alpha);
		constexpr double f1 = -1 * oneDividedByAlpha1PlusAlpha2_1 * Atom1::orbital_1s.primitiveGaussians[0].alpha * Atom2::orbital_2px.primitiveGaussians[0].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_2 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[0].alpha + Atom2::orbital_2px.primitiveGaussians[1].alpha);
		constexpr double f2 = -1 * oneDividedByAlpha1PlusAlpha2_2 * Atom1::orbital_1s.primitiveGaussians[0].alpha * Atom2::orbital_2px.primitiveGaussians[1].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_3 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[0].alpha + Atom2::orbital_2px.primitiveGaussians[2].alpha);
		constexpr double f3 = -1 * oneDividedByAlpha1PlusAlpha2_3 * Atom1::orbital_1s.primitiveGaussians[0].alpha * Atom2::orbital_2px.primitiveGaussians[2].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_4 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[1].alpha + Atom2::orbital_2px.primitiveGaussians[0].alpha);
		constexpr double f4 = -1 * oneDividedByAlpha1PlusAlpha2_4 * Atom1::orbital_1s.primitiveGaussians[1].alpha * Atom2::orbital_2px.primitiveGaussians[0].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_5 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[1].alpha + Atom2::orbital_2px.primitiveGaussians[1].alpha);
		constexpr double f5 = -1 * oneDividedByAlpha1PlusAlpha2_5 * Atom1::orbital_1s.primitiveGaussians[1].alpha * Atom2::orbital_2px.primitiveGaussians[1].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_6 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[1].alpha + Atom2::orbital_2px.primitiveGaussians[2].alpha);
		constexpr double f6 = -1 * oneDividedByAlpha1PlusAlpha2_6 * Atom1::orbital_1s.primitiveGaussians[1].alpha * Atom2::orbital_2px.primitiveGaussians[2].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_7 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[2].alpha + Atom2::orbital_2px.primitiveGaussians[0].alpha);
		constexpr double f7 = -1 * oneDividedByAlpha1PlusAlpha2_7 * Atom1::orbital_1s.primitiveGaussians[2].alpha * Atom2::orbital_2px.primitiveGaussians[0].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_8 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[2].alpha + Atom2::orbital_2px.primitiveGaussians[1].alpha);
		constexpr double f8 = -1 * oneDividedByAlpha1PlusAlpha2_8 * Atom1::orbital_1s.primitiveGaussians[2].alpha * Atom2::orbital_2px.primitiveGaussians[1].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_9 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[2].alpha + Atom2::orbital_2px.primitiveGaussians[2].alpha);
		constexpr double f9 = -1 * oneDividedByAlpha1PlusAlpha2_9 * Atom1::orbital_1s.primitiveGaussians[2].alpha * Atom2::orbital_2px.primitiveGaussians[2].alpha;

		return { f1, f2, f3, f4, f5, f6, f7, f8, f9 };
	}
	template<typename Atom1, typename Atom2>
	consteval std::array<double, 9> GenerateOneDividedByAlpha1PlusAlpha2_1s2px() noexcept
	{
		constexpr double d1 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[0].alpha + Atom2::orbital_2px.primitiveGaussians[0].alpha);
		constexpr double d2 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[0].alpha + Atom2::orbital_2px.primitiveGaussians[1].alpha);
		constexpr double d3 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[0].alpha + Atom2::orbital_2px.primitiveGaussians[2].alpha);
		constexpr double d4 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[1].alpha + Atom2::orbital_2px.primitiveGaussians[0].alpha);
		constexpr double d5 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[1].alpha + Atom2::orbital_2px.primitiveGaussians[1].alpha);
		constexpr double d6 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[1].alpha + Atom2::orbital_2px.primitiveGaussians[2].alpha);
		constexpr double d7 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[2].alpha + Atom2::orbital_2px.primitiveGaussians[0].alpha);
		constexpr double d8 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[2].alpha + Atom2::orbital_2px.primitiveGaussians[1].alpha);
		constexpr double d9 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[2].alpha + Atom2::orbital_2px.primitiveGaussians[2].alpha);

		return { d1, d2, d3, d4, d5, d6, d7, d8, d9 };
	}

	// 1s - 2py
	template<typename Atom1, typename Atom2>
	consteval std::array<double, 9> GenerateOverlapFactors_1s2py() noexcept
	{
		constexpr double oneDividedByAlpha1PlusAlpha2_1 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[0].alpha + Atom2::orbital_2py.primitiveGaussians[0].alpha);
		constexpr double d1 =
			Atom1::orbital_1s.primitiveGaussians[0].normalizationFactor *
			Atom1::orbital_1s.primitiveGaussians[0].coefficient *
			Atom2::orbital_2py.primitiveGaussians[0].normalizationFactor *
			Atom2::orbital_2py.primitiveGaussians[0].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_1, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_2 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[0].alpha + Atom2::orbital_2py.primitiveGaussians[1].alpha);
		constexpr double d2 =
			Atom1::orbital_1s.primitiveGaussians[0].normalizationFactor *
			Atom1::orbital_1s.primitiveGaussians[0].coefficient *
			Atom2::orbital_2py.primitiveGaussians[1].normalizationFactor *
			Atom2::orbital_2py.primitiveGaussians[1].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_2, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_3 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[0].alpha + Atom2::orbital_2py.primitiveGaussians[2].alpha);
		constexpr double d3 =
			Atom1::orbital_1s.primitiveGaussians[0].normalizationFactor *
			Atom1::orbital_1s.primitiveGaussians[0].coefficient *
			Atom2::orbital_2py.primitiveGaussians[2].normalizationFactor *
			Atom2::orbital_2py.primitiveGaussians[2].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_3, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_4 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[1].alpha + Atom2::orbital_2py.primitiveGaussians[0].alpha);
		constexpr double d4 =
			Atom1::orbital_1s.primitiveGaussians[1].normalizationFactor *
			Atom1::orbital_1s.primitiveGaussians[1].coefficient *
			Atom2::orbital_2py.primitiveGaussians[0].normalizationFactor *
			Atom2::orbital_2py.primitiveGaussians[0].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_4, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_5 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[1].alpha + Atom2::orbital_2py.primitiveGaussians[1].alpha);
		constexpr double d5 =
			Atom1::orbital_1s.primitiveGaussians[1].normalizationFactor *
			Atom1::orbital_1s.primitiveGaussians[1].coefficient *
			Atom2::orbital_2py.primitiveGaussians[1].normalizationFactor *
			Atom2::orbital_2py.primitiveGaussians[1].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_5, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_6 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[1].alpha + Atom2::orbital_2py.primitiveGaussians[2].alpha);
		constexpr double d6 =
			Atom1::orbital_1s.primitiveGaussians[1].normalizationFactor *
			Atom1::orbital_1s.primitiveGaussians[1].coefficient *
			Atom2::orbital_2py.primitiveGaussians[2].normalizationFactor *
			Atom2::orbital_2py.primitiveGaussians[2].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_6, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_7 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[2].alpha + Atom2::orbital_2py.primitiveGaussians[0].alpha);
		constexpr double d7 =
			Atom1::orbital_1s.primitiveGaussians[2].normalizationFactor *
			Atom1::orbital_1s.primitiveGaussians[2].coefficient *
			Atom2::orbital_2py.primitiveGaussians[0].normalizationFactor *
			Atom2::orbital_2py.primitiveGaussians[0].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_7, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_8 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[2].alpha + Atom2::orbital_2py.primitiveGaussians[1].alpha);
		constexpr double d8 =
			Atom1::orbital_1s.primitiveGaussians[2].normalizationFactor *
			Atom1::orbital_1s.primitiveGaussians[2].coefficient *
			Atom2::orbital_2py.primitiveGaussians[1].normalizationFactor *
			Atom2::orbital_2py.primitiveGaussians[1].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_8, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_9 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[2].alpha + Atom2::orbital_2py.primitiveGaussians[2].alpha);
		constexpr double d9 =
			Atom1::orbital_1s.primitiveGaussians[2].normalizationFactor *
			Atom1::orbital_1s.primitiveGaussians[2].coefficient *
			Atom2::orbital_2py.primitiveGaussians[2].normalizationFactor *
			Atom2::orbital_2py.primitiveGaussians[2].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_9, 1.5);

		return { d1, d2, d3, d4, d5, d6, d7, d8, d9 };
	}
	template<typename Atom1, typename Atom2>
	consteval std::array<double, 9> GenerateOverlapExponentialFactors_1s2py() noexcept
	{
		constexpr double oneDividedByAlpha1PlusAlpha2_1 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[0].alpha + Atom2::orbital_2py.primitiveGaussians[0].alpha);
		constexpr double f1 = -1 * oneDividedByAlpha1PlusAlpha2_1 * Atom1::orbital_1s.primitiveGaussians[0].alpha * Atom2::orbital_2py.primitiveGaussians[0].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_2 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[0].alpha + Atom2::orbital_2py.primitiveGaussians[1].alpha);
		constexpr double f2 = -1 * oneDividedByAlpha1PlusAlpha2_2 * Atom1::orbital_1s.primitiveGaussians[0].alpha * Atom2::orbital_2py.primitiveGaussians[1].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_3 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[0].alpha + Atom2::orbital_2py.primitiveGaussians[2].alpha);
		constexpr double f3 = -1 * oneDividedByAlpha1PlusAlpha2_3 * Atom1::orbital_1s.primitiveGaussians[0].alpha * Atom2::orbital_2py.primitiveGaussians[2].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_4 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[1].alpha + Atom2::orbital_2py.primitiveGaussians[0].alpha);
		constexpr double f4 = -1 * oneDividedByAlpha1PlusAlpha2_4 * Atom1::orbital_1s.primitiveGaussians[1].alpha * Atom2::orbital_2py.primitiveGaussians[0].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_5 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[1].alpha + Atom2::orbital_2py.primitiveGaussians[1].alpha);
		constexpr double f5 = -1 * oneDividedByAlpha1PlusAlpha2_5 * Atom1::orbital_1s.primitiveGaussians[1].alpha * Atom2::orbital_2py.primitiveGaussians[1].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_6 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[1].alpha + Atom2::orbital_2py.primitiveGaussians[2].alpha);
		constexpr double f6 = -1 * oneDividedByAlpha1PlusAlpha2_6 * Atom1::orbital_1s.primitiveGaussians[1].alpha * Atom2::orbital_2py.primitiveGaussians[2].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_7 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[2].alpha + Atom2::orbital_2py.primitiveGaussians[0].alpha);
		constexpr double f7 = -1 * oneDividedByAlpha1PlusAlpha2_7 * Atom1::orbital_1s.primitiveGaussians[2].alpha * Atom2::orbital_2py.primitiveGaussians[0].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_8 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[2].alpha + Atom2::orbital_2py.primitiveGaussians[1].alpha);
		constexpr double f8 = -1 * oneDividedByAlpha1PlusAlpha2_8 * Atom1::orbital_1s.primitiveGaussians[2].alpha * Atom2::orbital_2py.primitiveGaussians[1].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_9 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[2].alpha + Atom2::orbital_2py.primitiveGaussians[2].alpha);
		constexpr double f9 = -1 * oneDividedByAlpha1PlusAlpha2_9 * Atom1::orbital_1s.primitiveGaussians[2].alpha * Atom2::orbital_2py.primitiveGaussians[2].alpha;

		return { f1, f2, f3, f4, f5, f6, f7, f8, f9 };
	}
	template<typename Atom1, typename Atom2>
	consteval std::array<double, 9> GenerateOneDividedByAlpha1PlusAlpha2_1s2py() noexcept
	{
		constexpr double d1 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[0].alpha + Atom2::orbital_2py.primitiveGaussians[0].alpha);
		constexpr double d2 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[0].alpha + Atom2::orbital_2py.primitiveGaussians[1].alpha);
		constexpr double d3 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[0].alpha + Atom2::orbital_2py.primitiveGaussians[2].alpha);
		constexpr double d4 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[1].alpha + Atom2::orbital_2py.primitiveGaussians[0].alpha);
		constexpr double d5 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[1].alpha + Atom2::orbital_2py.primitiveGaussians[1].alpha);
		constexpr double d6 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[1].alpha + Atom2::orbital_2py.primitiveGaussians[2].alpha);
		constexpr double d7 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[2].alpha + Atom2::orbital_2py.primitiveGaussians[0].alpha);
		constexpr double d8 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[2].alpha + Atom2::orbital_2py.primitiveGaussians[1].alpha);
		constexpr double d9 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[2].alpha + Atom2::orbital_2py.primitiveGaussians[2].alpha);

		return { d1, d2, d3, d4, d5, d6, d7, d8, d9 };
	}

	// 1s - 2pz
	template<typename Atom1, typename Atom2>
	consteval std::array<double, 9> GenerateOverlapFactors_1s2pz() noexcept
	{
		constexpr double oneDividedByAlpha1PlusAlpha2_1 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[0].alpha + Atom2::orbital_2pz.primitiveGaussians[0].alpha);
		constexpr double d1 =
			Atom1::orbital_1s.primitiveGaussians[0].normalizationFactor *
			Atom1::orbital_1s.primitiveGaussians[0].coefficient *
			Atom2::orbital_2pz.primitiveGaussians[0].normalizationFactor *
			Atom2::orbital_2pz.primitiveGaussians[0].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_1, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_2 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[0].alpha + Atom2::orbital_2pz.primitiveGaussians[1].alpha);
		constexpr double d2 =
			Atom1::orbital_1s.primitiveGaussians[0].normalizationFactor *
			Atom1::orbital_1s.primitiveGaussians[0].coefficient *
			Atom2::orbital_2pz.primitiveGaussians[1].normalizationFactor *
			Atom2::orbital_2pz.primitiveGaussians[1].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_2, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_3 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[0].alpha + Atom2::orbital_2pz.primitiveGaussians[2].alpha);
		constexpr double d3 =
			Atom1::orbital_1s.primitiveGaussians[0].normalizationFactor *
			Atom1::orbital_1s.primitiveGaussians[0].coefficient *
			Atom2::orbital_2pz.primitiveGaussians[2].normalizationFactor *
			Atom2::orbital_2pz.primitiveGaussians[2].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_3, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_4 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[1].alpha + Atom2::orbital_2pz.primitiveGaussians[0].alpha);
		constexpr double d4 =
			Atom1::orbital_1s.primitiveGaussians[1].normalizationFactor *
			Atom1::orbital_1s.primitiveGaussians[1].coefficient *
			Atom2::orbital_2pz.primitiveGaussians[0].normalizationFactor *
			Atom2::orbital_2pz.primitiveGaussians[0].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_4, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_5 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[1].alpha + Atom2::orbital_2pz.primitiveGaussians[1].alpha);
		constexpr double d5 =
			Atom1::orbital_1s.primitiveGaussians[1].normalizationFactor *
			Atom1::orbital_1s.primitiveGaussians[1].coefficient *
			Atom2::orbital_2pz.primitiveGaussians[1].normalizationFactor *
			Atom2::orbital_2pz.primitiveGaussians[1].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_5, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_6 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[1].alpha + Atom2::orbital_2pz.primitiveGaussians[2].alpha);
		constexpr double d6 =
			Atom1::orbital_1s.primitiveGaussians[1].normalizationFactor *
			Atom1::orbital_1s.primitiveGaussians[1].coefficient *
			Atom2::orbital_2pz.primitiveGaussians[2].normalizationFactor *
			Atom2::orbital_2pz.primitiveGaussians[2].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_6, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_7 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[2].alpha + Atom2::orbital_2pz.primitiveGaussians[0].alpha);
		constexpr double d7 =
			Atom1::orbital_1s.primitiveGaussians[2].normalizationFactor *
			Atom1::orbital_1s.primitiveGaussians[2].coefficient *
			Atom2::orbital_2pz.primitiveGaussians[0].normalizationFactor *
			Atom2::orbital_2pz.primitiveGaussians[0].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_7, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_8 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[2].alpha + Atom2::orbital_2pz.primitiveGaussians[1].alpha);
		constexpr double d8 =
			Atom1::orbital_1s.primitiveGaussians[2].normalizationFactor *
			Atom1::orbital_1s.primitiveGaussians[2].coefficient *
			Atom2::orbital_2pz.primitiveGaussians[1].normalizationFactor *
			Atom2::orbital_2pz.primitiveGaussians[1].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_8, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_9 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[2].alpha + Atom2::orbital_2pz.primitiveGaussians[2].alpha);
		constexpr double d9 =
			Atom1::orbital_1s.primitiveGaussians[2].normalizationFactor *
			Atom1::orbital_1s.primitiveGaussians[2].coefficient *
			Atom2::orbital_2pz.primitiveGaussians[2].normalizationFactor *
			Atom2::orbital_2pz.primitiveGaussians[2].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_9, 1.5);

		return { d1, d2, d3, d4, d5, d6, d7, d8, d9 };
	}
	template<typename Atom1, typename Atom2>
	consteval std::array<double, 9> GenerateOverlapExponentialFactors_1s2pz() noexcept
	{
		constexpr double oneDividedByAlpha1PlusAlpha2_1 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[0].alpha + Atom2::orbital_2pz.primitiveGaussians[0].alpha);
		constexpr double f1 = -1 * oneDividedByAlpha1PlusAlpha2_1 * Atom1::orbital_1s.primitiveGaussians[0].alpha * Atom2::orbital_2pz.primitiveGaussians[0].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_2 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[0].alpha + Atom2::orbital_2pz.primitiveGaussians[1].alpha);
		constexpr double f2 = -1 * oneDividedByAlpha1PlusAlpha2_2 * Atom1::orbital_1s.primitiveGaussians[0].alpha * Atom2::orbital_2pz.primitiveGaussians[1].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_3 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[0].alpha + Atom2::orbital_2pz.primitiveGaussians[2].alpha);
		constexpr double f3 = -1 * oneDividedByAlpha1PlusAlpha2_3 * Atom1::orbital_1s.primitiveGaussians[0].alpha * Atom2::orbital_2pz.primitiveGaussians[2].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_4 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[1].alpha + Atom2::orbital_2pz.primitiveGaussians[0].alpha);
		constexpr double f4 = -1 * oneDividedByAlpha1PlusAlpha2_4 * Atom1::orbital_1s.primitiveGaussians[1].alpha * Atom2::orbital_2pz.primitiveGaussians[0].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_5 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[1].alpha + Atom2::orbital_2pz.primitiveGaussians[1].alpha);
		constexpr double f5 = -1 * oneDividedByAlpha1PlusAlpha2_5 * Atom1::orbital_1s.primitiveGaussians[1].alpha * Atom2::orbital_2pz.primitiveGaussians[1].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_6 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[1].alpha + Atom2::orbital_2pz.primitiveGaussians[2].alpha);
		constexpr double f6 = -1 * oneDividedByAlpha1PlusAlpha2_6 * Atom1::orbital_1s.primitiveGaussians[1].alpha * Atom2::orbital_2pz.primitiveGaussians[2].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_7 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[2].alpha + Atom2::orbital_2pz.primitiveGaussians[0].alpha);
		constexpr double f7 = -1 * oneDividedByAlpha1PlusAlpha2_7 * Atom1::orbital_1s.primitiveGaussians[2].alpha * Atom2::orbital_2pz.primitiveGaussians[0].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_8 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[2].alpha + Atom2::orbital_2pz.primitiveGaussians[1].alpha);
		constexpr double f8 = -1 * oneDividedByAlpha1PlusAlpha2_8 * Atom1::orbital_1s.primitiveGaussians[2].alpha * Atom2::orbital_2pz.primitiveGaussians[1].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_9 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[2].alpha + Atom2::orbital_2pz.primitiveGaussians[2].alpha);
		constexpr double f9 = -1 * oneDividedByAlpha1PlusAlpha2_9 * Atom1::orbital_1s.primitiveGaussians[2].alpha * Atom2::orbital_2pz.primitiveGaussians[2].alpha;

		return { f1, f2, f3, f4, f5, f6, f7, f8, f9 };
	}
	template<typename Atom1, typename Atom2>
	consteval std::array<double, 9> GenerateOneDividedByAlpha1PlusAlpha2_1s2pz() noexcept
	{
		constexpr double d1 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[0].alpha + Atom2::orbital_2pz.primitiveGaussians[0].alpha);
		constexpr double d2 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[0].alpha + Atom2::orbital_2pz.primitiveGaussians[1].alpha);
		constexpr double d3 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[0].alpha + Atom2::orbital_2pz.primitiveGaussians[2].alpha);
		constexpr double d4 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[1].alpha + Atom2::orbital_2pz.primitiveGaussians[0].alpha);
		constexpr double d5 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[1].alpha + Atom2::orbital_2pz.primitiveGaussians[1].alpha);
		constexpr double d6 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[1].alpha + Atom2::orbital_2pz.primitiveGaussians[2].alpha);
		constexpr double d7 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[2].alpha + Atom2::orbital_2pz.primitiveGaussians[0].alpha);
		constexpr double d8 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[2].alpha + Atom2::orbital_2pz.primitiveGaussians[1].alpha);
		constexpr double d9 = 1.0 / (Atom1::orbital_1s.primitiveGaussians[2].alpha + Atom2::orbital_2pz.primitiveGaussians[2].alpha);

		return { d1, d2, d3, d4, d5, d6, d7, d8, d9 };
	}

	// 2s - 2s
	template<typename Atom1, typename Atom2>
	consteval std::array<double, 9> GenerateOverlapFactors_2s2s() noexcept
	{
		constexpr double oneDividedByAlpha1PlusAlpha2_1 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[0].alpha + Atom2::orbital_2s.primitiveGaussians[0].alpha);
		constexpr double d1 =
			Atom1::orbital_2s.primitiveGaussians[0].normalizationFactor *
			Atom1::orbital_2s.primitiveGaussians[0].coefficient *
			Atom2::orbital_2s.primitiveGaussians[0].normalizationFactor *
			Atom2::orbital_2s.primitiveGaussians[0].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_1, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_2 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[0].alpha + Atom2::orbital_2s.primitiveGaussians[1].alpha);
		constexpr double d2 =
			Atom1::orbital_2s.primitiveGaussians[0].normalizationFactor *
			Atom1::orbital_2s.primitiveGaussians[0].coefficient *
			Atom2::orbital_2s.primitiveGaussians[1].normalizationFactor *
			Atom2::orbital_2s.primitiveGaussians[1].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_2, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_3 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[0].alpha + Atom2::orbital_2s.primitiveGaussians[2].alpha);
		constexpr double d3 =
			Atom1::orbital_2s.primitiveGaussians[0].normalizationFactor *
			Atom1::orbital_2s.primitiveGaussians[0].coefficient *
			Atom2::orbital_2s.primitiveGaussians[2].normalizationFactor *
			Atom2::orbital_2s.primitiveGaussians[2].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_3, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_4 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[1].alpha + Atom2::orbital_2s.primitiveGaussians[0].alpha);
		constexpr double d4 =
			Atom1::orbital_2s.primitiveGaussians[1].normalizationFactor *
			Atom1::orbital_2s.primitiveGaussians[1].coefficient *
			Atom2::orbital_2s.primitiveGaussians[0].normalizationFactor *
			Atom2::orbital_2s.primitiveGaussians[0].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_4, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_5 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[1].alpha + Atom2::orbital_2s.primitiveGaussians[1].alpha);
		constexpr double d5 =
			Atom1::orbital_2s.primitiveGaussians[1].normalizationFactor *
			Atom1::orbital_2s.primitiveGaussians[1].coefficient *
			Atom2::orbital_2s.primitiveGaussians[1].normalizationFactor *
			Atom2::orbital_2s.primitiveGaussians[1].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_5, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_6 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[1].alpha + Atom2::orbital_2s.primitiveGaussians[2].alpha);
		constexpr double d6 =
			Atom1::orbital_2s.primitiveGaussians[1].normalizationFactor *
			Atom1::orbital_2s.primitiveGaussians[1].coefficient *
			Atom2::orbital_2s.primitiveGaussians[2].normalizationFactor *
			Atom2::orbital_2s.primitiveGaussians[2].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_6, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_7 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[2].alpha + Atom2::orbital_2s.primitiveGaussians[0].alpha);
		constexpr double d7 =
			Atom1::orbital_2s.primitiveGaussians[2].normalizationFactor *
			Atom1::orbital_2s.primitiveGaussians[2].coefficient *
			Atom2::orbital_2s.primitiveGaussians[0].normalizationFactor *
			Atom2::orbital_2s.primitiveGaussians[0].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_7, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_8 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[2].alpha + Atom2::orbital_2s.primitiveGaussians[1].alpha);
		constexpr double d8 =
			Atom1::orbital_2s.primitiveGaussians[2].normalizationFactor *
			Atom1::orbital_2s.primitiveGaussians[2].coefficient *
			Atom2::orbital_2s.primitiveGaussians[1].normalizationFactor *
			Atom2::orbital_2s.primitiveGaussians[1].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_8, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_9 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[2].alpha + Atom2::orbital_2s.primitiveGaussians[2].alpha);
		constexpr double d9 =
			Atom1::orbital_2s.primitiveGaussians[2].normalizationFactor *
			Atom1::orbital_2s.primitiveGaussians[2].coefficient *
			Atom2::orbital_2s.primitiveGaussians[2].normalizationFactor *
			Atom2::orbital_2s.primitiveGaussians[2].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_9, 1.5);

		return { d1, d2, d3, d4, d5, d6, d7, d8, d9 };
	}
	template<typename Atom1, typename Atom2>
	consteval std::array<double, 9> GenerateOverlapExponentialFactors_2s2s() noexcept
	{
		constexpr double oneDividedByAlpha1PlusAlpha2_1 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[0].alpha + Atom2::orbital_2s.primitiveGaussians[0].alpha);
		constexpr double f1 = -1 * oneDividedByAlpha1PlusAlpha2_1 * Atom1::orbital_2s.primitiveGaussians[0].alpha * Atom2::orbital_2s.primitiveGaussians[0].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_2 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[0].alpha + Atom2::orbital_2s.primitiveGaussians[1].alpha);
		constexpr double f2 = -1 * oneDividedByAlpha1PlusAlpha2_2 * Atom1::orbital_2s.primitiveGaussians[0].alpha * Atom2::orbital_2s.primitiveGaussians[1].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_3 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[0].alpha + Atom2::orbital_2s.primitiveGaussians[2].alpha);
		constexpr double f3 = -1 * oneDividedByAlpha1PlusAlpha2_3 * Atom1::orbital_2s.primitiveGaussians[0].alpha * Atom2::orbital_2s.primitiveGaussians[2].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_4 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[1].alpha + Atom2::orbital_2s.primitiveGaussians[0].alpha);
		constexpr double f4 = -1 * oneDividedByAlpha1PlusAlpha2_4 * Atom1::orbital_2s.primitiveGaussians[1].alpha * Atom2::orbital_2s.primitiveGaussians[0].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_5 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[1].alpha + Atom2::orbital_2s.primitiveGaussians[1].alpha);
		constexpr double f5 = -1 * oneDividedByAlpha1PlusAlpha2_5 * Atom1::orbital_2s.primitiveGaussians[1].alpha * Atom2::orbital_2s.primitiveGaussians[1].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_6 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[1].alpha + Atom2::orbital_2s.primitiveGaussians[2].alpha);
		constexpr double f6 = -1 * oneDividedByAlpha1PlusAlpha2_6 * Atom1::orbital_2s.primitiveGaussians[1].alpha * Atom2::orbital_2s.primitiveGaussians[2].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_7 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[2].alpha + Atom2::orbital_2s.primitiveGaussians[0].alpha);
		constexpr double f7 = -1 * oneDividedByAlpha1PlusAlpha2_7 * Atom1::orbital_2s.primitiveGaussians[2].alpha * Atom2::orbital_2s.primitiveGaussians[0].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_8 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[2].alpha + Atom2::orbital_2s.primitiveGaussians[1].alpha);
		constexpr double f8 = -1 * oneDividedByAlpha1PlusAlpha2_8 * Atom1::orbital_2s.primitiveGaussians[2].alpha * Atom2::orbital_2s.primitiveGaussians[1].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_9 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[2].alpha + Atom2::orbital_2s.primitiveGaussians[2].alpha);
		constexpr double f9 = -1 * oneDividedByAlpha1PlusAlpha2_9 * Atom1::orbital_2s.primitiveGaussians[2].alpha * Atom2::orbital_2s.primitiveGaussians[2].alpha;

		return { f1, f2, f3, f4, f5, f6, f7, f8, f9 };
	}

	// 2s - 2px
	template<typename Atom1, typename Atom2>
	consteval std::array<double, 9> GenerateOverlapFactors_2s2px() noexcept
	{
		constexpr double oneDividedByAlpha1PlusAlpha2_1 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[0].alpha + Atom2::orbital_2px.primitiveGaussians[0].alpha);
		constexpr double d1 =
			Atom1::orbital_2s.primitiveGaussians[0].normalizationFactor *
			Atom1::orbital_2s.primitiveGaussians[0].coefficient *
			Atom2::orbital_2px.primitiveGaussians[0].normalizationFactor *
			Atom2::orbital_2px.primitiveGaussians[0].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_1, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_2 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[0].alpha + Atom2::orbital_2px.primitiveGaussians[1].alpha);
		constexpr double d2 =
			Atom1::orbital_2s.primitiveGaussians[0].normalizationFactor *
			Atom1::orbital_2s.primitiveGaussians[0].coefficient *
			Atom2::orbital_2px.primitiveGaussians[1].normalizationFactor *
			Atom2::orbital_2px.primitiveGaussians[1].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_2, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_3 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[0].alpha + Atom2::orbital_2px.primitiveGaussians[2].alpha);
		constexpr double d3 =
			Atom1::orbital_2s.primitiveGaussians[0].normalizationFactor *
			Atom1::orbital_2s.primitiveGaussians[0].coefficient *
			Atom2::orbital_2px.primitiveGaussians[2].normalizationFactor *
			Atom2::orbital_2px.primitiveGaussians[2].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_3, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_4 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[1].alpha + Atom2::orbital_2px.primitiveGaussians[0].alpha);
		constexpr double d4 =
			Atom1::orbital_2s.primitiveGaussians[1].normalizationFactor *
			Atom1::orbital_2s.primitiveGaussians[1].coefficient *
			Atom2::orbital_2px.primitiveGaussians[0].normalizationFactor *
			Atom2::orbital_2px.primitiveGaussians[0].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_4, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_5 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[1].alpha + Atom2::orbital_2px.primitiveGaussians[1].alpha);
		constexpr double d5 =
			Atom1::orbital_2s.primitiveGaussians[1].normalizationFactor *
			Atom1::orbital_2s.primitiveGaussians[1].coefficient *
			Atom2::orbital_2px.primitiveGaussians[1].normalizationFactor *
			Atom2::orbital_2px.primitiveGaussians[1].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_5, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_6 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[1].alpha + Atom2::orbital_2px.primitiveGaussians[2].alpha);
		constexpr double d6 =
			Atom1::orbital_2s.primitiveGaussians[1].normalizationFactor *
			Atom1::orbital_2s.primitiveGaussians[1].coefficient *
			Atom2::orbital_2px.primitiveGaussians[2].normalizationFactor *
			Atom2::orbital_2px.primitiveGaussians[2].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_6, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_7 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[2].alpha + Atom2::orbital_2px.primitiveGaussians[0].alpha);
		constexpr double d7 =
			Atom1::orbital_2s.primitiveGaussians[2].normalizationFactor *
			Atom1::orbital_2s.primitiveGaussians[2].coefficient *
			Atom2::orbital_2px.primitiveGaussians[0].normalizationFactor *
			Atom2::orbital_2px.primitiveGaussians[0].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_7, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_8 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[2].alpha + Atom2::orbital_2px.primitiveGaussians[1].alpha);
		constexpr double d8 =
			Atom1::orbital_2s.primitiveGaussians[2].normalizationFactor *
			Atom1::orbital_2s.primitiveGaussians[2].coefficient *
			Atom2::orbital_2px.primitiveGaussians[1].normalizationFactor *
			Atom2::orbital_2px.primitiveGaussians[1].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_8, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_9 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[2].alpha + Atom2::orbital_2px.primitiveGaussians[2].alpha);
		constexpr double d9 =
			Atom1::orbital_2s.primitiveGaussians[2].normalizationFactor *
			Atom1::orbital_2s.primitiveGaussians[2].coefficient *
			Atom2::orbital_2px.primitiveGaussians[2].normalizationFactor *
			Atom2::orbital_2px.primitiveGaussians[2].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_9, 1.5);

		return { d1, d2, d3, d4, d5, d6, d7, d8, d9 };
	}
	template<typename Atom1, typename Atom2>
	consteval std::array<double, 9> GenerateOverlapExponentialFactors_2s2px() noexcept
	{
		constexpr double oneDividedByAlpha1PlusAlpha2_1 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[0].alpha + Atom2::orbital_2px.primitiveGaussians[0].alpha);
		constexpr double f1 = -1 * oneDividedByAlpha1PlusAlpha2_1 * Atom1::orbital_2s.primitiveGaussians[0].alpha * Atom2::orbital_2px.primitiveGaussians[0].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_2 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[0].alpha + Atom2::orbital_2px.primitiveGaussians[1].alpha);
		constexpr double f2 = -1 * oneDividedByAlpha1PlusAlpha2_2 * Atom1::orbital_2s.primitiveGaussians[0].alpha * Atom2::orbital_2px.primitiveGaussians[1].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_3 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[0].alpha + Atom2::orbital_2px.primitiveGaussians[2].alpha);
		constexpr double f3 = -1 * oneDividedByAlpha1PlusAlpha2_3 * Atom1::orbital_2s.primitiveGaussians[0].alpha * Atom2::orbital_2px.primitiveGaussians[2].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_4 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[1].alpha + Atom2::orbital_2px.primitiveGaussians[0].alpha);
		constexpr double f4 = -1 * oneDividedByAlpha1PlusAlpha2_4 * Atom1::orbital_2s.primitiveGaussians[1].alpha * Atom2::orbital_2px.primitiveGaussians[0].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_5 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[1].alpha + Atom2::orbital_2px.primitiveGaussians[1].alpha);
		constexpr double f5 = -1 * oneDividedByAlpha1PlusAlpha2_5 * Atom1::orbital_2s.primitiveGaussians[1].alpha * Atom2::orbital_2px.primitiveGaussians[1].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_6 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[1].alpha + Atom2::orbital_2px.primitiveGaussians[2].alpha);
		constexpr double f6 = -1 * oneDividedByAlpha1PlusAlpha2_6 * Atom1::orbital_2s.primitiveGaussians[1].alpha * Atom2::orbital_2px.primitiveGaussians[2].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_7 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[2].alpha + Atom2::orbital_2px.primitiveGaussians[0].alpha);
		constexpr double f7 = -1 * oneDividedByAlpha1PlusAlpha2_7 * Atom1::orbital_2s.primitiveGaussians[2].alpha * Atom2::orbital_2px.primitiveGaussians[0].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_8 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[2].alpha + Atom2::orbital_2px.primitiveGaussians[1].alpha);
		constexpr double f8 = -1 * oneDividedByAlpha1PlusAlpha2_8 * Atom1::orbital_2s.primitiveGaussians[2].alpha * Atom2::orbital_2px.primitiveGaussians[1].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_9 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[2].alpha + Atom2::orbital_2px.primitiveGaussians[2].alpha);
		constexpr double f9 = -1 * oneDividedByAlpha1PlusAlpha2_9 * Atom1::orbital_2s.primitiveGaussians[2].alpha * Atom2::orbital_2px.primitiveGaussians[2].alpha;

		return { f1, f2, f3, f4, f5, f6, f7, f8, f9 };
	}
	template<typename Atom1, typename Atom2>
	consteval std::array<double, 9> GenerateOneDividedByAlpha1PlusAlpha2_2s2px() noexcept
	{
		constexpr double d1 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[0].alpha + Atom2::orbital_2px.primitiveGaussians[0].alpha);
		constexpr double d2 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[0].alpha + Atom2::orbital_2px.primitiveGaussians[1].alpha);
		constexpr double d3 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[0].alpha + Atom2::orbital_2px.primitiveGaussians[2].alpha);
		constexpr double d4 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[1].alpha + Atom2::orbital_2px.primitiveGaussians[0].alpha);
		constexpr double d5 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[1].alpha + Atom2::orbital_2px.primitiveGaussians[1].alpha);
		constexpr double d6 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[1].alpha + Atom2::orbital_2px.primitiveGaussians[2].alpha);
		constexpr double d7 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[2].alpha + Atom2::orbital_2px.primitiveGaussians[0].alpha);
		constexpr double d8 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[2].alpha + Atom2::orbital_2px.primitiveGaussians[1].alpha);
		constexpr double d9 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[2].alpha + Atom2::orbital_2px.primitiveGaussians[2].alpha);

		return { d1, d2, d3, d4, d5, d6, d7, d8, d9 };
	}

	// 2s - 2py
	template<typename Atom1, typename Atom2>
	consteval std::array<double, 9> GenerateOverlapFactors_2s2py() noexcept
	{
		constexpr double oneDividedByAlpha1PlusAlpha2_1 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[0].alpha + Atom2::orbital_2py.primitiveGaussians[0].alpha);
		constexpr double d1 =
			Atom1::orbital_2s.primitiveGaussians[0].normalizationFactor *
			Atom1::orbital_2s.primitiveGaussians[0].coefficient *
			Atom2::orbital_2py.primitiveGaussians[0].normalizationFactor *
			Atom2::orbital_2py.primitiveGaussians[0].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_1, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_2 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[0].alpha + Atom2::orbital_2py.primitiveGaussians[1].alpha);
		constexpr double d2 =
			Atom1::orbital_2s.primitiveGaussians[0].normalizationFactor *
			Atom1::orbital_2s.primitiveGaussians[0].coefficient *
			Atom2::orbital_2py.primitiveGaussians[1].normalizationFactor *
			Atom2::orbital_2py.primitiveGaussians[1].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_2, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_3 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[0].alpha + Atom2::orbital_2py.primitiveGaussians[2].alpha);
		constexpr double d3 =
			Atom1::orbital_2s.primitiveGaussians[0].normalizationFactor *
			Atom1::orbital_2s.primitiveGaussians[0].coefficient *
			Atom2::orbital_2py.primitiveGaussians[2].normalizationFactor *
			Atom2::orbital_2py.primitiveGaussians[2].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_3, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_4 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[1].alpha + Atom2::orbital_2py.primitiveGaussians[0].alpha);
		constexpr double d4 =
			Atom1::orbital_2s.primitiveGaussians[1].normalizationFactor *
			Atom1::orbital_2s.primitiveGaussians[1].coefficient *
			Atom2::orbital_2py.primitiveGaussians[0].normalizationFactor *
			Atom2::orbital_2py.primitiveGaussians[0].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_4, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_5 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[1].alpha + Atom2::orbital_2py.primitiveGaussians[1].alpha);
		constexpr double d5 =
			Atom1::orbital_2s.primitiveGaussians[1].normalizationFactor *
			Atom1::orbital_2s.primitiveGaussians[1].coefficient *
			Atom2::orbital_2py.primitiveGaussians[1].normalizationFactor *
			Atom2::orbital_2py.primitiveGaussians[1].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_5, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_6 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[1].alpha + Atom2::orbital_2py.primitiveGaussians[2].alpha);
		constexpr double d6 =
			Atom1::orbital_2s.primitiveGaussians[1].normalizationFactor *
			Atom1::orbital_2s.primitiveGaussians[1].coefficient *
			Atom2::orbital_2py.primitiveGaussians[2].normalizationFactor *
			Atom2::orbital_2py.primitiveGaussians[2].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_6, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_7 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[2].alpha + Atom2::orbital_2py.primitiveGaussians[0].alpha);
		constexpr double d7 =
			Atom1::orbital_2s.primitiveGaussians[2].normalizationFactor *
			Atom1::orbital_2s.primitiveGaussians[2].coefficient *
			Atom2::orbital_2py.primitiveGaussians[0].normalizationFactor *
			Atom2::orbital_2py.primitiveGaussians[0].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_7, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_8 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[2].alpha + Atom2::orbital_2py.primitiveGaussians[1].alpha);
		constexpr double d8 =
			Atom1::orbital_2s.primitiveGaussians[2].normalizationFactor *
			Atom1::orbital_2s.primitiveGaussians[2].coefficient *
			Atom2::orbital_2py.primitiveGaussians[1].normalizationFactor *
			Atom2::orbital_2py.primitiveGaussians[1].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_8, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_9 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[2].alpha + Atom2::orbital_2py.primitiveGaussians[2].alpha);
		constexpr double d9 =
			Atom1::orbital_2s.primitiveGaussians[2].normalizationFactor *
			Atom1::orbital_2s.primitiveGaussians[2].coefficient *
			Atom2::orbital_2py.primitiveGaussians[2].normalizationFactor *
			Atom2::orbital_2py.primitiveGaussians[2].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_9, 1.5);

		return { d1, d2, d3, d4, d5, d6, d7, d8, d9 };
	}
	template<typename Atom1, typename Atom2>
	consteval std::array<double, 9> GenerateOverlapExponentialFactors_2s2py() noexcept
	{
		constexpr double oneDividedByAlpha1PlusAlpha2_1 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[0].alpha + Atom2::orbital_2py.primitiveGaussians[0].alpha);
		constexpr double f1 = -1 * oneDividedByAlpha1PlusAlpha2_1 * Atom1::orbital_2s.primitiveGaussians[0].alpha * Atom2::orbital_2py.primitiveGaussians[0].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_2 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[0].alpha + Atom2::orbital_2py.primitiveGaussians[1].alpha);
		constexpr double f2 = -1 * oneDividedByAlpha1PlusAlpha2_2 * Atom1::orbital_2s.primitiveGaussians[0].alpha * Atom2::orbital_2py.primitiveGaussians[1].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_3 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[0].alpha + Atom2::orbital_2py.primitiveGaussians[2].alpha);
		constexpr double f3 = -1 * oneDividedByAlpha1PlusAlpha2_3 * Atom1::orbital_2s.primitiveGaussians[0].alpha * Atom2::orbital_2py.primitiveGaussians[2].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_4 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[1].alpha + Atom2::orbital_2py.primitiveGaussians[0].alpha);
		constexpr double f4 = -1 * oneDividedByAlpha1PlusAlpha2_4 * Atom1::orbital_2s.primitiveGaussians[1].alpha * Atom2::orbital_2py.primitiveGaussians[0].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_5 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[1].alpha + Atom2::orbital_2py.primitiveGaussians[1].alpha);
		constexpr double f5 = -1 * oneDividedByAlpha1PlusAlpha2_5 * Atom1::orbital_2s.primitiveGaussians[1].alpha * Atom2::orbital_2py.primitiveGaussians[1].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_6 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[1].alpha + Atom2::orbital_2py.primitiveGaussians[2].alpha);
		constexpr double f6 = -1 * oneDividedByAlpha1PlusAlpha2_6 * Atom1::orbital_2s.primitiveGaussians[1].alpha * Atom2::orbital_2py.primitiveGaussians[2].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_7 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[2].alpha + Atom2::orbital_2py.primitiveGaussians[0].alpha);
		constexpr double f7 = -1 * oneDividedByAlpha1PlusAlpha2_7 * Atom1::orbital_2s.primitiveGaussians[2].alpha * Atom2::orbital_2py.primitiveGaussians[0].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_8 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[2].alpha + Atom2::orbital_2py.primitiveGaussians[1].alpha);
		constexpr double f8 = -1 * oneDividedByAlpha1PlusAlpha2_8 * Atom1::orbital_2s.primitiveGaussians[2].alpha * Atom2::orbital_2py.primitiveGaussians[1].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_9 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[2].alpha + Atom2::orbital_2py.primitiveGaussians[2].alpha);
		constexpr double f9 = -1 * oneDividedByAlpha1PlusAlpha2_9 * Atom1::orbital_2s.primitiveGaussians[2].alpha * Atom2::orbital_2py.primitiveGaussians[2].alpha;

		return { f1, f2, f3, f4, f5, f6, f7, f8, f9 };
	}
	template<typename Atom1, typename Atom2>
	consteval std::array<double, 9> GenerateOneDividedByAlpha1PlusAlpha2_2s2py() noexcept
	{
		constexpr double d1 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[0].alpha + Atom2::orbital_2py.primitiveGaussians[0].alpha);
		constexpr double d2 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[0].alpha + Atom2::orbital_2py.primitiveGaussians[1].alpha);
		constexpr double d3 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[0].alpha + Atom2::orbital_2py.primitiveGaussians[2].alpha);
		constexpr double d4 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[1].alpha + Atom2::orbital_2py.primitiveGaussians[0].alpha);
		constexpr double d5 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[1].alpha + Atom2::orbital_2py.primitiveGaussians[1].alpha);
		constexpr double d6 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[1].alpha + Atom2::orbital_2py.primitiveGaussians[2].alpha);
		constexpr double d7 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[2].alpha + Atom2::orbital_2py.primitiveGaussians[0].alpha);
		constexpr double d8 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[2].alpha + Atom2::orbital_2py.primitiveGaussians[1].alpha);
		constexpr double d9 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[2].alpha + Atom2::orbital_2py.primitiveGaussians[2].alpha);

		return { d1, d2, d3, d4, d5, d6, d7, d8, d9 };
	}

	// 2s - 2pz
	template<typename Atom1, typename Atom2>
	consteval std::array<double, 9> GenerateOverlapFactors_2s2pz() noexcept
	{
		constexpr double oneDividedByAlpha1PlusAlpha2_1 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[0].alpha + Atom2::orbital_2pz.primitiveGaussians[0].alpha);
		constexpr double d1 =
			Atom1::orbital_2s.primitiveGaussians[0].normalizationFactor *
			Atom1::orbital_2s.primitiveGaussians[0].coefficient *
			Atom2::orbital_2pz.primitiveGaussians[0].normalizationFactor *
			Atom2::orbital_2pz.primitiveGaussians[0].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_1, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_2 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[0].alpha + Atom2::orbital_2pz.primitiveGaussians[1].alpha);
		constexpr double d2 =
			Atom1::orbital_2s.primitiveGaussians[0].normalizationFactor *
			Atom1::orbital_2s.primitiveGaussians[0].coefficient *
			Atom2::orbital_2pz.primitiveGaussians[1].normalizationFactor *
			Atom2::orbital_2pz.primitiveGaussians[1].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_2, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_3 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[0].alpha + Atom2::orbital_2pz.primitiveGaussians[2].alpha);
		constexpr double d3 =
			Atom1::orbital_2s.primitiveGaussians[0].normalizationFactor *
			Atom1::orbital_2s.primitiveGaussians[0].coefficient *
			Atom2::orbital_2pz.primitiveGaussians[2].normalizationFactor *
			Atom2::orbital_2pz.primitiveGaussians[2].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_3, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_4 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[1].alpha + Atom2::orbital_2pz.primitiveGaussians[0].alpha);
		constexpr double d4 =
			Atom1::orbital_2s.primitiveGaussians[1].normalizationFactor *
			Atom1::orbital_2s.primitiveGaussians[1].coefficient *
			Atom2::orbital_2pz.primitiveGaussians[0].normalizationFactor *
			Atom2::orbital_2pz.primitiveGaussians[0].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_4, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_5 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[1].alpha + Atom2::orbital_2pz.primitiveGaussians[1].alpha);
		constexpr double d5 =
			Atom1::orbital_2s.primitiveGaussians[1].normalizationFactor *
			Atom1::orbital_2s.primitiveGaussians[1].coefficient *
			Atom2::orbital_2pz.primitiveGaussians[1].normalizationFactor *
			Atom2::orbital_2pz.primitiveGaussians[1].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_5, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_6 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[1].alpha + Atom2::orbital_2pz.primitiveGaussians[2].alpha);
		constexpr double d6 =
			Atom1::orbital_2s.primitiveGaussians[1].normalizationFactor *
			Atom1::orbital_2s.primitiveGaussians[1].coefficient *
			Atom2::orbital_2pz.primitiveGaussians[2].normalizationFactor *
			Atom2::orbital_2pz.primitiveGaussians[2].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_6, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_7 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[2].alpha + Atom2::orbital_2pz.primitiveGaussians[0].alpha);
		constexpr double d7 =
			Atom1::orbital_2s.primitiveGaussians[2].normalizationFactor *
			Atom1::orbital_2s.primitiveGaussians[2].coefficient *
			Atom2::orbital_2pz.primitiveGaussians[0].normalizationFactor *
			Atom2::orbital_2pz.primitiveGaussians[0].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_7, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_8 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[2].alpha + Atom2::orbital_2pz.primitiveGaussians[1].alpha);
		constexpr double d8 =
			Atom1::orbital_2s.primitiveGaussians[2].normalizationFactor *
			Atom1::orbital_2s.primitiveGaussians[2].coefficient *
			Atom2::orbital_2pz.primitiveGaussians[1].normalizationFactor *
			Atom2::orbital_2pz.primitiveGaussians[1].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_8, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_9 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[2].alpha + Atom2::orbital_2pz.primitiveGaussians[2].alpha);
		constexpr double d9 =
			Atom1::orbital_2s.primitiveGaussians[2].normalizationFactor *
			Atom1::orbital_2s.primitiveGaussians[2].coefficient *
			Atom2::orbital_2pz.primitiveGaussians[2].normalizationFactor *
			Atom2::orbital_2pz.primitiveGaussians[2].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_9, 1.5);

		return { d1, d2, d3, d4, d5, d6, d7, d8, d9 };
	}
	template<typename Atom1, typename Atom2>
	consteval std::array<double, 9> GenerateOverlapExponentialFactors_2s2pz() noexcept
	{
		constexpr double oneDividedByAlpha1PlusAlpha2_1 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[0].alpha + Atom2::orbital_2pz.primitiveGaussians[0].alpha);
		constexpr double f1 = -1 * oneDividedByAlpha1PlusAlpha2_1 * Atom1::orbital_2s.primitiveGaussians[0].alpha * Atom2::orbital_2pz.primitiveGaussians[0].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_2 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[0].alpha + Atom2::orbital_2pz.primitiveGaussians[1].alpha);
		constexpr double f2 = -1 * oneDividedByAlpha1PlusAlpha2_2 * Atom1::orbital_2s.primitiveGaussians[0].alpha * Atom2::orbital_2pz.primitiveGaussians[1].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_3 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[0].alpha + Atom2::orbital_2pz.primitiveGaussians[2].alpha);
		constexpr double f3 = -1 * oneDividedByAlpha1PlusAlpha2_3 * Atom1::orbital_2s.primitiveGaussians[0].alpha * Atom2::orbital_2pz.primitiveGaussians[2].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_4 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[1].alpha + Atom2::orbital_2pz.primitiveGaussians[0].alpha);
		constexpr double f4 = -1 * oneDividedByAlpha1PlusAlpha2_4 * Atom1::orbital_2s.primitiveGaussians[1].alpha * Atom2::orbital_2pz.primitiveGaussians[0].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_5 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[1].alpha + Atom2::orbital_2pz.primitiveGaussians[1].alpha);
		constexpr double f5 = -1 * oneDividedByAlpha1PlusAlpha2_5 * Atom1::orbital_2s.primitiveGaussians[1].alpha * Atom2::orbital_2pz.primitiveGaussians[1].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_6 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[1].alpha + Atom2::orbital_2pz.primitiveGaussians[2].alpha);
		constexpr double f6 = -1 * oneDividedByAlpha1PlusAlpha2_6 * Atom1::orbital_2s.primitiveGaussians[1].alpha * Atom2::orbital_2pz.primitiveGaussians[2].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_7 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[2].alpha + Atom2::orbital_2pz.primitiveGaussians[0].alpha);
		constexpr double f7 = -1 * oneDividedByAlpha1PlusAlpha2_7 * Atom1::orbital_2s.primitiveGaussians[2].alpha * Atom2::orbital_2pz.primitiveGaussians[0].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_8 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[2].alpha + Atom2::orbital_2pz.primitiveGaussians[1].alpha);
		constexpr double f8 = -1 * oneDividedByAlpha1PlusAlpha2_8 * Atom1::orbital_2s.primitiveGaussians[2].alpha * Atom2::orbital_2pz.primitiveGaussians[1].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_9 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[2].alpha + Atom2::orbital_2pz.primitiveGaussians[2].alpha);
		constexpr double f9 = -1 * oneDividedByAlpha1PlusAlpha2_9 * Atom1::orbital_2s.primitiveGaussians[2].alpha * Atom2::orbital_2pz.primitiveGaussians[2].alpha;

		return { f1, f2, f3, f4, f5, f6, f7, f8, f9 };
	}
	template<typename Atom1, typename Atom2>
	consteval std::array<double, 9> GenerateOneDividedByAlpha1PlusAlpha2_2s2pz() noexcept
	{
		constexpr double d1 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[0].alpha + Atom2::orbital_2pz.primitiveGaussians[0].alpha);
		constexpr double d2 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[0].alpha + Atom2::orbital_2pz.primitiveGaussians[1].alpha);
		constexpr double d3 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[0].alpha + Atom2::orbital_2pz.primitiveGaussians[2].alpha);
		constexpr double d4 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[1].alpha + Atom2::orbital_2pz.primitiveGaussians[0].alpha);
		constexpr double d5 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[1].alpha + Atom2::orbital_2pz.primitiveGaussians[1].alpha);
		constexpr double d6 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[1].alpha + Atom2::orbital_2pz.primitiveGaussians[2].alpha);
		constexpr double d7 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[2].alpha + Atom2::orbital_2pz.primitiveGaussians[0].alpha);
		constexpr double d8 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[2].alpha + Atom2::orbital_2pz.primitiveGaussians[1].alpha);
		constexpr double d9 = 1.0 / (Atom1::orbital_2s.primitiveGaussians[2].alpha + Atom2::orbital_2pz.primitiveGaussians[2].alpha);

		return { d1, d2, d3, d4, d5, d6, d7, d8, d9 };
	}

	// 2px - 2px
	template<typename Atom1, typename Atom2>
	consteval std::array<double, 9> GenerateOverlapFactors_2px2px() noexcept
	{
		constexpr double oneDividedByAlpha1PlusAlpha2_1 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[0].alpha + Atom2::orbital_2px.primitiveGaussians[0].alpha);
		constexpr double d1 =
			Atom1::orbital_2px.primitiveGaussians[0].normalizationFactor *
			Atom1::orbital_2px.primitiveGaussians[0].coefficient *
			Atom2::orbital_2px.primitiveGaussians[0].normalizationFactor *
			Atom2::orbital_2px.primitiveGaussians[0].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_1, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_2 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[0].alpha + Atom2::orbital_2px.primitiveGaussians[1].alpha);
		constexpr double d2 =
			Atom1::orbital_2px.primitiveGaussians[0].normalizationFactor *
			Atom1::orbital_2px.primitiveGaussians[0].coefficient *
			Atom2::orbital_2px.primitiveGaussians[1].normalizationFactor *
			Atom2::orbital_2px.primitiveGaussians[1].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_2, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_3 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[0].alpha + Atom2::orbital_2px.primitiveGaussians[2].alpha);
		constexpr double d3 =
			Atom1::orbital_2px.primitiveGaussians[0].normalizationFactor *
			Atom1::orbital_2px.primitiveGaussians[0].coefficient *
			Atom2::orbital_2px.primitiveGaussians[2].normalizationFactor *
			Atom2::orbital_2px.primitiveGaussians[2].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_3, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_4 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[1].alpha + Atom2::orbital_2px.primitiveGaussians[0].alpha);
		constexpr double d4 =
			Atom1::orbital_2px.primitiveGaussians[1].normalizationFactor *
			Atom1::orbital_2px.primitiveGaussians[1].coefficient *
			Atom2::orbital_2px.primitiveGaussians[0].normalizationFactor *
			Atom2::orbital_2px.primitiveGaussians[0].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_4, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_5 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[1].alpha + Atom2::orbital_2px.primitiveGaussians[1].alpha);
		constexpr double d5 =
			Atom1::orbital_2px.primitiveGaussians[1].normalizationFactor *
			Atom1::orbital_2px.primitiveGaussians[1].coefficient *
			Atom2::orbital_2px.primitiveGaussians[1].normalizationFactor *
			Atom2::orbital_2px.primitiveGaussians[1].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_5, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_6 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[1].alpha + Atom2::orbital_2px.primitiveGaussians[2].alpha);
		constexpr double d6 =
			Atom1::orbital_2px.primitiveGaussians[1].normalizationFactor *
			Atom1::orbital_2px.primitiveGaussians[1].coefficient *
			Atom2::orbital_2px.primitiveGaussians[2].normalizationFactor *
			Atom2::orbital_2px.primitiveGaussians[2].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_6, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_7 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[2].alpha + Atom2::orbital_2px.primitiveGaussians[0].alpha);
		constexpr double d7 =
			Atom1::orbital_2px.primitiveGaussians[2].normalizationFactor *
			Atom1::orbital_2px.primitiveGaussians[2].coefficient *
			Atom2::orbital_2px.primitiveGaussians[0].normalizationFactor *
			Atom2::orbital_2px.primitiveGaussians[0].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_7, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_8 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[2].alpha + Atom2::orbital_2px.primitiveGaussians[1].alpha);
		constexpr double d8 =
			Atom1::orbital_2px.primitiveGaussians[2].normalizationFactor *
			Atom1::orbital_2px.primitiveGaussians[2].coefficient *
			Atom2::orbital_2px.primitiveGaussians[1].normalizationFactor *
			Atom2::orbital_2px.primitiveGaussians[1].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_8, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_9 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[2].alpha + Atom2::orbital_2px.primitiveGaussians[2].alpha);
		constexpr double d9 =
			Atom1::orbital_2px.primitiveGaussians[2].normalizationFactor *
			Atom1::orbital_2px.primitiveGaussians[2].coefficient *
			Atom2::orbital_2px.primitiveGaussians[2].normalizationFactor *
			Atom2::orbital_2px.primitiveGaussians[2].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_9, 1.5);

		return { d1, d2, d3, d4, d5, d6, d7, d8, d9 };
	}
	template<typename Atom1, typename Atom2>
	consteval std::array<double, 9> GenerateOverlapExponentialFactors_2px2px() noexcept
	{
		constexpr double oneDividedByAlpha1PlusAlpha2_1 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[0].alpha + Atom2::orbital_2px.primitiveGaussians[0].alpha);
		constexpr double f1 = -1 * oneDividedByAlpha1PlusAlpha2_1 * Atom1::orbital_2px.primitiveGaussians[0].alpha * Atom2::orbital_2px.primitiveGaussians[0].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_2 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[0].alpha + Atom2::orbital_2px.primitiveGaussians[1].alpha);
		constexpr double f2 = -1 * oneDividedByAlpha1PlusAlpha2_2 * Atom1::orbital_2px.primitiveGaussians[0].alpha * Atom2::orbital_2px.primitiveGaussians[1].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_3 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[0].alpha + Atom2::orbital_2px.primitiveGaussians[2].alpha);
		constexpr double f3 = -1 * oneDividedByAlpha1PlusAlpha2_3 * Atom1::orbital_2px.primitiveGaussians[0].alpha * Atom2::orbital_2px.primitiveGaussians[2].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_4 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[1].alpha + Atom2::orbital_2px.primitiveGaussians[0].alpha);
		constexpr double f4 = -1 * oneDividedByAlpha1PlusAlpha2_4 * Atom1::orbital_2px.primitiveGaussians[1].alpha * Atom2::orbital_2px.primitiveGaussians[0].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_5 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[1].alpha + Atom2::orbital_2px.primitiveGaussians[1].alpha);
		constexpr double f5 = -1 * oneDividedByAlpha1PlusAlpha2_5 * Atom1::orbital_2px.primitiveGaussians[1].alpha * Atom2::orbital_2px.primitiveGaussians[1].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_6 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[1].alpha + Atom2::orbital_2px.primitiveGaussians[2].alpha);
		constexpr double f6 = -1 * oneDividedByAlpha1PlusAlpha2_6 * Atom1::orbital_2px.primitiveGaussians[1].alpha * Atom2::orbital_2px.primitiveGaussians[2].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_7 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[2].alpha + Atom2::orbital_2px.primitiveGaussians[0].alpha);
		constexpr double f7 = -1 * oneDividedByAlpha1PlusAlpha2_7 * Atom1::orbital_2px.primitiveGaussians[2].alpha * Atom2::orbital_2px.primitiveGaussians[0].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_8 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[2].alpha + Atom2::orbital_2px.primitiveGaussians[1].alpha);
		constexpr double f8 = -1 * oneDividedByAlpha1PlusAlpha2_8 * Atom1::orbital_2px.primitiveGaussians[2].alpha * Atom2::orbital_2px.primitiveGaussians[1].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_9 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[2].alpha + Atom2::orbital_2px.primitiveGaussians[2].alpha);
		constexpr double f9 = -1 * oneDividedByAlpha1PlusAlpha2_9 * Atom1::orbital_2px.primitiveGaussians[2].alpha * Atom2::orbital_2px.primitiveGaussians[2].alpha;

		return { f1, f2, f3, f4, f5, f6, f7, f8, f9 };
	}
	template<typename Atom1, typename Atom2>
	consteval std::array<double, 9> GenerateOneDividedByAlpha1PlusAlpha2_2px2px() noexcept
	{
		constexpr double d1 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[0].alpha + Atom2::orbital_2px.primitiveGaussians[0].alpha);
		constexpr double d2 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[0].alpha + Atom2::orbital_2px.primitiveGaussians[1].alpha);
		constexpr double d3 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[0].alpha + Atom2::orbital_2px.primitiveGaussians[2].alpha);
		constexpr double d4 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[1].alpha + Atom2::orbital_2px.primitiveGaussians[0].alpha);
		constexpr double d5 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[1].alpha + Atom2::orbital_2px.primitiveGaussians[1].alpha);
		constexpr double d6 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[1].alpha + Atom2::orbital_2px.primitiveGaussians[2].alpha);
		constexpr double d7 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[2].alpha + Atom2::orbital_2px.primitiveGaussians[0].alpha);
		constexpr double d8 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[2].alpha + Atom2::orbital_2px.primitiveGaussians[1].alpha);
		constexpr double d9 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[2].alpha + Atom2::orbital_2px.primitiveGaussians[2].alpha);

		return { d1, d2, d3, d4, d5, d6, d7, d8, d9 };
	}

	// 2py - 2py
	template<typename Atom1, typename Atom2>
	consteval std::array<double, 9> GenerateOverlapFactors_2py2py() noexcept
	{
		constexpr double oneDividedByAlpha1PlusAlpha2_1 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[0].alpha + Atom2::orbital_2py.primitiveGaussians[0].alpha);
		constexpr double d1 =
			Atom1::orbital_2py.primitiveGaussians[0].normalizationFactor *
			Atom1::orbital_2py.primitiveGaussians[0].coefficient *
			Atom2::orbital_2py.primitiveGaussians[0].normalizationFactor *
			Atom2::orbital_2py.primitiveGaussians[0].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_1, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_2 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[0].alpha + Atom2::orbital_2py.primitiveGaussians[1].alpha);
		constexpr double d2 =
			Atom1::orbital_2py.primitiveGaussians[0].normalizationFactor *
			Atom1::orbital_2py.primitiveGaussians[0].coefficient *
			Atom2::orbital_2py.primitiveGaussians[1].normalizationFactor *
			Atom2::orbital_2py.primitiveGaussians[1].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_2, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_3 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[0].alpha + Atom2::orbital_2py.primitiveGaussians[2].alpha);
		constexpr double d3 =
			Atom1::orbital_2py.primitiveGaussians[0].normalizationFactor *
			Atom1::orbital_2py.primitiveGaussians[0].coefficient *
			Atom2::orbital_2py.primitiveGaussians[2].normalizationFactor *
			Atom2::orbital_2py.primitiveGaussians[2].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_3, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_4 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[1].alpha + Atom2::orbital_2py.primitiveGaussians[0].alpha);
		constexpr double d4 =
			Atom1::orbital_2py.primitiveGaussians[1].normalizationFactor *
			Atom1::orbital_2py.primitiveGaussians[1].coefficient *
			Atom2::orbital_2py.primitiveGaussians[0].normalizationFactor *
			Atom2::orbital_2py.primitiveGaussians[0].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_4, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_5 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[1].alpha + Atom2::orbital_2py.primitiveGaussians[1].alpha);
		constexpr double d5 =
			Atom1::orbital_2py.primitiveGaussians[1].normalizationFactor *
			Atom1::orbital_2py.primitiveGaussians[1].coefficient *
			Atom2::orbital_2py.primitiveGaussians[1].normalizationFactor *
			Atom2::orbital_2py.primitiveGaussians[1].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_5, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_6 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[1].alpha + Atom2::orbital_2py.primitiveGaussians[2].alpha);
		constexpr double d6 =
			Atom1::orbital_2py.primitiveGaussians[1].normalizationFactor *
			Atom1::orbital_2py.primitiveGaussians[1].coefficient *
			Atom2::orbital_2py.primitiveGaussians[2].normalizationFactor *
			Atom2::orbital_2py.primitiveGaussians[2].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_6, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_7 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[2].alpha + Atom2::orbital_2py.primitiveGaussians[0].alpha);
		constexpr double d7 =
			Atom1::orbital_2py.primitiveGaussians[2].normalizationFactor *
			Atom1::orbital_2py.primitiveGaussians[2].coefficient *
			Atom2::orbital_2py.primitiveGaussians[0].normalizationFactor *
			Atom2::orbital_2py.primitiveGaussians[0].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_7, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_8 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[2].alpha + Atom2::orbital_2py.primitiveGaussians[1].alpha);
		constexpr double d8 =
			Atom1::orbital_2py.primitiveGaussians[2].normalizationFactor *
			Atom1::orbital_2py.primitiveGaussians[2].coefficient *
			Atom2::orbital_2py.primitiveGaussians[1].normalizationFactor *
			Atom2::orbital_2py.primitiveGaussians[1].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_8, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_9 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[2].alpha + Atom2::orbital_2py.primitiveGaussians[2].alpha);
		constexpr double d9 =
			Atom1::orbital_2py.primitiveGaussians[2].normalizationFactor *
			Atom1::orbital_2py.primitiveGaussians[2].coefficient *
			Atom2::orbital_2py.primitiveGaussians[2].normalizationFactor *
			Atom2::orbital_2py.primitiveGaussians[2].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_9, 1.5);

		return { d1, d2, d3, d4, d5, d6, d7, d8, d9 };
	}
	template<typename Atom1, typename Atom2>
	consteval std::array<double, 9> GenerateOverlapExponentialFactors_2py2py() noexcept
	{
		constexpr double oneDividedByAlpha1PlusAlpha2_1 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[0].alpha + Atom2::orbital_2py.primitiveGaussians[0].alpha);
		constexpr double f1 = -1 * oneDividedByAlpha1PlusAlpha2_1 * Atom1::orbital_2py.primitiveGaussians[0].alpha * Atom2::orbital_2py.primitiveGaussians[0].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_2 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[0].alpha + Atom2::orbital_2py.primitiveGaussians[1].alpha);
		constexpr double f2 = -1 * oneDividedByAlpha1PlusAlpha2_2 * Atom1::orbital_2py.primitiveGaussians[0].alpha * Atom2::orbital_2py.primitiveGaussians[1].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_3 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[0].alpha + Atom2::orbital_2py.primitiveGaussians[2].alpha);
		constexpr double f3 = -1 * oneDividedByAlpha1PlusAlpha2_3 * Atom1::orbital_2py.primitiveGaussians[0].alpha * Atom2::orbital_2py.primitiveGaussians[2].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_4 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[1].alpha + Atom2::orbital_2py.primitiveGaussians[0].alpha);
		constexpr double f4 = -1 * oneDividedByAlpha1PlusAlpha2_4 * Atom1::orbital_2py.primitiveGaussians[1].alpha * Atom2::orbital_2py.primitiveGaussians[0].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_5 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[1].alpha + Atom2::orbital_2py.primitiveGaussians[1].alpha);
		constexpr double f5 = -1 * oneDividedByAlpha1PlusAlpha2_5 * Atom1::orbital_2py.primitiveGaussians[1].alpha * Atom2::orbital_2py.primitiveGaussians[1].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_6 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[1].alpha + Atom2::orbital_2py.primitiveGaussians[2].alpha);
		constexpr double f6 = -1 * oneDividedByAlpha1PlusAlpha2_6 * Atom1::orbital_2py.primitiveGaussians[1].alpha * Atom2::orbital_2py.primitiveGaussians[2].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_7 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[2].alpha + Atom2::orbital_2py.primitiveGaussians[0].alpha);
		constexpr double f7 = -1 * oneDividedByAlpha1PlusAlpha2_7 * Atom1::orbital_2py.primitiveGaussians[2].alpha * Atom2::orbital_2py.primitiveGaussians[0].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_8 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[2].alpha + Atom2::orbital_2py.primitiveGaussians[1].alpha);
		constexpr double f8 = -1 * oneDividedByAlpha1PlusAlpha2_8 * Atom1::orbital_2py.primitiveGaussians[2].alpha * Atom2::orbital_2py.primitiveGaussians[1].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_9 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[2].alpha + Atom2::orbital_2py.primitiveGaussians[2].alpha);
		constexpr double f9 = -1 * oneDividedByAlpha1PlusAlpha2_9 * Atom1::orbital_2py.primitiveGaussians[2].alpha * Atom2::orbital_2py.primitiveGaussians[2].alpha;

		return { f1, f2, f3, f4, f5, f6, f7, f8, f9 };
	}
	template<typename Atom1, typename Atom2>
	consteval std::array<double, 9> GenerateOneDividedByAlpha1PlusAlpha2_2py2py() noexcept
	{
		constexpr double d1 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[0].alpha + Atom2::orbital_2py.primitiveGaussians[0].alpha);
		constexpr double d2 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[0].alpha + Atom2::orbital_2py.primitiveGaussians[1].alpha);
		constexpr double d3 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[0].alpha + Atom2::orbital_2py.primitiveGaussians[2].alpha);
		constexpr double d4 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[1].alpha + Atom2::orbital_2py.primitiveGaussians[0].alpha);
		constexpr double d5 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[1].alpha + Atom2::orbital_2py.primitiveGaussians[1].alpha);
		constexpr double d6 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[1].alpha + Atom2::orbital_2py.primitiveGaussians[2].alpha);
		constexpr double d7 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[2].alpha + Atom2::orbital_2py.primitiveGaussians[0].alpha);
		constexpr double d8 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[2].alpha + Atom2::orbital_2py.primitiveGaussians[1].alpha);
		constexpr double d9 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[2].alpha + Atom2::orbital_2py.primitiveGaussians[2].alpha);

		return { d1, d2, d3, d4, d5, d6, d7, d8, d9 };
	}

	// 2pz - 2pz
	template<typename Atom1, typename Atom2>
	consteval std::array<double, 9> GenerateOverlapFactors_2pz2pz() noexcept
	{
		constexpr double oneDividedByAlpha1PlusAlpha2_1 = 1.0 / (Atom1::orbital_2pz.primitiveGaussians[0].alpha + Atom2::orbital_2pz.primitiveGaussians[0].alpha);
		constexpr double d1 =
			Atom1::orbital_2pz.primitiveGaussians[0].normalizationFactor *
			Atom1::orbital_2pz.primitiveGaussians[0].coefficient *
			Atom2::orbital_2pz.primitiveGaussians[0].normalizationFactor *
			Atom2::orbital_2pz.primitiveGaussians[0].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_1, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_2 = 1.0 / (Atom1::orbital_2pz.primitiveGaussians[0].alpha + Atom2::orbital_2pz.primitiveGaussians[1].alpha);
		constexpr double d2 =
			Atom1::orbital_2pz.primitiveGaussians[0].normalizationFactor *
			Atom1::orbital_2pz.primitiveGaussians[0].coefficient *
			Atom2::orbital_2pz.primitiveGaussians[1].normalizationFactor *
			Atom2::orbital_2pz.primitiveGaussians[1].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_2, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_3 = 1.0 / (Atom1::orbital_2pz.primitiveGaussians[0].alpha + Atom2::orbital_2pz.primitiveGaussians[2].alpha);
		constexpr double d3 =
			Atom1::orbital_2pz.primitiveGaussians[0].normalizationFactor *
			Atom1::orbital_2pz.primitiveGaussians[0].coefficient *
			Atom2::orbital_2pz.primitiveGaussians[2].normalizationFactor *
			Atom2::orbital_2pz.primitiveGaussians[2].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_3, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_4 = 1.0 / (Atom1::orbital_2pz.primitiveGaussians[1].alpha + Atom2::orbital_2pz.primitiveGaussians[0].alpha);
		constexpr double d4 =
			Atom1::orbital_2pz.primitiveGaussians[1].normalizationFactor *
			Atom1::orbital_2pz.primitiveGaussians[1].coefficient *
			Atom2::orbital_2pz.primitiveGaussians[0].normalizationFactor *
			Atom2::orbital_2pz.primitiveGaussians[0].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_4, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_5 = 1.0 / (Atom1::orbital_2pz.primitiveGaussians[1].alpha + Atom2::orbital_2pz.primitiveGaussians[1].alpha);
		constexpr double d5 =
			Atom1::orbital_2pz.primitiveGaussians[1].normalizationFactor *
			Atom1::orbital_2pz.primitiveGaussians[1].coefficient *
			Atom2::orbital_2pz.primitiveGaussians[1].normalizationFactor *
			Atom2::orbital_2pz.primitiveGaussians[1].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_5, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_6 = 1.0 / (Atom1::orbital_2pz.primitiveGaussians[1].alpha + Atom2::orbital_2pz.primitiveGaussians[2].alpha);
		constexpr double d6 =
			Atom1::orbital_2pz.primitiveGaussians[1].normalizationFactor *
			Atom1::orbital_2pz.primitiveGaussians[1].coefficient *
			Atom2::orbital_2pz.primitiveGaussians[2].normalizationFactor *
			Atom2::orbital_2pz.primitiveGaussians[2].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_6, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_7 = 1.0 / (Atom1::orbital_2pz.primitiveGaussians[2].alpha + Atom2::orbital_2pz.primitiveGaussians[0].alpha);
		constexpr double d7 =
			Atom1::orbital_2pz.primitiveGaussians[2].normalizationFactor *
			Atom1::orbital_2pz.primitiveGaussians[2].coefficient *
			Atom2::orbital_2pz.primitiveGaussians[0].normalizationFactor *
			Atom2::orbital_2pz.primitiveGaussians[0].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_7, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_8 = 1.0 / (Atom1::orbital_2pz.primitiveGaussians[2].alpha + Atom2::orbital_2pz.primitiveGaussians[1].alpha);
		constexpr double d8 =
			Atom1::orbital_2pz.primitiveGaussians[2].normalizationFactor *
			Atom1::orbital_2pz.primitiveGaussians[2].coefficient *
			Atom2::orbital_2pz.primitiveGaussians[1].normalizationFactor *
			Atom2::orbital_2pz.primitiveGaussians[1].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_8, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_9 = 1.0 / (Atom1::orbital_2pz.primitiveGaussians[2].alpha + Atom2::orbital_2pz.primitiveGaussians[2].alpha);
		constexpr double d9 =
			Atom1::orbital_2pz.primitiveGaussians[2].normalizationFactor *
			Atom1::orbital_2pz.primitiveGaussians[2].coefficient *
			Atom2::orbital_2pz.primitiveGaussians[2].normalizationFactor *
			Atom2::orbital_2pz.primitiveGaussians[2].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_9, 1.5);

		return { d1, d2, d3, d4, d5, d6, d7, d8, d9 };
	}
	template<typename Atom1, typename Atom2>
	consteval std::array<double, 9> GenerateOverlapExponentialFactors_2pz2pz() noexcept
	{
		constexpr double oneDividedByAlpha1PlusAlpha2_1 = 1.0 / (Atom1::orbital_2pz.primitiveGaussians[0].alpha + Atom2::orbital_2pz.primitiveGaussians[0].alpha);
		constexpr double f1 = -1 * oneDividedByAlpha1PlusAlpha2_1 * Atom1::orbital_2pz.primitiveGaussians[0].alpha * Atom2::orbital_2pz.primitiveGaussians[0].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_2 = 1.0 / (Atom1::orbital_2pz.primitiveGaussians[0].alpha + Atom2::orbital_2pz.primitiveGaussians[1].alpha);
		constexpr double f2 = -1 * oneDividedByAlpha1PlusAlpha2_2 * Atom1::orbital_2pz.primitiveGaussians[0].alpha * Atom2::orbital_2pz.primitiveGaussians[1].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_3 = 1.0 / (Atom1::orbital_2pz.primitiveGaussians[0].alpha + Atom2::orbital_2pz.primitiveGaussians[2].alpha);
		constexpr double f3 = -1 * oneDividedByAlpha1PlusAlpha2_3 * Atom1::orbital_2pz.primitiveGaussians[0].alpha * Atom2::orbital_2pz.primitiveGaussians[2].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_4 = 1.0 / (Atom1::orbital_2pz.primitiveGaussians[1].alpha + Atom2::orbital_2pz.primitiveGaussians[0].alpha);
		constexpr double f4 = -1 * oneDividedByAlpha1PlusAlpha2_4 * Atom1::orbital_2pz.primitiveGaussians[1].alpha * Atom2::orbital_2pz.primitiveGaussians[0].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_5 = 1.0 / (Atom1::orbital_2pz.primitiveGaussians[1].alpha + Atom2::orbital_2pz.primitiveGaussians[1].alpha);
		constexpr double f5 = -1 * oneDividedByAlpha1PlusAlpha2_5 * Atom1::orbital_2pz.primitiveGaussians[1].alpha * Atom2::orbital_2pz.primitiveGaussians[1].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_6 = 1.0 / (Atom1::orbital_2pz.primitiveGaussians[1].alpha + Atom2::orbital_2pz.primitiveGaussians[2].alpha);
		constexpr double f6 = -1 * oneDividedByAlpha1PlusAlpha2_6 * Atom1::orbital_2pz.primitiveGaussians[1].alpha * Atom2::orbital_2pz.primitiveGaussians[2].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_7 = 1.0 / (Atom1::orbital_2pz.primitiveGaussians[2].alpha + Atom2::orbital_2pz.primitiveGaussians[0].alpha);
		constexpr double f7 = -1 * oneDividedByAlpha1PlusAlpha2_7 * Atom1::orbital_2pz.primitiveGaussians[2].alpha * Atom2::orbital_2pz.primitiveGaussians[0].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_8 = 1.0 / (Atom1::orbital_2pz.primitiveGaussians[2].alpha + Atom2::orbital_2pz.primitiveGaussians[1].alpha);
		constexpr double f8 = -1 * oneDividedByAlpha1PlusAlpha2_8 * Atom1::orbital_2pz.primitiveGaussians[2].alpha * Atom2::orbital_2pz.primitiveGaussians[1].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_9 = 1.0 / (Atom1::orbital_2pz.primitiveGaussians[2].alpha + Atom2::orbital_2pz.primitiveGaussians[2].alpha);
		constexpr double f9 = -1 * oneDividedByAlpha1PlusAlpha2_9 * Atom1::orbital_2pz.primitiveGaussians[2].alpha * Atom2::orbital_2pz.primitiveGaussians[2].alpha;

		return { f1, f2, f3, f4, f5, f6, f7, f8, f9 };
	}
	template<typename Atom1, typename Atom2>
	consteval std::array<double, 9> GenerateOneDividedByAlpha1PlusAlpha2_2pz2pz() noexcept
	{
		constexpr double d1 = 1.0 / (Atom1::orbital_2pz.primitiveGaussians[0].alpha + Atom2::orbital_2pz.primitiveGaussians[0].alpha);
		constexpr double d2 = 1.0 / (Atom1::orbital_2pz.primitiveGaussians[0].alpha + Atom2::orbital_2pz.primitiveGaussians[1].alpha);
		constexpr double d3 = 1.0 / (Atom1::orbital_2pz.primitiveGaussians[0].alpha + Atom2::orbital_2pz.primitiveGaussians[2].alpha);
		constexpr double d4 = 1.0 / (Atom1::orbital_2pz.primitiveGaussians[1].alpha + Atom2::orbital_2pz.primitiveGaussians[0].alpha);
		constexpr double d5 = 1.0 / (Atom1::orbital_2pz.primitiveGaussians[1].alpha + Atom2::orbital_2pz.primitiveGaussians[1].alpha);
		constexpr double d6 = 1.0 / (Atom1::orbital_2pz.primitiveGaussians[1].alpha + Atom2::orbital_2pz.primitiveGaussians[2].alpha);
		constexpr double d7 = 1.0 / (Atom1::orbital_2pz.primitiveGaussians[2].alpha + Atom2::orbital_2pz.primitiveGaussians[0].alpha);
		constexpr double d8 = 1.0 / (Atom1::orbital_2pz.primitiveGaussians[2].alpha + Atom2::orbital_2pz.primitiveGaussians[1].alpha);
		constexpr double d9 = 1.0 / (Atom1::orbital_2pz.primitiveGaussians[2].alpha + Atom2::orbital_2pz.primitiveGaussians[2].alpha);

		return { d1, d2, d3, d4, d5, d6, d7, d8, d9 };
	}

	// 2px - 2py
	template<typename Atom1, typename Atom2>
	consteval std::array<double, 9> GenerateOverlapFactors_2px2py() noexcept
	{
		constexpr double oneDividedByAlpha1PlusAlpha2_1 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[0].alpha + Atom2::orbital_2py.primitiveGaussians[0].alpha);
		constexpr double d1 =
			Atom1::orbital_2px.primitiveGaussians[0].normalizationFactor *
			Atom1::orbital_2px.primitiveGaussians[0].coefficient *
			Atom2::orbital_2py.primitiveGaussians[0].normalizationFactor *
			Atom2::orbital_2py.primitiveGaussians[0].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_1, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_2 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[0].alpha + Atom2::orbital_2py.primitiveGaussians[1].alpha);
		constexpr double d2 =
			Atom1::orbital_2px.primitiveGaussians[0].normalizationFactor *
			Atom1::orbital_2px.primitiveGaussians[0].coefficient *
			Atom2::orbital_2py.primitiveGaussians[1].normalizationFactor *
			Atom2::orbital_2py.primitiveGaussians[1].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_2, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_3 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[0].alpha + Atom2::orbital_2py.primitiveGaussians[2].alpha);
		constexpr double d3 =
			Atom1::orbital_2px.primitiveGaussians[0].normalizationFactor *
			Atom1::orbital_2px.primitiveGaussians[0].coefficient *
			Atom2::orbital_2py.primitiveGaussians[2].normalizationFactor *
			Atom2::orbital_2py.primitiveGaussians[2].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_3, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_4 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[1].alpha + Atom2::orbital_2py.primitiveGaussians[0].alpha);
		constexpr double d4 =
			Atom1::orbital_2px.primitiveGaussians[1].normalizationFactor *
			Atom1::orbital_2px.primitiveGaussians[1].coefficient *
			Atom2::orbital_2py.primitiveGaussians[0].normalizationFactor *
			Atom2::orbital_2py.primitiveGaussians[0].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_4, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_5 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[1].alpha + Atom2::orbital_2py.primitiveGaussians[1].alpha);
		constexpr double d5 =
			Atom1::orbital_2px.primitiveGaussians[1].normalizationFactor *
			Atom1::orbital_2px.primitiveGaussians[1].coefficient *
			Atom2::orbital_2py.primitiveGaussians[1].normalizationFactor *
			Atom2::orbital_2py.primitiveGaussians[1].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_5, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_6 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[1].alpha + Atom2::orbital_2py.primitiveGaussians[2].alpha);
		constexpr double d6 =
			Atom1::orbital_2px.primitiveGaussians[1].normalizationFactor *
			Atom1::orbital_2px.primitiveGaussians[1].coefficient *
			Atom2::orbital_2py.primitiveGaussians[2].normalizationFactor *
			Atom2::orbital_2py.primitiveGaussians[2].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_6, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_7 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[2].alpha + Atom2::orbital_2py.primitiveGaussians[0].alpha);
		constexpr double d7 =
			Atom1::orbital_2px.primitiveGaussians[2].normalizationFactor *
			Atom1::orbital_2px.primitiveGaussians[2].coefficient *
			Atom2::orbital_2py.primitiveGaussians[0].normalizationFactor *
			Atom2::orbital_2py.primitiveGaussians[0].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_7, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_8 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[2].alpha + Atom2::orbital_2py.primitiveGaussians[1].alpha);
		constexpr double d8 =
			Atom1::orbital_2px.primitiveGaussians[2].normalizationFactor *
			Atom1::orbital_2px.primitiveGaussians[2].coefficient *
			Atom2::orbital_2py.primitiveGaussians[1].normalizationFactor *
			Atom2::orbital_2py.primitiveGaussians[1].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_8, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_9 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[2].alpha + Atom2::orbital_2py.primitiveGaussians[2].alpha);
		constexpr double d9 =
			Atom1::orbital_2px.primitiveGaussians[2].normalizationFactor *
			Atom1::orbital_2px.primitiveGaussians[2].coefficient *
			Atom2::orbital_2py.primitiveGaussians[2].normalizationFactor *
			Atom2::orbital_2py.primitiveGaussians[2].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_9, 1.5);

		return { d1, d2, d3, d4, d5, d6, d7, d8, d9 };
	}
	template<typename Atom1, typename Atom2>
	consteval std::array<double, 9> GenerateOverlapExponentialFactors_2px2py() noexcept
	{
		constexpr double oneDividedByAlpha1PlusAlpha2_1 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[0].alpha + Atom2::orbital_2py.primitiveGaussians[0].alpha);
		constexpr double f1 = -1 * oneDividedByAlpha1PlusAlpha2_1 * Atom1::orbital_2px.primitiveGaussians[0].alpha * Atom2::orbital_2py.primitiveGaussians[0].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_2 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[0].alpha + Atom2::orbital_2py.primitiveGaussians[1].alpha);
		constexpr double f2 = -1 * oneDividedByAlpha1PlusAlpha2_2 * Atom1::orbital_2px.primitiveGaussians[0].alpha * Atom2::orbital_2py.primitiveGaussians[1].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_3 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[0].alpha + Atom2::orbital_2py.primitiveGaussians[2].alpha);
		constexpr double f3 = -1 * oneDividedByAlpha1PlusAlpha2_3 * Atom1::orbital_2px.primitiveGaussians[0].alpha * Atom2::orbital_2py.primitiveGaussians[2].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_4 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[1].alpha + Atom2::orbital_2py.primitiveGaussians[0].alpha);
		constexpr double f4 = -1 * oneDividedByAlpha1PlusAlpha2_4 * Atom1::orbital_2px.primitiveGaussians[1].alpha * Atom2::orbital_2py.primitiveGaussians[0].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_5 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[1].alpha + Atom2::orbital_2py.primitiveGaussians[1].alpha);
		constexpr double f5 = -1 * oneDividedByAlpha1PlusAlpha2_5 * Atom1::orbital_2px.primitiveGaussians[1].alpha * Atom2::orbital_2py.primitiveGaussians[1].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_6 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[1].alpha + Atom2::orbital_2py.primitiveGaussians[2].alpha);
		constexpr double f6 = -1 * oneDividedByAlpha1PlusAlpha2_6 * Atom1::orbital_2px.primitiveGaussians[1].alpha * Atom2::orbital_2py.primitiveGaussians[2].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_7 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[2].alpha + Atom2::orbital_2py.primitiveGaussians[0].alpha);
		constexpr double f7 = -1 * oneDividedByAlpha1PlusAlpha2_7 * Atom1::orbital_2px.primitiveGaussians[2].alpha * Atom2::orbital_2py.primitiveGaussians[0].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_8 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[2].alpha + Atom2::orbital_2py.primitiveGaussians[1].alpha);
		constexpr double f8 = -1 * oneDividedByAlpha1PlusAlpha2_8 * Atom1::orbital_2px.primitiveGaussians[2].alpha * Atom2::orbital_2py.primitiveGaussians[1].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_9 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[2].alpha + Atom2::orbital_2py.primitiveGaussians[2].alpha);
		constexpr double f9 = -1 * oneDividedByAlpha1PlusAlpha2_9 * Atom1::orbital_2px.primitiveGaussians[2].alpha * Atom2::orbital_2py.primitiveGaussians[2].alpha;

		return { f1, f2, f3, f4, f5, f6, f7, f8, f9 };
	}
	template<typename Atom1, typename Atom2>
	consteval std::array<double, 9> GenerateOneDividedByAlpha1PlusAlpha2_2px2py() noexcept
	{
		constexpr double d1 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[0].alpha + Atom2::orbital_2py.primitiveGaussians[0].alpha);
		constexpr double d2 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[0].alpha + Atom2::orbital_2py.primitiveGaussians[1].alpha);
		constexpr double d3 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[0].alpha + Atom2::orbital_2py.primitiveGaussians[2].alpha);
		constexpr double d4 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[1].alpha + Atom2::orbital_2py.primitiveGaussians[0].alpha);
		constexpr double d5 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[1].alpha + Atom2::orbital_2py.primitiveGaussians[1].alpha);
		constexpr double d6 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[1].alpha + Atom2::orbital_2py.primitiveGaussians[2].alpha);
		constexpr double d7 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[2].alpha + Atom2::orbital_2py.primitiveGaussians[0].alpha);
		constexpr double d8 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[2].alpha + Atom2::orbital_2py.primitiveGaussians[1].alpha);
		constexpr double d9 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[2].alpha + Atom2::orbital_2py.primitiveGaussians[2].alpha);

		return { d1, d2, d3, d4, d5, d6, d7, d8, d9 };
	}

	// 2px - 2pz
	template<typename Atom1, typename Atom2>
	consteval std::array<double, 9> GenerateOverlapFactors_2px2pz() noexcept
	{
		constexpr double oneDividedByAlpha1PlusAlpha2_1 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[0].alpha + Atom2::orbital_2pz.primitiveGaussians[0].alpha);
		constexpr double d1 =
			Atom1::orbital_2px.primitiveGaussians[0].normalizationFactor *
			Atom1::orbital_2px.primitiveGaussians[0].coefficient *
			Atom2::orbital_2pz.primitiveGaussians[0].normalizationFactor *
			Atom2::orbital_2pz.primitiveGaussians[0].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_1, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_2 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[0].alpha + Atom2::orbital_2pz.primitiveGaussians[1].alpha);
		constexpr double d2 =
			Atom1::orbital_2px.primitiveGaussians[0].normalizationFactor *
			Atom1::orbital_2px.primitiveGaussians[0].coefficient *
			Atom2::orbital_2pz.primitiveGaussians[1].normalizationFactor *
			Atom2::orbital_2pz.primitiveGaussians[1].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_2, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_3 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[0].alpha + Atom2::orbital_2pz.primitiveGaussians[2].alpha);
		constexpr double d3 =
			Atom1::orbital_2px.primitiveGaussians[0].normalizationFactor *
			Atom1::orbital_2px.primitiveGaussians[0].coefficient *
			Atom2::orbital_2pz.primitiveGaussians[2].normalizationFactor *
			Atom2::orbital_2pz.primitiveGaussians[2].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_3, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_4 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[1].alpha + Atom2::orbital_2pz.primitiveGaussians[0].alpha);
		constexpr double d4 =
			Atom1::orbital_2px.primitiveGaussians[1].normalizationFactor *
			Atom1::orbital_2px.primitiveGaussians[1].coefficient *
			Atom2::orbital_2pz.primitiveGaussians[0].normalizationFactor *
			Atom2::orbital_2pz.primitiveGaussians[0].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_4, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_5 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[1].alpha + Atom2::orbital_2pz.primitiveGaussians[1].alpha);
		constexpr double d5 =
			Atom1::orbital_2px.primitiveGaussians[1].normalizationFactor *
			Atom1::orbital_2px.primitiveGaussians[1].coefficient *
			Atom2::orbital_2pz.primitiveGaussians[1].normalizationFactor *
			Atom2::orbital_2pz.primitiveGaussians[1].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_5, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_6 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[1].alpha + Atom2::orbital_2pz.primitiveGaussians[2].alpha);
		constexpr double d6 =
			Atom1::orbital_2px.primitiveGaussians[1].normalizationFactor *
			Atom1::orbital_2px.primitiveGaussians[1].coefficient *
			Atom2::orbital_2pz.primitiveGaussians[2].normalizationFactor *
			Atom2::orbital_2pz.primitiveGaussians[2].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_6, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_7 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[2].alpha + Atom2::orbital_2pz.primitiveGaussians[0].alpha);
		constexpr double d7 =
			Atom1::orbital_2px.primitiveGaussians[2].normalizationFactor *
			Atom1::orbital_2px.primitiveGaussians[2].coefficient *
			Atom2::orbital_2pz.primitiveGaussians[0].normalizationFactor *
			Atom2::orbital_2pz.primitiveGaussians[0].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_7, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_8 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[2].alpha + Atom2::orbital_2pz.primitiveGaussians[1].alpha);
		constexpr double d8 =
			Atom1::orbital_2px.primitiveGaussians[2].normalizationFactor *
			Atom1::orbital_2px.primitiveGaussians[2].coefficient *
			Atom2::orbital_2pz.primitiveGaussians[1].normalizationFactor *
			Atom2::orbital_2pz.primitiveGaussians[1].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_8, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_9 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[2].alpha + Atom2::orbital_2pz.primitiveGaussians[2].alpha);
		constexpr double d9 =
			Atom1::orbital_2px.primitiveGaussians[2].normalizationFactor *
			Atom1::orbital_2px.primitiveGaussians[2].coefficient *
			Atom2::orbital_2pz.primitiveGaussians[2].normalizationFactor *
			Atom2::orbital_2pz.primitiveGaussians[2].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_9, 1.5);

		return { d1, d2, d3, d4, d5, d6, d7, d8, d9 };
	}
	template<typename Atom1, typename Atom2>
	consteval std::array<double, 9> GenerateOverlapExponentialFactors_2px2pz() noexcept
	{
		constexpr double oneDividedByAlpha1PlusAlpha2_1 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[0].alpha + Atom2::orbital_2pz.primitiveGaussians[0].alpha);
		constexpr double f1 = -1 * oneDividedByAlpha1PlusAlpha2_1 * Atom1::orbital_2px.primitiveGaussians[0].alpha * Atom2::orbital_2pz.primitiveGaussians[0].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_2 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[0].alpha + Atom2::orbital_2pz.primitiveGaussians[1].alpha);
		constexpr double f2 = -1 * oneDividedByAlpha1PlusAlpha2_2 * Atom1::orbital_2px.primitiveGaussians[0].alpha * Atom2::orbital_2pz.primitiveGaussians[1].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_3 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[0].alpha + Atom2::orbital_2pz.primitiveGaussians[2].alpha);
		constexpr double f3 = -1 * oneDividedByAlpha1PlusAlpha2_3 * Atom1::orbital_2px.primitiveGaussians[0].alpha * Atom2::orbital_2pz.primitiveGaussians[2].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_4 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[1].alpha + Atom2::orbital_2pz.primitiveGaussians[0].alpha);
		constexpr double f4 = -1 * oneDividedByAlpha1PlusAlpha2_4 * Atom1::orbital_2px.primitiveGaussians[1].alpha * Atom2::orbital_2pz.primitiveGaussians[0].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_5 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[1].alpha + Atom2::orbital_2pz.primitiveGaussians[1].alpha);
		constexpr double f5 = -1 * oneDividedByAlpha1PlusAlpha2_5 * Atom1::orbital_2px.primitiveGaussians[1].alpha * Atom2::orbital_2pz.primitiveGaussians[1].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_6 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[1].alpha + Atom2::orbital_2pz.primitiveGaussians[2].alpha);
		constexpr double f6 = -1 * oneDividedByAlpha1PlusAlpha2_6 * Atom1::orbital_2px.primitiveGaussians[1].alpha * Atom2::orbital_2pz.primitiveGaussians[2].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_7 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[2].alpha + Atom2::orbital_2pz.primitiveGaussians[0].alpha);
		constexpr double f7 = -1 * oneDividedByAlpha1PlusAlpha2_7 * Atom1::orbital_2px.primitiveGaussians[2].alpha * Atom2::orbital_2pz.primitiveGaussians[0].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_8 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[2].alpha + Atom2::orbital_2pz.primitiveGaussians[1].alpha);
		constexpr double f8 = -1 * oneDividedByAlpha1PlusAlpha2_8 * Atom1::orbital_2px.primitiveGaussians[2].alpha * Atom2::orbital_2pz.primitiveGaussians[1].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_9 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[2].alpha + Atom2::orbital_2pz.primitiveGaussians[2].alpha);
		constexpr double f9 = -1 * oneDividedByAlpha1PlusAlpha2_9 * Atom1::orbital_2px.primitiveGaussians[2].alpha * Atom2::orbital_2pz.primitiveGaussians[2].alpha;

		return { f1, f2, f3, f4, f5, f6, f7, f8, f9 };
	}
	template<typename Atom1, typename Atom2>
	consteval std::array<double, 9> GenerateOneDividedByAlpha1PlusAlpha2_2px2pz() noexcept
	{
		constexpr double d1 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[0].alpha + Atom2::orbital_2pz.primitiveGaussians[0].alpha);
		constexpr double d2 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[0].alpha + Atom2::orbital_2pz.primitiveGaussians[1].alpha);
		constexpr double d3 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[0].alpha + Atom2::orbital_2pz.primitiveGaussians[2].alpha);
		constexpr double d4 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[1].alpha + Atom2::orbital_2pz.primitiveGaussians[0].alpha);
		constexpr double d5 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[1].alpha + Atom2::orbital_2pz.primitiveGaussians[1].alpha);
		constexpr double d6 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[1].alpha + Atom2::orbital_2pz.primitiveGaussians[2].alpha);
		constexpr double d7 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[2].alpha + Atom2::orbital_2pz.primitiveGaussians[0].alpha);
		constexpr double d8 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[2].alpha + Atom2::orbital_2pz.primitiveGaussians[1].alpha);
		constexpr double d9 = 1.0 / (Atom1::orbital_2px.primitiveGaussians[2].alpha + Atom2::orbital_2pz.primitiveGaussians[2].alpha);

		return { d1, d2, d3, d4, d5, d6, d7, d8, d9 };
	}

	// 2py - 2pz
	template<typename Atom1, typename Atom2>
	consteval std::array<double, 9> GenerateOverlapFactors_2py2pz() noexcept
	{
		constexpr double oneDividedByAlpha1PlusAlpha2_1 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[0].alpha + Atom2::orbital_2pz.primitiveGaussians[0].alpha);
		constexpr double d1 =
			Atom1::orbital_2py.primitiveGaussians[0].normalizationFactor *
			Atom1::orbital_2py.primitiveGaussians[0].coefficient *
			Atom2::orbital_2pz.primitiveGaussians[0].normalizationFactor *
			Atom2::orbital_2pz.primitiveGaussians[0].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_1, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_2 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[0].alpha + Atom2::orbital_2pz.primitiveGaussians[1].alpha);
		constexpr double d2 =
			Atom1::orbital_2py.primitiveGaussians[0].normalizationFactor *
			Atom1::orbital_2py.primitiveGaussians[0].coefficient *
			Atom2::orbital_2pz.primitiveGaussians[1].normalizationFactor *
			Atom2::orbital_2pz.primitiveGaussians[1].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_2, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_3 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[0].alpha + Atom2::orbital_2pz.primitiveGaussians[2].alpha);
		constexpr double d3 =
			Atom1::orbital_2py.primitiveGaussians[0].normalizationFactor *
			Atom1::orbital_2py.primitiveGaussians[0].coefficient *
			Atom2::orbital_2pz.primitiveGaussians[2].normalizationFactor *
			Atom2::orbital_2pz.primitiveGaussians[2].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_3, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_4 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[1].alpha + Atom2::orbital_2pz.primitiveGaussians[0].alpha);
		constexpr double d4 =
			Atom1::orbital_2py.primitiveGaussians[1].normalizationFactor *
			Atom1::orbital_2py.primitiveGaussians[1].coefficient *
			Atom2::orbital_2pz.primitiveGaussians[0].normalizationFactor *
			Atom2::orbital_2pz.primitiveGaussians[0].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_4, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_5 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[1].alpha + Atom2::orbital_2pz.primitiveGaussians[1].alpha);
		constexpr double d5 =
			Atom1::orbital_2py.primitiveGaussians[1].normalizationFactor *
			Atom1::orbital_2py.primitiveGaussians[1].coefficient *
			Atom2::orbital_2pz.primitiveGaussians[1].normalizationFactor *
			Atom2::orbital_2pz.primitiveGaussians[1].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_5, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_6 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[1].alpha + Atom2::orbital_2pz.primitiveGaussians[2].alpha);
		constexpr double d6 =
			Atom1::orbital_2py.primitiveGaussians[1].normalizationFactor *
			Atom1::orbital_2py.primitiveGaussians[1].coefficient *
			Atom2::orbital_2pz.primitiveGaussians[2].normalizationFactor *
			Atom2::orbital_2pz.primitiveGaussians[2].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_6, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_7 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[2].alpha + Atom2::orbital_2pz.primitiveGaussians[0].alpha);
		constexpr double d7 =
			Atom1::orbital_2py.primitiveGaussians[2].normalizationFactor *
			Atom1::orbital_2py.primitiveGaussians[2].coefficient *
			Atom2::orbital_2pz.primitiveGaussians[0].normalizationFactor *
			Atom2::orbital_2pz.primitiveGaussians[0].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_7, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_8 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[2].alpha + Atom2::orbital_2pz.primitiveGaussians[1].alpha);
		constexpr double d8 =
			Atom1::orbital_2py.primitiveGaussians[2].normalizationFactor *
			Atom1::orbital_2py.primitiveGaussians[2].coefficient *
			Atom2::orbital_2pz.primitiveGaussians[1].normalizationFactor *
			Atom2::orbital_2pz.primitiveGaussians[1].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_8, 1.5);

		constexpr double oneDividedByAlpha1PlusAlpha2_9 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[2].alpha + Atom2::orbital_2pz.primitiveGaussians[2].alpha);
		constexpr double d9 =
			Atom1::orbital_2py.primitiveGaussians[2].normalizationFactor *
			Atom1::orbital_2py.primitiveGaussians[2].coefficient *
			Atom2::orbital_2pz.primitiveGaussians[2].normalizationFactor *
			Atom2::orbital_2pz.primitiveGaussians[2].coefficient *
			gcem::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2_9, 1.5);

		return { d1, d2, d3, d4, d5, d6, d7, d8, d9 };
	}
	template<typename Atom1, typename Atom2>
	consteval std::array<double, 9> GenerateOverlapExponentialFactors_2py2pz() noexcept
	{
		constexpr double oneDividedByAlpha1PlusAlpha2_1 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[0].alpha + Atom2::orbital_2pz.primitiveGaussians[0].alpha);
		constexpr double f1 = -1 * oneDividedByAlpha1PlusAlpha2_1 * Atom1::orbital_2py.primitiveGaussians[0].alpha * Atom2::orbital_2pz.primitiveGaussians[0].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_2 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[0].alpha + Atom2::orbital_2pz.primitiveGaussians[1].alpha);
		constexpr double f2 = -1 * oneDividedByAlpha1PlusAlpha2_2 * Atom1::orbital_2py.primitiveGaussians[0].alpha * Atom2::orbital_2pz.primitiveGaussians[1].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_3 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[0].alpha + Atom2::orbital_2pz.primitiveGaussians[2].alpha);
		constexpr double f3 = -1 * oneDividedByAlpha1PlusAlpha2_3 * Atom1::orbital_2py.primitiveGaussians[0].alpha * Atom2::orbital_2pz.primitiveGaussians[2].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_4 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[1].alpha + Atom2::orbital_2pz.primitiveGaussians[0].alpha);
		constexpr double f4 = -1 * oneDividedByAlpha1PlusAlpha2_4 * Atom1::orbital_2py.primitiveGaussians[1].alpha * Atom2::orbital_2pz.primitiveGaussians[0].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_5 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[1].alpha + Atom2::orbital_2pz.primitiveGaussians[1].alpha);
		constexpr double f5 = -1 * oneDividedByAlpha1PlusAlpha2_5 * Atom1::orbital_2py.primitiveGaussians[1].alpha * Atom2::orbital_2pz.primitiveGaussians[1].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_6 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[1].alpha + Atom2::orbital_2pz.primitiveGaussians[2].alpha);
		constexpr double f6 = -1 * oneDividedByAlpha1PlusAlpha2_6 * Atom1::orbital_2py.primitiveGaussians[1].alpha * Atom2::orbital_2pz.primitiveGaussians[2].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_7 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[2].alpha + Atom2::orbital_2pz.primitiveGaussians[0].alpha);
		constexpr double f7 = -1 * oneDividedByAlpha1PlusAlpha2_7 * Atom1::orbital_2py.primitiveGaussians[2].alpha * Atom2::orbital_2pz.primitiveGaussians[0].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_8 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[2].alpha + Atom2::orbital_2pz.primitiveGaussians[1].alpha);
		constexpr double f8 = -1 * oneDividedByAlpha1PlusAlpha2_8 * Atom1::orbital_2py.primitiveGaussians[2].alpha * Atom2::orbital_2pz.primitiveGaussians[1].alpha;

		constexpr double oneDividedByAlpha1PlusAlpha2_9 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[2].alpha + Atom2::orbital_2pz.primitiveGaussians[2].alpha);
		constexpr double f9 = -1 * oneDividedByAlpha1PlusAlpha2_9 * Atom1::orbital_2py.primitiveGaussians[2].alpha * Atom2::orbital_2pz.primitiveGaussians[2].alpha;

		return { f1, f2, f3, f4, f5, f6, f7, f8, f9 };
	}
	template<typename Atom1, typename Atom2>
	consteval std::array<double, 9> GenerateOneDividedByAlpha1PlusAlpha2_2py2pz() noexcept
	{
		constexpr double d1 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[0].alpha + Atom2::orbital_2pz.primitiveGaussians[0].alpha);
		constexpr double d2 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[0].alpha + Atom2::orbital_2pz.primitiveGaussians[1].alpha);
		constexpr double d3 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[0].alpha + Atom2::orbital_2pz.primitiveGaussians[2].alpha);
		constexpr double d4 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[1].alpha + Atom2::orbital_2pz.primitiveGaussians[0].alpha);
		constexpr double d5 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[1].alpha + Atom2::orbital_2pz.primitiveGaussians[1].alpha);
		constexpr double d6 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[1].alpha + Atom2::orbital_2pz.primitiveGaussians[2].alpha);
		constexpr double d7 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[2].alpha + Atom2::orbital_2pz.primitiveGaussians[0].alpha);
		constexpr double d8 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[2].alpha + Atom2::orbital_2pz.primitiveGaussians[1].alpha);
		constexpr double d9 = 1.0 / (Atom1::orbital_2py.primitiveGaussians[2].alpha + Atom2::orbital_2pz.primitiveGaussians[2].alpha);

		return { d1, d2, d3, d4, d5, d6, d7, d8, d9 };
	}



}
}