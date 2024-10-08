#pragma once

#include "Atom.h"
#include "DataDrivenTest.h"


namespace myhf
{
namespace test
{

	[[nodiscard]] Eigen::MatrixXd OverlapMatrix(std::span<Atom> atoms) noexcept;
//	[[nodiscard]] Eigen::MatrixXd KineticEnergyMatrix(std::span<Atom> atoms, const Basis& basis) noexcept;
//	[[nodiscard]] void OverlapAndKineticEnergyMatrix(std::span<Atom> atoms, const Basis& basis, Eigen::MatrixXd& overlapMatrix, Eigen::MatrixXd& kineticEnergyMatrix) noexcept;
//	[[nodiscard]] Eigen::MatrixXd NuclearElectronAttractionEnergyMatrix(std::span<Atom> atoms, const Basis& basis) noexcept;
//	[[nodiscard]] Eigen::MatrixXd NuclearElectronAttractionEnergyMatrix_Par(std::span<Atom> atoms, const Basis& basis) noexcept;
//
//
//	[[nodiscard]] Eigen::MatrixXd OverlapMatrix_2(std::span<Atom> atoms, const Basis& basis) noexcept;


//	template<typename Atom1, typename Atom2>
//	Eigen::MatrixXd Overlap(Atom1 atom1, const Vec3d& position1, Atom2 atom2, const Vec3d& position2) noexcept
//	{
//		unsigned int numberOfContractedGaussians = Atom1::NumberOfContractedGaussians + Atom2::NumberOfContractedGaussians;
//		Eigen::MatrixXd overlap = Eigen::MatrixXd::Identity(numberOfContractedGaussians, numberOfContractedGaussians);
//		return overlap;
//	}
//
//	template<>
//	Eigen::MatrixXd Overlap<Hydrogen, Hydrogen>(Hydrogen atom1, const Vec3d& position1, Hydrogen atom2, const Vec3d& position2) noexcept;


	template<typename Atom1, typename Atom2>
	void OverlapAtomAtomImpl(const Vec3d& position1, const Vec3d& position2, unsigned int row, unsigned int col, Eigen::MatrixXd& overlapMatrix) noexcept
	{
		static constexpr unsigned int numberOfContractedGaussians = Atom1::NumberOfContractedGaussians + Atom2::NumberOfContractedGaussians;
		Eigen::MatrixXd overlap = Eigen::MatrixXd::Identity(numberOfContractedGaussians, numberOfContractedGaussians);
		return;
	}

	template<>
	void OverlapAtomAtomImpl<Hydrogen, Hydrogen>(const Vec3d& position1, const Vec3d& position2, unsigned int row, unsigned int col, Eigen::MatrixXd& overlapMatrix) noexcept;
	template<>
	void OverlapAtomAtomImpl<Hydrogen, Oxygen>(const Vec3d& position1, const Vec3d& position2, unsigned int row, unsigned int col, Eigen::MatrixXd& overlapMatrix) noexcept;


	template<typename Atom1>
	void OverlapWithSelfImpl(unsigned int topLeftIndex, Eigen::MatrixXd& overlapMatrix) noexcept
	{
		return;
	}

	template<>
	void OverlapWithSelfImpl<Hydrogen>(unsigned int topLeftIndex, Eigen::MatrixXd& overlapMatrix) noexcept;


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


}
}