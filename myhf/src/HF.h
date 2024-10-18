#pragma once

#include "Atom.h"
#include "Basis.h"

namespace myhf
{

[[nodiscard]] Eigen::MatrixXd OverlapMatrix(std::span<Atom> atoms, const Basis& basis) noexcept;
[[nodiscard]] Eigen::MatrixXd KineticEnergyMatrix(std::span<Atom> atoms, const Basis& basis) noexcept;
[[nodiscard]] void OverlapAndKineticEnergyMatrix(std::span<Atom> atoms, const Basis& basis, Eigen::MatrixXd& overlapMatrix, Eigen::MatrixXd& kineticEnergyMatrix) noexcept;


[[nodiscard]] Eigen::MatrixXd NuclearElectronAttractionEnergyMatrix(std::span<Atom> atoms, const Basis& basis) noexcept;
[[nodiscard]] Eigen::MatrixXd NuclearElectronAttractionEnergyMatrix_2(std::span<Atom> atoms, const Basis& basis) noexcept;
[[nodiscard]] Eigen::MatrixXd NuclearElectronAttractionEnergyMatrix_Par(std::span<Atom> atoms, const Basis& basis) noexcept;


[[nodiscard]] Eigen::MatrixXd OverlapMatrix_2(std::span<Atom> atoms, const Basis& basis) noexcept;

}