#pragma once

#include "Atom.h"
#include "Basis.h"

namespace myhf
{

[[nodiscard]] Eigen::MatrixXd OverlapMatrix(std::span<Atom> atoms, const Basis& basis) noexcept;
[[nodiscard]] Eigen::MatrixXd KineticEnergyMatrix(std::span<Atom> atoms, const Basis& basis) noexcept;
[[nodiscard]] void OverlapAndKineticEnergyMatrix(std::span<Atom> atoms, const Basis& basis, Eigen::MatrixXd& overlapMatrix, Eigen::MatrixXd& kineticEnergyMatrix) noexcept;


[[nodiscard]] double NuclearElectronAttractionOfTwoPrimitiveGaussians(double alpha1, double alpha2, const Vec3d& position1, const Vec3d& position2, const Vec3d& nuclearCenter, const QuantumNumbers& angularMomentum1, const QuantumNumbers& angularMomentum2) noexcept;


}