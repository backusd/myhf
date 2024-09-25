#include "pch.h"
#include "HF.h"

#include <execution>
#include <print>
#include <iostream>

using Eigen::MatrixXd;

namespace myhf
{
// Helper class for keeping track of the recurrence relation values with data that resides on the stack
class StackMatrix2D
{
public:
	StackMatrix2D(unsigned int rows, unsigned int columns) noexcept :
		rowCount(rows), colCount(columns), values(nullptr)
	{
		assert(rowCount > 0);
		assert(colCount > 0);
		values = static_cast<double*>(alloca(sizeof(double) * rowCount * colCount));
	}
	// delete copies/moves because we do not need to support that use case
	StackMatrix2D(const StackMatrix2D&) = delete;
	StackMatrix2D(StackMatrix2D&&) = delete;
	StackMatrix2D& operator=(const StackMatrix2D&) = delete;
	StackMatrix2D& operator=(StackMatrix2D&&) = delete;

	constexpr double& operator()(unsigned int row, unsigned int col) noexcept
	{
		assert(row >= 0 && row < rowCount);
		assert(col >= 0 && col < colCount);
		return values[row * colCount + col];
	}
	constexpr unsigned int RowCount() const noexcept { return rowCount; }
	constexpr unsigned int ColumnCount() const noexcept { return colCount; }

private:
	unsigned int rowCount;
	unsigned int colCount;
	double* values;
};

// Overlap Functions
static double OverlapOfTwoPrimitiveGaussiansAlongAxis(double oneDividedByAlpha1PlusAlpha2, double alpha1, double alpha2, double position1, double position2, unsigned int angularMomentum1, unsigned int angularMomentum2) noexcept
{
	if (angularMomentum1 == 0 && angularMomentum2 == 0)
		return 1.0;

	// Supposed we need to return the value 'matrix[a, b]'. To compute that value, you must
	// have already computed 'matrix[a+1, b-1]' - this is the value located down 1 and left 1 spot
	// in the matrix. If you continuing this logic, you will eventually get to the first column of the
	// matrix, which will tell us how many rows we need. Notice that along that diagonal in the matrix, 
	// the row index plus the column index have a constant value (a + b). To adjust for the off-by-one 
	// error, we also add one.
	StackMatrix2D matrix(angularMomentum1 + angularMomentum2 + 1, angularMomentum2 + 1); 

	// Initial conditions
	matrix(0, 0) = 1.0;
	double startingValue = -(position1 - ((alpha1 * position1 + alpha2 * position2) * oneDividedByAlpha1PlusAlpha2));
	matrix(1, 0) = startingValue;

	// Recurrence Index
	for (unsigned int iii = 2; iii < matrix.RowCount(); ++iii)
		matrix(iii, 0) = startingValue * matrix(iii - 1, 0) + ((iii - 1) * 0.5 * oneDividedByAlpha1PlusAlpha2) * matrix(iii - 2, 0);

	// Transfer Equation
	//
	// The recurrence relation requires filling out each column before moving on to the next one.
	// Therefore, the outer loop must be the loop over the columns
	for (unsigned int col = 1; col < matrix.ColumnCount(); ++col)
	{
		// The conditional here is a little weird. Basically, we have set up a matrix to store the recurrence values in. 
		// However, the bottom right half of values will never get used (and cannot be calculated without expanding the matrix). 
		// For example, consider trying to compute matrix[0, 2]. The matrix will look something like this:
		//		v1	v4	v6
		//		v2	v5	-
		//		v3	-	-
		// We only need to compute up to v6 and don't need the '-' values and we wouldn't even be able to compute them unless we added more rows.

		for (unsigned int row = 0; row + col <= angularMomentum1 + angularMomentum2; ++row)
			matrix(row, col) = matrix(row + 1, col - 1) + (position1 - position2) * matrix(row, col - 1);
	}

	return matrix(angularMomentum1, angularMomentum2);
}
static double OverlapOfTwoPrimitiveGaussians(double alpha1, double alpha2, const Vec3d& position1, const Vec3d& position2, const QuantumNumbers& angularMomentum1, const QuantumNumbers& angularMomentum2) noexcept
{
	assert(std::max(angularMomentum1.l, angularMomentum2.l) + 1 < 4);
	assert(std::max(angularMomentum1.m, angularMomentum2.m) + 1 < 4);
	assert(std::max(angularMomentum1.n, angularMomentum2.n) + 1 < 4);

	const double oneDividedByAlpha1PlusAlpha2 = 1 / (alpha1 + alpha2);
	double s_x_value = OverlapOfTwoPrimitiveGaussiansAlongAxis(oneDividedByAlpha1PlusAlpha2, alpha1, alpha2, position1.x, position2.x, angularMomentum1.l, angularMomentum2.l);
	double s_y_value = OverlapOfTwoPrimitiveGaussiansAlongAxis(oneDividedByAlpha1PlusAlpha2, alpha1, alpha2, position1.y, position2.y, angularMomentum1.m, angularMomentum2.m);
	double s_z_value = OverlapOfTwoPrimitiveGaussiansAlongAxis(oneDividedByAlpha1PlusAlpha2, alpha1, alpha2, position1.z, position2.z, angularMomentum1.n, angularMomentum2.n);

	Vec3d diff = position2 - position1;
	return std::exp(-1 * alpha1 * alpha2 * oneDividedByAlpha1PlusAlpha2 * diff.Dot(diff)) *
		std::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2, 1.5) *
		s_x_value * s_y_value * s_z_value;
}
static double OverlapOfTwoOrbitals(const ContractedGaussianOrbital& orbital1, const Vec3d& position1, const ContractedGaussianOrbital& orbital2, const Vec3d& position2) noexcept
{
	double res = 0;

	for (auto& gaussian1 : orbital1.gaussianOrbitals)
	{
		for (auto& gaussian2 : orbital2.gaussianOrbitals)
		{
			double overlap = OverlapOfTwoPrimitiveGaussians(gaussian1.alpha, gaussian2.alpha, position1, position2, orbital1.angularMomentum, orbital2.angularMomentum);
			res += gaussian1.normalizationFactor * gaussian2.normalizationFactor * gaussian1.coefficient * gaussian2.coefficient * overlap;
		}
	}
	return res;
}
MatrixXd OverlapMatrix(std::span<Atom> atoms, const Basis& basis) noexcept
{
	MatrixXd overlapMatrix;

	assert(atoms.size() > 1);
	unsigned int numberOfContractedGaussians = 0;
	for (const Atom& atom : atoms)
		numberOfContractedGaussians += basis.GetAtom(atom.type).NumberOfContractedGaussians();

	overlapMatrix = MatrixXd::Identity(numberOfContractedGaussians, numberOfContractedGaussians);

	unsigned int i = 0;
	for (const auto& atom1 : atoms)
	{
		for (const auto& shell1 : basis.GetAtom(atom1.type).shells)
		{
			for (const auto& orbital1 : shell1.basisFunctions)
			{
				unsigned int j = 0;

				for (const auto& atom2 : atoms)
				{
					for (const auto& shell2 : basis.GetAtom(atom2.type).shells)
					{
						for (const auto& orbital2 : shell2.basisFunctions)
						{
							if (j > i)
							{
								overlapMatrix(i, j) = OverlapOfTwoOrbitals(orbital1, atom1.position, orbital2, atom2.position);
								overlapMatrix(j, i) = overlapMatrix(i, j);
							}

							++j;
						}
					}
				}
				++i;
			}
		}
	}

	return overlapMatrix;
}


// Kinetic Functions
static std::pair<double, double> KineticEnergyOfTwoPrimitiveGaussiansAlongAxis(double oneDividedByAlpha1PlusAlpha2, double alpha1, double alpha2, double position1, double position2, unsigned int angularMomentum1, unsigned int angularMomentum2)
{
	// NOTE: We do not actually need a k-matrix because the values of the k-matrix is are soley based
	// on the s-matrix, so we only need to compute the necessary values of the s-matrix to compute the
	// necessary k-matrix value.

	// Suppose we need to return the value 'matrix[a, b]'. To compute that value, you must
	// have already computed 'matrix[a+1, b-1]' - this is the value located down 1 and left 1 spot
	// in the matrix. If you continuing this logic, you will eventually get to the first column of the
	// matrix, which will tell us how many rows we need. Notice that along that diagonal in the matrix, 
	// the row index plus the column index have a constant value (a + b). To adjust for the off-by-one 
	// error, we also add one.
	StackMatrix2D s(angularMomentum1 + angularMomentum2 + 3, angularMomentum2 + 2);

	// Initial conditions
	s(0, 0) = 1.0;
	double startingValue = -(position1 - ((alpha1 * position1 + alpha2 * position2) * oneDividedByAlpha1PlusAlpha2));
	s(1, 0) = startingValue;

	// Recurrence Index
	for (unsigned int iii = 2; iii < s.RowCount(); ++iii)
		s(iii, 0) = startingValue * s(iii - 1, 0) + ((iii - 1) * 0.5 * oneDividedByAlpha1PlusAlpha2) * s(iii - 2, 0);

	// Transfer Equation
	//
	// The recurrence relation requires filling out each column before moving on to the next one.
	// Therefore, the outer loop must be the loop over the columns
	for (unsigned int col = 1; col < s.ColumnCount(); ++col)
	{
		// The conditional here is a little weird. Basically, we have set up a matrix to store the recurrence values in. 
		// However, the bottom right half of values will never get used (and cannot be calculated without expanding the matrix). 
		// For example, consider trying to compute matrix[0, 2]. The matrix will look something like this:
		//		v1	v4	v6
		//		v2	v5	-
		//		v3	-	-
		// We only need to compute up to v6 and don't need the '-' values and we wouldn't even be able to compute them unless we added more rows.

		for (unsigned int row = 0; row + col <= angularMomentum1 + angularMomentum2 + 2; ++row)
			s(row, col) = s(row + 1, col - 1) + (position1 - position2) * s(row, col - 1);
	}

	double k = 0.0;
	if (angularMomentum1 == 0 && angularMomentum2 == 0)
		k = 2 * alpha1 * alpha2 * s(1, 1);
	else if (angularMomentum2 == 0)
		k = -1 * static_cast<double>(angularMomentum1) * alpha2 * s(angularMomentum1 - 1, 1) + 2 * alpha1 * alpha2 * s(angularMomentum1 + 1, 1);
	else if (angularMomentum1 == 0)
		k = -1 * static_cast<int>(angularMomentum2) * alpha1 * s(1, angularMomentum2 - 1) + 2 * alpha1 * alpha2 * s(1, angularMomentum2 + 1);
	else
		k = 0.5 * (
			angularMomentum1 * angularMomentum2 * s(angularMomentum1 - 1, angularMomentum2 - 1) -
			2 * angularMomentum1 * alpha2 * s(angularMomentum1 - 1, angularMomentum2 + 1) -
			2 * angularMomentum2 * alpha1 * s(angularMomentum1 + 1, angularMomentum2 - 1) +
			4 * alpha1 * alpha2 * s(angularMomentum1 + 1, angularMomentum2 + 1)
		);

	return { s(angularMomentum1, angularMomentum2), k }; 
}
static double KineticEnergyOfTwoPrimitiveGaussians(double alpha1, double alpha2, const Vec3d& position1, const Vec3d& position2, const QuantumNumbers& angularMomentum1, const QuantumNumbers& angularMomentum2) noexcept
{
	const double oneDividedByAlpha1PlusAlpha2 = 1 / (alpha1 + alpha2);
	auto [s_x_value, k_x_value] = KineticEnergyOfTwoPrimitiveGaussiansAlongAxis(oneDividedByAlpha1PlusAlpha2, alpha1, alpha2, position1.x, position2.x, angularMomentum1.l, angularMomentum2.l);
	auto [s_y_value, k_y_value] = KineticEnergyOfTwoPrimitiveGaussiansAlongAxis(oneDividedByAlpha1PlusAlpha2, alpha1, alpha2, position1.y, position2.y, angularMomentum1.m, angularMomentum2.m);
	auto [s_z_value, k_z_value] = KineticEnergyOfTwoPrimitiveGaussiansAlongAxis(oneDividedByAlpha1PlusAlpha2, alpha1, alpha2, position1.z, position2.z, angularMomentum1.n, angularMomentum2.n);	

	Vec3d diff = position2 - position1;	
	return std::exp(-1 * alpha1 * alpha2 * oneDividedByAlpha1PlusAlpha2 * diff.Dot(diff)) * 
		   std::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2, 1.5) *
		   (k_x_value * s_y_value * s_z_value +
		    s_x_value * k_y_value * s_z_value +
		    s_x_value * s_y_value * k_z_value);
}
static double KineticEnergyOfTwoOrbitals(const ContractedGaussianOrbital& orbital1, const Vec3d& position1, const ContractedGaussianOrbital& orbital2, const Vec3d& position2) noexcept
{
	double res = 0;

	for (auto& gaussian1 : orbital1.gaussianOrbitals)
	{
		for (auto& gaussian2 : orbital2.gaussianOrbitals)
		{
			double kineticEnergy = KineticEnergyOfTwoPrimitiveGaussians(gaussian1.alpha, gaussian2.alpha, position1, position2, orbital1.angularMomentum, orbital2.angularMomentum);
			res += gaussian1.normalizationFactor * gaussian2.normalizationFactor * gaussian1.coefficient * gaussian2.coefficient * kineticEnergy;
		}
	}
	return res;
}
static void KineticEnergyMatrixFillDiagonal(std::span<Atom> atoms, const Basis& basis, Eigen::MatrixXd& kineticEnergyMatrix)
{
	// Small optimization: The kinetic matrix will very likely have repeated values along the diagonal. We only need to calculate
	// the value for each atom/orbital type once, so just iterate down the diagonal and copy values where possible.
	unsigned int i = 0;
	for (size_t atomIndex1 = 0; atomIndex1 < atoms.size(); ++atomIndex1)
	{
		const auto& atom1 = atoms[atomIndex1];
		const auto& shells1 = basis.GetAtom(atom1.type).shells;

		for (size_t shellIndex1 = 0; shellIndex1 < shells1.size(); ++shellIndex1)
		{
			const auto& shell1 = shells1[shellIndex1];

			for (size_t orbitalIndex1 = 0; orbitalIndex1 < shell1.basisFunctions.size(); ++orbitalIndex1)
			{
				const auto& orbital1 = shell1.basisFunctions[orbitalIndex1];
				unsigned int totalAngularMomentum1 = orbital1.angularMomentum.AngularMomentum();

				// If the diagonal value is 0, then it has not yet been calculated
				if (kineticEnergyMatrix(i, i) == 0)
				{
					kineticEnergyMatrix(i, i) = KineticEnergyOfTwoOrbitals(orbital1, atom1.position, orbital1, atom1.position);

					// Now that we have calculated the diagonal value, search forward to fill in that value whereever else we can

					// NOTE: First, we must finish iterating over the shells/orbitals of the current atom. For example, going
					//       from 2px to 2py
					unsigned int j = i + 1;

					// Start 1 passed the current orbital and iterate until you finish iterating over the entire shell
					for (size_t orbitalIndex2 = orbitalIndex1 + 1; orbitalIndex2 < shell1.basisFunctions.size(); ++orbitalIndex2)
					{
						if (shell1.basisFunctions[orbitalIndex2].angularMomentum.AngularMomentum() == totalAngularMomentum1)
							kineticEnergyMatrix(j, j) = kineticEnergyMatrix(i, i);

						++j;
					}

					// Iterate over the remaining shells in the atom to advance the j value
					for (size_t shellIndex2 = shellIndex1 + 1; shellIndex2 < shells1.size(); ++shellIndex2)
						j += static_cast<unsigned int>(shells1[shellIndex2].basisFunctions.size());

					// Now that we have completed the current atom, iterate over all remaining atoms
					for (size_t atomIndex2 = atomIndex1 + 1; atomIndex2 < atoms.size(); ++atomIndex2)
					{
						const auto& atom2 = atoms[atomIndex2];
						const auto& shells2 = basis.GetAtom(atom2.type).shells;

						// Check to see if the atom types match
						if (atom2.type == atom1.type)
						{
							// Atom types match, so now iterate over the shell values
							for (size_t shellIndex2 = 0; shellIndex2 < shells2.size(); ++shellIndex2)
							{
								const auto& shell2 = shells2[shellIndex2];

								// Check if shells match
								if (shellIndex1 == shellIndex2)
								{
									// Shells index match, so iterate over the orbitals of the shell
									for (size_t orbitalIndex2 = 0; orbitalIndex2 < shell2.basisFunctions.size(); ++orbitalIndex2)
									{
										if (shell2.basisFunctions[orbitalIndex2].angularMomentum.AngularMomentum() == totalAngularMomentum1)
											kineticEnergyMatrix(j, j) = kineticEnergyMatrix(i, i);

										++j;
									}
								}
								else
								{
									// Shell index does not match, so just fast forward the j value
									j += static_cast<unsigned int>(shell2.basisFunctions.size());
								}
							}
						}
						else
						{
							// atom types do not match, so just fast forward the j value
							for (const auto& shell2 : shells2)
								j += static_cast<unsigned int>(shell2.basisFunctions.size());
						}
					}
				}

				++i;
			}
		}
	}
}
MatrixXd KineticEnergyMatrix(std::span<Atom> atoms, const Basis& basis) noexcept
{
	MatrixXd kineticEnergyMatrix;

	assert(atoms.size() > 1);
	unsigned int numberOfContractedGaussians = 0;
	for (const Atom& atom : atoms)
		numberOfContractedGaussians += basis.GetAtom(atom.type).NumberOfContractedGaussians();

	kineticEnergyMatrix = MatrixXd::Zero(numberOfContractedGaussians, numberOfContractedGaussians);

	// Special function that computes the diagonal elements in a more optimized way than if we were to do it 
	// in the section below
	KineticEnergyMatrixFillDiagonal(atoms, basis, kineticEnergyMatrix);

	unsigned int i = 0;
	for (const auto& atom1 : atoms)
	{
		for (const auto& shell1 : basis.GetAtom(atom1.type).shells)
		{
			for (const auto& orbital1 : shell1.basisFunctions)
			{
				unsigned int j = 0;

				for (const auto& atom2 : atoms)
				{
					for (const auto& shell2 : basis.GetAtom(atom2.type).shells)
					{
						for (const auto& orbital2 : shell2.basisFunctions)
						{
							if (j > i)
							{
								kineticEnergyMatrix(i, j) = KineticEnergyOfTwoOrbitals(orbital1, atom1.position, orbital2, atom2.position);								
								kineticEnergyMatrix(j, i) = kineticEnergyMatrix(i, j);
							}

							++j;
						}
					}
				}
				++i;
			}
		}
	}

	return kineticEnergyMatrix;
}


static std::pair<double, double> OverlapAndKineticEnergyOfTwoPrimitiveGaussians(double alpha1, double alpha2, const Vec3d& position1, const Vec3d& position2, const QuantumNumbers& angularMomentum1, const QuantumNumbers& angularMomentum2) noexcept
{
	const double oneDividedByAlpha1PlusAlpha2 = 1 / (alpha1 + alpha2);
	auto [s_x_value, k_x_value] = KineticEnergyOfTwoPrimitiveGaussiansAlongAxis(oneDividedByAlpha1PlusAlpha2, alpha1, alpha2, position1.x, position2.x, angularMomentum1.l, angularMomentum2.l);
	auto [s_y_value, k_y_value] = KineticEnergyOfTwoPrimitiveGaussiansAlongAxis(oneDividedByAlpha1PlusAlpha2, alpha1, alpha2, position1.y, position2.y, angularMomentum1.m, angularMomentum2.m);
	auto [s_z_value, k_z_value] = KineticEnergyOfTwoPrimitiveGaussiansAlongAxis(oneDividedByAlpha1PlusAlpha2, alpha1, alpha2, position1.z, position2.z, angularMomentum1.n, angularMomentum2.n);

	Vec3d diff = position2 - position1;
	double factor = std::exp(-1 * alpha1 * alpha2 * oneDividedByAlpha1PlusAlpha2 * diff.Dot(diff)) *
		std::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2, 1.5);

	return { factor * s_x_value * s_y_value * s_z_value, 
		   factor *
		   (k_x_value * s_y_value * s_z_value +
			s_x_value * k_y_value * s_z_value +
			s_x_value * s_y_value * k_z_value) };
}
static std::pair<double, double> OverlapAndKineticEnergyOfTwoOrbitals(const ContractedGaussianOrbital& orbital1, const Vec3d& position1, const ContractedGaussianOrbital& orbital2, const Vec3d& position2) noexcept
{
	double overlapResult = 0.0;
	double kineticResult = 0.0;

	for (auto& gaussian1 : orbital1.gaussianOrbitals)
	{
		for (auto& gaussian2 : orbital2.gaussianOrbitals)
		{
			double factor = gaussian1.normalizationFactor * gaussian2.normalizationFactor * gaussian1.coefficient * gaussian2.coefficient;
			auto [overlap, kinetic] = OverlapAndKineticEnergyOfTwoPrimitiveGaussians(gaussian1.alpha, gaussian2.alpha, position1, position2, orbital1.angularMomentum, orbital2.angularMomentum);
			
			overlapResult += factor * overlap;
			kineticResult += factor * kinetic;
		}
	}
	return { overlapResult, kineticResult };
}
void OverlapAndKineticEnergyMatrix(std::span<Atom> atoms, const Basis& basis, Eigen::MatrixXd& overlapMatrix, Eigen::MatrixXd& kineticEnergyMatrix) noexcept
{
	assert(atoms.size() > 1);
	unsigned int numberOfContractedGaussians = 0;
	for (const Atom& atom : atoms)
		numberOfContractedGaussians += basis.GetAtom(atom.type).NumberOfContractedGaussians();

	overlapMatrix = MatrixXd::Identity(numberOfContractedGaussians, numberOfContractedGaussians);
	kineticEnergyMatrix = MatrixXd::Zero(numberOfContractedGaussians, numberOfContractedGaussians);
	
	// Special function that computes the diagonal elements in a more optimized way than if we were to do it 
	// in the section below
	KineticEnergyMatrixFillDiagonal(atoms, basis, kineticEnergyMatrix);

	unsigned int i = 0;
	for (const auto& atom1 : atoms)
	{
		for (const auto& shell1 : basis.GetAtom(atom1.type).shells)
		{
			for (const auto& orbital1 : shell1.basisFunctions)
			{
				unsigned int j = 0;

				for (const auto& atom2 : atoms)
				{
					for (const auto& shell2 : basis.GetAtom(atom2.type).shells)
					{
						for (const auto& orbital2 : shell2.basisFunctions)
						{
							// We have already compute the diagonal values, so we only need to compute off diagonal values
							if (j > i)
							{
								auto [overlap, kinetic] = OverlapAndKineticEnergyOfTwoOrbitals(orbital1, atom1.position, orbital2, atom2.position);
								
								overlapMatrix(i, j) = overlap;
								overlapMatrix(j, i) = overlap;

								kineticEnergyMatrix(i, j) = kinetic;
								kineticEnergyMatrix(j, i) = kinetic; 
							}

							++j;
						}
					}
				}
				++i;
			}
		}
	}
}


}