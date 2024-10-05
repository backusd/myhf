#include "pch.h"
#include "HF.h"

#include <execution>
#include <print>
#include <iostream>
#include <thread>

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

class StackMatrix2D_Functions
{
public:
	StackMatrix2D_Functions(unsigned int rows, unsigned int columns) noexcept :
		rowCount(rows), 
		colCount(columns)
	{
		assert(rowCount > 0);
		assert(colCount > 0);
		assert(rowCount * colCount <= 16);
	}
	// delete copies/moves because we do not need to support that use case
	StackMatrix2D_Functions(const StackMatrix2D_Functions&) = delete;
	StackMatrix2D_Functions(StackMatrix2D_Functions&&) = delete;
	StackMatrix2D_Functions& operator=(const StackMatrix2D_Functions&) = delete;
	StackMatrix2D_Functions& operator=(StackMatrix2D_Functions&&) = delete;

	constexpr std::function<double(double)>& operator()(unsigned int row, unsigned int col) noexcept
	{
		assert(row >= 0 && row < rowCount);
		assert(col >= 0 && col < colCount);
		assert(row * colCount + col < 16);
		return values[row * colCount + col];
	}
	constexpr unsigned int RowCount() const noexcept { return rowCount; }
	constexpr unsigned int ColumnCount() const noexcept { return colCount; }

private:
	unsigned int rowCount;
	unsigned int colCount;
	std::array<std::function<double(double)>, 16> values;
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

	struct ComputeItem
	{
		unsigned int index_i;
		unsigned int index_j;
		const ContractedGaussianOrbital& orbital_1;
		const ContractedGaussianOrbital& orbital_2;
		const Vec3d& position_1;
		const Vec3d& position_2;
	};
	std::vector<ComputeItem> computeItems;
	computeItems.reserve(numberOfContractedGaussians * numberOfContractedGaussians); // Note, we are over reserving here, but thats fine to waste the extra memory

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
								computeItems.emplace_back(i, j, orbital1, orbital2, atom1.position, atom2.position);
							}

							++j;
						}
					}
				}
				++i;
			}
		}
	}

	std::for_each(std::execution::par_unseq, computeItems.begin(), computeItems.end(),
		[&overlapMatrix](const ComputeItem& item)
		{
			double result = OverlapOfTwoOrbitals(item.orbital_1, item.position_1, item.orbital_2, item.position_2);
			overlapMatrix(item.index_i, item.index_j) = result;
			overlapMatrix(item.index_j, item.index_i) = result;
		});

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

	struct ComputeItem
	{
		unsigned int index_i;
		unsigned int index_j;
		const ContractedGaussianOrbital& orbital_1;
		const ContractedGaussianOrbital& orbital_2;
		const Vec3d& position_1;
		const Vec3d& position_2;
	};
	std::vector<ComputeItem> computeItems;
	computeItems.reserve(numberOfContractedGaussians * numberOfContractedGaussians); // Note, we are over reserving here, but thats fine to waste the extra memory

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
								computeItems.emplace_back(i, j, orbital1, orbital2, atom1.position, atom2.position);
							}

							++j;
						}
					}
				}
				++i;
			}
		}
	}

	std::for_each(std::execution::par_unseq, computeItems.begin(), computeItems.end(),
		[&kineticEnergyMatrix](const ComputeItem& item)
		{
			double result = KineticEnergyOfTwoOrbitals(item.orbital_1, item.position_1, item.orbital_2, item.position_2);
			kineticEnergyMatrix(item.index_i, item.index_j) = result;
			kineticEnergyMatrix(item.index_j, item.index_i) = result;
		});

	return kineticEnergyMatrix;
}

// Overlap & Kinetic Functions
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

	struct ComputeItem
	{
		unsigned int index_i;
		unsigned int index_j;
		const ContractedGaussianOrbital& orbital_1;
		const ContractedGaussianOrbital& orbital_2;
		const Vec3d& position_1;
		const Vec3d& position_2;
	};
	std::vector<ComputeItem> computeItems;
	computeItems.reserve(numberOfContractedGaussians * numberOfContractedGaussians); // Note, we are over reserving here, but thats fine to waste the extra memory

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
								computeItems.emplace_back(i, j, orbital1, orbital2, atom1.position, atom2.position);
							}

							++j;
						}
					}
				}
				++i;
			}
		}
	}

	std::for_each(std::execution::par_unseq, computeItems.begin(), computeItems.end(), 
		[&overlapMatrix, &kineticEnergyMatrix](const ComputeItem& item) 
		{
			auto [overlap, kinetic] = OverlapAndKineticEnergyOfTwoOrbitals(item.orbital_1, item.position_1, item.orbital_2, item.position_2);

			overlapMatrix(item.index_i, item.index_j) = overlap;
			overlapMatrix(item.index_j, item.index_i) = overlap;

			kineticEnergyMatrix(item.index_i, item.index_j) = kinetic;
			kineticEnergyMatrix(item.index_j, item.index_i) = kinetic;
		});
}

// Nuclear-Electron Attraction Functions
static double Abscissa(int n, int i) noexcept
{
	const double val = i * std::numbers::pi / (1. + n);
	const double sinVal = std::sin(val);
	const double cosVal = std::cos(val);
	return (n + 1. - 2. * i) / (n + 1.) + 2. / std::numbers::pi * (1. + 2. / 3. * sinVal * sinVal) * cosVal * sinVal;
}
static double Omega(int n, int i) noexcept
{
	const double sinVal = std::sin(i * std::numbers::pi / (n + 1.));
	return 16. / (3 * (1. + n)) * sinVal * sinVal * sinVal * sinVal;
}
static double ChebyshevIntegral(double eps, int M, const std::function<double(double)>& F) noexcept
{
	double c0 = std::cos(std::numbers::pi / 6.);
	double s0 = 0.5;
	double c1, s1, q, p, chp, c, s, x;

	double err = 10.;

	int n = 3;
	c1 = s0;
	s1 = c0;

	const double res = Abscissa(2, 1);

	q = (F(res) + F(-res)) * Omega(2, 1);
	p = F(0);

	chp = q + p;

	int j = 0;

	while (err > eps && 2. * n * (1. - j) + j * 4. * n / 3. - 1. <= M)
	{
		j = 1 - j;

		c1 = j * c1 + (1. - j) * c0;
		s1 = j * s1 + (1. - j) * s0;
		c0 = j * c0 + (1. - j) * std::sqrt((1. + c0) / 2.);
		s0 = j * s0 + (1. - j) * s0 / (c0 + c0);

		c = c0;
		s = s0;

		for (int i = 1; i < n; i += 2)
		{
			x = 1. + 2. / (3. * std::numbers::pi) * s * c * (3. + 2. * s * s) - ((double)i) / n;

			const int div = (int)((i + 2. * j) / 3.);
			if (3 * div >= i + j)
				chp += (F(-x) + F(x)) * s * s * s * s;

			x = s;
			s = s * c1 + c * s1;
			c = c * c1 - x * s1;
		}

		n *= 1 + j;
		p += (1. - j) * (chp - q);

		err = 16. * std::abs((1. - j) * (q - 3. * p / 2.) + j * (chp - q - q)) / (3. * n);
		q = (1. - j) * q + j * chp;
	}

	chp = 16. * q / (3. * n);

	return chp;
}
static double NuclearElectronAttractionOfTwoPrimitiveGaussians(double alpha1, double alpha2, const Vec3d& position1, const Vec3d& position2, const Vec3d& nuclearCenter, const QuantumNumbers& angularMomentum1, const QuantumNumbers& angularMomentum2) noexcept
{
	double oneDividedByAlpha1PlusAlpha2 = 1 / (alpha1 + alpha2);

	// NOTE: Because the functions within the each matrix capture a reference to the matrix itself, we cannot easily
	//       separate this into a separate function, otherwise we would have to worry about dangling references. For 
	//       simplicity, we create each matrix on the stack here, and that guarantees the lifetime of each function we create
	StackMatrix2D_Functions eta_x(angularMomentum1.l + angularMomentum2.l + 1, angularMomentum2.l + 1);
	StackMatrix2D_Functions eta_y(angularMomentum1.m + angularMomentum2.m + 1, angularMomentum2.m + 1);
	StackMatrix2D_Functions eta_z(angularMomentum1.n + angularMomentum2.n + 1, angularMomentum2.n + 1);

	eta_x(0, 0) = [](double) -> double { return 1.0; };
	eta_y(0, 0) = [](double) -> double { return 1.0; };
	eta_z(0, 0) = [](double) -> double { return 1.0; };

	// eta_x
	if (eta_x.RowCount() > 1)
	{
		double p = alpha1 * position1.x + alpha2 * position2.x;
		double pTimesOneDividedByAlpha1PlusAlpha2 = p * oneDividedByAlpha1PlusAlpha2;

		// Initial conditions
		eta_x(1, 0) = 
			[a = -1 * (position1.x - pTimesOneDividedByAlpha1PlusAlpha2), 
			 b = -1 * (pTimesOneDividedByAlpha1PlusAlpha2 - nuclearCenter.x)](double t) -> double
			{
				return a + b * t * t;
			};

		// Recurrence Index
		for (unsigned int row = 2; row < eta_x.RowCount(); ++row)
		{
			eta_x(row, 0) =
				[&eta_x, row, oneDividedByAlpha1PlusAlpha2,
				a = -1 * (position1.x - pTimesOneDividedByAlpha1PlusAlpha2),
				b = -1 * (pTimesOneDividedByAlpha1PlusAlpha2 - nuclearCenter.x),
				c = (row - 1) * 0.5 * oneDividedByAlpha1PlusAlpha2](double t) -> double
				{
					return a + b * t * t * eta_x(row - 1, 0)(t) + c * (1 - t * t) * eta_x(row - 2, 0)(t);
				};
		}

		// Transfer Equation
		//     The recurrence relation requires filling out each column before moving on to the next one.
		//     Therefore, the outer loop must be the loop over the columns
		for (unsigned int col = 1; col < eta_x.ColumnCount(); ++col)
		{
			for (unsigned int row = 0; row + col <= angularMomentum1.l + angularMomentum2.l; ++row)
			{
				eta_x(row, col) = [&eta_x, row, col, a = position1.x - position2.x](double t) -> double
					{
						return eta_x(row + 1, col - 1)(t) + a * eta_x(row, col - 1)(t);
					};
			}
		}
	}

	// eta_y
	if (eta_y.RowCount() > 1)
	{
		double p = alpha1 * position1.y + alpha2 * position2.y;
		double pTimesOneDividedByAlpha1PlusAlpha2 = p * oneDividedByAlpha1PlusAlpha2; 

		// Initial conditions
		eta_y(1, 0) = 
			[a = -1 * (position1.y - pTimesOneDividedByAlpha1PlusAlpha2), 
			 b = -1 * (pTimesOneDividedByAlpha1PlusAlpha2 - nuclearCenter.y)](double t) -> double
			{
				return a + b * t * t;
			};

		// Recurrence Index
		for (unsigned int row = 2; row < eta_y.RowCount(); ++row)
		{
			eta_y(row, 0) =
				[&eta_y, row, oneDividedByAlpha1PlusAlpha2,
				a = -1 * (position1.y - pTimesOneDividedByAlpha1PlusAlpha2),
				b = -1 * (pTimesOneDividedByAlpha1PlusAlpha2 - nuclearCenter.y),
				c = (row - 1) * 0.5 * oneDividedByAlpha1PlusAlpha2](double t) -> double
				{
					return a + b * t * t * eta_y(row - 1, 0)(t) + c * (1 - t * t) * eta_y(row - 2, 0)(t);
				};
		}

		// Transfer Equation
		//     The recurrence relation requires filling out each column before moving on to the next one.
		//     Therefore, the outer loop must be the loop over the columns
		for (unsigned int col = 1; col < eta_y.ColumnCount(); ++col)
		{
			for (unsigned int row = 0; row + col <= angularMomentum1.m + angularMomentum2.m; ++row)
			{
				eta_y(row, col) = [&eta_y, row, col, a = position1.y - position2.y](double t) -> double
					{
						return eta_y(row + 1, col - 1)(t) + a * eta_y(row, col - 1)(t);
					};
			}
		}
	}

	// eta_z
	if (eta_z.RowCount() > 1)
	{
		double p = alpha1 * position1.z + alpha2 * position2.z;
		double pTimesOneDividedByAlpha1PlusAlpha2 = p * oneDividedByAlpha1PlusAlpha2;

		// Initial conditions
		eta_z(1, 0) = 
			[a = -1 * (position1.z - pTimesOneDividedByAlpha1PlusAlpha2), 
			 b = -1 * (pTimesOneDividedByAlpha1PlusAlpha2 - nuclearCenter.z)](double t) -> double
			{
				return a + b * t * t;
			};

		// Recurrence Index
		for (unsigned int row = 2; row < eta_z.RowCount(); ++row)
		{
			eta_z(row, 0) =
				[&eta_z, row, oneDividedByAlpha1PlusAlpha2,
				a = -1 * (position1.z - pTimesOneDividedByAlpha1PlusAlpha2),
				b = -1 * (pTimesOneDividedByAlpha1PlusAlpha2 - nuclearCenter.z),
				c = (row - 1) * 0.5 * oneDividedByAlpha1PlusAlpha2](double t) -> double
				{
					return a + b * t * t * eta_z(row - 1, 0)(t) + c * (1 - t * t) * eta_z(row - 2, 0)(t);
				};
		}

		// Transfer Equation
		//     The recurrence relation requires filling out each column before moving on to the next one.
		//     Therefore, the outer loop must be the loop over the columns
		for (unsigned int col = 1; col < eta_z.ColumnCount(); ++col)
		{
			for (unsigned int row = 0; row + col <= angularMomentum1.n + angularMomentum2.n; ++row)
			{
				eta_z(row, col) = [&eta_z, row, col, a = position1.z - position2.z](double t) -> double
					{
						return eta_z(row + 1, col - 1)(t) + a * eta_z(row, col - 1)(t);
					};
			}
		}
	}


	Vec3d diff = position2 - position1; 
	double factor = std::exp(-1 * alpha1 * alpha2 * oneDividedByAlpha1PlusAlpha2 * diff.Dot(diff)) * 2 * std::numbers::pi * oneDividedByAlpha1PlusAlpha2;

	Vec3d v = (oneDividedByAlpha1PlusAlpha2 * (alpha1 * position1 + alpha2 * position2)) - nuclearCenter;
	auto func = [&eta_x, &eta_y, &eta_z,
				 &angularMomentum1, &angularMomentum2,
				 a = alpha1 + alpha2,
	             b = v.Dot(v)](double t) -> double
		{
			t = (t + 1) / 2; 
			return 0.5 * std::exp(-1 * a * t * t * b) *
				eta_x(angularMomentum1.l, angularMomentum2.l)(t) *
				eta_y(angularMomentum1.m, angularMomentum2.m)(t) *
				eta_z(angularMomentum1.n, angularMomentum2.n)(t);			
		};

	return factor * ChebyshevIntegral(1E-10, 50000, func);
}
static double NuclearElectronAttractionEnergyOfTwoOrbitals(const Vec3d& nuclearCenter, const ContractedGaussianOrbital& orbital1, const Vec3d& position1, const ContractedGaussianOrbital& orbital2, const Vec3d& position2) noexcept
{
	double res = 0;

	for (auto& gaussian1 : orbital1.gaussianOrbitals)
	{
		for (auto& gaussian2 : orbital2.gaussianOrbitals)
		{
			double nuclear = NuclearElectronAttractionOfTwoPrimitiveGaussians(gaussian1.alpha, gaussian2.alpha, position1, position2, nuclearCenter, orbital1.angularMomentum, orbital2.angularMomentum);
			res += gaussian1.normalizationFactor * gaussian2.normalizationFactor * gaussian1.coefficient * gaussian2.coefficient * nuclear;
		}
	}
	return res;
}
MatrixXd NuclearElectronAttractionEnergyMatrix(std::span<Atom> atoms, const Basis& basis) noexcept
{
	MatrixXd nuclearMatrix;

	assert(atoms.size() > 1);
	unsigned int numberOfContractedGaussians = 0;
	for (const Atom& atom : atoms)
		numberOfContractedGaussians += basis.GetAtom(atom.type).NumberOfContractedGaussians();

	nuclearMatrix = MatrixXd::Zero(numberOfContractedGaussians, numberOfContractedGaussians);

	struct ComputeItem
	{
		unsigned int index_i;
		unsigned int index_j;
		const ContractedGaussianOrbital& orbital_1;
		const ContractedGaussianOrbital& orbital_2;
		const Vec3d& position_1;
		const Vec3d& position_2;
	};
	std::vector<ComputeItem> computeItems;
	computeItems.reserve(numberOfContractedGaussians * numberOfContractedGaussians); // Note, we are over reserving here, but thats fine to waste the extra memory

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
							if (j >= i)
							{
								computeItems.emplace_back(i, j, orbital1, orbital2, atom1.position, atom2.position);
							}

							++j;
						}
					}
				}
				++i;
			}
		}
	}

	std::for_each(std::execution::par_unseq, computeItems.begin(), computeItems.end(), 
		[&nuclearMatrix, &atoms](const ComputeItem& item)
		{
			for (const auto& nuclearCenterAtom : atoms)
				nuclearMatrix(item.index_i, item.index_j) -= nuclearCenterAtom.Z() * NuclearElectronAttractionEnergyOfTwoOrbitals(nuclearCenterAtom.position, item.orbital_1, item.position_1, item.orbital_2, item.position_2);

			if (item.index_i != item.index_j)
				nuclearMatrix(item.index_j, item.index_i) = nuclearMatrix(item.index_i, item.index_j);
		});

	return nuclearMatrix;
}



}