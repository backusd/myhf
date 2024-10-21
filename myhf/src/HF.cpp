#include "pch.h"
#include "HF.h"
#include "Profiling.h"

#include <execution>
#include <print>
#include <iostream>
#include <thread>
#include <ranges>

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
static double OverlapOfTwoPrimitiveGaussiansAlongAxis_2(double oneDividedByAlpha1PlusAlpha2, double alpha1, double alpha2, double position1, double position2, unsigned int angularMomentum1, unsigned int angularMomentum2) noexcept
{
	if (angularMomentum1 == 0 && angularMomentum2 == 0)
		return 1.0;

	unsigned int angularMomentumSum = angularMomentum1 + angularMomentum2;

	// Supposed we need to return the value 'matrix[a, b]'. To compute that value, you must
	// have already computed 'matrix[a+1, b-1]' - this is the value located down 1 and left 1 spot
	// in the matrix. If you continuing this logic, you will eventually get to the first column of the
	// matrix, which will tell us how many rows we need. Notice that along that diagonal in the matrix, 
	// the row index plus the column index have a constant value (a + b). To adjust for the off-by-one 
	// error, we also add one.
	StackMatrix2D matrix(angularMomentumSum + 1, angularMomentum2 + 1);

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
	double positionDelta = position1 - position2;
	for (unsigned int col = 1; col < matrix.ColumnCount(); ++col)
	{
		// The conditional here is a little weird. Basically, we have set up a matrix to store the recurrence values in. 
		// However, the bottom right half of values will never get used (and cannot be calculated without expanding the matrix). 
		// For example, consider trying to compute matrix[0, 2]. The matrix will look something like this:
		//		v1	v4	v6
		//		v2	v5	-
		//		v3	-	-
		// We only need to compute up to v6 and don't need the '-' values and we wouldn't even be able to compute them unless we added more rows.

		for (unsigned int row = 0; row + col <= angularMomentumSum; ++row)
			matrix(row, col) = matrix(row + 1, col - 1) + positionDelta * matrix(row, col - 1);
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
static double OverlapOfTwoPrimitiveGaussians_2(double alpha1, double alpha2, const Vec3d& position1, const Vec3d& position2, const QuantumNumbers& angularMomentum1, const QuantumNumbers& angularMomentum2) noexcept
{
	assert(std::max(angularMomentum1.l, angularMomentum2.l) + 1 < 4);
	assert(std::max(angularMomentum1.m, angularMomentum2.m) + 1 < 4);
	assert(std::max(angularMomentum1.n, angularMomentum2.n) + 1 < 4);

	const double oneDividedByAlpha1PlusAlpha2 = 1 / (alpha1 + alpha2);
	double s_x_value = OverlapOfTwoPrimitiveGaussiansAlongAxis_2(oneDividedByAlpha1PlusAlpha2, alpha1, alpha2, position1.x, position2.x, angularMomentum1.l, angularMomentum2.l);
	double s_y_value = OverlapOfTwoPrimitiveGaussiansAlongAxis_2(oneDividedByAlpha1PlusAlpha2, alpha1, alpha2, position1.y, position2.y, angularMomentum1.m, angularMomentum2.m);
	double s_z_value = OverlapOfTwoPrimitiveGaussiansAlongAxis_2(oneDividedByAlpha1PlusAlpha2, alpha1, alpha2, position1.z, position2.z, angularMomentum1.n, angularMomentum2.n);

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
static double OverlapOfTwoOrbitals_2(const ContractedGaussianOrbital& orbital1, const Vec3d& position1, const ContractedGaussianOrbital& orbital2, const Vec3d& position2) noexcept
{
	double res = 0;

	for (auto& gaussian1 : orbital1.gaussianOrbitals)
	{
		for (auto& gaussian2 : orbital2.gaussianOrbitals)
		{
			double overlap = OverlapOfTwoPrimitiveGaussians_2(gaussian1.alpha, gaussian2.alpha, position1, position2, orbital1.angularMomentum, orbital2.angularMomentum);
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
MatrixXd OverlapMatrix_2(std::span<Atom> atoms, const Basis& basis) noexcept
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
			double result = OverlapOfTwoOrbitals_2(item.orbital_1, item.position_1, item.orbital_2, item.position_2);
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
struct QuantumNumbersWithMetadata
{
	constexpr QuantumNumbersWithMetadata(unsigned int l, unsigned int m, unsigned int n) noexcept :
		quantumNumbers(l, m, n),
		totalIndexOfPrev(0),
		totalIndexOfPrevPrev(0),
		offset(0),
		N(0)
	{
		canonicalIndex = quantumNumbers.GetCanonicalIndex();
		totalCanonicalIndex = quantumNumbers.GetTotalCanonicalIndex();

		// Compute the total canonical index of each of the next quantum numbers
		nextQuantumNumbers[0] = QuantumNumbers(l + 1, m, n).GetTotalCanonicalIndex();
		nextQuantumNumbers[1] = QuantumNumbers(l, m + 1, n).GetTotalCanonicalIndex();
		nextQuantumNumbers[2] = QuantumNumbers(l, m, n + 1).GetTotalCanonicalIndex();

		const unsigned int maxIndex = quantumNumbers.MaxComponentVal();
		if (maxIndex > 0)
		{
			QuantumNumbers prev = quantumNumbers;
			QuantumNumbers next = quantumNumbers;
			if (quantumNumbers.l == maxIndex)
			{
				// An offset of 0, indicates the 'l' quantum number is the one that changes when computing
				// the prev value. This also has the impact of selecting the .x value from a Vec3d, when used later on
				offset = 0;

				--prev.l;
				totalIndexOfPrev = prev.GetTotalCanonicalIndex();

				N = prev.l;

				if (prev.l > 0)
				{
					--prev.l;
					totalIndexOfPrevPrev = prev.GetTotalCanonicalIndex();
				}
			}
			else if (quantumNumbers.m == maxIndex)
			{
				// An offset of 1, indicates the 'm' quantum number is the one that changes when computing
				// the prev value. This also has the impact of selecting the .y value from a Vec3d, when used later on
				offset = 1;

				--prev.m;
				totalIndexOfPrev = prev.GetTotalCanonicalIndex();

				N = prev.m;

				if (prev.m > 0)
				{
					--prev.m;
					totalIndexOfPrevPrev = prev.GetTotalCanonicalIndex();
				}
			}
			else
			{
				// An offset of 2, indicates the 'n' quantum number is the one that changes when computing
				// the prev value. This also has the impact of selecting the .z value from a Vec3d, when used later on
				offset = 2;

				--prev.n;
				totalIndexOfPrev = prev.GetTotalCanonicalIndex();

				N = prev.n;

				if (prev.n > 0)
				{
					--prev.n;
					totalIndexOfPrevPrev = prev.GetTotalCanonicalIndex();
				}
			}
		}
	}
	QuantumNumbers quantumNumbers;
	unsigned int canonicalIndex;
	unsigned int totalCanonicalIndex;

	unsigned int totalIndexOfPrev;
	unsigned int totalIndexOfPrevPrev;

	std::array<unsigned int, 3> nextQuantumNumbers;

	unsigned int offset;
	unsigned int N; // This value is used in the vertical recurrence relationship
};
static constexpr std::array<QuantumNumbersWithMetadata, 35> allQuantumNumbersInCanonicalOrder = {
		QuantumNumbersWithMetadata{0, 0, 0}, // s
		QuantumNumbersWithMetadata{1, 0, 0}, // p
		QuantumNumbersWithMetadata{0, 1, 0},
		QuantumNumbersWithMetadata{0, 0, 1},
		QuantumNumbersWithMetadata{2, 0, 0}, // d
		QuantumNumbersWithMetadata{1, 1, 0},
		QuantumNumbersWithMetadata{1, 0, 1},
		QuantumNumbersWithMetadata{0, 2, 0},
		QuantumNumbersWithMetadata{0, 1, 1},
		QuantumNumbersWithMetadata{0, 0, 2},
		QuantumNumbersWithMetadata{3, 0, 0}, // f
		QuantumNumbersWithMetadata{2, 1, 0},
		QuantumNumbersWithMetadata{2, 0, 1},
		QuantumNumbersWithMetadata{1, 2, 0},
		QuantumNumbersWithMetadata{1, 1, 1},
		QuantumNumbersWithMetadata{1, 0, 2},
		QuantumNumbersWithMetadata{0, 3, 0},
		QuantumNumbersWithMetadata{0, 2, 1},
		QuantumNumbersWithMetadata{0, 1, 2},
		QuantumNumbersWithMetadata{0, 0, 3},
		QuantumNumbersWithMetadata{4, 0, 0}, // g
		QuantumNumbersWithMetadata{3, 1, 0},
		QuantumNumbersWithMetadata{3, 0, 1},
		QuantumNumbersWithMetadata{2, 2, 0},
		QuantumNumbersWithMetadata{2, 1, 1},
		QuantumNumbersWithMetadata{2, 0, 2},
		QuantumNumbersWithMetadata{1, 3, 0},
		QuantumNumbersWithMetadata{1, 2, 1},
		QuantumNumbersWithMetadata{1, 1, 2},
		QuantumNumbersWithMetadata{1, 0, 3},
		QuantumNumbersWithMetadata{0, 4, 0},
		QuantumNumbersWithMetadata{0, 3, 1},
		QuantumNumbersWithMetadata{0, 2, 2},
		QuantumNumbersWithMetadata{0, 1, 3},
		QuantumNumbersWithMetadata{0, 0, 4},
};
class NuclearMatrix
{
public:
	NuclearMatrix() noexcept = default;
	// delete copies/moves because we do not need to support that use case
	NuclearMatrix(const NuclearMatrix&) = delete;
	NuclearMatrix(NuclearMatrix&&) = delete;
	NuclearMatrix& operator=(const NuclearMatrix&) = delete;
	NuclearMatrix& operator=(NuclearMatrix&&) = delete;

	constexpr double& operator()(unsigned int row, unsigned int col) noexcept
	{
		assert(row >= 0 && row < m_rowCount);
		assert(col >= 0 && col < m_colCount);
		return values[row * m_colCount + col];
	}
	constexpr unsigned int RowCount() const noexcept { return m_rowCount; }
	constexpr unsigned int ColumnCount() const noexcept { return m_colCount; }

private:
	// The number of rows required is the max quantum number's total canonical index + 1. So for now, let's at least support up to f orbitals
	static constexpr unsigned int m_rowCount = QuantumNumbers(0, 0, 4).GetTotalCanonicalIndex() + 1;

	// The number of columns required is the angular momentum of both orbitals + 1. Again, lets support up to f orbitals
	static constexpr unsigned int m_colCount = 2 * QuantumNumbers(0, 0, 4).AngularMomentum() + 1;

	std::array<double, m_rowCount * m_colCount> values = {}; // The braces are required so the array is 0-initialized
};

static double BoysFunction(double m, double x)
{
	PROFILE_SCOPE("BoysFunction");

	if (0.0 == m && 0.0 == x)
		return 1.0;

	double res = 0;

	const double t = 1E-10;
	if (std::abs(x) < t)
		x = t;

	if (x > 160)
	{
		return static_cast<double>(
			MathUtils::DoubleFactorial(static_cast<unsigned int>(2 * m - 1)) /
			std::pow(2, m + 1) *
			std::sqrt(std::numbers::pi / std::pow(x, 2 * m + 1))
			);
	}

	MathUtils::IncompleteGamma(m + 1. / 2., x, res);

	return res / (2. * std::pow(x, m + 1. / 2.));
}
static inline auto GenerateBoysFunctions(int maxM, double T)
{
	PROFILE_SCOPE("GenerateBoysFunctions");

	// The maximum number of boys functions required is the angular momentum of both orbitals + 1. So here we support up to f orbitals
	std::array<double, 2 * QuantumNumbers(0, 0, 4).AngularMomentum() + 1> results{};

	results[maxM] = BoysFunction(maxM, T);

	for (int m = maxM - 1; m >= 0; --m)
		results[m] = (2. * T * results[m + 1ULL] + std::exp(-T)) / (2. * m + 1.);

	return results;
}
static void VerticalRecursion(NuclearMatrix& mat, double alpha, const Vec3d& Rp, const Vec3d& center1, const Vec3d& difN, unsigned int maxL)
{
	PROFILE_SCOPE("VerticalRecursion");

	if (maxL == 0)
		return;

	double difScalar, difNScalar;

	const unsigned int size = maxL + 1;
	const auto difRp = Rp - center1;

	// Iterate over the quantum numbers with total angular momentum = 1. These are just the 3 p quantum numbers
	// We separate these out because the calculation is slightly different given that they do not have a prevprev quantum number
	for (unsigned int qnIndex = 1; qnIndex <= 3; ++qnIndex)
	{
		const QuantumNumbersWithMetadata& currentQNMd = allQuantumNumbersInCanonicalOrder[qnIndex];
		const QuantumNumbers& currentQN = currentQNMd.quantumNumbers;
		unsigned int curIndex = currentQNMd.totalCanonicalIndex;

		const QuantumNumbersWithMetadata& prevQNMd = allQuantumNumbersInCanonicalOrder[currentQNMd.totalIndexOfPrev];
		const QuantumNumbers& prevQN = prevQNMd.quantumNumbers;
		unsigned int prevIndex = prevQNMd.totalCanonicalIndex;

		unsigned int stoppingPoint = size - currentQN.AngularMomentum();
		for (unsigned int m = 0; m < stoppingPoint; ++m)
		{
			difScalar = difRp[currentQNMd.offset];
			difNScalar = difN[currentQNMd.offset];
			mat(curIndex, m) = difScalar * mat(prevIndex, m) + difNScalar * mat(prevIndex, m + 1ULL);
		}
	}

	// Now continue a similar calculation, but start with the d quantum numbers with total angular momentum 2 and continue
	// until you've covered all sets of quantum numbers with total angular momentum <= maxL
	for (unsigned int qnIndex = 4; allQuantumNumbersInCanonicalOrder[qnIndex].quantumNumbers.AngularMomentum() < size; ++qnIndex)
	{
		const QuantumNumbersWithMetadata& currentQNMd = allQuantumNumbersInCanonicalOrder[qnIndex];
		const QuantumNumbers& currentQN = currentQNMd.quantumNumbers;
		unsigned int curIndex = currentQNMd.totalCanonicalIndex;

		const QuantumNumbersWithMetadata& prevQNMd = allQuantumNumbersInCanonicalOrder[currentQNMd.totalIndexOfPrev];
		const QuantumNumbers& prevQN = prevQNMd.quantumNumbers;
		unsigned int prevIndex = prevQNMd.totalCanonicalIndex;

		const QuantumNumbersWithMetadata& prevPrevQNMd = allQuantumNumbersInCanonicalOrder[currentQNMd.totalIndexOfPrevPrev];
		const QuantumNumbers& prevPrevQN = prevPrevQNMd.quantumNumbers;
		unsigned int prevPrevIndex = prevPrevQNMd.totalCanonicalIndex;

		unsigned int stoppingPoint = size - currentQN.AngularMomentum();
		for (unsigned int m = 0; m < stoppingPoint; ++m)
		{
			difScalar = difRp[currentQNMd.offset];
			difNScalar = difN[currentQNMd.offset];

			// NOTE: currentQNMd.N will be 0 for certain quantum numbers. In order for it to not be 0, the current QuantumNumber needs to 
			//       have a value >= 2 for at least one of its components. When this is not the case, it is not possible to compute the
			//       prevprev QuantumNumber. So in theory, this could probably be optimized to not attempt to execute that part of the 
			//       computation when it knows N will be 0, but for simplicity, we lleave it like this for now
			mat(curIndex, m) = difScalar * mat(prevIndex, m) + difNScalar * mat(prevIndex, m + 1ULL) +
				currentQNMd.N / (2. * alpha) * (mat(prevPrevIndex, m) - mat(prevPrevIndex, m + 1ULL));
		}
	}
}
static void GetNuclearVertical(NuclearMatrix& mat, const Vec3d& nuclearCenter, const PrimitiveGaussian& gaussian1, const Vec3d& position1, const PrimitiveGaussian& gaussian2, const Vec3d& position2, unsigned int maxL1, unsigned int maxL2)
{
	PROFILE_SCOPE("GetNuclearVertical");

	const double alpha = gaussian1.alpha + gaussian2.alpha;
	const Vec3d Rp = (gaussian1.alpha * position1 + gaussian2.alpha * position2) / alpha;
	const Vec3d difN = nuclearCenter - Rp;
	const Vec3d dif = position1 - position2;

	const unsigned int maxL = maxL1 + maxL2;
	const unsigned int size = maxL + 1;

	// auxiliary integrals
	auto boys = GenerateBoysFunctions(maxL, alpha * (difN * difN));

	QuantumNumbers maxQN(0, 0, maxL);

	const double factor = 2. * std::numbers::pi / alpha * std::exp(-gaussian1.alpha * gaussian2.alpha / alpha * dif * dif);

	for (unsigned int i = 0; i < size; ++i)
		mat(0, i) = factor * boys[i];

	VerticalRecursion(mat, alpha, Rp, position1, difN, maxL);
}
static void HorizontalRecursion(NuclearMatrix& mat, const Vec3d& dif, unsigned int maxL1, unsigned int maxL2)
{
	PROFILE_SCOPE("HorizontalRecursion");

	unsigned int maxL = maxL1 + maxL2;

	for (unsigned int qnIndex_y = 1; allQuantumNumbersInCanonicalOrder[qnIndex_y].quantumNumbers.AngularMomentum() <= maxL2; ++qnIndex_y)
	{
		for (unsigned int qnIndex_x = 0; allQuantumNumbersInCanonicalOrder[qnIndex_x].quantumNumbers.AngularMomentum() < maxL; ++qnIndex_x)
		{
			unsigned int nextIndex_x = allQuantumNumbersInCanonicalOrder[qnIndex_x].nextQuantumNumbers[allQuantumNumbersInCanonicalOrder[qnIndex_y].offset];	
			unsigned int prevIndex_y = allQuantumNumbersInCanonicalOrder[qnIndex_y].totalIndexOfPrev;
			double difScalar = dif[allQuantumNumbersInCanonicalOrder[qnIndex_y].offset];

			mat(qnIndex_x, qnIndex_y) = mat(nextIndex_x, prevIndex_y) + difScalar * mat(qnIndex_x, prevIndex_y);
		}

		// Decrement maxL by 1 whenever we go up in overall angular momentum (s -> p -> d -> f)
		maxL -= allQuantumNumbersInCanonicalOrder[qnIndex_y + 1].quantumNumbers.AngularMomentum() - allQuantumNumbersInCanonicalOrder[qnIndex_y].quantumNumbers.AngularMomentum();
	}
}
static double NuclearElectronAttractionEnergyOfTwoOrbitals(const Vec3d& nuclearCenter, const ContractedGaussianOrbital& orbital1, const Vec3d& position1, const ContractedGaussianOrbital& orbital2, const Vec3d& position2) noexcept
{
	PROFILE_SCOPE("NuclearElectronAttractionEnergyOfTwoOrbitals");

	// The algorithm assumes that the angular momentum of orbital1 is <= that of orbital2. So we need to do a swap here if that isn't the case
	const ContractedGaussianOrbital& orbital_1 = (orbital1.angularMomentum >= orbital2.angularMomentum) ? orbital1 : orbital2;
	const ContractedGaussianOrbital& orbital_2 = (orbital1.angularMomentum >= orbital2.angularMomentum) ? orbital2 : orbital1;

	const Vec3d& position_1 = (orbital1.angularMomentum >= orbital2.angularMomentum) ? position1 : position2;
	const Vec3d& position_2 = (orbital1.angularMomentum >= orbital2.angularMomentum) ? position2 : position1;

	QuantumNumbers maxQN(0, 0, orbital_1.angularMomentum.AngularMomentum() + orbital_2.angularMomentum.AngularMomentum());
	unsigned int horizNuclearRowCount = maxQN.GetTotalCanonicalIndex() + 1ULL;

	NuclearMatrix horizNuclear;
	NuclearMatrix nuclear;

	for (const auto& gaussian_1 : orbital_1.gaussianOrbitals)
	{
		for (const auto& gaussian_2 : orbital_2.gaussianOrbitals)
		{
			GetNuclearVertical(nuclear, nuclearCenter, gaussian_1, position_1, gaussian_2, position_2, orbital_1.angularMomentum.AngularMomentum(), orbital_2.angularMomentum.AngularMomentum());

			double factor = gaussian_1.normalizationFactor * gaussian_2.normalizationFactor * gaussian_1.coefficient * gaussian_2.coefficient;

			for (unsigned int row = 0; row < horizNuclearRowCount; ++row)
				horizNuclear(row, 0) += factor * nuclear(row, 0);
		}
	}

	HorizontalRecursion(horizNuclear, position_1 - position_2, orbital_1.angularMomentum.AngularMomentum(), orbital_2.angularMomentum.AngularMomentum());
	return horizNuclear(orbital_1.angularMomentum.GetTotalCanonicalIndex(), orbital_2.angularMomentum.GetTotalCanonicalIndex());
}
MatrixXd NuclearElectronAttractionEnergyMatrix(std::span<Atom> atoms, const Basis& basis) noexcept
{
	MatrixXd nuclearMatrix;

	PROFILE_SCOPE("NuclearElectronAttractionEnergyMatrix");

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

	{
		PROFILE_SCOPE("Make compute items");

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
	}
	std::for_each(computeItems.begin(), computeItems.end(),
		[&nuclearMatrix, &atoms](const ComputeItem& item)
		{
			PROFILE_SCOPE("Compute Item");

			for (const auto& nuclearCenterAtom : atoms)
			{
				PROFILE_SCOPE("Iterating on atom");
				nuclearMatrix(item.index_i, item.index_j) -= nuclearCenterAtom.Z() * NuclearElectronAttractionEnergyOfTwoOrbitals(nuclearCenterAtom.position, item.orbital_1, item.position_1, item.orbital_2, item.position_2);
			}
			if (item.index_i != item.index_j)
				nuclearMatrix(item.index_j, item.index_i) = nuclearMatrix(item.index_i, item.index_j);
		});

	return nuclearMatrix;
}
MatrixXd NuclearElectronAttractionEnergyMatrix_Par(std::span<Atom> atoms, const Basis& basis) noexcept
{
	MatrixXd nuclearMatrix;

	PROFILE_SCOPE("NuclearElectronAttractionEnergyMatrix");

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

	{
		PROFILE_SCOPE("Make compute items");

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
	}
	std::for_each(std::execution::par_unseq, computeItems.begin(), computeItems.end(),
		[&nuclearMatrix, &atoms](const ComputeItem& item)
		{
			PROFILE_SCOPE("Compute Item");

			for (const auto& nuclearCenterAtom : atoms)
			{
				PROFILE_SCOPE("Iterating on atom");
				nuclearMatrix(item.index_i, item.index_j) -= nuclearCenterAtom.Z() * NuclearElectronAttractionEnergyOfTwoOrbitals(nuclearCenterAtom.position, item.orbital_1, item.position_1, item.orbital_2, item.position_2);
			}
			if (item.index_i != item.index_j)
				nuclearMatrix(item.index_j, item.index_i) = nuclearMatrix(item.index_i, item.index_j);
		});

	return nuclearMatrix;
}







}