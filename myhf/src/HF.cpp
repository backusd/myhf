#include "pch.h"
#include "HF.h"

#include <execution>

using Eigen::MatrixXd;

namespace myhf
{
static double OverlapOfTwoPrimitiveGaussiansAlongAxis(double oneDividedByAlpha1PlusAlpha2, double alpha1, double alpha2, double position1, double position2, unsigned int angularMomentum1, unsigned int angularMomentum2) noexcept
{
	// Helper class for keeping track of the recurrence relation values with data that resides on the stack
	class StackMatrix2D
	{
	public:
		// Example:
		//    Supposed we need to return the value 'matrix[a, b]'. To compute that value, you must
		//    have already computed 'matrix[a+1, b-1]' - this is the value located down 1 and left 1 spot
		//	  in the matrix. If you continuing this logic, you will eventually get to the first column of the
		//	  matrix, which will tell us how many rows we need. Notice that along that diagonal in the matrix, 
		//    the row index plus the column index have a constant value (a + b). To adjust for the off-by-one 
		//    error, we also add one.
		StackMatrix2D(unsigned int angularMomentum1, unsigned int angularMomentum2) noexcept :
			rowCount(angularMomentum1 + angularMomentum2 + 1), colCount(angularMomentum2 + 1), values(nullptr)
		{
			values = static_cast<double*>(alloca(sizeof(double) * rowCount * colCount));
		}
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

	if (angularMomentum1 == 0 && angularMomentum2 == 0)
		return 1.0;

	StackMatrix2D matrix(angularMomentum1, angularMomentum2); 

	// X: Initial conditions
	matrix(0, 0) = 1.0;
	double startingValue = -(position1 - ((alpha1 * position1 + alpha2 * position2) * oneDividedByAlpha1PlusAlpha2));
	matrix(1, 0) = startingValue;

	// X: Recurrence Index
	for (unsigned int iii = 2; iii < matrix.RowCount(); ++iii)
		matrix(iii, 0) = startingValue * matrix(iii - 1, 0) + ((iii - 1) * 0.5 * oneDividedByAlpha1PlusAlpha2) * matrix(iii - 2, 0);

	// X: Transfer Equation
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

	int i = 0;
	for (const auto& atom1 : atoms)
	{
		for (const auto& shell1 : basis.GetAtom(atom1.type).shells)
		{
			for (const auto& orbital1 : shell1.basisFunctions)
			{
				int j = 0;

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

}