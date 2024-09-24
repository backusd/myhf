#include "pch.h"
#include "HF.h"

#include <execution>
#include <print>

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

//	std::println("==================================================");
//	double factor = std::exp(-1 * alpha1 * alpha2 * oneDividedByAlpha1PlusAlpha2 * diff.Dot(diff)) * std::pow(std::numbers::pi * oneDividedByAlpha1PlusAlpha2, 1.5);
//	std::println("Factor: {0}", factor);
//	std::println("k_x   : {0}", k_x_value);
//	std::println("k_y   : {0}", k_y_value);
//	std::println("k_z   : {0}", k_z_value);
//	std::println("s_x   : {0}", s_x_value);
//	std::println("s_y   : {0}", s_y_value);
//	std::println("s_z   : {0}", s_z_value);
	
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
MatrixXd KineticEnergyMatrix(std::span<Atom> atoms, const Basis& basis) noexcept
{
	MatrixXd kineticMatrix;

	assert(atoms.size() > 1);
	unsigned int numberOfContractedGaussians = 0;
	for (const Atom& atom : atoms)
		numberOfContractedGaussians += basis.GetAtom(atom.type).NumberOfContractedGaussians();

	kineticMatrix = MatrixXd::Zero(numberOfContractedGaussians, numberOfContractedGaussians);

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
							if (j >= i)
							{
								kineticMatrix(i, j) = KineticEnergyOfTwoOrbitals(orbital1, atom1.position, orbital2, atom2.position);
								
								if (j > i)
									kineticMatrix(j, i) = kineticMatrix(i, j);
							}

							++j;
						}
					}
				}
				++i;
			}
		}
	}

	return kineticMatrix;
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
							if (j >= i)
							{
								auto [overlap, kinetic] = OverlapAndKineticEnergyOfTwoOrbitals(orbital1, atom1.position, orbital2, atom2.position);
								overlapMatrix(i, j) = overlap;
								kineticEnergyMatrix(i, j) = kinetic;

								if (j > i)
								{
									overlapMatrix(j, i) = overlapMatrix(i, j);
									kineticEnergyMatrix(j, i) = kineticEnergyMatrix(i, j);
								}
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