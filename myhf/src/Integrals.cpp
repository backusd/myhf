#include "pch.h"
#include "Integrals.h"



namespace myhf
{

//double GetOverlap(const Atom& atom1, const PrimitiveGaussian& gaussian1, const Atom& atom2, const PrimitiveGaussian& gaussian2)
//{
////	const std::tuple<unsigned int, unsigned int, double, double > params(gaussian1.shellID, gaussian2.shellID, gaussian1.alpha, gaussian2.alpha);
////
////	//auto it = overlapIntegralsMap.find(params);
////	//if (overlapIntegralsMap.end() != it) return it->second.getOverlap(gaussian1.angularMomentum, gaussian2.angularMomentum);
////
////	auto it = momentIntegralsMap.find(params);
////	if (momentIntegralsMap.end() != it) return it->second.getOverlap(gaussian1.angularMomentum, gaussian2.angularMomentum);
////
////	// unfortunately it's not yet calculated
////	//GaussianOverlap overlap;
////	//auto result = overlapIntegralsMap.insert(std::make_pair(params, overlap));
////	GaussianMoment moment;
////	auto result = momentIntegralsMap.insert(std::make_pair(params, moment));
//
//	QuantumNumbers maxQN1 = atom1.GetMaxQN(gaussian1.alpha);
//	QuantumNumbers maxQN2 = atom2.GetMaxQN(gaussian2.alpha);
//
////	if (extendForKinetic)
////	{
////		// calculating the kinetic integrals needs +1 quantum numbers for overlap integrals
////		++maxQN1.l;
////		++maxQN1.m;
////		++maxQN1.n;
////
////		++maxQN2.n;
////		++maxQN2.l;
////		++maxQN2.m;
////	}
//
//	// calculate the integrals and that's about it
//
//	GaussianMoment moment;
//	moment.Reset(gaussian1.alpha, gaussian2.alpha, atom1.position, atom2.position, maxQN1, maxQN2);
//
//	return moment.getOverlap(gaussian1.angularMomentum, gaussian2.angularMomentum);
//}


GaussianMoment::GaussianMoment(double alpha1, double alpha2, const Vec3d& center1, const Vec3d& center2, const Vec3d& center3, const QuantumNumbers& maxQN1, const QuantumNumbers& maxQN2)
	: factor(0)
{
	Reset(alpha1, alpha2, center1, center2, center3, maxQN1, maxQN2);
}


void GaussianMoment::Reset(double alpha1, double alpha2, const Vec3d& center1, const Vec3d& center2, const Vec3d& center3, const QuantumNumbers& maxQN1, const QuantumNumbers& maxQN2)
{
	matrixX = Eigen::MatrixXd::Zero(3ULL + maxQN1.l + maxQN2.l, maxQN2.l + 2ULL);
	matrixY = Eigen::MatrixXd::Zero(3ULL + maxQN1.m + maxQN2.m, maxQN2.m + 2ULL);
	matrixZ = Eigen::MatrixXd::Zero(3ULL + maxQN1.n + maxQN2.n, maxQN2.n + 2ULL);

	matrixX1 = Eigen::MatrixXd::Zero(2ULL + maxQN1.l + maxQN2.l, maxQN2.l + 2ULL);
	matrixY1 = Eigen::MatrixXd::Zero(2ULL + maxQN1.m + maxQN2.m, maxQN2.m + 2ULL);
	matrixZ1 = Eigen::MatrixXd::Zero(2ULL + maxQN1.n + maxQN2.n, maxQN2.n + 2ULL);

	CalculateMoment(matrixX, matrixX1, alpha1, alpha2, center1.x, center2.x, center3.x, maxQN1.l, maxQN2.l);
	CalculateMoment(matrixY, matrixY1, alpha1, alpha2, center1.y, center2.y, center3.y, maxQN1.m, maxQN2.m);
	CalculateMoment(matrixZ, matrixZ1, alpha1, alpha2, center1.z, center2.z, center3.z, maxQN1.n, maxQN2.n);

	const Vec3d dif = center1 - center2;
	factor = exp(-alpha1 * alpha2 / (alpha1 + alpha2) * dif * dif) * pow(std::numbers::pi / (alpha1 + alpha2), 3. / 2.);
}

double GaussianMoment::getMoment(const QuantumNumbers& QN1, const QuantumNumbers& QN2, bool momentX, bool momentY, bool momentZ) const
{
	return factor * (momentX ? matrixX1(QN1.l, QN2.l) : matrixX(QN1.l, QN2.l)) * (momentY ? matrixY1(QN1.m, QN2.m) : matrixY(QN1.m, QN2.m)) * (momentZ ? matrixZ1(QN1.n, QN2.n) : matrixZ(QN1.n, QN2.n));
}

double GaussianMoment::getOverlap(const QuantumNumbers& QN1, const QuantumNumbers& QN2) const
{
	return factor * matrixX(QN1.l, QN2.l) * matrixY(QN1.m, QN2.m) * matrixZ(QN1.n, QN2.n); // or getMoment(QN1, QN2, false, false, false) but this is slightly faster
}

void GaussianMoment::CalculateMoment(Eigen::MatrixXd& matrix, Eigen::MatrixXd& matrix1, double alpha1, double alpha2, double center1, double center2, double center3, unsigned int maxQN1, unsigned int maxQN2)
{
	const double alpha = alpha1 + alpha2;
	const double productCenter = (alpha1 * center1 + alpha2 * center2) / alpha;
	const double dif = center1 - center2;
	const double dif1 = center1 - center3;
	const double difCenter = productCenter - center1;

	matrix(0, 0) = 1;
	matrix(1, 0) = difCenter;

	// recurrence index
	unsigned int limit = maxQN1 + maxQN2 + 2;

	// vertical recurrence relation - the same as for overlap
	for (unsigned int i = 2; i <= limit; ++i)
		matrix(i, 0) = difCenter * matrix(i - 1, 0) + (i - 1.) / (2. * alpha) * matrix(i - 2, 0);


	--limit;

	// for first column of matrix1 as well
	for (unsigned int i = 0; i < limit; ++i)
		matrix1(i, 0) = matrix(i + 1ULL, 0) + dif1 * matrix(i, 0);


	// transfer equation - the horizontal recurrence relation
	const unsigned int horzLimit = maxQN2 + 1;
	for (unsigned int j = 1; j <= horzLimit; ++j, --limit)
	{
		const unsigned int jminus1 = j - 1ULL;

		// *** this is correct, is exactly as for overlap ********************************

		for (unsigned int i = 0; i <= limit; ++i)
			matrix(i, j) = matrix(i + 1ULL, jminus1) + dif * matrix(i, jminus1);

		// *******************************************************************************

		for (unsigned int i = 0; i < limit; ++i)
			matrix1(i, j) = matrix(i + 1ULL, j) + dif1 * matrix(i, j);
	}

	matrix = matrix.block(0, 0, maxQN1 + 1ULL, maxQN2 + 1ULL).eval();
	matrix1 = matrix1.block(0, 0, maxQN1 + 1ULL, maxQN2 + 1ULL).eval();
}
}