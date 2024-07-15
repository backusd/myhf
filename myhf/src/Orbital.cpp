#include "pch.h"
#include "Orbital.h"


namespace myhf
{
	double PrimitiveGaussian::operator()(const Vec3d& r, const Vec3d& center, const QuantumNumbers& angularMomentum) const noexcept
	{
		const Vec3d R = r - center;
		return coeffProdNorm * std::pow(R.x, angularMomentum.l) * std::pow(R.y, angularMomentum.m) * std::pow(R.z, angularMomentum.n) * std::exp(-alpha * R * R); 
	}

	Vec3d PrimitiveGaussian::GetGradient(const Vec3d& r, const Vec3d& center, const QuantumNumbers& angularMomentum) const noexcept
	{
		const Vec3d R = r - center;

		const double pRX = std::pow(R.x, angularMomentum.l);
		const double pRY = std::pow(R.y, angularMomentum.m);
		const double pRZ = std::pow(R.z, angularMomentum.n);
		const double powProd = pRX * pRY * pRZ;
		const double expaR2 = std::exp(-alpha * R * R);

		const Vec3d eD = -alpha * 2. * R;

		double valX = eD.x;
		if (angularMomentum.l > 0)
			valX += static_cast<double>(angularMomentum.l) / R.x;

		double valY = eD.y;
		if (angularMomentum.m > 0)
			valY += static_cast<double>(angularMomentum.m) / R.y;

		double valZ = eD.z;
		if (angularMomentum.n > 0)
			valZ += static_cast<double>(angularMomentum.n) / R.z;

		return coeffProdNorm * expaR2 * powProd * Vec3d(valX, valY, valZ);

		// Numerical, for tests:
		// 
		// double h = 0.0001;
		// 
		// double val = operator()(r);
		// 
		// double dx = (operator()(Vector3D<double>(r.X + h, r.Y, r.Z)) - val) / h;
		// double dy = (operator()(Vector3D<double>(r.X, r.Y + h, r.Z)) - val) / h;
		// double dz = (operator()(Vector3D<double>(r.X, r.Y, r.Z + h)) - val) / h;
		// 
		// return Vector3D<double>(dx, dy, dz);
		// 
	}

	double PrimitiveGaussian::GetLaplacian(const Vec3d& r, const Vec3d& center, const QuantumNumbers& angularMomentum) const noexcept
	{
		const Vec3d R = r - center;

		const double pRX = std::pow(R.x, angularMomentum.l);
		const double pRY = std::pow(R.y, angularMomentum.m);
		const double pRZ = std::pow(R.z, angularMomentum.n);
		const double powProd = pRX * pRY * pRZ;
		const double expaR2 = std::exp(-alpha * R * R);

		const double twoalpha = 2. * alpha;
		const Vec3d eD = -twoalpha * R;

		double valX = eD.x;
		if (angularMomentum.l > 0)
			valX += static_cast<double>(angularMomentum.l) / R.x;
		valX *= -R.x * twoalpha;

		double valY = eD.y;
		if (angularMomentum.m > 0)
			valY += static_cast<double>(angularMomentum.m) / R.y;
		valY *= -R.y * twoalpha;

		double valZ = eD.z;
		if (angularMomentum.n > 0)
			valZ += static_cast<double>(angularMomentum.n) / R.z;
		valZ *= -R.z * twoalpha;

		double valX2 = -twoalpha * (angularMomentum.l + 1.);
		if (angularMomentum.l > 1)
			valX2 += static_cast<double>(angularMomentum.l) * (angularMomentum.l - 1.) / (R.x * R.x);

		double valY2 = -twoalpha * (angularMomentum.m + 1.);
		if (angularMomentum.m > 1)
			valY2 += static_cast<double>(angularMomentum.m) * (angularMomentum.m - 1.) / (R.y * R.y);

		double valZ2 = -twoalpha * (angularMomentum.n + 1.);
		if (angularMomentum.n > 1)
			valZ2 += static_cast<double>(angularMomentum.n) * (angularMomentum.n - 1.) / (R.z * R.z);

		return (valX + valX2 + valY + valY2 + valZ + valZ2) * coeffProdNorm * powProd * expaR2;

		// 
		// double h = 0.0001;
		// double h2 = h * h;
		// 
		// double val = operator()(r);
		// 
		// double d2x = operator()(Vector3D<double>(r.X + h, r.Y, r.Z));
		// double d2y = operator()(Vector3D<double>(r.X, r.Y + h, r.Z));
		// double d2z = operator()(Vector3D<double>(r.X, r.Y, r.Z + h));
		// 
		// d2x += operator()(Vector3D<double>(r.X - h, r.Y, r.Z)) - 2. * val;
		// d2y += operator()(Vector3D<double>(r.X, r.Y - h, r.Z)) - 2. * val;
		// d2z += operator()(Vector3D<double>(r.X, r.Y, r.Z - h)) - 2. * val;
		// 
		// d2x /= h2;
		// d2y /= h2;
		// d2z /= h2;
		// 
		// return d2x + d2y + d2z;
		// 
	}








double ContractedGaussianOrbital::operator()(const Vec3d& r, const Vec3d& center) const noexcept
{
	double res = 0;

	for (const auto& orbital : gaussianOrbitals)
		res += orbital(r, center, angularMomentum);

	return res;
}

Vec3d ContractedGaussianOrbital::GetGradient(const Vec3d& r, const Vec3d& center) const noexcept
{
	Vec3d res(0, 0, 0);

	for (const auto& orbital : gaussianOrbitals)
		res += orbital.GetGradient(r, center, angularMomentum);

	return res;
}

double ContractedGaussianOrbital::GetLaplacian(const Vec3d& r, const Vec3d& center) const noexcept
{
	double res = 0;

	for (const auto& orbital : gaussianOrbitals)
		res += orbital.GetLaplacian(r, center, angularMomentum);

	return res;
}
}
