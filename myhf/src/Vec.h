#pragma once

namespace myhf
{
template<typename T>
struct Vec3
{
	T x;
	T y;
	T z;

	constexpr Vec3() noexcept : x(0), y(0), z(0) {}
	constexpr Vec3(T all) noexcept : x(all), y(all), z(all) {}
	constexpr Vec3(T _x, T _y, T _z) noexcept : x(_x), y(_y), z(_z) {}
	constexpr Vec3(const Vec3<T>& rhs) noexcept : x(rhs.x), y(rhs.y), z(rhs.z) {}
	constexpr Vec3(Vec3<T>&& rhs) noexcept : x(rhs.x), y(rhs.y), z(rhs.z) {}
	constexpr Vec3<T>& operator=(const Vec3<T>& rhs) noexcept
	{
		x = rhs.x;
		y = rhs.y;
		z = rhs.z;
		return *this;
	}
	constexpr Vec3<T>& operator=(Vec3<T>&& rhs) noexcept
	{
		x = rhs.x;
		y = rhs.y;
		z = rhs.z;
		return *this;
	}

	constexpr const Vec3<T>& operator+() const noexcept { return *this; }
	constexpr Vec3<T> operator-() const noexcept { return { -x, -y, -z }; }

	constexpr Vec3<T> operator+(const Vec3<T>& rhs) const noexcept { return { x + rhs.x, y + rhs.y, z + rhs.z }; }
	constexpr Vec3<T> operator-(const Vec3<T>& rhs) const noexcept { return { x - rhs.x, y - rhs.y, z - rhs.z }; }
	constexpr T operator*(const Vec3<T>& rhs) const noexcept { return x * rhs.x + y * rhs.y + z * rhs.z; }
	constexpr Vec3<T> operator%(const Vec3<T>& rhs) const noexcept { return { y * rhs.z - z * rhs.y, z * rhs.x - x * rhs.z, x * rhs.y - y * rhs.x }; }

	constexpr T Dot(const Vec3<T>& rhs) const noexcept { return this->operator*(rhs); }

	constexpr Vec3<T> operator*(T s) const noexcept { return { x * s, y * s, z * s }; }
	constexpr Vec3<T> operator/(T s) const noexcept { return { x / s, y / s, z / s }; }

	constexpr Vec3<T>& operator+=(const Vec3<T>& rhs) noexcept { x += rhs.x; y += rhs.y; z += rhs.z; return *this; }
	constexpr Vec3<T>& operator-=(const Vec3<T>& rhs) noexcept { x -= rhs.x; y -= rhs.y; z -= rhs.z; return *this; }

	constexpr Vec3<T>& operator*=(T s) noexcept { x *= s; y *= s; z *= s; return *this; }
	constexpr Vec3<T>& operator/=(T s) noexcept { assert(s != 0); x /= s; y /= s; z /= s; return *this; }

	constexpr double Length() const noexcept { return std::sqrt(this->Dot(*this)); }
	constexpr Vec3<T> Normalize() const noexcept
	{
		const T len = Length();
		return len == 0 ? *this : *this / len;
	}

	template<typename A>
	Vec3<T> RotateAround(const Vec3<T>& other, A angle) const noexcept
	{
		Vec3<T> result{};
		const Vec3<T> u = other.Normalize();

		const T sinAngle = std::sin(angle);
		const T cosAngle = std::cos(angle);
		const T oneMinus = 1. - cosAngle;

		result.x = (cosAngle + oneMinus * u.x * u.x) * x + (oneMinus * u.x * u.y - u.z * sinAngle) * y + (oneMinus * u.x * u.z + u.y * sinAngle) * z;
		result.y = (oneMinus * u.x * u.y + u.z * sinAngle) * x + (cosAngle + oneMinus * u.y * u.y) * y + (oneMinus * u.y * u.z - u.x * sinAngle) * z;
		result.z = (oneMinus * u.x * u.z - u.y * sinAngle) * x + (oneMinus * u.y * u.z + u.x * sinAngle) * y + (cosAngle + oneMinus * u.z * u.z) * z;

		return result;
	}
	template<typename A>
	Vec3<T> RotateTowards(const Vec3<T>& other, A angle) const noexcept
	{
		const Vec3<T> a = this->operator%(other);
		return RotateAround(a, angle);
	}

	constexpr double GetTheta() const noexcept
	{
		double cosTheta = z / Length();
		if (std::isnan(cosTheta) || std::isinf(cosTheta) || cosTheta > 1. || cosTheta < -1)
			cosTheta = (cosTheta < 0 ? -1 : 1);

		return std::acos(cosTheta);
	}
	constexpr double GetPhi() const noexcept { return std::atan2(y, x); }

	friend std::ostream& operator<<(std::ostream& os, const Vec3<T>& v)
	{
		os << '[' << v.x << ", " << v.y << ", " << v.z << ']';
		return os;
	}
};

template<typename T> Vec3<T> operator*(T o, const Vec3<T>& t) { return t * o; }
template<typename T> bool operator==(const Vec3<T>& f, const Vec3<T>& t) { return f.X == t.X && f.Y == t.Y && f.Z == t.Z; }

using Vec3d = Vec3<double>;
using Vec3f = Vec3<float>;

}

template <typename T> 
struct std::formatter<myhf::Vec3<T>>
{
	constexpr auto parse(std::format_parse_context& ctx) noexcept
	{
		return ctx.begin();
	}
	auto format(const myhf::Vec3<T>& vec, std::format_context& ctx) const
	{
		return std::format_to(ctx.out(), "({}, {}, {})", vec.x, vec.y, vec.z);
	}
};