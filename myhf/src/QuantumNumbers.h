#pragma once
#include "pch.h"

namespace myhf
{
struct QuantumNumbers
{
	unsigned int l, m, n;

	constexpr QuantumNumbers(unsigned int L = 0, unsigned int M = 0, unsigned int N = 0) noexcept : 
		l(L), m(M), n(N)
	{}

	constexpr char AtomicOrbital() const noexcept
	{
		switch (AngularMomentum())
		{
		case 0: return 's';
		case 1: return 'p';
		case 2: return 'd';
		case 3: return 'f';
		case 4: return 'g';
		case 5: return 'h';
		}

		return -1;
	}

	constexpr unsigned int AngularMomentum() const noexcept
	{
		return l + m + n;
	}

	// I don't like that an int is used to differentiate which value to return
	// Hold off on implementing this until it gets used
//	inline unsigned int N(unsigned int ind)
//	{
//		if (0 == ind) return l;
//		else if (1 == ind) return m;
//
//		assert(2 == ind);
//
//		return n;
//	}

	constexpr unsigned int MaxComponentVal() const noexcept { return std::max(std::max(l, m), n); }
	constexpr unsigned int MinComponentVal() const noexcept { return std::min(std::min(l, m), n); }

	// we're using Cartesian Gaussians, there are (L + 1) * (L + 2) of them in a shell with angular momentum L
	// don't confuse it with the usual 2L+1 number
	static constexpr unsigned int NumOrbitals(unsigned int L) noexcept
	{
		return (L + 1) * (L + 2) / 2;
	}
	constexpr unsigned int NumOrbitals() const noexcept
	{
		return NumOrbitals(AngularMomentum());
	}

	constexpr unsigned int GetCanonicalIndex() const noexcept
	{
		// Example of canonical indices for L = 3 (L = l + m + n)
		//		Quantum Numbers		Canonical Index
		//		(3, 0, 0)			0
		//		(2, 1, 0)			1
		//		(2, 0, 1)			2	
		//		(1, 2, 0)			3
		//		(1, 1, 1)			4
		//		(1, 0, 2)			5
		//		(0, 3, 0)			6
		//		(0, 2, 1)			7
		//		(0, 1, 2)			8
		//		(0, 0, 3)			9
		return (m + n + 1) * (m + n + 2) / 2 - m - 1;
	}

	constexpr unsigned int GetTotalCanonicalIndex() const noexcept
	{
		const unsigned int L = AngularMomentum();

		const unsigned int val = L * (L * (L + 3) + 2) / 6;

		return val + GetCanonicalIndex();
	}

	constexpr void CanonicalIncrement() noexcept
	{
		const unsigned int L = AngularMomentum();

		if (n == L)
		{
			l = L + 1;
			m = n = 0;
		}
		else if (m > 0)
		{
			--m;
			++n;
		}
		else if (l > 0)
		{
			--l;
			m = L - l;
			n = 0;
		}
	}



	// I don't like using a call operator that isn't 100% clear what it is doing without having
	// to look at its implementation
//	inline operator unsigned int() const { return AngularMomentum(); }

	constexpr QuantumNumbers& operator++() noexcept
	{
		CanonicalIncrement();
		return *this;
	}

	constexpr QuantumNumbers operator++(int) noexcept
	{
		QuantumNumbers temp(*this);
		operator++();
		return temp;
	}
};

constexpr bool operator<(const QuantumNumbers& lhs, const QuantumNumbers& rhs) noexcept { return lhs.AngularMomentum() < rhs.AngularMomentum(); }
constexpr bool operator<(const QuantumNumbers& lhs, unsigned int rhs) noexcept { return lhs.AngularMomentum() < rhs; }
constexpr bool operator<(const QuantumNumbers& lhs, int rhs) noexcept { return lhs.AngularMomentum() < static_cast<unsigned int>(rhs); }
constexpr bool operator<(unsigned int lhs, const QuantumNumbers& rhs) noexcept { return lhs < rhs.AngularMomentum(); }


constexpr bool operator>(const QuantumNumbers& lhs, const QuantumNumbers& rhs) noexcept { return lhs.AngularMomentum() > rhs.AngularMomentum(); }
constexpr bool operator>(const QuantumNumbers& lhs, unsigned int rhs) noexcept { return lhs.AngularMomentum() > rhs; }
constexpr bool operator>(const QuantumNumbers& lhs, int rhs) noexcept { return lhs.AngularMomentum() > static_cast<unsigned int>(rhs); }
constexpr bool operator>(unsigned int lhs, const QuantumNumbers& rhs) noexcept { return lhs > rhs.AngularMomentum(); }


constexpr bool operator<=(const QuantumNumbers& lhs, const QuantumNumbers& rhs) noexcept { return lhs.AngularMomentum() <= rhs.AngularMomentum(); }
constexpr bool operator<=(const QuantumNumbers& lhs, unsigned int rhs) noexcept { return lhs.AngularMomentum() <= rhs; }
constexpr bool operator<=(const QuantumNumbers& lhs, int rhs) noexcept { return lhs.AngularMomentum() <= static_cast<unsigned int>(rhs); }
constexpr bool operator<=(unsigned int lhs, const QuantumNumbers& rhs) noexcept { return lhs <= rhs.AngularMomentum(); }


constexpr bool operator>=(const QuantumNumbers& lhs, const QuantumNumbers& rhs) noexcept { return lhs.AngularMomentum() >= rhs.AngularMomentum(); }
constexpr bool operator>=(const QuantumNumbers& lhs, unsigned int rhs) noexcept { return lhs.AngularMomentum() >= rhs; }
constexpr bool operator>=(const QuantumNumbers& lhs, int rhs) noexcept { return lhs.AngularMomentum() >= static_cast<unsigned int>(rhs); }
constexpr bool operator>=(unsigned int lhs, const QuantumNumbers& rhs) noexcept { return lhs >= rhs.AngularMomentum(); }
}