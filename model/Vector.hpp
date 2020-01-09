// No Copyright. Vladislav Aleinik 2019
#ifndef POTTS_MODEL_VECTOR_HPP_INCLUDED
#define POTTS_MODEL_VECTOR_HPP_INCLUDED

#include <cmath>
#include <limits>
#include <immintrin.h>

// 3D-vector with basic arithmetic
union Vector
{
	struct
	{
		double x, y, z;
	};
	__m256d reg256;
	int64_t conds[4];

	Vector() = default;

	Vector(double newX, double newY, double newZ);

	inline Vector& operator+=(const Vector& vec);
	inline Vector& operator-=(const Vector& vec);

	inline Vector operator+(const Vector& vec) const;
	inline Vector operator-(const Vector& vec) const;
	
	inline Vector& operator*=(double k);
	inline Vector& operator/=(double k);

	inline Vector operator*(double k) const;
	inline Vector operator/(double k) const;

	double lenSqr() const;
	double length() const;
	
	inline void setLength(double newLen);

	inline double scalar(const Vector& v) const;
};

Vector::Vector(double newX, double newY, double newZ) : 
	x (newX), y (newY), z (newZ)
{}

inline Vector& Vector::operator+=(const Vector& vec)
{
	reg256 = _mm256_add_pd(reg256, vec.reg256);

	return *this;
}

inline Vector Vector::operator+(const Vector& vec) const
{
	Vector toReturn;
	toReturn.reg256 = _mm256_add_pd(reg256, vec.reg256);
	
	return toReturn;
}

inline Vector& Vector::operator-=(const Vector& vec)
{
	reg256 = _mm256_sub_pd(reg256, vec.reg256);

	return *this;
}

inline Vector Vector::operator-(const Vector& vec) const 
{
	Vector toReturn;
	toReturn.reg256 = _mm256_sub_pd(reg256, vec.reg256);

	return toReturn;
}

inline Vector& Vector::operator*=(double k)
{
	reg256 = _mm256_mul_pd(reg256, _mm256_set1_pd(k));

	return *this;
}

inline Vector Vector::operator*(double k) const
{
	Vector toReturn;
	toReturn.reg256 = _mm256_mul_pd(reg256, _mm256_set1_pd(k));

	return toReturn;
}

inline Vector& Vector::operator/=(double k)
{
	reg256 = _mm256_mul_pd(reg256, _mm256_set1_pd(1/k));

	return *this;
}

inline Vector Vector::operator/(double k) const
{
	Vector toReturn;
	toReturn.reg256 = _mm256_mul_pd(reg256, _mm256_set1_pd(1/k));

	return toReturn;
}

double Vector::lenSqr() const
{
	Vector len;
	len.reg256 = _mm256_mul_pd(reg256, reg256);

	return len.x + len.y + len.z;
}

double Vector::length() const
{
	return std::sqrt(lenSqr());
}

inline void Vector::setLength(double newLen)
{
	double curLen = length();

	if (std::abs(curLen) < 10 * std::numeric_limits<double>::epsilon()) return;

	operator*=(newLen / curLen);
}

inline double Vector::scalar(const Vector& v) const
{
	return v.x*x + v.y*y + v.z*z;
}

#endif  // POTTS_MODEL_VECTOR_HPP_INCLUDED