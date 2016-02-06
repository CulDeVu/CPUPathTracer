#pragma once

#include <math.h>

using namespace std;

struct vec3
{
	vec3() : vec3(0, 0, 0) {}
	vec3(float a, float b, float c) : x(a), y(b), z(c) {}

	float x, y, z;

	float length() { return sqrt(x*x + y*y + z*z); }
	float lengthsqr() { return x*x + y*y + z*z; }
	vec3 normalized() {
		float len = length();
		vec3 me = vec3(x / len, y / len, z / len);
		return me;
	}
	vec3 operator+(vec3 other) const
	{
		vec3 t = vec3(other.x + x, other.y + y, other.z + z);
		return t;
	}
	vec3 operator-(vec3 other) const
	{
		vec3 t = vec3(x - other.x, y - other.y, z - other.z);
		return t;
	}
	vec3 operator*(float f) const
	{
		vec3 t = vec3(x * f, y * f, z * f);
		return t;
	}
	vec3 operator/(float f) const
	{
		vec3 t = vec3(x / f, y / f, z / f);
		return t;
	}
	bool operator==(const vec3& rhs)
	{
		return (x == rhs.x) && (y == rhs.y) && (z == rhs.z);
	}
	bool operator!=(const vec3& rhs)
	{
		return !(*this == rhs);
	}
};

float dot(vec3 v1, vec3 v2)
{
	return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}
vec3 cross(vec3 v1, vec3 v2)
{
	return vec3(v1.y*v2.z - v1.z*v2.y, v1.z*v2.x - v1.x*v2.z, v1.x*v2.y - v1.y*v2.x);
}