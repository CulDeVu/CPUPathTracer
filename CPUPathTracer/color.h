#pragma once

#include <iostream>
#include <math.h>

#include "vec3.h"
#include "useful.h"

using namespace std;

struct color
{
	color() : r(0), g(0), b(0) {}
	color(float r, float g, float b) : r(r), g(g), b(b) {}

	float r, g, b;

	color normalized()
	{
		if (r < 0 || g < 0 || b < 0 ||
			r > maxFloat || g > maxFloat || b > maxFloat)
			cout << "BAD COLOR VALUE: " << r << " " << g << " " << b << endl;

		//return color(r / (r + 1), g / (g + 1), b / (b + 1));
		return color(min(1, r), min(1, g), min(1, b));
	}
	color denormalized()
	{
		if (r > 255 || g > 255 || b > 255)
			cout << "BAD COLOR VALUE: DENORMALIZED" << endl;

		return color(r / (1 - r), g / (1 - g), b / (1 - b));
	}
	color clamp()
	{
		return color(min(1, r), min(1, g), min(1, b));
	}
	float greyscale()
	{
		return r * 0.21f + g * 0.72f + b * 0.07f;
	}

	color mul(color other)
	{
		return color(r * other.r, g * other.g, b * other.b);
	}

	color operator+(color other) const
	{
		return color(r + other.r, g + other.g, b + other.b);
	}
	color operator-(color other) const
	{
		return color(r - other.r, g - other.g, b - other.b);
	}
	color operator*(float f) const
	{
		return color(r * f, g * f, b * f);
	}
	color operator/(float f) const
	{
		return color(r / f, g / f, b / f);
	}
	bool operator==(const color& rhs)
	{
		return (r == rhs.r && g == rhs.g && b == rhs.b);
	}
	bool operator!=(const color& rhs)
	{
		return !(*this == rhs);
	}

	color fromVec3(vec3 other)
	{
		color ret;
		other = normalize(other);
		//ret.r = (other.x + 1) / 2;
		//ret.g = (other.y + 1) / 2;
		ret.b = (other.z + 1) / 2;
		return ret;
	}
};

float dot(color a, color b)
{
	return a.r*b.r + a.g*b.g + a.b*b.b;
}
float length(color a)
{
	return sqrt(dot(a, a));
}
color clamp(color c, float a, float b)
{
	color ret;
	ret.r = clamp(c.r, a, b);
	ret.g = clamp(c.g, a, b);
	ret.b = clamp(c.b, a, b);
	return ret;
}

// useful constants
color SKY_BLACK = color(0, 0, 0);
color SKY_BRIGHT = color(12, 12, 12);