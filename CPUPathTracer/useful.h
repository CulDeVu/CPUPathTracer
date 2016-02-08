#pragma once

#include <random>

const float maxFloat = 2147483647.0f;

std::random_device rd;
std::mt19937 mt(rd());
std::uniform_real_distribution<float> unitDist(0.f, 1.f);
float nrand()
{
	return unitDist(mt);
}

float min(float a, float b)
{
	if (a < b)
		return a;
	return b;
}
float max(float a, float b)
{
	if (a > b)
		return a;
	return b;
}
float min(float a, float b, float c)
{
	float x = min(a, b);
	if (c < x)
		return c;
	return x;
}
float max(float a, float b, float c)
{
	float x = max(a, b);
	if (c > x)
		return c;
	return x;
}

vec3 getTangent(vec3 norm)
{
	vec3 tangent;
	vec3 c1 = cross(norm, vec3(0, 0, 1));
	vec3 c2 = cross(norm, vec3(0, 1, 0));
	if (dot(c1, c1) > dot(c2, c2))
		tangent = c1;
	else
		tangent = c2;
	return tangent;
}

vec3 normalSpaceToWorldSpace(vec3 in, vec3 norm)
{
	vec3 tangent = getTangent(norm);
	vec3 bitangent = cross(norm, tangent);

	vec3 castRay = normalize(tangent*in.x + bitangent*in.z + norm*in.y);
	return castRay;
}
vec3 worldSpaceToNormalSpace(vec3 in, vec3 norm)
{
	vec3 tangent = getTangent(norm);
	vec3 bitangent = cross(norm, tangent);

	vec3 ret = vec3(dot(in, tangent), dot(in, norm), dot(in, bitangent));
	return ret;
}