#pragma once

#include "color.h"
#include "texture.h"
#include "vec3.h"

const float PI = 3.14159f;

vec3 randCosineWeightedRay(vec3 norm)
{
	float rx = 1, rz = 1;
	while (rx*rx + rz*rz >= 1)
	{
		rx = 2 * nrand() - 1.0f;
		rz = 2 * nrand() - 1.0f;
	}
	float ry = sqrt(1 - rx*rx - rz*rz);

	vec3 tangent = getTangent(norm);
	vec3 bitangent = cross(norm, tangent);

	vec3 castRay = (tangent*rx + bitangent*rz + norm*ry).normalized();
	return castRay;
}

class material
{
public:
	material()
		: Albedo(color(0.4f, 0.4f, 0.4f)), 
		Fresnel(color(0.0f, 0.0f, 0.0f)), 
		m(1),
		Emmision(color(0, 0, 0)) {}
	material(color a, color f, float r)
		: Albedo(a), Fresnel(f), m(r) {}

	color Albedo;
	color Fresnel;
	float m;
	color Emmision;
	int albedoTex;
};

float Fresnel(vec3 l, vec3 m, float rf0)
{
	//float LDotM = dot(l, m);
	//return rf0 + (1 - rf0) * pow(1 - LDotM, 5);
	//return rf0;

	float n = (1 + sqrt(rf0)) / max(1 - sqrt(rf0), 0.001);
	float c = dot(l, m);
	float g_temp = n * n - 1 + c * c;
	if (g_temp < 0)
		return 1;
	float g = sqrt(g_temp);

	float first = pow((g - c) / (g + c), 2) / 2;
	float second = 1 + pow((c * (g + c) - 1) / (c * (g - c) + 1), 2);
	float result = first * second;

	if (result != result)
	{
		printf("Fresnel became invalid: n: %f, c: %f, g: %f\n", n, c, g);
	}

	return result;
}
float BeckmanD(vec3 m, vec3 norm, float alpha)
{
	float NDotM = dot(m, norm);
	
	float positivity = (NDotM > 0) ? 1.0f : 0.0f;

	float theta_m = acos(NDotM);
	float coef = -pow(tan(theta_m), 2) / (alpha * alpha);
	float denom = pow(alpha, 2) * pow(NDotM, 4);
	if (denom < 0.001)
		denom = 0.001;
	float total = positivity * max(0.001, exp(coef)) / (PI * denom);
	if (total > 500)
		printf("asdfasdfe");
	return total;
}
float SmithG1Approx(vec3 v, vec3 m, vec3 norm, float alpha)
{	
	float VDotM = dot(v, m);
	float VDotN = dot(v, norm);

	float theta_v = acos(VDotN);
	float a = 1 / (alpha * tan(theta_v));

	float positivity = (VDotM / VDotN > 0) ? 1.0f : 0.0f;

	if (a < 1.6)
	{
		return (3.535 * a + 2.181 * a * a) / (1 + 2.276 * a + 2.577 * a * a);
	}
	else
	{
		return 1;
	}
}

vec3 sampleMicrofacet(material mat, vec3 norm, float u, float v, vec3 iDir, float* pdf)
{
	float roughness = mat.m;
	float alpha = roughness;
	//float alpha_t = (1.2 - 0.2 * sqrt(dot(iDir, norm))) * roughness;
	//float alpha_t = 1.2 * roughness;
	float alpha_t = alpha;
	if (alpha_t != alpha_t ||
		alpha_t < 0.001)
	{
		printf("%f\n", (double)dot(iDir, norm));
	}

	float d1 = nrand();
	float d2 = nrand();

	float theta_m = atan(sqrt(-alpha_t * alpha_t * log(1 - d1)));
	float phi_m = 2 * PI * d2;

	float mY = cos(theta_m);
	float mX = cos(phi_m) * sin(theta_m);
	float mZ = sin(phi_m) * sin(theta_m);

	vec3 tangent = getTangent(norm);
	vec3 bitangent = cross(norm, tangent);
	vec3 m = (tangent * mX + norm * mY + bitangent * mZ).normalized();

	float NDotM = dot(norm, m);
	float IDotM = dot(iDir, m);

	float pdf_m = BeckmanD(m, norm, alpha_t) * NDotM;

	vec3 oDir = m * IDotM * 2.0f - iDir;
	float ODotH = dot(oDir, m);
	if (ODotH < 0.001)
		ODotH = 0.001;

	float  pdf_o = pdf_m / (4 * ODotH);
	*pdf = max(pdf_o, 0.01);

	if (pdf_o < 0)
		printf("asfdassdf");
	if (IDotM < 0)
		*pdf = 123456;
	if (dot(oDir, norm) < 0)
		*pdf = 123456;

	return oDir;
}

/*vec3 sampleBRDF(material mat, vec3 norm, float u, float v, vec3 iDir, float* pdf)
{
	float roughness = mat.m;
	float alpha = roughness;
	float alpha_t = (1.2 - 0.2 * sqrt(dot(iDir, norm))) * roughness;
	//float alpha_t = 1.2 * roughness;
	//float alpha_t = alpha;
	if (alpha_t != alpha_t ||
		alpha_t < 0.001)
	{
		printf("%f\n", (double)dot(iDir, norm));
	}
	
	float d1 = nrand();
	float d2 = nrand();
	float d3 = nrand();

	float theta_m = atan(sqrt(-alpha_t * alpha_t * log(1 - d1)));
	float phi_m = 2 * PI * d2;

	float mY = cos(theta_m);
	float mX = cos(phi_m) * sin(theta_m);
	float mZ = sin(phi_m) * sin(theta_m);

	vec3 tangent = getTangent(norm);
	vec3 bitangent = cross(norm, tangent);
	vec3 m = (tangent * mX + norm * mY + bitangent * mZ).normalized();

	float NDotM = dot(norm, m);
	float IDotM = dot(iDir, m);

	float fWeight = Fresnel(iDir, m, mat.Fresnel.r);
	if (d3 > fWeight)
	{
		fallback:
		vec3 oDir = randCosineWeightedRay(norm);
		float cos_theta = max(dot(oDir, norm), 0.01);
		*pdf = cos_theta / PI;
		return oDir;
	}

	if (IDotM < 0)
		goto fallback; // i really dont know what to do here, but this is definately a possibility

	float pdf_m = BeckmanD(m, norm, alpha_t) * NDotM;

	vec3 oDir = m * IDotM * 2.0f - iDir;
	float ODotH = dot(oDir, m);
	if (ODotH < 0.001)
		ODotH = 0.001;

	if (dot(oDir, norm) < 0)
		goto fallback;

	float  pdf_o = pdf_m / (4 * ODotH);
	//pdf_o *= fWeight;
	*pdf = max(pdf_o, 0.01);
	if (pdf_o < 0)
		printf("asfdassdf");
	return oDir;
}*/

color microfacetBRDF(material mat, vec3 norm, float u, float v, vec3 iDir, vec3 oDir)
{
	vec3 lDir = iDir;
	vec3 vDir = oDir;
	vec3 hDir = (lDir + vDir).normalized();

	float LDotH = max(0, dot(lDir, hDir));
	float NDotH = max(0, dot(norm, hDir));

	float LDotN = dot(lDir, norm);
	float ODotN = dot(oDir, norm);

	float roughness = mat.m;
	float alpha = roughness;

	color rf_o = mat.Fresnel;
	color FresnelTerm = color(1, 1, 1);// *Fresnel(iDir, hDir, rf_o.r);

	float DistributionTerm = BeckmanD(hDir, norm, alpha);

	float GeometryTerm = SmithG1Approx(iDir, hDir, norm, alpha) * SmithG1Approx(oDir, hDir, norm, alpha);

	float denom = 4 * LDotN * ODotN;

	color f_microfacet = FresnelTerm * GeometryTerm * DistributionTerm / denom;

	if (f_microfacet.r != f_microfacet.r)
	{
		printf("asdfasfda");
	}
	if (f_microfacet.r < 0 || 1 < f_microfacet.r)
	{
		//printf("Energy conservation is being violated\n");
	}

	return f_microfacet;
}

color lambertBRDF(material mat, vec3 norm, float u, float v, vec3 iDir, vec3 oDir)
{
	color diff = fetch(mat.albedoTex, u, v) / PI;
	return diff;
}

/*color BRDF(material mat, vec3 norm, float u, float v, vec3 iDir, vec3 oDir)
{
	vec3 lDir = iDir;
	vec3 vDir = oDir;
	vec3 hDir = (lDir + vDir).normalized();

	float LDotH = max(0, dot(lDir, hDir));
	float NDotH = max(0, dot(norm, hDir));

	float LDotN = dot(lDir, norm);
	float ODorN = dot(oDir, norm);

	float roughness = mat.m;
	float alpha = roughness;

	color rf_o = mat.Fresnel;
	//color FresnelTerm = rf_o + (color(1, 1, 1) - rf_o) * pow(1 - LDotH, 5);
	//color FresnelTerm = rf_o;
	color FresnelTerm = color(1, 1, 1) * Fresnel(iDir, hDir, rf_o.r);
	//color FresnelTerm = color(1, 1, 1);

	float DistributionTerm = BeckmanD(hDir, norm, alpha);

	float GeometryTerm = SmithG1Approx(iDir, hDir, norm, alpha) * SmithG1Approx(oDir, hDir, norm, alpha);

	float denom = 4 * LDotN * ODorN;

	color f_microfacet = FresnelTerm * GeometryTerm * DistributionTerm / denom;

	color diffuse = fetch(mat.albedoTex, u, v) / 3.14159;
	color spec = f_microfacet;
	color total = diffuse.mul(color(1, 1, 1) - FresnelTerm) + spec;

	//if (total.r > 500)
		//printf("asdfsaf");

	return total;
}*/

