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

	vec3 castRay = normalize(tangent*rx + bitangent*rz + norm*ry);
	return castRay;
}

enum MATERIAL_LAYER_TYPE
{
	NONE, DIFFUSE, MICROFACET, EMMISION
};
class materialLayer
{
public:
	MATERIAL_LAYER_TYPE type;

	// top layers
	//int fresnelTexId;
	float fresnel;

	// microfacet
	int roughnessTexId;

	// color information
	int hueTexId;
	color hue;

	color getHue(vec2 uv)
	{
		if (hueTexId == -1)
			return hue;
		return fetch(hueTexId, uv.x, uv.y);
	}

	float getRoughness(vec2 uv)
	{
		if (roughnessTexId == -1)
			return 1;
		return 1.f - fetch(roughnessTexId, uv.x, uv.y).r;
	}

	materialLayer()
		: fresnel(0), roughnessTexId(-1), hueTexId(-1) {}
};

materialLayer* createDiffuseLayer(int texId)
{
	materialLayer* newLayer = new materialLayer();
	newLayer->type = DIFFUSE;
	newLayer->fresnel = 1;
	newLayer->hueTexId = texId;
	return newLayer;
}
materialLayer* createMicrofacetLayer(float fresnel, int texId)
{
	materialLayer* newLayer = new materialLayer();
	newLayer->type = MICROFACET;
	newLayer->fresnel = fresnel;
	newLayer->roughnessTexId = texId;
	//newLayer->hueTexId = createSolidTexture(color(1, 1, 1));
	return newLayer;
}
materialLayer* createEmmisionLayer(color c)
{
	materialLayer* newLayer = new materialLayer();
	newLayer->type = EMMISION;
	newLayer->fresnel = 1;
	newLayer->hueTexId = -1;
	newLayer->hue = c;
	return newLayer;
}

class layeredMaterial
{
public:
	vector<materialLayer*> layers;

	void addLayerTop(materialLayer* l)
	{
		layers.insert(layers.begin(), l);
	}
	void addLayerBottom(materialLayer* l)
	{
		layers.push_back(l);
	}
};

//-----------------------------------------------------------------------------
// Microfacet functions
//-----------------------------------------------------------------------------
float Fresnel(vec3 l, vec3 m, float rf0)
{
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

	float theta_m = glm::clamp(acos(NDotM), -PI / 2 + 0.001f, PI / 2 - 0.001f);
	float coef = -pow(tan(theta_m), 2) / (alpha * alpha);
	float denom = pow(alpha, 2) * pow(NDotM, 4);
	if (denom < 0.001)
		denom = 0.001;
	float total = positivity * max(0.001, exp(coef)) / (PI * denom);
	return total;
}
float SmithG1Approx(vec3 v, vec3 m, vec3 norm, float alpha)
{	
	float VDotM = dot(v, m);
	float VDotN = dot(v, norm);

	float theta_v = acos(VDotN);
	float a = 1 / max(0.01, alpha * tan(theta_v));

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

//-----------------------------------------------------------------------------
// Importance Sampling Functions
//-----------------------------------------------------------------------------
vec3 importanceSampleHalfVectorDiffuse(materialLayer* ml, vec3 norm, vec2 uv, vec3 oDir)
{
	vec3 iDir = randCosineWeightedRay(norm);
	return normalize(oDir + iDir);
}
vec3 importanceSampleHalfVectorMicrofacet(materialLayer* ml, vec3 norm, vec2 uv, vec3 oDir)
{
	float alpha_t = ml->getRoughness(uv);

	float d1 = nrand();
	float d2 = nrand();

	float theta_m = atan(sqrt(-alpha_t * alpha_t * log(1 - d1)));
	float phi_m = 2 * PI * d2;

	float mY = cos(theta_m);
	float mX = cos(phi_m) * sin(theta_m);
	float mZ = sin(phi_m) * sin(theta_m);

	vec3 tangent = getTangent(norm);
	vec3 bitangent = cross(norm, tangent);
	vec3 m = normalize(tangent * mX + norm * mY + bitangent * mZ);

	return m;
}
vec3 importanceSampleHalfVector(materialLayer* ml, vec3 norm, vec2 uv, vec3 oDir)
{
	if (ml->type == DIFFUSE)
		return importanceSampleHalfVectorDiffuse(ml, norm, uv, oDir);
	if (ml->type == MICROFACET)
		return importanceSampleHalfVectorMicrofacet(ml, norm, uv, oDir);
	return norm;
}

//-----------------------------------------------------------------------------
// Probability densities
//-----------------------------------------------------------------------------
float probabilityDensityDiffuse(materialLayer* ml, vec3 norm, vec2 uv, vec3 half, vec3 iDir)
{
	float cosAngle = dot(iDir, norm);
	return cosAngle / PI;
}
float probabilityDensityMicrofacet(materialLayer* ml, vec3 norm, vec2 uv, vec3 half, vec3 iDir)
{
	float alpha = ml->getRoughness(uv);
	
	float pdf = BeckmanD(half, norm, alpha) * dot(half, norm);
	pdf /= 4 * dot(iDir, half);
	return pdf;
}
float probabilityDensity(materialLayer* ml, vec3 norm, vec2 uv, vec3 half, vec3 iDir)
{
	if (ml->type == DIFFUSE)
		return probabilityDensityDiffuse(ml, norm, uv, half, iDir);
	if (ml->type == MICROFACET)
		return probabilityDensityMicrofacet(ml, norm, uv, half, iDir);
	return 1;
}

//-----------------------------------------------------------------------------
// BRDFs
//-----------------------------------------------------------------------------
color BRDFDiffuse(materialLayer* ml, vec3 norm, vec2 uv, vec3 oDir, vec3 iDir)
{
	return ml->getHue(uv) / PI;
}
color BRDFMicrofacet(materialLayer* ml, vec3 norm, vec2 uv, vec3 oDir, vec3 iDir)
{
	vec3 lDir = iDir;
	vec3 vDir = oDir;
	vec3 hDir = normalize(lDir + vDir);

	float LDotH = max(0.001, dot(lDir, hDir));
	float NDotH = max(0.001, dot(norm, hDir));

	float LDotN = dot(lDir, norm);
	float ODotN = dot(oDir, norm);

	float roughness = ml->getRoughness(uv);
	float alpha = roughness;

	//color rf_o = mat.Fresnel;
	color FresnelTerm = color(1, 1, 1);// *Fresnel(iDir, hDir, rf_o.r);

	float DistributionTerm = BeckmanD(hDir, norm, alpha);

	float GeometryTerm = SmithG1Approx(iDir, hDir, norm, alpha) * SmithG1Approx(oDir, hDir, norm, alpha);

	float denom = max(0.001, 4 * LDotN * ODotN);

	color f_microfacet = FresnelTerm * GeometryTerm * DistributionTerm / denom;

	if (f_microfacet.r != f_microfacet.r)
	{
		printf("Microfacet BRDF experiencing discontinuities\n");
	}

	return clamp(f_microfacet, 0, 10);
}
color BRDF(materialLayer* ml, vec3 norm, vec2 uv, vec3 oDir, vec3 iDir)
{
	if (ml->type == DIFFUSE)
		return BRDFDiffuse(ml, norm, uv, oDir, iDir);
	if (ml->type == MICROFACET)
		return BRDFMicrofacet(ml, norm, uv, oDir, iDir);
	return color(1, 0, 1);
}

//-----------------------------------------------------------------------------
// Weights
//-----------------------------------------------------------------------------
color sampleWeightDiffuse(materialLayer* ml, vec3 norm, vec2 uv, vec3 oDir, vec3 half, vec3 iDir)
{
	return ml->getHue(uv);
}
color sampleWeightMicrofacet(materialLayer* ml, vec3 norm, vec2 uv, vec3 oDir, vec3 half, vec3 iDir)
{
	float alpha = ml->getRoughness(uv);
	float  G = SmithG1Approx(iDir, half, norm, alpha) * SmithG1Approx(oDir, half, norm, alpha);
	color final = color(1, 1, 1) * G * dot(oDir, half) / max(0.001, dot(iDir, norm) * dot(half, norm));
	return clamp(final, 0, 10);
}
color sampleWeight(materialLayer* ml, vec3 norm, vec2 uv, vec3 oDir, vec3 half, vec3 iDir)
{
	if (ml->type == DIFFUSE)
		return sampleWeightDiffuse(ml, norm, uv, oDir, half, iDir);
	if (ml->type == MICROFACET)
		return sampleWeightMicrofacet(ml, norm, uv, oDir, half, iDir);
	return color(1, 0, 1);
}
