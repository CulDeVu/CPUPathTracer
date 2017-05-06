#define _CRT_SECURE_NO_DEPRECATE

#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <climits>
#include <time.h>
#include <chrono>
#include <random>

using namespace std;
using std::cout;
using std::cin;

#include "glm\glm.hpp"
using glm::vec3;
using namespace glm;

#include "triangle.h"
#include "color.h"
#include "material.h"
#include "model.h"
#include "useful.h"
#include "vec3.h"

#include "embree2\rtcore.h"
#include "embree2\rtcore_ray.h"

const int imageWidth = 500,
		  imageHeight = 500;
const int NUM_SAMPLES = 50,
		  NUM_BOUNCES = 1;

const color SKY_ILLUMINATION = SKY_BLACK;

color buffer[imageWidth][imageHeight];
#include "GLRender.h"
//#include "fastBilateral.h"

#pragma comment(lib, "embree.lib")

vec3 randRayInSphere() {
	float rx = 0, ry = 0, rz = 0;

	// technically, this is faster, sooooo :P
	while (rx*rx + ry*ry + rz*rz <= 1.0f)
	{
		rx = 2 * nrand() - 1.0f;
		ry = 2 * nrand() - 1.0f;
		rz = 2 * nrand() - 1.0f;
		if (rx*rx + ry*ry + rz*rz >= 0.98f)
			continue;
		break;
	}
	return normalize(vec3(rx, ry, rz));
}

vec3 randHemisphereRay(vec3 norm)
{
	vec3 castRay = randRayInSphere();
	if (dot(castRay, norm) < 0)
		castRay = -castRay;
	return normalize(castRay);
}

RTCScene scene;

RTCRay makeRay(vec3 o, vec3 dir)
{
	RTCRay ray;
	ray.org[0] = o.x; ray.org[1] = o.y; ray.org[2] = o.z;
	ray.dir[0] = dir.x; ray.dir[1] = dir.y; ray.dir[2] = dir.z;
	ray.tnear = 0;
	ray.tfar = maxFloat;
	ray.geomID = RTC_INVALID_GEOMETRY_ID;
	return ray;
}

float roundDown(float x)
{
	float ret = int(x);
	if (x < 0)
		--ret;
	return ret;
}

struct pathNode
{
	vec3 pt;
	vec3 iDir;
	vec3 norm;
	color f_r;
	float pdf;
};

pathNode importanceSample(vec3 pt, layeredMaterial mat, vec3 norm, vec2 uv, vec3 oDir)
{
	vec3 half = importanceSampleHalfVector(mat.layers[0], norm, uv.x, uv.y, oDir);

	vec3 iDir = half * dot(oDir, half) * 2.0f - oDir;

	color weight = sampleWeight(mat.layers[0], norm, uv.x, uv.y, oDir, half, iDir);

	pathNode ret;
	ret.pt = pt;
	ret.iDir = iDir;
	ret.norm = norm;
	ret.f_r = weight;
	ret.pdf = abs(dot(norm, iDir));

	return ret;
}

void importanceSampleBRDF(vec3 pt, materialLayer* ml, vec3 norm, vec2 uv, vec3 oDir, vec3* ret_iDir, color* ret_weight)
{
	vec3 half = importanceSampleHalfVector(ml, norm, uv.x, uv.y, oDir);

	vec3 iDir = half * dot(oDir, half) * 2.0f - oDir;
	
	*ret_iDir = iDir;
	//*ret_pdf = probabilityDensity(ml, norm, uv.x, uv.y, oDir, half, iDir);
	*ret_weight = sampleWeight(ml, norm, uv.x, uv.y, oDir, half, iDir);
}
void importanceSampleLight(vec3 pt, materialLayer* ml, vec3 norm, vec2 uv, vec3 oDir, vec3* ret_iDir, color* ret_weight)
{
	float r1, r2;
	r1 = nrand();
	r2 = nrand();

	vec3 x2 = vec3(3 * r1 - 1.5, 6, 12 * r2 - 6);
	vec3 dirXtoX2 = x2 - pt;
	vec3 iDir = normalize(x2 - pt);

	vec3 norm2 = vec3(0, -1, 0);
	float prob = 1.f / (3.f * 12.f);
	float cosAngle = dot(norm, iDir);
	float cosAngle2 = dot(norm2, -iDir);
	float G = max(0.001, cosAngle2) / dot(dirXtoX2, dirXtoX2);

	color brdf = BRDF(ml, norm, uv.x, uv.y, iDir, oDir);
	int chi = dot(norm, norm2) < 0 ? 1 : 0;

	*ret_iDir = iDir;
	//*ret_pdf = max(0.001f, prob / G);
	*ret_weight = brdf * chi * abs(cosAngle * cosAngle2) / (dot(dirXtoX2, dirXtoX2) * prob);
}


color radiance2(vec3 o, vec3 ray, int bounces)
{
	vector<pathNode> path;

	pathNode camera;
	camera.pt = o;
	camera.iDir = ray;
	camera.f_r = color(1, 1, 1);
	//camera.prob_dA = 1;

	vec3 pos = o;
	vec3 oDir = -ray;
	for (int i = 0; i < bounces; ++i)
	{
		RTCRay rtcODir = makeRay(pos, oDir);
		rtcIntersect(scene, rtcODir);
		if (rtcODir.tfar == maxFloat)
			return color(0, 0, 0);

		model* curModel = models[rtcODir.geomID];

		// calculate the UV coords
		float u = 1, v = 1;
		if (curModel->uv.size() != 0)
		{
			int primId = rtcODir.primID;
			float u0, u1, u2;
			float v0, v1, v2;

			u0 = curModel->uv[2 * curModel->indices[3 * primId + 0]];
			u1 = curModel->uv[2 * curModel->indices[3 * primId + 1]];
			u2 = curModel->uv[2 * curModel->indices[3 * primId + 2]];

			v0 = curModel->uv[2 * curModel->indices[3 * primId + 0] + 1];
			v1 = curModel->uv[2 * curModel->indices[3 * primId + 1] + 1];
			v2 = curModel->uv[2 * curModel->indices[3 * primId + 2] + 1];

			u = rtcODir.u * u1 + rtcODir.v * u2 + (1 - rtcODir.u - rtcODir.v) * u0;
			v = rtcODir.u * v1 + rtcODir.v * v2 + (1 - rtcODir.u - rtcODir.v) * v0;
			u = u - roundDown(u);
			v = v - roundDown(v);
		}

		float t = rtcODir.tfar;
		t -= 0.001f;
		if (t < 0.001f)
			return color(0, 0, 0);

		vec3 pos = o + ray * t;
		vec3 normal = normalize(vec3(-rtcODir.Ng[0], -rtcODir.Ng[1], -rtcODir.Ng[2]));
		if (dot(normal, ray * -1.0f) < 0)
			normal = normal * -1.0f;

		layeredMaterial mat = curModel->mat;

		color total = color();

		// simulates black body lights
		if (mat.layers[0]->type == EMMISION)
			return mat.layers[0]->getHue(0, 0);

		float d1 = nrand();
		for (int i = 0; i < mat.layers.size(); ++i)
		{

		}
	}
}

color radiance(vec3 o, vec3 ray, float bounces)
{
	vec3 oDir = -ray;
	
	RTCRay rtcNegODir = makeRay(o, ray);
	rtcIntersect(scene, rtcNegODir);
	if (rtcNegODir.tfar == maxFloat)
	{
		return SKY_ILLUMINATION;
	}

	model* curModel = models[rtcNegODir.geomID];

	// calculate the UV coords
	float u = 1, v = 1;
	if (curModel->uv.size() != 0)
	{
		int primId = rtcNegODir.primID;
		float u0, u1, u2;
		float v0, v1, v2;

		u0 = curModel->uv[2 * curModel->indices[3 * primId + 0]];
		u1 = curModel->uv[2 * curModel->indices[3 * primId + 1]];
		u2 = curModel->uv[2 * curModel->indices[3 * primId + 2]];

		v0 = curModel->uv[2 * curModel->indices[3 * primId + 0] + 1];
		v1 = curModel->uv[2 * curModel->indices[3 * primId + 1] + 1];
		v2 = curModel->uv[2 * curModel->indices[3 * primId + 2] + 1];

		u = rtcNegODir.u * u1 + rtcNegODir.v * u2 + (1 - rtcNegODir.u - rtcNegODir.v) * u0;
		v = rtcNegODir.u * v1 + rtcNegODir.v * v2 + (1 - rtcNegODir.u - rtcNegODir.v) * v0;
		u = u - roundDown(u);
		v = v - roundDown(v);
	}

	float t = rtcNegODir.tfar;
	t -= 0.001f;
	if (t < 0.001f)
		return color(0, 0, 0);

	vec3 pos = o + ray * t;
	vec3 normal = normalize(vec3(-rtcNegODir.Ng[0], -rtcNegODir.Ng[1], -rtcNegODir.Ng[2]));
	if (dot(normal, ray * -1.0f) < 0)
		normal = normal * -1.0f;

	layeredMaterial mat = curModel->mat;

	color total = color();

	// simulates black body lights
	if (mat.layers[0]->type == EMMISION)
		return mat.layers[0]->getHue(0, 0);

	if (bounces == 0)
		return total;

	// microfacet importance sampling
	float d1 = nrand();
	for (int i = 0; i < mat.layers.size(); ++i)
	{
		/*vec3 half = importanceSampleHalfVector(mat.layers[i], normal, u, v, oDir);

		float f = Fresnel(oDir, half, mat.layers[i]->Fresnel);
		if (d1 > f)
		{
			d1 -= f;
			continue;
		}

		vec3 iDir = half * dot(oDir, half) * 2.0f - oDir;

		color L_i = radiance(pos, iDir, bounces - 1);
		color weight = sampleWeight(mat.layers[i], normal, u, v, oDir, half, iDir);

		total = L_i.mul(weight);

		break;*/

		vec3 iDir;
		color weight;

		enum sampleStrategy
		{
			ss_BRDF, ss_LIGHT
		};
		sampleStrategy curStrategy;
		float mul = 1.0f;

		/*if (bounces == 1)
		{
			curStrategy = ss_LIGHT;
			mul = 1.0f;
		}
		else
		{
			float r1 = nrand();
			if (r1 < 0.5f && mat.layers[i]->type == DIFFUSE)
			{
				curStrategy = ss_LIGHT;
				mul = 0.0f;
			}
			else
			{
				curStrategy = ss_BRDF;
				mul = 1.0f;
			}
		}*/
		curStrategy = ss_BRDF;

		if (curStrategy == ss_LIGHT)
		{
			importanceSampleLight(pos, mat.layers[i], normal, vec2(u, v), oDir, &iDir, &weight);
		}
		else
			importanceSampleBRDF(pos, mat.layers[i], normal, vec2(u, v), oDir, &iDir, &weight);

		vec3 half = normalize(iDir + oDir);

		float f = Fresnel(oDir, half, mat.layers[i]->getFresnel(vec2(u, v)));
		if (d1 > f)
		{
			d1 -= f;
			continue;
		}

		color L_i;
		if (curStrategy == ss_LIGHT)
			L_i = radiance(pos, iDir, 0);
		else
			L_i = radiance(pos, iDir, bounces - 1);

		total = L_i.mul(weight);

		break;
	}

	return total;
}

float w(vec2 d, color c, float var_d, float var_c)
{
	float s = 1, r = 1;
	if (var_d != 0)
		s = exp(-dot(d, d) / (2 * var_d));
	if (var_c != 0)
		r = exp(-dot(c, c) / (2 * var_c));
	return s*r;
}
void bilateral()
{
	const int kernelSize = 21;
	const int hsize = (kernelSize - 1) / 2;

	color* computedBuffer = new color[imageWidth*imageHeight];

	// d = position difference
	// c = color difference

	// for each pixel in image
	for (int x = 0; x < imageWidth; ++x)
	{
		for (int y = 0; y < imageHeight; ++y)
		{
			// point p = (x,y)
			// point q = (i,j)
			
			color c_p = buffer[x][y];

			// compute the limits of the kernel
			int min_i = max(x - hsize, 0);
			int max_i = min(x + hsize, imageWidth - 1);
			int min_j = max(y - hsize, 0);
			int max_j = min(y + hsize, imageHeight - 1);

			// add up the sum and sum squared for mean and variance calculation
			float sum_sqr_d = 0,
				sum_d = 0,
				sum_sqr_c = 0,
				sum_c = 0;
			for (int i = min_i; i <= max_i; ++i)
			{
				for (int j = min_j; j <= max_j; ++j)
				{
					color c_q = buffer[i][j];
					
					vec2 d = vec2(i, j) - vec2(x, y);
					color c = c_q - c_p;

					sum_sqr_d += dot(d, d);
					sum_d += length(d);

					sum_sqr_c += dot(c, c);
					sum_c += length(c);
				}
			}

			float mean_d = sum_d / (kernelSize*kernelSize);
			float mean_c = sum_c / (kernelSize*kernelSize);

			float var_d = sum_sqr_d / (kernelSize*kernelSize) - mean_d*mean_d;
			float var_c = sum_sqr_c / (kernelSize*kernelSize) - mean_c*mean_c;

			color num = color(0, 0, 0);
			float denom = 0;

			for (int i = min_i; i <= max_i; ++i)
			{
				for (int j = min_j; j <= max_j; ++j)
				{
					color c_q = buffer[i][j];

					float w_pq = w(vec2(i, j) - vec2(x, y), c_p - c_q, var_d, var_c);

					num = num + c_q * w_pq;
					denom += w_pq;
				}
			}
			
			computedBuffer[x*imageHeight + y] = num /  denom;
		}
	}

	memcpy(buffer, computedBuffer, imageWidth*imageHeight * sizeof(color));
	delete[] computedBuffer;
}

int main()
{
	RTCDevice device = rtcNewDevice(NULL);

	threadedRenderWindow();

	scene = rtcDeviceNewScene(device, RTC_SCENE_STATIC, RTC_INTERSECT1);
	addObj(scene, "models/sponza_light.obj");
	//addObj(scene, "C:/Users/Daniel/Downloads/sanMiguel/sanMiguel.obj");
	//addObj(scene, "models/a.obj");
	//addObj(scene, "models/teapot.obj", vec3(0.8, 0, -2), 1.0);
	//addObj(scene, "models/lenin.obj", vec3(0.0, -0.1, -1), 1.0);
	//addObj(scene, "models/cube.obj", vec3(0, -0.9, -2), 1.0f);
	//addObj(scene, "models/dragon_simple.obj", vec3(0, 0, -1), 2.5f);
	//addObj(scene, "models/CornellBox-Original.obj", vec3(0, 0, 2));
	//models[models.size() - 1]->mat.albedoTex = createSolidTexture(color(0.2, 0.2, 0.2));
	//models[models.size() - 1]->mat.m = 0.02;
	//models[models.size() - 1]->mat.Fresnel = color(0.8, 0.8, 0.8);

	//models[models.size() - 1]->mat.addLayerTop(createEmmisionLayer(color(17, 12, 4)));

	for (int i = 0; i < models.size(); ++i)
	{
		if (models[i]->mat.layers[0]->hueTexId == 21)
		{
			//models[i]->mat.Fresnel = color(0.6, 0.6, 0.6);
			//models[i]->mat.m = 0.001;
			//models[i]->mat.addLayerTop(createMicrofacetLayer(0.6, createSolidTexture(color(0.001, 0.001, 0.001))));
			//printf("%d\n", i);
		}
		/*else
		{
			models[i]->mat.Fresnel = color(0, 0, 0);
			models[i]->mat.m = 1;
		}*/
	}

	rtcCommit(scene);

	auto start = chrono::high_resolution_clock::now();

	//#pragma omp parallel for schedule(dynamic)
	for (int i = 1; i <= NUM_SAMPLES; ++i)
	{

		cout << "sample " << i << " of " << NUM_SAMPLES << endl;
		for (int y = 0; y < imageHeight; ++y)
		{
			for (int x = 0; x < imageWidth; ++x)
			{
				color totalColor = buffer[x][y];

				float rx = 1.0f * (nrand() - 0.5f) / (imageWidth),
				ry = 1.0f * (nrand() - 0.5f) / (imageHeight);
				//float rx = 0, ry = 0;

				float aspect = (float)imageWidth / imageHeight;
				float fovW = 3.14159 / 4;
				float zNear = 0.5;
				float nearhW = tan(fovW) * zNear;
				float nearhH = nearhW / aspect;

				vec3 ray = vec3(((float)x / imageWidth - 1.0f / 2) * nearhW + rx, ((float)y / imageHeight - 1.0f / 2) * nearhH + ry, -zNear);

				//vec3 o = vec3(-2.77, -10, 3.5f);
				vec3 o = vec3(0, 1, 5);
				color sampleColor = radiance(o, normalize(ray), NUM_BOUNCES);
				//color sampleColor = asdfjkl(vec3(0, 1, 3.5f), normalize(ray));

				totalColor = (totalColor * (i - 1.0f) + sampleColor) / i;

				buffer[x][y] = totalColor;
			}
		}
	}

	bilateral();

	for (int y = imageHeight - 1; y >= 0; --y)
	{
		for (int x = 0; x < imageWidth; ++x)
		{
			buffer[x][y] = buffer[x][y].normalized();
		}
	}

	auto end = chrono::high_resolution_clock::now();

	cout << "Elapsed time: " << chrono::duration_cast<chrono::hours>(end - start).count() << " hours" << endl;
	cout << "Elapsed time: " << chrono::duration_cast<chrono::minutes>(end - start).count() << " minutes" << endl;
	cout << "Elapsed time: " << chrono::duration_cast<chrono::seconds>(end - start).count() << " seconds" << endl;
	cout << "Elapsed time: " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " milliseconds" << endl;
	cout << "     Seconds per pixel: " << (float)chrono::duration_cast<chrono::seconds>(end - start).count() / (imageWidth*imageHeight) << endl;
	cout << "Milliseconds per pixel: " << (float)chrono::duration_cast<chrono::milliseconds>(end - start).count() / (imageWidth*imageHeight) << endl;
	cout << "Microseconds per pixel: " << (float)chrono::duration_cast<chrono::microseconds>(end - start).count() / (imageWidth*imageHeight) << endl;

	//denoise();

	ofstream file("image2.ppm");
	file << "P3 " << imageWidth << " " << imageHeight << " 255" << endl;

	for (int y = imageHeight - 1; y >= 0; --y)
	{
		for (int x = 0; x < imageWidth; ++x)
		{
			if ((int)(buffer[x][y].r * 255) == 0x80000000 ||
				(int)(buffer[x][y].g * 255) == 0x80000000 ||
				(int)(buffer[x][y].b * 255) == 0x80000000)
			{
				//cout << "a divide by zero happened" << endl;
				/*for (int i = 0; i < 50; ++i)
					file << "\n";*/
				file << 0 << " " << 0 << " " << 0 << " ";
			}
			else
				file << (int)(buffer[x][y].r * 255) << " " << (int)(buffer[x][y].g * 255) << " " << (int)(buffer[x][y].b * 255) << " ";
		}
	}

	file.close();

	cout << "Finished" << endl;
	getchar();

	//rtcDeleteScene(scene);
	//rtcExit();

	exitThreadedRenderWindow();

	return 0;
}
