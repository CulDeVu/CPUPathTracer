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
		  NUM_BOUNCES = 2;

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

void importanceSampleBRDF(materialLayer* ml, vec3 norm, vec2 uv, vec3 oDir, vec3* ret_iDir, color* ret_weight)
{
	vec3 half = importanceSampleHalfVector(ml, norm, uv, oDir);

	vec3 iDir = half * dot(oDir, half) * 2.0f - oDir;
	
	*ret_iDir = iDir;
	//*ret_pdf = probabilityDensity(ml, norm, uv.x, uv.y, oDir, half, iDir);
	*ret_weight = sampleWeight(ml, norm, uv, oDir, half, iDir);
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

	color brdf = BRDF(ml, norm, uv, oDir, iDir);
	int chi = 1.f;// dot(norm, norm2) < 0 ? 1 : 0;

	*ret_iDir = iDir;
	//*ret_pdf = max(0.001f, prob / G);
	*ret_weight = clamp(brdf * chi * abs(cosAngle * cosAngle2) / (dot(dirXtoX2, dirXtoX2) * prob), 0, 1);
}

materialLayer* randomlyPickMaterialLayer(layeredMaterial* mat)
{
	float r = nrand();

	for (int i = 0; i < mat->layers.size(); ++i)
	{
		float s = mat->layers[i]->fresnel;
		if (r <= s)
		{
			return mat->layers[i];
		}
		r -= s;
	}

	// should never happen
	printf("randomlyPickMaterialLayer did something very unexpected. recovering gracefully?\n");
	return mat->layers[mat->layers.size() - 1];
}

// NOTE: this works due to BRDF reciprocity. If this changes, change this
void pickMaterialLayerAndBRDFRay(layeredMaterial* mat, vec3 inDir, vec3 norm, vec2 uv, materialLayer** outMat, vec3* outDir, float* outInvPdf)
{
	float d1 = nrand();
	for (int i = 0; i < mat->layers.size(); ++i)
	{
		materialLayer* ml = mat->layers[i];

		vec3 half = importanceSampleHalfVector(ml, norm, uv, inDir);
		vec3 iDir = half * dot(inDir, half) * 2.0f - inDir;

		float f = Fresnel(inDir, half, ml->fresnel);
		if (d1 <= f)
		{
			*outDir = iDir;
			*outMat = ml;
			*outInvPdf = 1.f / probabilityDensity(ml, norm, uv, half, iDir);
			return;
		}

		d1 -= f;
	}

	printf("pickMaterialLayerAndBRDFRay: THIS SHOULD NEVER HAPPEN. THIS IS BAD.\n");
}

struct intersectionInfo
{
	float t;
	vec2 uv;
	vec3 pos;
	vec3 normal;
	layeredMaterial* mat;
};
void getIntersectionInfo(vec3 o, vec3 ray, intersectionInfo* ret)
{
	RTCRay rtcNegODir = makeRay(o, ray);
	rtcIntersect(scene, rtcNegODir);

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
		ret->uv.x = u - roundDown(u);
		ret->uv.y = v - roundDown(v);
	}

	ret->t = rtcNegODir.tfar - 0.001f;

	ret->pos = o + ray * ret->t;
	ret->normal = normalize(vec3(-rtcNegODir.Ng[0], -rtcNegODir.Ng[1], -rtcNegODir.Ng[2]));
	if (dot(ret->normal, ray * -1.0f) < 0)
		ret->normal *= -1.0f;

	ret->mat = &(curModel->mat);
}

color radiance_BDPT(vec3 o, vec3 ray, float bounces)
{
	vec3 pos[4];
	vec2 uv[4];
	vec3 normal[4];
	materialLayer* mat[4];
	float invPdf[4];
	vec3 traceRay[4];
	
	// light path
	{
		float r1, r2;
		r1 = nrand();
		r2 = nrand();

		pos[0] = vec3(3 * r1 - 1.5, 6, 12 * r2 - 6);
		uv[0] = vec2();
		normal[0] = vec3(0, -1, 0);
		mat[0] = lights[0]->mat.layers[0];
		invPdf[0] = 3 * 12;
	}

	// camera path
	{
		pos[3] = o;
		uv[3] = vec2();
		normal[3] = ray;
		mat[3] = nullptr;
		invPdf[3] = 1;
		traceRay[3] = ray;
	}
	{
		intersectionInfo trace;
		getIntersectionInfo(o, ray, &trace);
		pos[2] = trace.pos;
		uv[2] = trace.uv;
		normal[2] = trace.normal;
		//randomlyPickMaterialLayer(trace.mat); // trace.mat->layers[trace.mat->layers.size() - 1];
		pickMaterialLayerAndBRDFRay(trace.mat, -traceRay[3], normal[2], uv[2], &mat[2], &traceRay[2], &invPdf[1]);
		invPdf[2] = 1;
	}
	{
		vec3 iDir = randCosineWeightedRay(normal[2]);

		intersectionInfo trace;
		getIntersectionInfo(pos[2], iDir, &trace);

		pos[1] = trace.pos;
		uv[1] = trace.uv;
		normal[1] = trace.normal;
		mat[1] = randomlyPickMaterialLayer(trace.mat);

		invPdf[1] = PI; // abs(dot(iDir, normal[2]));
	}

	int totalNumberOfSamples = 0;
	color totalColor = color();

	// direct lighting
	{
		vec3 dir = pos[0] - pos[2];
		vec3 ndir = normalize(dir);
		float dirLen = length(dir);
		
		RTCRay ray;
		ray.org[0] = pos[2].x; ray.org[1] = pos[2].y; ray.org[2] = pos[2].z;
		ray.dir[0] = ndir.x; ray.dir[1] = ndir.y; ray.dir[2] = ndir.z;
		ray.tnear = 0;
		ray.tfar = dirLen;
		ray.geomID = RTC_INVALID_GEOMETRY_ID;

		rtcOccluded(scene, ray);
		if (ray.geomID != 0)
		{
			float G = abs(dot(normal[2], ndir) * dot(normal[0], ndir)) / dot(dir, dir);
			G = clamp(G, 0.f, 1.f);
			totalColor = totalColor + BRDF(mat[2], normal[2], uv[2], normalize(pos[3] - pos[2]), ndir).mul(mat[0]->getHue(vec2())) * G * invPdf[0];
			++totalNumberOfSamples;
		}
	}

	// secondary lighting
	{
		vec3 dir = pos[0] - pos[1];
		vec3 ndir = normalize(dir);
		float dirLen = length(dir);

		RTCRay ray;
		ray.org[0] = pos[1].x; ray.org[1] = pos[1].y; ray.org[2] = pos[1].z;
		ray.dir[0] = ndir.x; ray.dir[1] = ndir.y; ray.dir[2] = ndir.z;
		ray.tnear = 0;
		ray.tfar = dirLen;
		ray.geomID = RTC_INVALID_GEOMETRY_ID;

		rtcOccluded(scene, ray);
		if (ray.geomID != 0 && mat[1]->type != EMMISION)
		{
			float G0 = abs(dot(normal[1], ndir) * dot(normal[0], ndir)) / max(0.01f, dot(dir, dir));
			//G0 = clamp(G0, 0.f, 0.5f);
			color L_i = BRDF(mat[1], normal[1], uv[1], normalize(pos[2] - pos[1]), ndir) * G0 * invPdf[0];
			//L_i = color(0.3, 0.3, 0.3);

			dir = pos[1] - pos[2];
			ndir = normalize(dir);

			float G1 = 1;// abs(dot(normal[2], ndir));
			color t = BRDF(mat[2], normal[2], uv[2], normalize(pos[3] - pos[2]), ndir).mul(L_i) * G1 * invPdf[1];
			//t = clamp(t, 0, 1);
			t = t.mul(mat[0]->getHue(vec2()));

			totalColor = totalColor + t;
			++totalNumberOfSamples;
		}
	}

	if (totalNumberOfSamples == 0)
		return color();
	return totalColor / totalNumberOfSamples;
}

color radiance(vec3 o, vec3 ray, float bounces)
{
	vec3 oDir = -ray;
	
	intersectionInfo trace;
	getIntersectionInfo(o, ray, &trace);

	if (trace.t == maxFloat)
		return SKY_ILLUMINATION;
	if (trace.t < 0.001f)
		return color(0, 0, 0);

	color total = color();

	// simulates black body lights
	if (trace.mat->layers[0]->type == EMMISION)
		return trace.mat->layers[0]->getHue(vec2());

	if (bounces == 0)
		return total;

	// microfacet importance sampling
	float d1 = nrand();
	for (int i = 0; i < trace.mat->layers.size(); ++i)
	{
		vec3 iDir;
		color weight;

		enum sampleStrategy
		{
			ss_BRDF, ss_LIGHT
		};
		sampleStrategy curStrategy;

		if (bounces == 1)
			curStrategy = ss_LIGHT;
		else
			curStrategy = ss_BRDF;

		if (curStrategy == ss_LIGHT)
		{
			importanceSampleLight(trace.pos, trace.mat->layers[i], trace.normal, trace.uv, oDir, &iDir, &weight);
		}
		else
			importanceSampleBRDF(trace.mat->layers[i], trace.normal, trace.uv, oDir, &iDir, &weight);

		vec3 half = normalize(iDir + oDir);

		float f = Fresnel(oDir, half, trace.mat->layers[i]->fresnel);
		if (d1 > f)
		{
			d1 -= f;
			continue;
		}

		color L_i;
		if (curStrategy == ss_LIGHT)
			L_i = radiance(trace.pos, iDir, 0);
		else
			L_i = radiance(trace.pos, iDir, bounces - 1);

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
	const int kernelSize = 9;
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
	addObj(scene, "models/dragon_simple.obj", vec3(0, 0, 0), 2.5f);
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

				//float rx = 1.0f * (nrand() - 0.5f) / (imageWidth),
				//ry = 1.0f * (nrand() - 0.5f) / (imageHeight);
				float rx = 0, ry = 0;

				float aspect = (float)imageWidth / imageHeight;
				float fovW = 3.14159 / 4;
				float zNear = 0.5;
				float nearhW = tan(fovW) * zNear;
				float nearhH = nearhW / aspect;

				vec3 ray = vec3(((float)x / imageWidth - 1.0f / 2) * nearhW + rx, ((float)y / imageHeight - 1.0f / 2) * nearhH + ry, -zNear);

				//vec3 o = vec3(-2.77, -10, 3.5f);
				vec3 o = vec3(0, 1, 3.5);
				color sampleColor = radiance_BDPT(o, normalize(ray), NUM_BOUNCES);
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
	GLRender_isNormalized = true;

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
