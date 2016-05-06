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

#include <glm\glm.hpp>
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
const int NUM_SAMPLES = 10,
		  NUM_BOUNCES = 2;

const color SKY_ILLUMINATION = color(12, 12, 12); // color(17, 12, 4);

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

//vector<shape*> shapes;
//vector<triangle> tris;

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
	//total = total + mat.Emmision;

	if (bounces == 0)
		return total;

	// sphere sampling
	/*{
		float invPDF = 1;

		vec3 d = vec3(0.0, 3, 1) - pos;
		float rad = 1;
		float dLen = d.length();
		if (dLen <= 0.01f)
			dLen = 0.01f;

		float cos_maxangle;
		if (1 - rad*rad / (dLen*dLen) <= 0)
			cos_maxangle = 0.0f;
		else
			cos_maxangle = sqrt(1 - rad*rad / (dLen*dLen));
		//float cos_maxangle = cos(asin(rad*rad / (dLen*dLen)));

		vec3 castRay = d.normalized();

		invPDF = (2 * 3.14159 * (1 - cos_maxangle));

		float cosAngle = dot(castRay.normalized(), normal);
		if (cosAngle < 0)
			cosAngle = 0;

		color f_r = BRDF(mat, normal, castRay, ray * -1);
		color L_i = radiance(pos, castRay, 0);

		total = total + f_r.mul(L_i) * cosAngle * invPDF / 2;
	}*/

	float d1 = nrand();
	for (int i = 0; i < mat.layers.size(); ++i)
	{
		vec3 half = importanceSampleHalfVector(mat.layers[i], normal, u, v, oDir);

		float f = Fresnel(oDir, half, mat.layers[i]->Fresnel);
		if (d1 > f)
		{
			d1 -= f;
			continue;
		}

		vec3 iDir = half * dot(oDir, half) * 2.0f - oDir;
		/*float pdf = probabilityDensity(mat.layers[i], normal, u, v, oDir, half, iDir);

		color f_r = BRDF(mat.layers[i], normal, u, v, oDir, iDir);
		color L_i = radiance(pos, iDir, bounces - 1);
		float cosAngle = max(dot(iDir, normal), 0.01);

		total = f_r.mul(L_i) * cosAngle / pdf;*/

		color L_i = radiance(pos, iDir, bounces - 1);
		color weight = sampleWeight(mat.layers[i], normal, u, v, oDir, half, iDir);

		total = L_i.mul(weight);

		break;
	}

	return total;
}

void denoise()
{
	// DOES NOT WORK AT ALL
	// TODO: THIS
	
	int neighborhoodHW = 0;
	
	float std_pos = (neighborhoodHW + 1) / 2.f;
	float std_color = 2.f / NUM_SAMPLES;

	color* buffer_copy = new color[imageWidth * imageHeight];
	memcpy(buffer_copy, buffer, imageWidth * imageHeight * sizeof(color));

	for (int x = neighborhoodHW; x < imageWidth - neighborhoodHW; ++x)
	{
		for (int y = neighborhoodHW; y < imageHeight - neighborhoodHW; ++y)
		{
			int minX = min(x - neighborhoodHW, 0);
			int maxX = max(x + neighborhoodHW + 1, imageWidth);
			int minY = min(y - neighborhoodHW, 0);
			int maxY = max(y + neighborhoodHW + 1, imageHeight);

			color numSum = color();
			float denomSum = 0;

			float I_x_y = buffer_copy[x * imageHeight + y].greyscale();
			
			for (int k = minX; k < maxX; ++k)
			{
				for (int l = minY; l < maxY; ++l)
				{
					//float I_k_l_grey = buffer_copy[k * imageHeight + l].greyscale();
					color I_k_l = buffer_copy[k * imageHeight + l];
					
					/*float exponent = -(
						((x - k)*(x - k) + (y - l)*(y - l)) / (2 * std_pos) +
						((I_x_y - I_k_l_grey) * (I_x_y - I_k_l_grey)) / (2 * std_color)
						);
					float w = exp(exponent);*/

					//numSum = numSum + I_k_l * w;
					//denomSum += w;
					//numSum = numSum + I_k_l;
					numSum = color(1, 1, 1);
					denomSum = 1;
				}
			}

			color g = numSum / denomSum;
			buffer[x][y] = g;
		}
	}
}

int main()
{
	RTCDevice device = rtcNewDevice(NULL);

	threadedRenderWindow();

	scene = rtcDeviceNewScene(device, RTC_SCENE_STATIC, RTC_INTERSECT1);
	addObj(scene, "models/sponza.obj");
	//addObj(scene, "models/teapot.obj", vec3(0.8, 0, -2), 1.0);
	//addObj(scene, "models/lenin.obj", vec3(0.0, -0.1, -1), 1.0);
	//addObj(scene, "models/cube.obj", vec3(0, -0.9, -2), 1.0f);
	//addObj(scene, "models/dragon_simple.obj", vec3(0, 0, 0), 3.0f);
	//models[models.size() - 1]->mat.albedoTex = createSolidTexture(color(0.2, 0.2, 0.2));
	//models[models.size() - 1]->mat.m = 0.02;
	//models[models.size() - 1]->mat.Fresnel = color(0.8, 0.8, 0.8);

	for (int i = 0; i < models.size(); ++i)
	{
		if (models[i]->mat.layers[0]->hueTexId == 21)
		{
			//models[i]->mat.Fresnel = color(0.6, 0.6, 0.6);
			//models[i]->mat.m = 0.001;
			models[i]->mat.addLayerTop(createMicrofacetLayer(0.6, 0.001));
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

				color sampleColor = radiance(vec3(0, 1, 5), normalize(ray), NUM_BOUNCES);

				totalColor = (totalColor * (i - 1.0f) + sampleColor) / i;

				buffer[x][y] = totalColor;
			}
		}
	}

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
	cin.get();

	//rtcDeleteScene(scene);
	//rtcExit();

	exitThreadedRenderWindow();

	return 0;
}
