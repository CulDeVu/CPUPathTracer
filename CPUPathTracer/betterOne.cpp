#include <iostream>
#include <fstream>
#include <string>
#include <stdlib.h>
#include <math.h>
#include <vector>
#include <climits>
#include <time.h>
#include <chrono>

using namespace std;
using std::cout;
using std::cin;

#include "useful.h"

const int imageWidth = 500,
imageHeight = 500;
int NUM_SAMPLES = 50,
NUM_BOUNCES = 2;

vec3 randRayInSphere() {
	float rx, ry, rz;

	// technically, this is faster, sooooo :P
	while (true)
	{
		rx = 2 * (float)rand() / RAND_MAX - 1.0f;
		ry = 2 * (float)rand() / RAND_MAX - 1.0f;
		rz = 2 * (float)rand() / RAND_MAX - 1.0f;
		if (rx*rx + ry*ry + rz*rz >= 0.98f)
			continue;
		break;
	}
	return newray(rx, ry, rz);
}

vec3 getTangent(vec3 norm)
{
	vec3 tangent;
	vec3 c1 = cross(norm, vec3(0, 0, 1));
	vec3 c2 = cross(norm, vec3(0, 1, 0));
	if (c1.lengthsqr() > c2.lengthsqr())
		tangent = c1;
	else
		tangent = c2;
	return tangent;
}

float curPDF = 0;
vec3 randHemisphereRay(vec3 norm)
{
	vec3 castRay = randRayInSphere();
	if (dot(castRay, norm) < 0)
		castRay = castRay * -1;
	curPDF = 1 / 3.14159;
	return castRay.normalized();
}
vec3 randHemisphereWeightedRay(vec3 norm)
{
	float u1 = (float)rand() / RAND_MAX;
	float u2 = (float)rand() / RAND_MAX;

	float rx, rz;
	while (true)
	{
		rx = 2 * (float)rand() / RAND_MAX - 1.0f;
		rz = 2 * (float)rand() / RAND_MAX - 1.0f;
		if (rx*rx + rz*rz >= 1)
			continue;
		break;
	}
	float ry = sqrt(1 - rx*rx - rz*rz);

	vec3 tangent = getTangent(norm);
	vec3 bitangent = cross(norm, tangent);

	vec3 castRay = (tangent*rx + bitangent*rz + norm*ry).normalized();
	curPDF = max(0.0, dot(norm, castRay) / 3.14159);
	return castRay;
}

color brdf(shape* s, vec3 hitPos, vec3 inDir, vec3 outDir)
{
	vec3 norm = s->getNorm(hitPos);
	color diffColor = s->getDiffuseColor(hitPos);
	color specColor = s->getSpecularColor(hitPos);
	float specN = s->getSpecularTerm(hitPos);

	float cosAngle = dot(outDir, norm);
	if (cosAngle < 0)
		cosAngle = 0;

	vec3 half_angle = outDir + inDir * -1;
	half_angle = half_angle.normalized();
	float nDotHPow = pow(dot(half_angle, norm), specN);
	if (nDotHPow < 0)
		nDotHPow = 0;

	//color f_r = diffColor / 3.14159f + specColor * (specN + 8) / (8 * 3.14159) * nDotHPow;
	color f_r = diffColor / 3.14159f;
	return f_r;
}

color buffer[imageWidth][imageHeight];
color diffBuffer[imageWidth][imageHeight];
int globalX, globalY;
vector<shape*> shapes;
int globalSampleNum = 0;

vec3 lightSamplePos = newray(0, 0, 0);
float lightSampleRad = 0;

color sample(vec3 o, vec3 ray, float bounces)
{
	float closestT = maxFloat;
	shape* closestShape = 0;

	for (int i = 0; i < shapes.size(); ++i)
	{
		float t = shapes[i]->intersect(o, ray);

		if (t > 0 && t < closestT)
		{
			closestT = t;
			closestShape = shapes[i];
		}
	}

	if (closestShape == 0)
		return color(0.0f);

	// i hate floating point numbers >:(
	closestT -= 0.01f;

	vec3 intersected = o + ray * closestT;
	vec3 norm = closestShape->getNorm(intersected);
	color diffColor = closestShape->getDiffuseColor(intersected);
	color specColor = closestShape->getSpecularColor(intersected);
	float specN = closestShape->getSpecularTerm(intersected);

	color emmisionTerm = closestShape->getEmmision(intersected);

	color reflectTerm = color(0);
	color directIllum = color(0);
	color specIllum = color(0);

	if (bounces == 0)
		return emmisionTerm;
	/*if (diffColor == color() && specColor == color())
	return emmisionTerm;*/

	int N = 0;

	/*float k_roulette = (float)rand() / RAND_MAX;
	float k_diffuse = (diffColor.r + diffColor.g + diffColor.b) / 3.0f;
	float k_specular = (specColor.r + specColor.g + specColor.b) / 3.0f;*/

	/*if (k_roulette < k_specular)
	{
	++N;

	ray = ray.normalized();

	vec3 reflectDir = ray - norm * dot(ray, norm) * 2.0f;
	reflectDir = reflectDir.normalized();

	vec3 castRay = randRayInHemisphereOf(reflectDir);
	//vec3 castRay = reflectDir;
	color f_r = brdf(closestShape, intersected, ray.normalized(), castRay);

	color L_i = sample(intersected, castRay, bounces - 1);

	float PDF = 1 / 3.14159; // hemisphere weighted sampling
	reflectTerm = f_r.mul(L_i) / PDF / k_specular;
	}
	else if (k_roulette < k_diffuse + k_specular)*/
	{
		++N;

		vec3 castRay = randHemisphereRay(norm);
		color f_r = brdf(closestShape, intersected, ray.normalized(), castRay);

		color L_i = sample(intersected, castRay, bounces - 1);

		reflectTerm = f_r.mul(L_i) * dot(norm, castRay) / curPDF;
		if (bounces == NUM_BOUNCES)
		{
			color celLight = L_i;
			celLight.r = (int)(celLight.r * 2) / 2.f;
			celLight.g = (int)(celLight.g * 2) / 2.f;
			celLight.b = (int)(celLight.b * 2) / 2.f;
			//reflectTerm = diffColor.mul(celLight);
			reflectTerm = L_i;

			diffBuffer[globalX][globalY] = diffColor + emmisionTerm.normalized();
		}
	}
	// light ray
	/*{
	++N;

	float PDF = 1;

	float rad = sqrt(3);
	vec3 d = vec3(-3.5f, -2, -10) + randRayInSphere()*rad - intersected;
	float dLen = d.length();
	if (dLen <= 0.01f)
	dLen = 0.01f;

	float cos_maxangle;
	if (1 - rad*rad / (dLen*dLen) <= 0)
	cos_maxangle = 0.0f;
	else
	cos_maxangle = sqrt(1 - rad*rad / (dLen*dLen));

	vec3 castRay = d.normalized();

	PDF = 1 / (2 * 3.14159 * (1 - cos_maxangle));

	float cosAngle = dot(castRay.normalized(), norm);
	if (cosAngle < 0)
	cosAngle = 0;

	color L_i = sample(intersected, castRay, 0);
	color f_r = brdf(closestShape, intersected, ray.normalized(), castRay);

	directIllum = directIllum + f_r.mul(L_i) * cosAngle / PDF;
	}
	// light ray
	/*{
	++N;

	float PDF = 1;

	float rad = 1.5f;
	vec3 d = vec3(3, -1.5f, -11) + randRayInSphere()*rad - intersected;
	float dLen = d.length();
	if (dLen <= 0.01f)
	dLen = 0.01f;

	float cos_maxangle;
	if (1 - rad*rad / (dLen*dLen) <= 0)
	cos_maxangle = 0.0f;
	else
	cos_maxangle = sqrt(1 - rad*rad / (dLen*dLen));

	vec3 castRay = d.normalized();

	PDF = 1 / (2 * 3.14159 * (1 - cos_maxangle));

	float cosAngle = dot(castRay.normalized(), norm);
	if (cosAngle < 0)
	cosAngle = 0;

	color L_i = sample(intersected, castRay, 0);
	color f_r = brdf(closestShape, intersected, ray.normalized(), castRay);

	directIllum = directIllum + f_r.mul(L_i) * cosAngle / PDF;
	}*/

	N = max(N, 1);

	reflectTerm = (directIllum + specIllum + reflectTerm) / N;

	color L_o = emmisionTerm + reflectTerm;

	return L_o;
}

int main()
{
	srand(time(0));

	cout << "Continue last render? (y/n)" << endl;
	char ans;
	int prevSamples = 0;
	cin >> ans;
	if (ans == 'y') {
		cout << "how many samples was it?" << endl;
		cin >> prevSamples;
		NUM_SAMPLES += prevSamples;

		string s;
		ifstream file("image2.ppm");
		file >> s; file >> s; file >> s; file >> s;
		for (int y = 0; y < imageHeight; ++y)
		{
			for (int x = 0; x < imageWidth; ++x)
			{
				color col = color();
				file >> s;
				col.r = (float)atoi(s.c_str()) / 255.0f;
				file >> s;
				col.g = (float)atoi(s.c_str()) / 255.0f;
				file >> s;
				col.b = (float)atoi(s.c_str()) / 255.0f;

				col = col.denormalized();
				col = col * prevSamples / (NUM_SAMPLES);
				buffer[x][imageHeight - y - 1] = col;

			}
		}
		file.close();
	}
	else {
		for (int y = 0; y < imageHeight; ++y)
			for (int x = 0; x < imageWidth; ++x)
				buffer[x][y] = color(0.0f);
	}

	plane* bottom = new plane(newray(0, -3, 0), newray(0, 1, 0), color(0.25f, 0.25f, 0.25f));
	bottom->setSpecularColor(color(0.75f, 0.75f, 0.75f));
	bottom->setSpecularTerm(50);
	shapes.push_back(bottom);

	plane* front = new plane(newray(0, 0, -20), newray(0, 0, 1), color(0.75f, 0.75f, 0.75f));
	shapes.push_back(front);

	plane* left = new plane(newray(-5, 0, 0), newray(1, 0, 0), color(0.75f, 0.25f, 0.25f));
	shapes.push_back(left);

	plane* right = new plane(newray(5, 0, 0), newray(-1, 0, 0), color(0.25f, 0.25f, 0.75f));
	shapes.push_back(right);

	plane* top = new plane(newray(0, 5, 0), newray(0, -1, 0), color(0.75f, 0.75f, 0.75f));
	shapes.push_back(top);

	plane* back = new plane(newray(0, 0, 5), newray(0, 0, -1), color(0.25f, 0.75f, 0.25f));
	shapes.push_back(back);

	circle* ball1 = new circle(newray(-1, -1.5, -15), 1.5f);
	ball1->setColor(color(0.15f, 0.15f, 0.4f));
	ball1->setSpecularColor(color(0.6f, 0.6f, 0.6f));
	ball1->setSpecularTerm(200);
	shapes.push_back(ball1);

	circle* ball3 = new circle(newray(3, -1.5f, -11), 1.5f);
	ball3->setColor(color(0.8f, 0.8f, 0.8f));
	ball3->setSpecularColor(color(0.2f, 0.2f, 0.2f));
	ball3->setSpecularTerm(100);
	shapes.push_back(ball3);

	cube* light = new cube(vec3(-3.5f, -2, -10), 1.0f, 1.0f, 1.0f);
	light->setColor(color(0));
	light->setEmmision(color(50.0f, 50.0f, 50.0f));
	shapes.push_back(light);

	/*cube* r = new cube(vec3(-4, 3, -10), 0.25f, 50, 0.25f);
	r->setEmmision(color(8, 0, 0));
	shapes.push_back(r);

	cube* g = new cube(vec3(0, 3, -10), 50, 0.25f, 0.25f);
	g->setEmmision(color(0, 8, 0));
	shapes.push_back(g);

	cube* b = new cube(vec3(4, 3, -10), 0.25f, 50, 0.25f);
	b->setEmmision(color(0, 0, 8));
	shapes.push_back(b);

	circle* ball1 = new circle(newray(0, -1.5, -11), 1.5f);
	ball1->setColor(color(0.8f, 0.8f, 0.8f));
	ball1->setSpecularColor(color(0.2f, 0.2f, 0.2f));
	ball1->setSpecularTerm(100);
	shapes.push_back(ball1);*/

	lightSamplePos = newray(0, 19.0f, -11);
	lightSampleRad = 14.2f;

	auto start = chrono::high_resolution_clock::now();

	//#pragma omp parallel for schedule(dynamic)
	for (int i = prevSamples + 1; i <= NUM_SAMPLES; ++i)
	{
		//lightSamplePos = randRayInSphere() * lightSampleRad + light->pos;

		globalSampleNum = i - 1;
		cout << "sample " << i << " of " << NUM_SAMPLES << endl;
		for (int y = -imageHeight / 2; y < imageHeight / 2; ++y)
		{
			for (float x = -imageWidth / 2; x < imageWidth / 2; ++x)
			{
				color totalColor = buffer[(int)(x + imageWidth / 2)][(int)(y + imageHeight / 2)];

				/*float rx = 1.0f * ((float)rand() / RAND_MAX - 0.5f) / (imageWidth),
				ry = 1.0f * ((float)rand() / RAND_MAX - 0.5f) / (imageHeight);*/
				float rx = 0, ry = 0;
				globalX = (int)(x + imageWidth / 2);
				globalY = (int)(y + imageHeight / 2);
				vec3 ray = newray(x / 2.0f / imageWidth + rx, y / 2.0f / imageWidth + ry - 0.055f, -0.5f);

				color sampleColor = sample(newray(0, 0, 1), ray, NUM_BOUNCES);

				totalColor = totalColor + sampleColor * (1.0f / NUM_SAMPLES);

				buffer[(int)(x + imageWidth / 2)][(int)(y + imageHeight / 2)] = totalColor;
			}
		}
	}

	for (int y = imageHeight - 1; y >= 0; --y)
		for (int x = 0; x < imageWidth; ++x)
		{
			buffer[x][y] = buffer[x][y].normalized();
		}

	auto end = chrono::high_resolution_clock::now();

	cout << "Elapsed time: " << chrono::duration_cast<chrono::hours>(end - start).count() << " hours" << endl;
	cout << "Elapsed time: " << chrono::duration_cast<chrono::minutes>(end - start).count() << " minutes" << endl;
	cout << "Elapsed time: " << chrono::duration_cast<chrono::seconds>(end - start).count() << " seconds" << endl;
	cout << "Elapsed time: " << chrono::duration_cast<chrono::milliseconds>(end - start).count() << " milliseconds" << endl;
	cout << "     Seconds per pixel: " << (float)chrono::duration_cast<chrono::seconds>(end - start).count() / (imageWidth*imageHeight) << endl;
	cout << "Milliseconds per pixel: " << (float)chrono::duration_cast<chrono::milliseconds>(end - start).count() / (imageWidth*imageHeight) << endl;
	cout << "Microseconds per pixel: " << (float)chrono::duration_cast<chrono::microseconds>(end - start).count() / (imageWidth*imageHeight) << endl;

	ofstream file("image2.ppm");
	file << "P3 " << imageWidth << " " << imageHeight << " 255" << endl;

	for (int y = imageHeight - 1; y >= 0; --y)
		for (int x = 0; x < imageWidth; ++x)
		{
			color celLight = buffer[x][y];
			celLight.r = (int)(celLight.r * 3 + 0.5) / 3.f;
			celLight.g = (int)(celLight.g * 3 + 0.5) / 3.f;
			celLight.b = (int)(celLight.b * 3 + 0.5) / 3.f;
			buffer[x][y] = diffBuffer[x][y].mul(celLight);
		}

	for (int y = imageHeight - 1; y >= 0; --y)
		for (int x = 0; x < imageWidth; ++x)
			file << (int)(buffer[x][y].r * 255) << " " << (int)(buffer[x][y].g * 255) << " " << (int)(buffer[x][y].b * 255) << " ";

	file.close();

	cout << "Finished" << endl;
	cin.get();

	return 0;
}
