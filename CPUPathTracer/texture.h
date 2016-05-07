#pragma once

#include <string>
#include <vector>

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

using std::string;
using std::vector;

struct texture
{
	float* data;
	int w, h;
	string filename;
};
vector<texture> textures;

color fetch(int texId, float u, float v)
{	
	//texture* tex = textures.begin() + texId;

	if (u < 0 || 1 < u ||
		v < 0 || 1 < v)
	{
		printf("THINGS AREN'T RIGHT OVER HERE\n");
		return color(1, 0, 1);
	}

	v = 1 - v; // textures upside down and shit :P

	if (v > 0.99)
		v = 0.99;
	if (u > 0.99)
		u = 0.99;

	if (texId >= textures.size())
	{
		printf("texId is larger that textures.size()\n");
		return color(0, 0, 0);
	}
	
	int y = v * textures[texId].h;
	int x = u * textures[texId].w;
	int i = 3 * (y * textures[texId].w + x);

	if (i + 2 >= textures[texId].w * textures[texId].h * 3)
	{
		printf("the fucking index is larger than the image again\n");
		printf("texId: %d, index: %d, w: %d, h: %d, w*h*3: %d, u: %f, v: %f", texId, i, textures[texId].w, textures[texId].h, textures[texId].w * textures[texId].h * 3,
			(double)u, (double)v);
		return color(0, 0, 0);
	}

	color ret;
	ret.r = textures[texId].data[i + 0];
	ret.g = textures[texId].data[i + 1];
	ret.b = textures[texId].data[i + 2];

	return ret;
}

int createSolidTexture(color c)
{
	texture t;
	t.w = 1;
	t.h = 1;
	t.filename = "";
	t.data = new float[3];
	t.data[0] = c.r;
	t.data[1] = c.g;
	t.data[2] = c.b;
	textures.push_back(t);
	return textures.size() - 1;
}

int loadTexture(string filename)
{
	if (textures.size() > 40)
		return -1;
	if (filename.size() == 0)
	{
		printf("No filename to load\n");
		return -1;
	}
	
	for (int i = 0; i < textures.size(); ++i)
	{
		if (textures[i].filename == filename &&
			textures[i].filename != "")
			return i;
	}

	int w, h, n;
	unsigned char *data = stbi_load(filename.c_str(), &w, &h, &n, 0);
	if (data == NULL)
	{
		printf("Could not find texture file: \"%s\"\n", filename.c_str());
		return -1;
	}
	printf("Loaded texture file: %s\n", filename.c_str());

	if (n < 3)
		printf("ERROR LOADING IMAGE: not enough components per pixel!\n");

	texture t;
	t.w = w;
	t.h = h;
	t.filename = filename;
	t.data = new float[w * h * 3];
	for (int i = 0; i < w * h; ++i)
	{
		float r = (float)(data[n * i + 0]) / 255.0f;
		float g = (float)(data[n * i + 1]) / 255.0f;
		float b = (float)(data[n * i + 2]) / 255.0f;

		t.data[3 * i + 0] = r;
		t.data[3 * i + 1] = g;
		t.data[3 * i + 2] = b;
	}

	textures.push_back(t);
	return textures.size() - 1;
}