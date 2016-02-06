#pragma once

#include "color.h"
#include "material.h"
#include "vec3.h"

class shape
{
public:
	vec3 pos;
	material mat;

	virtual float intersect(vec3, vec3) { return 0; }
	virtual vec3 getNorm(vec3) { return vec3(); }
	
	void setMaterial(material m) { mat = m; }
	virtual material getMaterial(vec3) { return mat; }
};

class plane : public shape
{
public:
	plane(vec3 p, vec3 n)
	{
		pos = p;
		norm = n;
	}

	virtual float intersect(vec3 o, vec3 ray)
	{
		if (dot(ray, norm) >= 0)
			return 0;
		return (dot(pos - o, norm) / dot(ray, norm));
	}

	virtual vec3 getNorm(vec3)
	{
		return norm;
	}

	/*virtual color getDiffuseColor(vec3 pos)
	{
		if (norm.y > 0.9f)
		{
			if ((int(floor(pos.z)) % 2 == 0) ^ (int(floor(pos.x)) % 2 == 0))
				return color(0.25f, 0.15f, 0.15f);
		}
		return col;
	}*/


protected:
	vec3 norm;
};

class circle : public shape
{
public:
	circle(vec3 p, float r)
	{
		pos = p;
		emm = color();
		radius = r;
		spec = color();
		shinyness = 0;
	}

	void setEmmision(color c)
	{
		emm = c;
	}

	void setColor(color c)
	{
		diffuse = c;
	}

	void setSpecularColor(color c)
	{
		spec = c;
	}

	virtual void setSpecularTerm(float t)
	{
		shinyness = t;
	}

	virtual float getSpecularTerm(vec3)
	{
		return shinyness;
	}

	virtual float intersect(vec3 o, vec3 ray)
	{
		vec3 nray = ray.normalized();
		float a = dot(ray, ray);
		float b = 2 * dot(ray, o - pos);
		float c = dot(o - pos, o - pos) - radius*radius;
		if (b * b - 4 * a * c < 0)
			return 0;
		float t0 = (-b - sqrt(b * b - 4 * a * c)) / (2 * a),
			t1 = (-b + sqrt(b * b - 4 * a * c)) / (2 * a);

		if (t0 <= 0 && t1 >= 0)
			return t1;
		if (t1 <= 0 && t0 >= 0)
			return t0;
		return min(t0, t1);
	}

	virtual vec3 getNorm(vec3 pt)
	{
		vec3 norm = pt - pos;
		norm = norm.normalized();
		return norm;
	}

	virtual color getDiffuseColor(vec3)
	{
		/*if (diffuse.r > 0.5f)
		{
		vec3 q = vec3(-2, -1, 2).normalized();
		if (dot(q, pos - this->pos) >= 1.25f)
		return color(0.05f, 0.05f, 0.05f);
		if (1.10f <= dot(q, pos - this->pos) && dot(q, pos - this->pos) <= 1.13f)
		return color(0.05f, 0.05f, 0.05f);

		if (abs(dot(q, pos - this->pos)) <= 0.4f)
		return color(0.05f, 0.05f, 0.05f);
		if (0.6f <= dot(q, pos - this->pos) && dot(q, pos - this->pos) <= 0.65f)
		return color(0.05f, 0.05f, 0.05f);
		}*/
		return diffuse;
	}

	virtual color getSpecularColor(vec3)
	{
		return spec;
	}

	virtual color getEmmision(vec3)
	{
		/*if (diffuse.r > 0.5f)
		{
		vec3 q = vec3(-2, -1, 2).normalized();
		if (1.10f <= abs(dot(q, pos - this->pos)) && abs(dot(q, pos - this->pos)) <= 1.13f)
		return color(20, 20, 20);

		if (0.6f <= abs(dot(q, pos - this->pos)) && abs(dot(q, pos - this->pos)) <= 0.65f)
		return color(20, 20, 20);
		}*/
		return emm;
	}

protected:
	color emm;
	float radius;
	color diffuse;
	color spec;
	float shinyness;
};

float rayPlaneIntersect(vec3 ray, vec3 pos, vec3 norm)
{
	if (dot(ray, norm) == 0)
		return 0;
	return dot(pos, norm) / dot(ray, norm);
}

/*class cube : public shape
{
public:
	cube(vec3 p, float rx, float ry, float rz)
	{
		pos = p;
		emm = color();
		spec = color();
		shinyness = 0;

		radX = rx; radY = ry; radZ = rz;
	}

	void setEmmision(color c)
	{
		emm = c;
	}

	void setColor(color c)
	{
		diffuse = c;
	}

	void setSpecularColor(color c)
	{
		spec = c;
	}

	virtual void setSpecularTerm(float t)
	{
		shinyness = t;
	}

	virtual float getSpecularTerm(vec3 pos)
	{
		return shinyness;
	}

	virtual float intersect(vec3 o, vec3 ray)
	{
		vec3 cent = pos - o;
		float t0, t1;

		t0 = rayPlaneIntersect(ray, cent + vec3(-radX, 0, 0), vec3(-1, 0, 0));
		t1 = rayPlaneIntersect(ray, cent + vec3(radX, 0, 0), vec3(1, 0, 0));

		float xMin = min(t0, t1);
		float xMax = max(t0, t1);

		t0 = rayPlaneIntersect(ray, cent + vec3(0, 0, -radZ), vec3(0, 0, -1));
		t1 = rayPlaneIntersect(ray, cent + vec3(0, 0, radZ), vec3(0, 0, 1));
		float zMin = min(t0, t1);
		float zMax = max(t0, t1);

		t0 = rayPlaneIntersect(ray, cent + vec3(0, -radY, 0), vec3(0, -1, 0));
		t1 = rayPlaneIntersect(ray, cent + vec3(0, radY, 0), vec3(0, 1, 0));
		float yMin = min(t0, t1);
		float yMax = max(t0, t1);

		if (min(xMax, yMax, zMax) < max(xMin, yMin, zMin))
			return maxFloat;

		return max(xMin, yMin, zMin);
	}

	virtual vec3 getNorm(vec3 pt)
	{
		/*vec3 norm = pt - pos;
		float t = max(abs(norm.x), abs(norm.y), abs(norm.z));
		if (abs(norm.x) == t)
		{
		return norm.x / radius;
		}
		if (abs(norm.y) == t)
		{
		return norm.y / radius;
		}
		if (abs(norm.z) == t)
		{
		return norm.z / radius;
		}
		return norm;*
		return vec3(0, 0, 0);
	}

	virtual color getDiffuseColor(vec3 pos)
	{
		return diffuse;
	}

	virtual color getSpecularColor(vec3 pos)
	{
		return spec;
	}

	virtual color getEmmision(vec3 pos)
	{
		return emm;
	}

protected:
	color emm;
	color diffuse;
	color spec;
	float shinyness;

	float radX, radY, radZ;
};*/