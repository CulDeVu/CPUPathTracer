#pragma once

#include <string>
#include <vector>

#include "material.h"
#include "useful.h"

#include "tiny_obj_loader.h"

#include "embree2\rtcore.h"
#include "embree2\rtcore_ray.h"

using std::string;
using std::vector;

using namespace tinyobj;

struct embVert
{
	float x, y, z, a;
};
struct embTriangle
{
	int v0, v1, v2;
};
struct uvTriangle
{
	float u0, v0,
		u1, v1,
		u2, v2;
};
struct model
{
	unsigned int geom_id;
	vector<unsigned int> indices;
	vector<float> uv;
	layeredMaterial mat;
};
vector<model*> models;
vector<model*> lights;

void addObj(RTCScene& scene, string filename, vec3 origin = vec3(), float scale = 1)
{
	printf("Loading .obj file: %s\n", filename.c_str());
	
	vector<shape_t> shapes;
	vector<material_t> materials;

	string err = LoadObj(shapes, materials, filename.c_str(), "models/");

	if (!err.empty())
	{
		printf("\n\nTINYOBJ ERROR: %s\n\n", err.c_str());
	}
	
	printf("Loaded .obj file. Transferring to Embree.\n");
	for (int i = 0; i < shapes.size(); ++i)
	{
		int mesh = rtcNewTriangleMesh(scene, 
			RTC_GEOMETRY_STATIC, 
			shapes[i].mesh.indices.size() / 3,
			shapes[i].mesh.positions.size() / 3);

		// setup vertex buffer
		embVert* verts = (embVert*)rtcMapBuffer(scene, mesh, RTC_VERTEX_BUFFER);
		for (int v = 0; v < shapes[i].mesh.positions.size() / 3; ++v)
		{
			verts[v].x = shapes[i].mesh.positions[3 * v + 0] * scale + origin.x;
			verts[v].y = shapes[i].mesh.positions[3 * v + 1] * scale + origin.y;
			verts[v].z = shapes[i].mesh.positions[3 * v + 2] * scale + origin.z;
		}
		rtcUnmapBuffer(scene, mesh, RTC_VERTEX_BUFFER);

		// setup index buffer
		embTriangle* tris = (embTriangle*)rtcMapBuffer(scene, mesh, RTC_INDEX_BUFFER);
		for (int v = 0; v < shapes[i].mesh.indices.size() / 3; ++v)
		{
			tris[v].v0 = shapes[i].mesh.indices[3 * v + 0];
			tris[v].v1 = shapes[i].mesh.indices[3 * v + 1];
			tris[v].v2 = shapes[i].mesh.indices[3 * v + 2];
		}
		rtcUnmapBuffer(scene, mesh, RTC_INDEX_BUFFER);

		float a0 = materials[shapes[i].mesh.material_ids[0]].diffuse[0];
		float a1 = materials[shapes[i].mesh.material_ids[0]].diffuse[1];
		float a2 = materials[shapes[i].mesh.material_ids[0]].diffuse[2];

		float e0 = materials[shapes[i].mesh.material_ids[0]].emission[0];
		float e1 = materials[shapes[i].mesh.material_ids[0]].emission[1];
		float e2 = materials[shapes[i].mesh.material_ids[0]].emission[2];
		
		float f0 = materials[shapes[i].mesh.material_ids[0]].specular[0];
		float f1 = materials[shapes[i].mesh.material_ids[0]].specular[1];
		float f2 = materials[shapes[i].mesh.material_ids[0]].specular[2];

		int tId = loadTexture(materials[shapes[i].mesh.material_ids[0]].diffuse_texname);
		if (tId == -1)
		{
			printf("instead, generated color texture: %f, %f, %f", (double)a0, (double)a1, (double)a2);
			tId = createSolidTexture(color(a0, a1, a2));
		}

		if (shapes[i].mesh.texcoords.size() == 0)
		{
			for (int j = 0; j < shapes[i].mesh.positions.size() / 3; ++j)
			{
				shapes[i].mesh.texcoords.push_back(0.5);
				shapes[i].mesh.texcoords.push_back(0.5);
			}
		}

		materialLayer* ml;
		//if (e0 == 0 && e1 == 0 && e2 == 0)
		{
			ml = createDiffuseLayer(tId);
		}

		// create model
		model* m = new model();
		m->geom_id = mesh;
		printf("%d\n", mesh);
		m->indices = shapes[i].mesh.indices;
		m->uv = shapes[i].mesh.texcoords;
		/*m->mat = material();
		m->mat.Albedo = color(a0, a1, a2);
		m->mat.Emmision = color(e0, e1, e2) * 1;
		m->mat.albedoTex = tId;*/
		m->mat.addLayerBottom(ml);
		models.push_back(m);
	}
}

/* adds a cube to the scene */
unsigned int addCube(RTCScene& scene_i)
{
	/* create a triangulated cube with 12 triangles and 8 vertices */
	unsigned int mesh = rtcNewTriangleMesh(scene_i, RTC_GEOMETRY_STATIC, 12, 8);

	/* set vertices */
	embVert* vertices = (embVert*)rtcMapBuffer(scene_i, mesh, RTC_VERTEX_BUFFER);
	vertices[0].x = -1; vertices[0].y = -1; vertices[0].z = -1;
	vertices[1].x = -1; vertices[1].y = -1; vertices[1].z = +1;
	vertices[2].x = -1; vertices[2].y = +1; vertices[2].z = -1;
	vertices[3].x = -1; vertices[3].y = +1; vertices[3].z = +1;
	vertices[4].x = +1; vertices[4].y = -1; vertices[4].z = -1;
	vertices[5].x = +1; vertices[5].y = -1; vertices[5].z = +1;
	vertices[6].x = +1; vertices[6].y = +1; vertices[6].z = -1;
	vertices[7].x = +1; vertices[7].y = +1; vertices[7].z = +1;
	rtcUnmapBuffer(scene_i, mesh, RTC_VERTEX_BUFFER);

	/* create triangle color array */
	//color* colors = (color*)alignedMalloc(12 * sizeof(color));
	color* colors = new color[12];

	/* set triangles and colors */
	int tri = 0;
	embTriangle* triangles = (embTriangle*)rtcMapBuffer(scene_i, mesh, RTC_INDEX_BUFFER);

	// left side
	colors[tri] = color(1, 0, 0); triangles[tri].v0 = 0; triangles[tri].v1 = 2; triangles[tri].v2 = 1; tri++;
	colors[tri] = color(1, 0, 0); triangles[tri].v0 = 1; triangles[tri].v1 = 2; triangles[tri].v2 = 3; tri++;

	// right side
	colors[tri] = color(0, 1, 0); triangles[tri].v0 = 4; triangles[tri].v1 = 5; triangles[tri].v2 = 6; tri++;
	colors[tri] = color(0, 1, 0); triangles[tri].v0 = 5; triangles[tri].v1 = 7; triangles[tri].v2 = 6; tri++;

	// bottom side
	colors[tri] = color(0.5f, 0, 0);  triangles[tri].v0 = 0; triangles[tri].v1 = 1; triangles[tri].v2 = 4; tri++;
	colors[tri] = color(0.5f, 0, 0);  triangles[tri].v0 = 1; triangles[tri].v1 = 5; triangles[tri].v2 = 4; tri++;

	// top side
	colors[tri] = color(1.0f, 0, 0);  triangles[tri].v0 = 2; triangles[tri].v1 = 6; triangles[tri].v2 = 3; tri++;
	colors[tri] = color(1.0f, 0, 0);  triangles[tri].v0 = 3; triangles[tri].v1 = 6; triangles[tri].v2 = 7; tri++;

	// front side
	colors[tri] = color(0, 0, 1); triangles[tri].v0 = 0; triangles[tri].v1 = 4; triangles[tri].v2 = 2; tri++;
	colors[tri] = color(0, 0, 1); triangles[tri].v0 = 2; triangles[tri].v1 = 4; triangles[tri].v2 = 6; tri++;

	// back side
	colors[tri] = color(1, 1, 0); triangles[tri].v0 = 1; triangles[tri].v1 = 3; triangles[tri].v2 = 5; tri++;
	colors[tri] = color(1, 1, 0); triangles[tri].v0 = 3; triangles[tri].v1 = 7; triangles[tri].v2 = 5; tri++;

	rtcUnmapBuffer(scene_i, mesh, RTC_INDEX_BUFFER);

	return mesh;
}