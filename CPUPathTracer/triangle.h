#pragma once

#include "vec3.h"
#include "color.h"

struct triangle
{
	int v0, v1, v2;

	triangle() {
		v0 = 0;
		v1 = 0;
		v2 = 0;
	}

	triangle(int p_0, int p_1, int p_2) {
		v0 = p_0;
		v1 = p_1;
		v2 = p_2;
	}
};