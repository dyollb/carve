// Copyright 2006-2015 Tobias Sargeant (tobias.sargeant@gmail.com).
//
// This file is part of the Carve CSG Library (http://carve-csg.com/)
//
// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation
// files (the "Software"), to deal in the Software without
// restriction, including without limitation the rights to use, copy,
// modify, merge, publish, distribute, sublicense, and/or sell copies
// of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
// BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
// ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
// CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.


#include <carve/csg.hpp>
#include <carve/csg_triangulator.hpp>

#include "geom_draw.hpp"
#include "geometry.hpp"
#include "scene.hpp"
#include <carve/input.hpp>

#include <gloop/gl.hpp>
#include <gloop/glu.hpp>
#include <gloop/glut.hpp>

#include <ctime>
#include <fstream>
#include <random>
#include <set>
#include <string>
#include <utility>

struct TestScene : public Scene
{
	GLuint draw_list_base;
	std::vector<bool> draw_flags;

	bool key(unsigned char k, int x, int y) override
	{
		const char* t;
		static const char* l = "1234567890!@#$%^&*()";
		t = strchr(l, k);
		if (t != nullptr)
		{
			int layer = t - l;
			if (layer < draw_flags.size())
			{
				draw_flags[layer] = !draw_flags[layer];
			}
		}
		return true;
	}

	GLvoid draw() override
	{
		for (int i = 0; i < draw_flags.size(); ++i)
		{
			if (draw_flags[i])
			{
				glCallList(draw_list_base + i);
			}
		}
	}

	TestScene(int argc, char** argv, int n_dlist) : Scene(argc, argv)
	{
		draw_list_base = glGenLists(n_dlist);

		draw_flags.resize(n_dlist, false);
	}

	~TestScene() override { glDeleteLists(draw_list_base, draw_flags.size()); }
};

std::mt19937 rng;
std::normal_distribution<double> norm;

carve::geom3d::Vector randomUnitVector()
{
	carve::geom3d::Vector vec =
			carve::geom::VECTOR(norm(rng), norm(rng), norm(rng));
	vec.normalize();
	return vec;
}

int main(int argc, char** argv)
{
	TestScene* scene = new TestScene(argc, argv, 1);

	glNewList(scene->draw_list_base + 0, GL_COMPILE);

	glEnable(GL_DEPTH_TEST);
	glDepthMask(GL_TRUE);
	glEnable(GL_CULL_FACE);
	glColor4f(.7, .7, .7, 1.0);
	glEnable(GL_LIGHTING);
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);

	carve::geom::tri<3> tri(randomUnitVector() * 10.0, randomUnitVector() * 10.0,
			randomUnitVector() * 10.0);

	carve::geom::vector<3> p = randomUnitVector() * 20.0;
	carve::geom::vector<3> tp = carve::geom::closestPoint(tri, p);
	double r = carve::geom::distance(p, tp);
	carve::geom::sphere<3> sphere(p, r);

	drawTri(tri);
	drawSphere(sphere);

	glEndList();

	scene->run();

	delete scene;

	return 0;
}
