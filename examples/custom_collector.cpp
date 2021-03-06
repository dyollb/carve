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

#include "geom_draw.hpp"
#include "geometry.hpp"
#include "scene.hpp"

#include <carve/input.hpp>

#if defined(__APPLE__)
#	include <GLUT/glut.h>
#	include <OpenGL/gl.h>
#	include <OpenGL/glu.h>
#else
#	include <GL/gl.h>
#	include <GL/glu.h>
#	include <GL/glut.h>
#endif

#include <ctime>
#include <fstream>
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

#define DIM 60

class Between : public carve::csg::CSG::Collector
{
	Between();
	Between(const Between&);
	Between& operator=(const Between&);

public:
	std::list<carve::mesh::MeshSet<3>::face_t*> faces;
	const carve::mesh::MeshSet<3>* src_a;
	const carve::mesh::MeshSet<3>* src_b;

	Between(const carve::mesh::MeshSet<3>* _src_a,
			const carve::mesh::MeshSet<3>* _src_b)
			: carve::csg::CSG::Collector(), src_a(_src_a), src_b(_src_b) {}

	~Between() override {}

	void collect(carve::csg::FaceLoopGroup* grp,
			carve::csg::CSG::Hooks& hooks) override
	{
		if (grp->face_loops.head->orig_face->mesh->meshset != src_a)
		{
			return;
		}
		if (grp->classificationAgainst(src_b->meshes[1]) != carve::csg::FACE_IN)
		{
			return;
		}
		if (grp->classificationAgainst(src_b->meshes[0]) != carve::csg::FACE_OUT)
		{
			return;
		}

		for (carve::csg::FaceLoop* f = grp->face_loops.head; f; f = f->next)
		{
			faces.push_back(
					f->orig_face->create(f->vertices.begin(), f->vertices.end(), false));
		}
	}

	carve::mesh::MeshSet<3>* done(carve::csg::CSG::Hooks& hooks) override
	{
		return new carve::mesh::MeshSet<3>(faces);
	}
};

int main(int argc, char** argv)
{
	carve::mesh::MeshSet<3>* a =
			makeTorus(30, 30, 2.0, 0.8, carve::math::Matrix::ROT(0.5, 1.0, 1.0, 1.0));

	carve::input::PolyhedronData data;

	for (int i = 0; i < DIM; i++)
	{
		double x = -3.0 + 6.0 * i / double(DIM - 1);
		for (int j = 0; j < DIM; j++)
		{
			double y = -3.0 + 6.0 * j / double(DIM - 1);
			double z =
					-1.0 +
					2.0 * cos(sqrt(x * x + y * y) * 2.0) / sqrt(1.0 + x * x + y * y);
			size_t n = data.addVertex(carve::geom::VECTOR(x, y, z));
			if (i && j)
			{
				data.addFace(n - DIM - 1, n - 1, n - DIM);
				data.addFace(n - 1, n, n - DIM);
			}
		}
	}

	for (int i = 0; i < DIM; i++)
	{
		double x = -3.0 + 6.0 * i / double(DIM - 1);
		for (int j = 0; j < DIM; j++)
		{
			double y = -3.0 + 6.0 * j / double(DIM - 1);
			double z =
					1.0 +
					2.0 * cos(sqrt(x * x + y * y) * 2.0) / sqrt(1.0 + x * x + y * y);
			size_t n = data.addVertex(carve::geom::VECTOR(x, y, z));
			if (i && j)
			{
				data.addFace(n - DIM - 1, n - 1, n - DIM);
				data.addFace(n - 1, n, n - DIM);
			}
		}
	}

	carve::mesh::MeshSet<3>* b = data.createMesh(carve::input::opts());
	CARVE_ASSERT(b->meshes.size() == 2);

	Between between_collector(a, b);
	carve::mesh::MeshSet<3>* c = carve::csg::CSG().compute(
			a, b, between_collector, nullptr, carve::csg::CSG::CLASSIFY_EDGE);

	TestScene* scene = new TestScene(argc, argv, 3);

	glNewList(scene->draw_list_base + 0, GL_COMPILE);
	drawMeshSet(a, .4, .6, .8, 1.0);
	glEndList();

	glNewList(scene->draw_list_base + 1, GL_COMPILE);
	drawMeshSet(b, .8, .6, .4, 1.0);
	glEndList();

	glNewList(scene->draw_list_base + 2, GL_COMPILE);
	drawMeshSet(c, .2, .2, .8, 1.0);
	drawMeshSetWireframe(c, -1, false, false);
	glEndList();

	scene->run();

	delete scene;

	delete a;
	delete b;
	delete c;

	return 0;
}
