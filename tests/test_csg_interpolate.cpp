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

#if defined(HAVE_CONFIG_H)
#	include <carve_config.h>
#endif

#include <carve/carve.hpp>
#include <carve/csg.hpp>
#include <carve/interpolator.hpp>

#include <ctime>
#include <fstream>
#include <set>
#include <string>
#include <utility>

#include "geom_draw.hpp"
#include "geometry.hpp"
#include "rgb.hpp"
#include "scene.hpp"

#if defined(__APPLE__)
#	include <GLUT/glut.h>
#	include <OpenGL/gl.h>
#	include <OpenGL/glu.h>
#else
#	include <GL/gl.h>
#	include <GL/glu.h>
#	include <GL/glut.h>
#endif

static inline void glVertex(const carve::geom3d::Vector& v)
{
	glVertex3f(g_scale * (v.x + g_translation.x),
			g_scale * (v.y + g_translation.y),
			g_scale * (v.z + g_translation.z));
}

static inline void glColor(const cRGBA& c)
{
	glColor4f(c.r, c.g, c.b, c.a);
}

carve::interpolate::FaceVertexAttr<cRGBA> fv_colours;

void drawColourPolyhedron(carve::mesh::MeshSet<3>* poly, float r, float g,
		float b, float a)
{
	cRGBA cdefault(r, g, b);
	glColor(cdefault);

	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glBegin(GL_TRIANGLES);

	for (carve::mesh::MeshSet<3>::face_iter i = poly->faceBegin();
			 i != poly->faceEnd(); ++i)
	{
		carve::mesh::MeshSet<3>::face_t* f = *i;
		if (f->nVertices() == 3)
		{
			glNormal3dv(f->plane.N.v);
			glColor(fv_colours.getAttribute(f, 0, cdefault));
			glVertex(f->edge->vert->v);
			glColor(fv_colours.getAttribute(f, 1, cdefault));
			glVertex(f->edge->next->vert->v);
			glColor(fv_colours.getAttribute(f, 2, cdefault));
			glVertex(f->edge->next->next->vert->v);
		}
	}
	glEnd();

	std::vector<std::pair<carve::geom3d::Vector, cRGBA>> v;
	for (carve::mesh::MeshSet<3>::face_iter i = poly->faceBegin();
			 i != poly->faceEnd(); ++i)
	{
		carve::mesh::MeshSet<3>::face_t* f = *i;
		if (f->nVertices() != 3)
		{
			v.resize(f->nVertices());
			for (carve::mesh::MeshSet<3>::face_t::edge_iter_t j = f->begin();
					 j != f->end(); ++j)
			{
				v[j.idx()] =
						std::make_pair(g_scale * (j->vert->v + g_translation),
								fv_colours.getAttribute(f, j.idx(), cdefault));
			}
			drawColourPoly(f->plane.N, v);
		}
	}
}

carve::mesh::MeshSet<3>* colourCube(
		const carve::math::Matrix& transform = carve::math::Matrix::IDENT())
{
	std::vector<carve::mesh::MeshSet<3>::vertex_t> v;
	v.push_back(carve::mesh::MeshSet<3>::vertex_t(
			transform * carve::geom::VECTOR(+1.0, +1.0, +1.0)));
	v.push_back(carve::mesh::MeshSet<3>::vertex_t(
			transform * carve::geom::VECTOR(-1.0, +1.0, +1.0)));
	v.push_back(carve::mesh::MeshSet<3>::vertex_t(
			transform * carve::geom::VECTOR(-1.0, -1.0, +1.0)));
	v.push_back(carve::mesh::MeshSet<3>::vertex_t(
			transform * carve::geom::VECTOR(+1.0, -1.0, +1.0)));
	v.push_back(carve::mesh::MeshSet<3>::vertex_t(
			transform * carve::geom::VECTOR(+1.0, +1.0, -1.0)));
	v.push_back(carve::mesh::MeshSet<3>::vertex_t(
			transform * carve::geom::VECTOR(-1.0, +1.0, -1.0)));
	v.push_back(carve::mesh::MeshSet<3>::vertex_t(
			transform * carve::geom::VECTOR(-1.0, -1.0, -1.0)));
	v.push_back(carve::mesh::MeshSet<3>::vertex_t(
			transform * carve::geom::VECTOR(+1.0, -1.0, -1.0)));

	std::vector<carve::mesh::MeshSet<3>::face_t*> faces;
	faces.reserve(6);

	faces.push_back(
			new carve::mesh::MeshSet<3>::face_t(&v[0], &v[1], &v[2], &v[3]));
	fv_colours.setAttribute(faces[0], 0, cRGBA(0, 0, 1));
	fv_colours.setAttribute(faces[0], 1, cRGBA(0, 0, 0));
	fv_colours.setAttribute(faces[0], 2, cRGBA(0, 1, 1));
	fv_colours.setAttribute(faces[0], 3, cRGBA(1, 0, 1));

	faces.push_back(
			new carve::mesh::MeshSet<3>::face_t(&v[7], &v[6], &v[5], &v[4]));
	fv_colours.setAttribute(faces[1], 0, cRGBA(0, 1, 0));
	fv_colours.setAttribute(faces[1], 1, cRGBA(0, 1, 1));
	fv_colours.setAttribute(faces[1], 2, cRGBA(0, 0, 0));
	fv_colours.setAttribute(faces[1], 3, cRGBA(1, 1, 0));

	faces.push_back(
			new carve::mesh::MeshSet<3>::face_t(&v[0], &v[4], &v[5], &v[1]));
	fv_colours.setAttribute(faces[2], 0, cRGBA(0, 1, 1));
	fv_colours.setAttribute(faces[2], 1, cRGBA(0, 1, 0));
	fv_colours.setAttribute(faces[2], 2, cRGBA(0, 0, 1));
	fv_colours.setAttribute(faces[2], 3, cRGBA(1, 1, 1));

	faces.push_back(
			new carve::mesh::MeshSet<3>::face_t(&v[1], &v[5], &v[6], &v[2]));
	fv_colours.setAttribute(faces[3], 0, cRGBA(1, 0, 0));
	fv_colours.setAttribute(faces[3], 1, cRGBA(1, 0, 1));
	fv_colours.setAttribute(faces[3], 2, cRGBA(1, 1, 0));
	fv_colours.setAttribute(faces[3], 3, cRGBA(0, 0, 0));

	faces.push_back(
			new carve::mesh::MeshSet<3>::face_t(&v[2], &v[6], &v[7], &v[3]));
	fv_colours.setAttribute(faces[4], 0, cRGBA(1, 0, 1));
	fv_colours.setAttribute(faces[4], 1, cRGBA(1, 0, 0));
	fv_colours.setAttribute(faces[4], 2, cRGBA(1, 1, 1));
	fv_colours.setAttribute(faces[4], 3, cRGBA(0, 0, 1));

	faces.push_back(
			new carve::mesh::MeshSet<3>::face_t(&v[3], &v[7], &v[4], &v[0]));
	fv_colours.setAttribute(faces[5], 0, cRGBA(1, 1, 0));
	fv_colours.setAttribute(faces[5], 1, cRGBA(1, 1, 1));
	fv_colours.setAttribute(faces[5], 2, cRGBA(1, 0, 0));
	fv_colours.setAttribute(faces[5], 3, cRGBA(0, 1, 0));

	return new carve::mesh::MeshSet<3>(faces);
}

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

int main(int argc, char** argv)
{
	installDebugHooks();

	TestScene* scene = new TestScene(argc, argv, 4);

	g_scale = 10.0;

	glNewList(scene->draw_list_base, GL_COMPILE);
	carve::mesh::MeshSet<3>* a =
			colourCube(carve::math::Matrix::ROT(.4, .2, .3, .4));
	carve::mesh::MeshSet<3>* b = makeTorus(20, 20, .9, .5);
	carve::csg::CSG csg;
	fv_colours.installHooks(csg);
	carve::mesh::MeshSet<3>* c = csg.compute(a, b, carve::csg::CSG::A_MINUS_B);
	glEndList();

	glNewList(scene->draw_list_base + 1, GL_COMPILE);
	drawColourPolyhedron(a, .6, .6, .6, 1.0);
	glEndList();

	glNewList(scene->draw_list_base + 2, GL_COMPILE);
	drawColourPolyhedron(b, .6, .6, .6, 1.0);
	glEndList();

	glNewList(scene->draw_list_base + 3, GL_COMPILE);
	drawColourPolyhedron(c, .6, .6, .6, 1.0);
	glEndList();

	scene->draw_flags[0] = false;
	scene->draw_flags[1] = false;
	scene->draw_flags[2] = false;
	scene->draw_flags[3] = true;

	scene->run();

	delete scene;

	return 0;
}
