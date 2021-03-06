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

#include <gtest/gtest.h>


#include <carve/carve.hpp>
#include <carve/mesh.hpp>
#include <carve/mesh_impl.hpp>
#include <carve/mesh_ops.hpp>

#include "write_ply.hpp"

#include "coords.h"

#include <memory>
#include <vector>

template<typename proj_t>
double area(carve::mesh::Edge<3>* edge, proj_t proj)
{
	double A = 0.0;
	carve::mesh::Edge<3>* e = edge;
	do
	{
		carve::geom2d::P2 p1 = proj(e->vert->v);
		carve::geom2d::P2 p2 = proj(e->next->vert->v);
		A += (p2.y + p1.y) * (p2.x - p1.x);
		e = e->next;
	} while (e != edge);
	return A / 2.0;
}

void triangulate(std::vector<carve::mesh::Vertex<3>>& vertices)
{
	std::vector<carve::mesh::Face<3>*> faces;
	std::vector<carve::mesh::Vertex<3>*> vptr;

	faces.resize(1);
	for (size_t i = 0; i < vertices.size(); ++i)
	{
		vptr.push_back(&vertices[i]);
	}
	faces[0] = new carve::mesh::Face<3>(vptr.begin(), vptr.end());
	//carve::mesh::MeshSet<3> mesh(faces);

	double a0 = area(faces[0]->edge, faces[0]->project);

	std::cerr << "AREA(LOOP): " << a0 << std::endl;

	std::vector<carve::mesh::Edge<3>*> triangles;
	//std::shared_ptr<int> triangles_guard(new int, [&triangles](int*)
	//{
	//	for (auto t: triangles)
	//		delete t;
	//});
	carve::mesh::triangulate(faces[0]->edge, faces[0]->project, std::back_inserter(triangles));
	ASSERT_EQ(triangles.size(), vertices.size() - 2);

	double a1 = 0.0;
	for (size_t i = 0; i < triangles.size(); ++i)
	{
		ASSERT_EQ(triangles[i]->loopLen(), 3);
		double a = area(triangles[i], faces[0]->project);
		// std::cerr << triangles[i]->face << " " << triangles[i]->next->face << " "
		// << triangles[i]->next->next->face << std::endl;
		ASSERT_LE(a, 0.0);
		a1 += a;
	}
	std::cerr << "AREA(TRI): " << a1 << std::endl;

	ASSERT_LE(fabs(a0 - a1), 1e-5);
}

TEST(MeshTriangulationTest, SimpleFace1)
{
	std::vector<carve::mesh::Vertex<3>> vertices;

	vertices.push_back(carve::geom::VECTOR(67.772, 49.906, 0.0));
	vertices.push_back(carve::geom::VECTOR(66.908, 48.229, 0.0));
	vertices.push_back(carve::geom::VECTOR(65.93, 46.44, 0.0));
	vertices.push_back(carve::geom::VECTOR(65.183, 45.64, 0.0));
	vertices.push_back(carve::geom::VECTOR(65.183, 41.324, 0.0));
	vertices.push_back(carve::geom::VECTOR(65.183, 42.239, 0.0));

	triangulate(vertices);
}

TEST(MeshTriangulationTest, SimpleFace2)
{
	std::vector<carve::mesh::Vertex<3>> vertices;

	vertices.push_back(carve::geom::VECTOR(67.902, 49.971, 0.0));
	vertices.push_back(carve::geom::VECTOR(63.584, 49.971, 0.0));
	vertices.push_back(carve::geom::VECTOR(62.772, 44.906, 0.0));
	vertices.push_back(carve::geom::VECTOR(61.908, 43.229, 0.0));
	vertices.push_back(carve::geom::VECTOR(60.93, 41.44, 0.0));
	vertices.push_back(carve::geom::VECTOR(60.183, 40.64, 0.0));
	vertices.push_back(carve::geom::VECTOR(60.183, 36.324, 0.0));
	vertices.push_back(carve::geom::VECTOR(60.183, 37.239, 0.0));
	vertices.push_back(carve::geom::VECTOR(61.908, 37.987, 0.0));
	vertices.push_back(carve::geom::VECTOR(63.584, 37.987, 0.0));
	vertices.push_back(carve::geom::VECTOR(65.197, 38.915, 0.0));
	vertices.push_back(carve::geom::VECTOR(67.902, 40.64, 0.0));

	triangulate(vertices);
}

TEST(MeshTriangulationTest, SimpleFace3)
{
	std::vector<carve::mesh::Vertex<3>> vertices;

	vertices.push_back(carve::geom::VECTOR(0, 0, 0));
	vertices.push_back(carve::geom::VECTOR(1, 0, 0));
	vertices.push_back(carve::geom::VECTOR(1, 1, 0));
	vertices.push_back(carve::geom::VECTOR(0, 1, 0));

	triangulate(vertices);
}

TEST(MeshTriangulationTest, SimpleFace4)
{
	std::vector<carve::mesh::Vertex<3>> vertices;

	vertices.push_back(carve::geom::VECTOR(0, 0, 0));
	vertices.push_back(carve::geom::VECTOR(1, 0, 0));
	vertices.push_back(carve::geom::VECTOR(1, 0.2, 0));
	vertices.push_back(carve::geom::VECTOR(0.2, 0.2, 0));
	vertices.push_back(carve::geom::VECTOR(0.2, 0.8, 0));
	vertices.push_back(carve::geom::VECTOR(1, 0.8, 0));
	vertices.push_back(carve::geom::VECTOR(1, 1, 0));
	vertices.push_back(carve::geom::VECTOR(0, 1, 0));

	triangulate(vertices);
}

TEST(MeshTriangulationTest, SimpleFace5)
{
	std::vector<carve::mesh::Vertex<3>> vertices;

	vertices.push_back(carve::geom::VECTOR(6.25561, 6.92795, 0.0));
	vertices.push_back(carve::geom::VECTOR(6.25561, 5.6227, 0.0));
	vertices.push_back(carve::geom::VECTOR(5, 5, 0.0));
	vertices.push_back(carve::geom::VECTOR(105, 40.4667, 0.0));
	vertices.push_back(carve::geom::VECTOR(55.6727, 69.3961, 0.0));

	triangulate(vertices);
}

TEST(MeshTriangulationTest, SimpleFace6)
{
	std::vector<carve::mesh::Vertex<3>> vertices;

	vertices.push_back(carve::geom::VECTOR(180.804, 180.005, 0.0));
	vertices.push_back(carve::geom::VECTOR(179.029, 180.005, 0.0));
	vertices.push_back(carve::geom::VECTOR(189.893, 180.005, 0.0));
	vertices.push_back(carve::geom::VECTOR(186.207, 181.794, 0.0));

	triangulate(vertices);
}

TEST(MeshTriangulationTest, Map)
{
	std::vector<carve::mesh::Vertex<3>> vertices;

	for (size_t i = 0; i < sizeof(map) / sizeof(map[0]); ++i)
	{
		vertices.push_back(carve::geom::VECTOR(map[i][0], map[i][1], 0.0));
	}

	triangulate(vertices);
}

TEST(MeshTriangulationTest, Floral)
{
	std::vector<carve::mesh::Vertex<3>> vertices;

	for (size_t i = 0; i < sizeof(floral) / sizeof(floral[0]); ++i)
	{
		vertices.push_back(carve::geom::VECTOR(floral[i][0], floral[i][1], 0.0));
	}

	triangulate(vertices);
}
