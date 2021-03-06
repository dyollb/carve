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


#include <carve/carve.hpp>
#include <carve/geom2d.hpp>
#include <carve/heap.hpp>
#include <carve/mesh_ops.hpp>
#include <carve/mesh_simplify.hpp>
#include <carve/pointset.hpp>
#include <carve/poly.hpp>
#include <carve/polyline.hpp>
#include <carve/rtree.hpp>

#include <algorithm>
#include <fstream>
#include <set>
#include <string>
#include <utility>

typedef carve::mesh::MeshSet<3> meshset_t;
using mesh_t = carve::mesh::Mesh<3>;
using vertex_t = mesh_t::vertex_t;
using edge_t = mesh_t::edge_t;
using face_t = mesh_t::face_t;

#include "read_ply.hpp"
#include "write_ply.hpp"

#include "opts.hpp"

int main(int argc, char** argv)
{
	try
	{
		carve::input::Input inputs;
		readPLY(std::string(argv[1]), inputs);
		carve::mesh::MeshSet<3>* p;
		p = carve::input::Input::create<carve::mesh::MeshSet<3>>(
				*inputs.input.begin());

		carve::mesh::MeshSimplifier simplifier;

		simplifier.removeFins(p);
		simplifier.removeLowVolumeManifolds(p, 1.0);

		// p->transform(carve::geom::quantize<10,3>());
		simplifier.simplify(p, 1e-2, 1.0, M_PI / 180.0, 2e-3);
		// std::cerr << "n_flips: " << simplifier.improveMesh_conservative(p) <<
		// std::endl;

		simplifier.removeFins(p);
		simplifier.removeLowVolumeManifolds(p, 1.0);

		writePLY(std::cout, p, true);
		return 0;
	}
	catch (carve::exception& e)
	{
		std::cerr << "exception: " << e.what() << std::endl;
	}
}
