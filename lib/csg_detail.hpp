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

#pragma once

#include <carve/carve.hpp>

#include <carve/collection_types.hpp>
#include <carve/faceloop.hpp>
#include <carve/mesh.hpp>
#include <carve/polyhedron_base.hpp>

namespace carve {
namespace csg {
namespace detail {
using EdgeIntInfo = std::map<
		carve::mesh::MeshSet<3>::vertex_t*,
		std::set<std::pair<carve::mesh::MeshSet<3>::face_t*, double>>>;

using VSet = std::unordered_set<carve::mesh::MeshSet<3>::vertex_t*>;
using FSet = std::unordered_set<carve::mesh::MeshSet<3>::face_t*>;

using VSetSmall = std::set<carve::mesh::MeshSet<3>::vertex_t*>;
using V2SetSmall = std::set<csg::V2>;
using FSetSmall = std::set<carve::mesh::MeshSet<3>::face_t*>;

using VVSMap = std::unordered_map<carve::mesh::MeshSet<3>::vertex_t*, VSetSmall>;
using EIntMap = std::unordered_map<carve::mesh::MeshSet<3>::edge_t*, EdgeIntInfo>;
using FVSMap = std::unordered_map<carve::mesh::MeshSet<3>::face_t*, VSetSmall>;

using VFSMap = std::unordered_map<carve::mesh::MeshSet<3>::vertex_t*, FSetSmall>;
using FV2SMap = std::unordered_map<carve::mesh::MeshSet<3>::face_t*, V2SetSmall>;

using EVVMap = std::unordered_map<carve::mesh::MeshSet<3>::edge_t*, std::vector<carve::mesh::MeshSet<3>::vertex_t*>>;

using VEVecMap = std::unordered_map<carve::mesh::MeshSet<3>::vertex_t*, std::vector<carve::mesh::MeshSet<3>::edge_t*>>;

class LoopEdges : public std::unordered_map<V2, std::list<FaceLoop*>, hash_pair>
{
	using super = std::unordered_map<V2, std::list<FaceLoop*>, hash_pair>;

public:
	void addFaceLoop(FaceLoop* fl);
	void sortFaceLoopLists();
	void removeFaceLoop(FaceLoop* fl);
};
}
}
} // namespace carve::csg::detail

inline std::ostream& operator<<(std::ostream& o, const carve::csg::detail::FSet& s)
{
	const char* sep = "";
	for (carve::csg::detail::FSet::const_iterator i = s.begin(); i != s.end();
			 ++i)
	{
		o << sep << *i;
		sep = ",";
	}
	return o;
}
