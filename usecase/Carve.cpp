#include "Carve.h"
#include "UndefMinMax.h"
#include "XCoreModelingApi.h"

//#include "FloatingPointExceptions.h"

#include "QGenTriangleMesh.h"

#pragma warning(push)
#pragma warning(disable : 4018)
#pragma warning(disable : 4101)
#pragma warning(disable : 4267)
#pragma warning(disable : 4244)
#include <carve/config.h>
#include <carve/csg.hpp>
#include <carve/csg_triangulator.hpp>
#include <carve/input.hpp>
#include <carve/mesh.hpp>
#include <carve/tree.hpp>

#include <carve/exact.hpp>
#include <carve/geom3d.hpp>
#include <carve/rtree.hpp>
#include <carve/triangle_intersection.hpp>
#pragma warning(pop)

namespace XCore {

carve::mesh::MeshSet<3>* ToCarve(const QGenTriangleMesh& mesh)
{
	using vec3 = carve::mesh::MeshSet<3>::vertex_t::vector_t;
	carve::mesh::MeshOptions opts;
	std::vector<vec3> verts;
	verts.reserve(mesh.verts_count());
	vec3 v;
	for (auto& p : mesh.get_verts_buffer())
	{
		v = p; //tr.Transform(p);
		verts.push_back(v);
	}

	std::vector<int> face_ids;
	face_ids.reserve(mesh.tris_count() * 4);
	for (auto& t : mesh.get_tris_buffer())
	{
		face_ids.push_back(3);
		face_ids.push_back(t.inds[0]);
		face_ids.push_back(t.inds[1]);
		face_ids.push_back(t.inds[2]);
	}

	return new carve::mesh::MeshSet<3>(verts, static_cast<size_t>(mesh.tris_count()), face_ids, opts);
}

QGenTriangleMesh_ptr FromCarve(const carve::mesh::MeshSet<3>& cmesh, bool with_improvement, bool compute_normals)
{
	QGenTriangleMesh::vec3_buffer verts;
	verts.reserve(cmesh.vertex_storage.size());
	for (auto& p : cmesh.vertex_storage)
	{
		verts.emplace_back(p.v.x, p.v.y, p.v.z);
	}

	QGenTriangleMesh::tri_buffer tris;
	const carve::mesh::MeshSet<3>::vertex_t* VBase = cmesh.vertex_storage.data();
	for (auto f = cmesh.faceBegin(); f != cmesh.faceEnd(); ++f)
	{
		const carve::mesh::MeshSet<3>::face_t* face = *f;

		if (face->nVertices() == 3)
		{
			QGenTriangleMesh::tri t;
			auto* e = face->edge;
			for (int c = 0; c < 3; c++, e = e->next)
			{
				t.inds[c] = static_cast<bit32>(e->v1() - VBase);
			}
			tris.push_back(t);
			continue;
		}

		std::vector<carve::triangulate::tri_idx> result;

		std::vector<carve::mesh::MeshSet<3>::vertex_t*> vloop;
		face->getVertices(vloop);

		carve::triangulate::triangulate(
				carve::mesh::MeshSet<3>::face_t::projection_mapping(face->project),
				vloop,
				result);

		if (with_improvement)
		{
			carve::triangulate::improve(
					carve::mesh::MeshSet<3>::face_t::projection_mapping(face->project),
					vloop,
					carve::mesh::vertex_distance(),
					result);
		}

		QGenTriangleMesh::tri t;
		for (size_t i = 0; i < result.size(); ++i)
		{
			t.inds[0] = static_cast<bit32>(vloop[result[i].a] - VBase);
			t.inds[1] = static_cast<bit32>(vloop[result[i].b] - VBase);
			t.inds[2] = static_cast<bit32>(vloop[result[i].c] - VBase);
			tris.push_back(t);
		}
	}

	QGenTriangleMesh_ptr mesh(new QGenTriangleMesh);
	mesh->set_verts_buffer(verts);
	mesh->set_tris_buffer(tris);

	// Build smooth normals clusters in temporary mesh to avoid
	// storing cached data.
	if (compute_normals)
	{
		QGenTriangleMesh tmp_mesh;
		tmp_mesh.set_verts_buffer(verts);
		tmp_mesh.set_tris_buffer(tris);

		// add some smoothing to reduce memory consumption of the mesh
		//MACOStmp_mesh.select_all_tris(true);
		//MACOStmp_mesh.auto_smooth(30.f);
		//MACOStmp_mesh.select_all_tris(false);

		// Copy smooth id back to mesh
		bit32 num_tris = tmp_mesh.tris_count();
		for (bit32 i = 0; i < num_tris; i++)
		{
			//MACOSmesh->set_tri_smooth_id(i, tmp_mesh.get_tris_buffer()[i].sm_group);
		}
	}

	return mesh;
}

carve::csg::CSG::OP ToCarve(eBoolean operation)
{
	carve::csg::CSG::OP op = carve::csg::CSG::UNION;
	switch (operation)
	{
	case eBoolean::kUnion:
		op = carve::csg::CSG::UNION;
		break;
	case eBoolean::kIntersection:
		op = carve::csg::CSG::INTERSECTION;
		break;
	case eBoolean::kA_Minus_B:
		op = carve::csg::CSG::A_MINUS_B;
		break;
	case eBoolean::kB_Minus_A:
		op = carve::csg::CSG::B_MINUS_A;
		break;
	case eBoolean::kSymmetricDifference:
		op = carve::csg::CSG::SYMMETRIC_DIFFERENCE;
		break;
	default:
		break;
	}
	return op;
}

QGenTriangleMesh_ptr Boolean(const QGenTriangleMesh& a, const QGenTriangleMesh& b, eBoolean operation)
{
	bool improve_triangulation = false;

	carve::csg::CSG::OP op = ToCarve(operation);

	QGenTriangleMesh_ptr output(new QGenTriangleMesh);

	carve::mesh::MeshSet<3>* A = ToCarve(a);
	carve::mesh::MeshSet<3>* B = ToCarve(b);

	carve::csg::CSG_TreeNode* lhs = new carve::csg::CSG_PolyNode(A, true);
	carve::csg::CSG_TreeNode* rhs = new carve::csg::CSG_PolyNode(B, true);
	lhs = new carve::csg::CSG_OPNode(lhs, rhs, op, false, carve::csg::CSG::CLASSIFY_EDGE);

	if (lhs != nullptr)
	{
		std::shared_ptr<carve::csg::CSG_TreeNode> lhs_raii(lhs);
		carve::mesh::MeshSet<3>* result_mesh = nullptr;
		bool is_temp;
		try
		{
			carve::csg::CSG csg;
			if (improve_triangulation)
			{
				csg.hooks.registerHook(new carve::csg::CarveTriangulatorWithImprovement, carve::csg::CSG::Hooks::PROCESS_OUTPUT_FACE_BIT);
			}
			else
			{
				csg.hooks.registerHook(new carve::csg::CarveTriangulator, carve::csg::CSG::Hooks::PROCESS_OUTPUT_FACE_BIT);
			}
			result_mesh = lhs->eval(is_temp, csg);
		}
		catch (carve::exception e)
		{
			throw std::runtime_error("CSG failed, exception: " + e.str());
		}

		if (result_mesh)
		{
			output = FromCarve(*result_mesh);
			if (is_temp)
				delete result_mesh;
		}
	}
	return output;
}

carve::mesh::MeshSet<3>* Boolean(carve::mesh::MeshSet<3>* a, carve::mesh::MeshSet<3>* b, eBoolean operation)
{
	carve::csg::CSG::OP op = ToCarve(operation);

	carve::csg::CSG_TreeNode* lhs = new carve::csg::CSG_PolyNode(a, false);
	carve::csg::CSG_TreeNode* rhs = new carve::csg::CSG_PolyNode(b, false);
	lhs = new carve::csg::CSG_OPNode(lhs, rhs, op, false, carve::csg::CSG::CLASSIFY_EDGE);

	if (lhs != nullptr)
	{
		std::shared_ptr<carve::csg::CSG_TreeNode> lhs_raii(lhs);
		carve::mesh::MeshSet<3>* result_mesh = nullptr;
		try
		{
			bool is_temp;
			carve::csg::CSG csg;
			result_mesh = lhs->eval(is_temp, csg);
		}
		catch (carve::exception e)
		{
			throw std::runtime_error("CSG failed, exception: " + e.str());
		}

		if (result_mesh)
		{
			return result_mesh;
		}
	}
	throw std::runtime_error("CSG failed");
}

std::vector<size_t> SelfIntersect(const QGenTriangleMesh& mesh)
{
	using vec3 = carve::geom::vector<3>;

	std::vector<size_t> tri_ids;
	if (mesh.tris_count() == 0)
		return tri_ids;

	carve::mesh::MeshSet<3>* poly = ToCarve(mesh);

	using face_rtree_t = carve::geom::RTreeNode<3, carve::mesh::Face<3>*>;
	face_rtree_t* tree = face_rtree_t::construct_STR(poly->faceBegin(), poly->faceEnd(), 4, 4);
	std::shared_ptr<face_rtree_t> tree_raii(tree);

	for (carve::mesh::MeshSet<3>::face_iter f = poly->faceBegin(); f != poly->faceEnd(); ++f)
	{
		carve::mesh::MeshSet<3>::face_t* fa = *f;
		if (fa->nVertices() != 3)
		{
			continue;
		}

		vec3 tri_a[3];
		tri_a[0] = fa->edge->vert->v;
		tri_a[1] = fa->edge->next->vert->v;
		tri_a[2] = fa->edge->next->next->vert->v;

		std::vector<const carve::mesh::MeshSet<3>::face_t*> near_faces;
		tree->search(fa->getAABB(), std::back_inserter(near_faces));

		for (size_t f2 = 0; f2 < near_faces.size(); ++f2)
		{
			const carve::mesh::MeshSet<3>::face_t* fb = near_faces[f2];
			if (fb->nVertices() != 3)
			{
				continue;
			}

			if (fa >= fb)
			{
				continue;
			}

			vec3 tri_b[3];
			tri_b[0] = fb->edge->vert->v;
			tri_b[1] = fb->edge->next->vert->v;
			tri_b[2] = fb->edge->next->next->vert->v;

			if (carve::geom::triangle_intersection_exact(tri_a, tri_b) == carve::geom::TR_TYPE_INT)
			{
				tri_ids.push_back(fa->id);
				tri_ids.push_back(fb->id);
				break;
			}
		}
	}

	if (tri_ids.size() > 1)
	{
		std::sort(tri_ids.begin(), tri_ids.end());
		auto it = std::unique(tri_ids.begin(), tri_ids.end());
		tri_ids.resize(std::distance(tri_ids.begin(), it));
	}
	return tri_ids;
}

} // namespace XCore
