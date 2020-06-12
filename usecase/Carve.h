/// \file Carve.h
///
/// \author Bryn Lloyd
/// \copyright 2017, IT'IS Foundation

#pragma once

#include "XCoreModeling.h"

//#include "../XCore/SharedPtr.h"

//#include "../XCoreMath/Vec3.h"

#include <vector>
#include <memory>

namespace carve { namespace mesh
{
	template<unsigned ndim> class MeshSet;
}}

//XCORE_DECL_CLASS_PTR(QGenTriangleMesh);
class QGenTriangleMesh;
using QGenTriangleMesh_ptr = std::shared_ptr<QGenTriangleMesh>;

namespace XCore
{
	/// Convert to a carve mesh
	XCOREMODELING_API carve::mesh::MeshSet<3>* ToCarve(const QGenTriangleMesh& mesh);
	
	/// Convert to QTech triangle mesh
	XCOREMODELING_API QGenTriangleMesh_ptr FromCarve(const carve::mesh::MeshSet<3>& cmesh, bool with_improvement = false, bool compute_normals = false);

	/// Compute between meshes
	enum eBoolean { kUnion, kIntersection, kA_Minus_B, kB_Minus_A, kSymmetricDifference };
	XCOREMODELING_API QGenTriangleMesh_ptr Boolean(const QGenTriangleMesh& a, const QGenTriangleMesh& b, eBoolean op);
	XCOREMODELING_API carve::mesh::MeshSet<3>* Boolean(carve::mesh::MeshSet<3>* a, carve::mesh::MeshSet<3>* b, eBoolean op);

	/// Test for self-intersections and return the intersecting triangle indices
	XCOREMODELING_API std::vector<size_t> SelfIntersect(const QGenTriangleMesh& mesh);
}
