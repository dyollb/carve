// ================================================================================ //
//						** Multi-domain preprocessing **							//
//																					//
// Multidomain preprocessing based on Carve lib.									//
// ================================================================================ //
/*!
*	@file	VkiMultiDomainPreprocessing.inl
*	@author	Alessandro Alaia (mail: alaia@itis.swiss.ch)
*	@brief	Implementation of template methods and functions
*/
# pragma once

namespace XCore { namespace Modeling
{
	// ============================================================================ //
	// IMPLEMENTATION OF EXTRACTION OPERATORS										//
	// ============================================================================ //

	// ---------------------------------------------------------------------------- //
	template< class RuleT >
	bool	ExtractFaceSubSet(
		const carve::mesh::MeshSet<3>	&mesh_in,
		carve::mesh::MeshSet<3>			*&mesh_out,
		const RuleT						&rule /*= [](const carve::mesh::MeshSet<3>::face_t *) { return true; }*/
	) {
		// Reset the output mesh.
		if (mesh_out)
		{
			delete mesh_out;
			mesh_out = nullptr;
		}
		
		// Typedef(s)
		typedef carve::mesh::MeshSet<3>::vertex_t::vector_t vec3;

		// Scope variables
		std::vector<int>	face_ids;
		std::vector<vec3>	verts;
		size_t				ncells = 0, 
							nverts = mesh_in.vertex_storage.size();

		// Initializations
		for ( const auto &m : mesh_in.meshes )
			ncells += m->faces.size();
		face_ids.reserve( 4*ncells );

		// Extract vertexes
		{
			verts.reserve( nverts );
			for ( const auto &v : mesh_in.vertex_storage )
				verts.push_back( v.v );
		}

		// Extract faces
		{
			// Scope variables
			auto	f = mesh_in.faceBegin(), fe = mesh_in.faceEnd();
			auto	VBase = mesh_in.vertex_storage.data();

			// Initializations
			ncells = 0;

			// Loop over faces in the input mesh
			for ( ; f != fe; ++f )
			{
				// Check if face has to be extracted
				if ( rule(*f) )
				{
					// Scope variables
					auto	edge = (*f)->edge;
					auto	nv = (*f)->nVertices();

					face_ids.push_back( (int) nv );
					for ( auto j = 0; j < nv; ++j )
					{
						face_ids.push_back( edge->v1() - VBase );
						edge = edge->next;
					} //next j
					++ncells;
				}
			} //next f
		}

		// Create a new carve mesh and populate with the extracted faces
		mesh_out = new carve::mesh::MeshSet<3>( verts, static_cast<size_t>( ncells ), face_ids );

		return true;
	}

	// ---------------------------------------------------------------------------- //
	template< class RuleT, class TagT >
	bool	ExtractTagSubSet(
		const carve::mesh::MeshSet<3>				&mesh_in,
		const carve::mesh::MeshSet<3>				&mesh_out,
		const carve::interpolate::FaceAttr< TagT >	&tags_in,
		carve::interpolate::FaceAttr< TagT >		*&tags_out,
		const RuleT									&rule /*= []( const carve::mesh::MeshSet<3>::face_t * ){ return true; }*/
	) {
		// Check inputs
		if ( tags_out )
		{
			delete tags_out;
			tags_out = nullptr;
		}
		if ( !tags_out )
			tags_out = new carve::interpolate::FaceAttr< TagT >();

		// Loop over faces of the first mesh amd extract tags of selected faces
		auto f = mesh_in.faceBegin(), fe = mesh_in.faceEnd(), g = mesh_out.faceBegin();
		for ( ; f != fe; ++f )
		{
			if ( rule( *f ) )
				tags_out->setAttribute( *(g++), const_cast< carve::interpolate::FaceAttr< TagT >& >( tags_in ).getAttribute( *f ) );
		} //next f

		return true;
	}

} }