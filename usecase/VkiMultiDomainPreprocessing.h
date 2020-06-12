// ================================================================================ //
//						** Multi-domain preprocessing **							//
//																					//
// Multidomain pre-processing based on Carve lib for multidomain tetrahedral mesh	//
// generation.																		//
// ================================================================================ //
/*!
 *	@file	VkiMultiDomainPreprocessing.h
 *	@author	Alessandro Alaia (mail: alaia@itis.swiss.ch)
 *	@brief	Multidomain pre-processing based on Carve lib for multidomain tetrahedral mesh generation (header to include in a project).
*/
# pragma once

// ================================================================================ //
// MACROS  																			//
// ================================================================================ //
#ifdef min
#undef min
#endif
#ifdef max
#undef max
#endif

// ================================================================================ //
// INCLUDES																			//
// ================================================================================ //

// XCore
# include "XCoreModeling.h"

// VTK
# include <vtkSmartPointer.h>
# include <vtkPolyData.h>
# include <vtkIntArray.h>
# include <vtkCellData.h>

//QT
# include "QGenTriangleMesh.h"


// ================================================================================ //
// FORWARD DECLARATIONS																//
// ================================================================================ //
namespace carve
{
	namespace mesh
	{
		template< unsigned >
		class MeshSet;
	}
	namespace interpolate
	{
		template< class >
		class FaceAttr;
	}
}
namespace XCore
{
	namespace Modeling
	{
		class CMultiDomainCollector;
	}
}

namespace XCore{ namespace Modeling
{

	// ============================================================================ //
	// DEFINITION OF MultiDomainPreprocessor										//
	// ============================================================================ //
	/*!
	 *	@brief Pre-processor for multi-domain geometries.
	 *
	 * 	Generate a conformal surface mesh from a set of intersecting and/or disconnected
	 * 	domains.
	*/
	class XCOREMODELING_API CMultiDomainPreprocessor
	{
		// Typedef(s) ------------------------------------------------------------- //
	public:
		// Indexing
		using DomainIDType = int;													///< domain ID type
		using PriorityScore = size_t;												///< priority score type

		// Tags & flags
		using TagType = std::pair< DomainIDType, DomainIDType >;					///< tag type
		using TagSysType = carve::interpolate::FaceAttr< TagType >;					///< type for tag container

		// Mesh types
		using MeshType = carve::mesh::MeshSet<3>;									///< mesh type
		using MeshPtrType = MeshType*;												///< mesh pointer type
		using MeshList = std::vector< MeshPtrType >;								///< list of meshes

		// Domain types
		using DomainType = MeshType;												///< domain type
		using DomainPtrType = MeshPtrType;											///< domain pointer type
		using DomainList = std::vector< DomainPtrType >;							///< domain list
		using DomainDescriptor = std::tuple< DomainPtrType, PriorityScore >;		///< domain descriptor
		using DomainDescriptorList = std::vector< DomainDescriptor >;				///< domain descriptor list

		// Flags
		using FlagType = char;														///< flag type
		using FlagSysType = carve::interpolate::FaceAttr< FlagType >;				///< flags container

		// Private
	private:
		using ErrorManager = std::string;											///< manager for error messages
		using WarningManager = std::string;											///< manager for warning messages

		// Friendships ------------------------------------------------------------ //
		friend class CMultiDomainCollector;

		// Static constant(s) ----------------------------------------------------- //
	public:
		static const DomainIDType	NULL_DOMAIN_ID;									///< null domain ID
		static const std::string	TAG_NAME;										///< name of the tag array which will be exported to vtk fileformat

		// Constructor(s) --------------------------------------------------------- //
	public:
		/*!
		 *	@brief Default constructor.
		 *
		 *	Initialize a empty multi-domain preprocessor.
		*/
		CMultiDomainPreprocessor();

		// Destructor(s) ---------------------------------------------------------- //
	public:
		/*!
		 *	@brief Release memory allocated for internal data structures.
		*/
		~CMultiDomainPreprocessor();
		
		// Configuration ---------------------------------------------------------- //
	public:
		/*!
		 *	@brief Add a new domain to the internal list of domains.
		 *
		 * 	@note Conversion between the input mesh container and carve mesh
		 * 	will be performed automatically provided that a conversion routine is made visible.
		 * 
		 *	@param [in]		domain		domain to be added to the internal list
		 *	@param [in]		priority	(default = 0) priority assigned to this domain.
		 *
		 *	@result Returns the ID assigned to the domain.
		*/
		DomainIDType		AddDomain( const vtkSmartPointer< vtkPolyData > &domain, PriorityScore priority );
		/*!
		 *	@brief Overloading of AddDomain for the case of QGenTriangleMesh.
		*/
		DomainIDType		AddDomain( const QGenTriangleMesh &domain, PriorityScore priority );
		/*!
		 *	@brief Set the priority of the domain with specified ID.
		 *
		 * 	@note If no domain with the specified ID exists, no action will be taken.	
		 *
		 *	This change will become effective after calling the Update method.
		 *
		 *	@param [in]		domainID	id of the domain
		 *	@param [in]		priority	priority to be assigned to the domain.
		*/
		void				SetPriority( DomainIDType domainID, PriorityScore priority );
		/*!
		 *	@brief Enforce user priorities even if this might results in domain loss.
		 *
		 * 	@param [in]		enforce		if (true) user's priorities will be enforced, if (false)
		 * 								a priorities will be changed (as little as possible) to
		 * 								avoid domain loss.
		*/
		void				EnforceUserPriority( bool enforce );

		// Getter(s) -------------------------------------------------------------- //
	public:
		/*!
		 *	@brief Returns (true) if the last call to Update triggered some warnings.
		*/
		bool				Warning() const;
		/*!
		 *	@brief Returns the warning messages generated after the last Update (if any).
		*/
		std::string			GetWarningMsg() const;
		/*!
		 *	@brief Returns (true) if the last call to Update triggered some errors
		*/
		bool				Error() const;
		/*!
		 *	@brief Returns the errors messages generated after the last Update (if any).
		*/
		std::string			GetErrorMsg() const;
		/*!
		 *	@brief Returns (true) if there are no domains added to the preprocessor.
		*/
		bool				Empty() const;
		/*!
		 *	@brief Returns pointer to the output domain with specified ID.
		 *
		 *	@note If Update has not been invoked yet, calling this method will trigger
		 *	an update.
		 *
		 *	@param [in]		domainID		id of the domain.
		*/
		const MeshType*		GetOutputDomain( DomainIDType domainID ) const;
		/*!
		 *	@brief Returns a pointer to the tags associated to each polygonal face of the output domain with specified ID.
		 *
		 *	@param [in]		domainID		id of the domain.
		*/
		const TagSysType*	GetOutputDomainTags( DomainIDType domainID ) const;
		/*!
		 *	@brief Returns all the domains merged into a single mesh container
		 *
		 * 	@note Duplicate nodes are not collapsed.
		 *
		 *	@note This method will trigger an automatic update, if Update has not been invoked previously.
		 *
		 *	@param [in,out] merged		on output stores all the domains merged into a single mesh entity.
		 *	@param [in,out]	tags		on output stores the tags associated with each original domain mapped onto the merged mesh
		*/
		void				GetMergedDomains( carve::mesh::MeshSet<3> *&merged, carve::interpolate::FaceAttr< TagType > *&tags ) const;
		/*!
		 *	@brief Returns the processed domains as a single vtkPolyData objects together with tags associated to each face.
		 *
		 *	@note Nodes are automatically merged during the conversion.
		 *	@note Tags are automatically added to the vtkPolyDataObject as cell data. The name of the cell dataset
		 *	is set in the static variable CMultiDomainPreprocessor::TAG_NAME
		 *
		 *	@param [in,out]	merged 		on output stores all the domains merged into a single vtkPolyData object.
		*/
		void				GetMergedDomains( vtkSmartPointer< vtkPolyData > &merged ) const;
		/*!
		 *	@brief Returns the processed domains as a single QGenTriangleMesh object together with the tags associated to each triangle.
		 *
		 *	@param [in,out]	merged 		on output stores all the domains merged into a single mesh entity.
		 *	@param [in,out]	tags		on output stores a pair of tags for each surface triangle indicating the ID of the input
		 *								domain on the left/right side of the triangle.
		*/
		void				GetMergedDomains( QGenTriangleMesh &merged, std::vector< std::array< int, 2 > > &tags ) const;

		// Modifier(s) ------------------------------------------------------------ //
	public:
		/*!
		 *	@brief Reset the preprocessor to the default state (clear the internal domain list, outputs, warnings and errors)
		*/
		void				Clear();
		/*!
		 *	@brief Process the input list of domains.
		 *
		 *	@result Returns (false) if some error occurred during the process.
		 *	To retrieve the error/warning message refer to GetErrorMsg ( and GetWarningMsg )
		*/
		bool				Update();

		// Private method(s) ------------------------------------------------------ //
	protected:
		/*!
		 *	@brief Check input for error(s).
		 *
		 * 	@result Returns (false) if errors/inconsistencies are found in the inputs.
		*/
		bool				CheckInput() const;
		/*!
		 *	@brief Check output for error(s).
		 *
		 * 	@result Returns (false) if errors/inconsistencies are found in the inputs.
		*/
		bool				CheckOutput() const;
		/*!
		 *	@brief Modify user-defined priorities based on topological properties of the domain
		 *	hierarchy.
		 *
		 * 	Domain priorities are re-assigned based on a partial inclusion hierarchy while trying
		 *	to change the priorities assigned by the user as little as possible.
		*/
		bool				AssignPrioritiesBaseOnTopology();

		// Member(s) -------------------------------------------------------------- //
	protected:
		DomainDescriptorList			m_domains;									///< list of domains to be processed
		DomainList						m_outputs;									///< list of output domains
		TagSysType						*m_tags;									///< tag for each face in the output domain
		mutable ErrorManager			m_errors;									///< list of errors
		mutable WarningManager			m_warnings;									///< list of warnings
		bool							m_is_up_to_date;							///< flag to update output domain after changes
		bool							m_enforce_user_priorities;					///< enforce user priority scheme
	};

	// ============================================================================ //
	// EXTRACTION OPERATORS                        									//
	// ============================================================================ //

	// ---------------------------------------------------------------------------- //
	/*!
	 *	@brief Populate mesh_out with a subset of mesh facets extracted from mesh_in according
	 *	to a user-defined rule.
	 *
	 *	@tparam			Rule		rule used to select faces from the input mesh
	 *								Rule can be any type ( functor, function pointer, lambda expression) implementing the following operator:
	 *								bool operator( const carve::mesh::MeshSet<3>::face_t * ) const
	 *								which returns (true) if the mesh facets pointed by the input pointer
	 *								must be selected.
	 *
	 *	@param [in]		mesh_in		source carve mesh
	 *	@param [in,out]	mesh_out	on output stores the mesh facets extracted from mesh_in.
	 *								(all the content previously stored in mesh_out will be overwritten).
	 *	@param [in]		rule		(default = extract all faces) rule to select mesh faces.
	*/
	template< class RuleT >
	bool	ExtractFaceSubSet(
		const carve::mesh::MeshSet<3>	&mesh_in, 
		carve::mesh::MeshSet<3>			*&mesh_out, 
		const RuleT						&rule = []( const carve::mesh::MeshSet<3>::face_t * ){ return true; }
	);

	// ---------------------------------------------------------------------------- //
	/*!
	 *	@brief Populate tags_out with a subset of tags extracted from tags_in according to a user 
	 *	specified rule.
	 *
	 *	@tparam			RuleT		rule user to select faces from the 1st carve mesh
	 *								(see ExtractFaceSubSet)
	 *	@tparam			TagT		tag type associated to each face
	 *
	 *	@param [in]		mesh_in		source carve mesh
	 *	@param [in]		mesh_out	carve mesh storing the facet extracted from mesh_in using "rule"
	 *	@param [in]		tags_in		facet tags associated to the input mesh
	 *	@param [in,out]	tags_out	on output stores the tags extracted from tags_in.
	 *	@param [in]		rule		(default = extract all tags) rule used to extract faces from the first mesh
	 *								(see ExtractFaceSubSet)
	*/
	template< class RuleT, class TagT >
	bool	ExtractTagSubSet(
		const carve::mesh::MeshSet<3>				&mesh_in,
		const carve::mesh::MeshSet<3>				&mesh_out,
		const carve::interpolate::FaceAttr< TagT >	&tags_in,
		carve::interpolate::FaceAttr< TagT >		*&tags_out,
		const RuleT									&rule = []( const carve::mesh::MeshSet<3>::face_t * ){ return true; }
	);

	// ============================================================================ //
	// CONVERSION OEPRATOR(S)                      									//
	// ============================================================================ //

	// vtk <-> carve ============================================================== //

	// ---------------------------------------------------------------------------- //
	/*!
	 *	@brief Create a carve mesh object from a vtkPolyData surface mesh.
	 *
	 *	@param [in]		mesh_in		input vtkPolyDataObject.
	 *	@param [in,out]	mesh_out	on output stores the vtkPolyData object converted into carve mesh format (previous content will be overwritten).
	*/
	bool	CopyVtkPolyDataToCarveMesh( const vtkSmartPointer<vtkPolyData> &mesh_in, carve::mesh::MeshSet<3> *&mesh_out );

	// ---------------------------------------------------------------------------- //
	/*!
	 *	@brief Create a vtkPolyData object from a carve mesh object.
	 *
	 *	@param [in]		mesh_in		input carve mesh
	 *	@param [in,out]	mesh_out	on output a vtkPolyData object populated with faces extracted from the carve mesh container 
	 *								(all the previously stored in the vtkPolyData object will be overwritten).
	 *	@param [in]		tags_in		tags associated to each carve mesh facet.
	*/
	bool	CopyCarveMeshToVtkPolyData( const carve::mesh::MeshSet<3> &mesh_int, vtkSmartPointer<vtkPolyData> &mesh_out, const carve::interpolate::FaceAttr< CMultiDomainPreprocessor::TagType > &tags_in );

	// QTech <-> carve ============================================================ //

	// ---------------------------------------------------------------------------- //
	/*!
	 *	@brief Convert a QGenTriangleMesh object into a carve mesh object.
	 *
	 *	@param [in]		mesh_in		QGenTriangleMesh object
	 *	@param [in,out]	mesh_out	on output stores the input mesh converted to carve mesh format (previous content will be overwritten).
	*/
	bool	CopyQGenTriangleMeshToCarveMesh( const QGenTriangleMesh &mesh_in, carve::mesh::MeshSet<3> *&mesh_out );

	// ---------------------------------------------------------------------------- //
	/*!
	 *	@brief Converts the input carve mesh into a QGenTriangleMesh object.
	 *
	 *	@note Non-triangular faces will be skipped.
	 *	@param [in]		mesh		carve mesh
	 *	@param [in,out]	trias		on output stores the triangular faces in the input carve mesh.
	*/
	bool CopyCarveMeshToQGenTriangleMesh( const carve::mesh::MeshSet<3> &mesh_in, QGenTriangleMesh &mesh_out );
	
	// ---------------------------------------------------------------------------- //
	/*!
	 *	@brief Converts the input carve mesh into a QGenTriangleMesh object.
	 *
	 *	@note Non-triangular faces will be skipped.
	 *	@param [in]		mesh		carve mesh
	 *	@param [in,out]	trias		on output stores the triangular faces in the input carve mesh
	 *	@param [in]		tags		(default = nullptr) pointer to tags to be exported in the QGenTriangleMesh object.
	 *	@param [in]		outTags		on output, stores the tags originally associated to the carve mesh object.
	*/
	bool	CopyCarveMeshToQGenTriangleMesh( const carve::mesh::MeshSet<3> &mesh_in, QGenTriangleMesh &mesh_out, const carve::interpolate::FaceAttr< CMultiDomainPreprocessor::TagType > &tags_in, std::vector< std::array< int, 2 > > &tags_out );
	
	// ============================================================================ //
	// UTILITIES                                   									//
	// ============================================================================ //
	
	// ---------------------------------------------------------------------------- //
	/*!
	 *	@brief Given two intersecting meshes, this routines cuts each mesh with the other.
	 *
	 *	@param [in]		mesh0		input #1
	 *	@param [in]		mesh1		input #2
	 *	@param [in,out]	out0		on output, stores mesh0 cut with mesh1
	 *	@param [in,out]	out1		on output, stores mesh1 cut with mesh0
	 *
	 *	@result Returns (true) if the symmetric cut is successful.
	*/
	XCOREMODELING_API bool ImprintMeshes( const QGenTriangleMesh &mesh0, const QGenTriangleMesh &mesh1, QGenTriangleMesh &out0, QGenTriangleMesh &out1, bool symmetric );

} }

// ================================================================================ //
// TEMPLATE IMPLEMENTATIONS															//
// ================================================================================ //
# include "VkiMultiDomainPreprocessing.inl"
