// ================================================================================ //
//						** Multi-domain preprocessing **							//
//																					//
// Implementation file for multi domain preprocessing.								//
// ================================================================================ //
/*!
 *	@file	VkiMultiDomainPreprocessing.cpp
 *	@author	Alessandro Alaia (mail: alaia@itis.swiss.ch)
 *	@brief	Multidomain pre-processing based on Carve lib for multidomain tetrahedral mesh generation (implementation file)
*/

// ================================================================================ //
// INCLUDES																			//
// ================================================================================ //
#include "XCoreModelingApi.h"
#include "UndefMinMax.h"

// Carve
# pragma warning( push )
# pragma warning( disable : 4267 )
# include <carve/config.h>
# include <carve/geom3d.hpp>
# include <carve/csg_triangulator.hpp>
# include <carve/aabb.hpp>
# include <carve/mesh.hpp>
# include <carve/csg.hpp>
# include <carve/interpolator.hpp>

#pragma warning( pop )

// XCore, keep this Order, carve stuff is needed below
# include "VkiMultiDomainPreprocessing.h"
# include "FrontMesher.h"

// Standard Template Library
# include <algorithm>
# include <sstream>
# include <exception>
# include <tuple>

namespace XCore { namespace Modeling
{
	// ============================================================================ //
	// IMPLEMENTATION OF CMultiDomainPreprocessor::CMultiDomainCollector			//
	// ============================================================================ //

	// Class definition =========================================================== //
	/*!
	 *	@brief Collector for multi-domain preprocessing.
	 *
	 *	Given two closed surfaces, this class imprint one mesh onto the other.
	*/
	class CMultiDomainCollector : public carve::csg::CSG::Collector 
	{
		// Typedef(s) ------------------------------------------------------------- //
	private:
		using face_t = carve::mesh::MeshSet<3>::face_t;								/// mesh face type
		using vertex_t = carve::mesh::MeshSet<3>::vertex_t;							///< mesh vertex type
		using MeshType = CMultiDomainPreprocessor::MeshType;						///< mesh type
		using MeshPtrType = CMultiDomainPreprocessor::MeshPtrType;					///< mesh pointer type
		using DomainIDType = CMultiDomainPreprocessor::DomainIDType;				///< domain ID type
		using TagType = CMultiDomainPreprocessor::TagType;							///< tag type
		using TagSysType = CMultiDomainPreprocessor::TagSysType;					///< tag management system
		using FlagType = CMultiDomainPreprocessor::FlagType;						///< flag type
		using FlagSysType = CMultiDomainPreprocessor::FlagSysType;					///< flag management system
		using DomainDescriptorList = CMultiDomainPreprocessor::DomainDescriptorList;///< domain decriptor type

		// Enum type(s) ----------------------------------------------------------- //
	private:
		enum ERASE_FACE
		{
			ERASE	= 0,															///< face will be removed from the face collection
			KEEP	= 1,															///< face will be kept in the face collection
			PENDING = 2																///< face will be processed later
		};

		// Friendships ------------------------------------------------------------ //
		friend class CMultiDomainPreprocessor;
		
		// Static member(s) ------------------------------------------------------- //
	private:
		/*!
		 *	@brief Set the default behavior of the collector in case of user-priorities are enforced.
		*/
		static void		UserPriorityEnforce( FlagType &new_flag );
		/*!
		 *	@brief Set the default behavior of the collector in case of user-priorities can be modified.
		*/
		static void		UserPriorityModify( FlagType &new_flag );		

		// Classes ---------------------------------------------------------------- //
	private:
		/*!
		 *	@brief Data structure holding the info for a mesh face.
		*/
		struct face_data_t
		{
			face_data_t( face_t* _face, const face_t* _orig_face, bool _flipped ) :
				face( _face ), 
				orig_face( _orig_face ),
				flipped( _flipped )
			{ };

			// Member(s)
			carve::mesh::MeshSet<3>::face_t* face;
			const carve::mesh::MeshSet<3>::face_t* orig_face;
			bool flipped;
		};

		// Constructor(s) --------------------------------------------------------- //
	private:
		/*!
		 *	@brief Default constructor.
		 *
		 *	Initialize a empty collector.
		*/
		CMultiDomainCollector( ) = delete;
		/*!
		 *	@brief Copy-constructor (deleted).
		 */
		CMultiDomainCollector( const CMultiDomainCollector& ) = delete;
	public:
		/*!
		 *	@brief Constructor #1.
		 *
		 *	Initialize a collector for the input meshes.
		 *
		 *	@param [in]		_src0		surface mesh enclosing domain #0
		 *	@param [in]		_src1 		surface mesh enclosing domain #1
		 *	@param [in]		_ID0		(default = CMultiDomainPreprocessor::NULL_DOMAIN_ID) id of the domain #0
		 *	@param [in]		_ID1		(default = CMultiDomainPreprocessor::NULL_DOMAIN_ID) id of the domain #1
		 *	@param [in]		_tags0		(default = nullptr) pointer to the tag manager for domain #0
		 *	@param [in]		_tags1		(default = nullptr) pointer to the tag manager for domain #1
		 *	@param [in]		_flag0		(default = nullptr) pointer to the flag manager for domain #0
		 *	@param [in]		_flag1		(default = nullptr) pointer to the flag manager for domain #1
		 *	@param [in]		_domains	pointer to domain descriptor list.
		*/
		CMultiDomainCollector( 
			const MeshType				*_src0, 
			const MeshType				*_src1, 
			DomainIDType				_ID0 = CMultiDomainPreprocessor::NULL_DOMAIN_ID,
			DomainIDType				_ID1 = CMultiDomainPreprocessor::NULL_DOMAIN_ID,
			TagSysType					*_tags0 = nullptr, 
			TagSysType					*_tags1 = nullptr,
			FlagSysType					*_flag0 = nullptr,
			FlagSysType					*_flag1 = nullptr,
			const DomainDescriptorList	*_domains = nullptr
		);

		// Destructor(s) ---------------------------------------------------------- //
	public:
		/*!
		 *	@brief Default destructor.
		*/
		~CMultiDomainCollector( ) override;

		// Operator(s) ------------------------------------------------------------ //
	private:
		/*!
		 *	@brief Copy assignment operator (deleted).
		*/
		CMultiDomainCollector& operator=( const CMultiDomainCollector& ) = delete;

		// Method(s) -------------------------------------------------------------- //
	private:
		/*!
		 *	@brief Forward a face from input #0 to output #0 together with its tag and flag.
		*/
		void			FWD0( const face_t *orig_face, const std::vector<vertex_t*> &vertices, carve::geom3d::Vector, bool, carve::csg::FaceClass face_class, carve::csg::CSG::Hooks &hooks );
		/*!
		 *	@brief Forward a face from input #1 to output #1 together with its tag and flag.
		*/
		void			FWD1( const face_t	*orig_face, const std::vector<vertex_t*> &vertices, carve::geom3d::Vector, bool, carve::csg::FaceClass face_class, carve::csg::CSG::Hooks &hooks );
		/*!
		 *	@brief Collect a processed face, forward it to the corresponding output, and adjust tags and flags.
		*/
		void			collect( const face_t *orig_face, const std::vector<vertex_t*>	&vertices, carve::geom3d::Vector normal, bool poly_a, carve::csg::FaceClass	face_class, carve::csg::CSG::Hooks &hooks );
		/*!
		 *	@brief Collect a group of faces and forward them to the corresponding output.
		*/
		void			collect( carve::csg::FaceLoopGroup* grp, carve::csg::CSG::Hooks& hooks ) override;
		/*!
		 *	@brief Finalize collection process (for debugging purposes).
		*/
		MeshType*		done( carve::csg::CSG::Hooks& hooks ) override;
		/*!
		 *	@brief Update flags for surface facets in a pending state.
		 *
		 *	@note This method is only required if the collector is allowed to modify user-priorities.
		 *	Otherwise no facets in a pending state will be found at the end of the processing step.
		*/
		void			UpdateFlags();
		
		// Setter(s) -------------------------------------------------------------- //
	public:
		/*!
		 *	@brief Set the pointer to the tag manager for domain #0.
		*/
		void			SetFlags0( FlagSysType *_flag0 );
		/*!
		 *	@brief Set the pointer to the flag manager for domain #1.
		*/
		void			SetFlags1( FlagSysType *_flag1 );
		/*!
		 *	@brief Set the pointer to the tag manager for domain #0.
		*/
		void			SetTags0( TagSysType *_tags0 );
		/*!
		 *	@brief Set the pointer to the tag manager for domain #1.
		*/
		void			SetTags1( TagSysType *_tags1 );
		/*!
		 *	@brief Set the ID for domain #0.
		*/
		void			SetID0( DomainIDType _ID0 );
		/*!
		 *	@brief Set the ID for domain #1.
		*/
		void			SetID1( DomainIDType _ID1 );
		/*!
		 *	@brief Set pointer to the domain descriptor list.
		*/
		void			SetDomains( const DomainDescriptorList *_domains );
		/*!
		 *	@brief Set pointer to the surface mesh enclosing domain #0.
		*/
		void			SetMesh0( const MeshType *_src0 );
		/*!
		 *	@brief Set pointer to the surface mesh enclosing domain #1.
		*/
		void			SetMesh1( const MeshType *_src1 );
		/*!
		 *	@brief Configure the collector depending on whether priorities assigned by the user must be enforced or not.
		 *
		 *	@param [in]	enforce		if (true) the collector will be configured in order not to alter user priorities.
		*/
		void			EnforceUserPriorities( bool enforce );

		// Getter(s) -------------------------------------------------------------- //
	public:
		/*!
		 *	@brief Returns the output for the domain #0.
		*/
		MeshType*		GetOutput0( carve::csg::CSG::Hooks hooks = carve::csg::CSG::Hooks() );
		/*!
		 *	@brief Returns the clean output for domain #0.
		 *
		 *	Only mesh facets marked as CMultiDomainCollector::ERASE_FACE::KEEP will be returned.
		*/
		MeshType*		GetCleanOutput0( carve::csg::CSG::Hooks hooks = carve::csg::CSG::Hooks() );
		/*!
		 *	@brief Returns the output for the domain #1.
		*/
		MeshType*		GetOutput1( carve::csg::CSG::Hooks hooks = carve::csg::CSG::Hooks() );
		/*!
		 *	@brief Returns the clean output for domain #1.
		 *
		 *	Only mesh facets marked as CMultiDomainCollector::ERASE_FACE::KEEP will be returned.
		*/
		MeshType*		GetCleanOutput1( carve::csg::CSG::Hooks hooks = carve::csg::CSG::Hooks() );
		/*!
		 *	@brief Finalize facet collection.
		*/
		void			Finalize();
	
		// Member(s) -------------------------------------------------------------- //
	private:
		std::list< face_data_t >	faces0, faces1;
		const MeshType				*src0, *src1;									///< pointers to input surf meshes
		DomainIDType				ID0, ID1;										///< IDs of input domains
		TagSysType					*tags0, *tags1;									///< pointers to tag managers for input domains
		FlagSysType					*flag0, *flag1;									///< pointers to flag managers for input domains
		void						(*tag_sys_cb)( FlagType& );						///< internal use only
		const DomainDescriptorList	*domains;										///< list of domain descriptors
		short						erase_new_marked_face0, erase_new_marked_face1;	///< internal use only (true, if a face in domain #0 or #1 must be erased at collection)
		bool						flag_are_up_to_date;							///< internal use only (true, if flags are up-to-date)
		bool						enforce_priorities;								///< internal use only (true, if user priorities are enforced)
	};

	// Static method(s) =========================================================== //

	// ---------------------------------------------------------------------------- //
	void CMultiDomainCollector::UserPriorityEnforce( FlagType &new_flag )
	{
		new_flag = ERASE_FACE::ERASE;
	}

	// ---------------------------------------------------------------------------- //
	void CMultiDomainCollector::UserPriorityModify( FlagType &new_flag )
	{
		new_flag = ERASE_FACE::PENDING;
	}

	// Constructor(s) ============================================================= //

	// ---------------------------------------------------------------------------- //
	CMultiDomainCollector::CMultiDomainCollector(
		const MeshType				*_src0,
		const MeshType				*_src1,
		DomainIDType				_ID0 		/*= NULL_DOMAIN_ID*/,
		DomainIDType				_ID1 		/*= NULL_DOMAIN_ID*/,
		TagSysType					*_tags0 	/*= nullptr*/,
		TagSysType					*_tags1 	/*= nullptr*/,
		FlagSysType					*_flag0 	/*= nullptr*/,
		FlagSysType					*_flag1 	/*= nullptr*/,
		const DomainDescriptorList	*_domains 	/*= nullptr*/
	) : carve::csg::CSG::Collector(),
		src0(_src0),
		src1(_src1),
		ID0(_ID0),
		ID1(_ID1),
		tags0(_tags0),
		tags1(_tags1),
		flag0(_flag0),
		flag1(_flag1),
		domains(_domains),
		erase_new_marked_face0( ERASE_FACE::KEEP ),
		erase_new_marked_face1( ERASE_FACE::KEEP ),
		flag_are_up_to_date(true),
		enforce_priorities(true)
	{ 
		if ( enforce_priorities )
			tag_sys_cb = &( this->UserPriorityEnforce );
		else
			tag_sys_cb = &( this->UserPriorityModify );
	}

	// Destructor(s) ============================================================== //

	// ---------------------------------------------------------------------------- //
	CMultiDomainCollector::~CMultiDomainCollector( ) 
	= default;

	// Method(s) ================================================================== //

	// ---------------------------------------------------------------------------- //
	void CMultiDomainCollector::FWD0( 
		const face_t					*orig_face, 
		const std::vector<vertex_t*>	&vertices, 
		carve::geom3d::Vector			/* normal */, 
		bool							/* poly_a */, 
		carve::csg::FaceClass			face_class, 
		carve::csg::CSG::Hooks			&hooks
	) {
		// Scope variables
		std::vector< face_t* >	new_faces;
		TagType					new_tag;
		FlagType				new_flag;

		// Initializations
		new_faces.reserve( 1 );
		new_faces.push_back(
			orig_face->create(
				vertices.begin(), 
				vertices.end(), false
			)
		);
		hooks.processOutputFace( new_faces, orig_face, false );
		for ( size_t i = 0; i < new_faces.size(); ++i )
		{
			faces0.push_back( face_data_t( new_faces[i], orig_face, false ) );
			new_tag = tags0->getAttribute( orig_face );
			new_flag = flag0->getAttribute( orig_face );
			if ( face_class == carve::csg::FACE_IN )
			{
				new_tag.second = ( ( new_tag.second == -1 ) || ( std::get<1>( domains->at( new_tag.second ) ) < std::get<1>( domains->at( ID1 ) ) ) ) ? ID1 : new_tag.second;
				if ( new_flag == ERASE_FACE::KEEP )
				{
					if ( std::get<1>( domains->at( ID0 ) ) > std::get<1>( domains->at( ID1 ) ) )
						new_flag = ERASE_FACE::KEEP;
					else
						(*tag_sys_cb)( new_flag );
				}
			}
			else if ( ( face_class == carve::csg::FACE_ON ) 
				   || ( face_class == carve::csg::FACE_ON_ORIENT_IN )
				   || ( face_class == carve::csg::FACE_ON_ORIENT_OUT ) )
			{
				if ( ( face_class == carve::csg::FACE_ON ) || ( face_class == carve::csg::FACE_ON_ORIENT_IN ) )
					new_tag.second = ( ( new_tag.second == -1 ) || ( std::get<1>(domains->at( new_tag.second ) ) < std::get<1>( domains->at(ID1) ) ) ) ? ID1 : new_tag.second;
				if (new_flag == ERASE_FACE::KEEP)
				{
					if ( std::get<1>(domains->at( ID0 ) ) < std::get<1>(domains->at( ID1 ) ) ) 
						new_flag = ERASE_FACE::ERASE;
				}
			}
			else if ( face_class == carve::csg::FACE_OUT )
				erase_new_marked_face0 = ERASE_FACE::ERASE;
			else
				throw std::runtime_error( "** ERROR domain #0 ** face classification " );
			tags0->setAttribute( faces0.back().face, new_tag );
			flag0->setAttribute( faces0.back().face, new_flag );
		}

		#if defined( CARVE_DEBUG ) && defined( DEBUG_PRINT_RESULT_FACES )
			std::cerr << "+" << ENUM(face_class) << " ";
			for ( unsigned i = 0; i < vertices.size(); ++i )
				std::cerr << " " << vertices[i] << ":" << *vertices[i];
			std::cerr << std::endl;
		#endif
	}

	// ---------------------------------------------------------------------------- //
	void CMultiDomainCollector::FWD1(
		const face_t					*orig_face, 
		const std::vector<vertex_t*>	&vertices, 
		carve::geom3d::Vector			/* normal */, 
		bool							/* poly_a */, 
		carve::csg::FaceClass			face_class, 
		carve::csg::CSG::Hooks			&hooks
	) {
		// Scope variables
		std::vector< face_t* >	new_faces;
		TagType					new_tag;
		FlagType				new_flag;

		// Initializations
		new_faces.reserve( 1 );
		new_faces.push_back(
			orig_face->create(
				vertices.begin(), 
				vertices.end(), false
			)
		);
		hooks.processOutputFace( new_faces, orig_face, false );
		for ( size_t i = 0; i < new_faces.size(); ++i )
		{
			faces1.push_back( face_data_t( new_faces[i], orig_face, false ) );
			new_tag = tags1->getAttribute( orig_face );
			new_flag = flag1->getAttribute( orig_face );
			if ( face_class == carve::csg::FACE_IN )
			{
				new_tag.second = ( ( new_tag.second == -1 ) || ( std::get<1>( domains->at( new_tag.second ) ) < std::get<1>( domains->at( ID0 ) ) ) ) ? ID0 : new_tag.second ;
				if ( new_flag == ERASE_FACE::KEEP )
				{
					if ( std::get<1>( domains->at( ID0 ) ) < std::get<1>( domains->at( ID1 ) ) )
						new_flag = ERASE_FACE::KEEP;
					(*tag_sys_cb )( new_flag );
				}
			}
			else if ( ( face_class == carve::csg::FACE_ON ) 
				   || ( face_class == carve::csg::FACE_ON_ORIENT_IN )
				   || ( face_class == carve::csg::FACE_ON_ORIENT_OUT ) )
			{
				if ( ( face_class == carve::csg::FACE_ON ) || ( face_class == carve::csg::FACE_ON_ORIENT_IN ) )
					new_tag.second = ( ( new_tag.second == -1 ) || ( std::get<1>( domains->at( new_tag.second ) ) < std::get<1>( domains->at( ID0 ) ) ) ) ? ID0 : new_tag.second ;
				if ( new_flag == ERASE_FACE::KEEP )
				{
					if ( std::get<1>( domains->at( ID1 ) ) < std::get<1>( domains->at( ID0 ) ) )
						new_flag = ERASE_FACE::ERASE;
				}
			}
			else if ( face_class == carve::csg::FACE_OUT )
				erase_new_marked_face1 = ERASE_FACE::ERASE;
			else
				std::cout << "** ERROR domain #1 ** face classification: " << face_class << std::endl;
			tags1->setAttribute( faces1.back().face, new_tag );
			flag1->setAttribute( faces1.back().face, new_flag );
		}

		#if defined( CARVE_DEBUG ) && defined( DEBUG_PRINT_RESULT_FACES )
			std::cerr << "+" << ENUM(face_class) << " ";
			for ( unsigned i = 0; i < vertices.size(); ++i )
				std::cerr << " " << vertices[i] << ":" << *vertices[i];
			std::cerr << std::endl;
		#endif
	}

	// ---------------------------------------------------------------------------- //
	void CMultiDomainCollector::collect(
		const face_t					*orig_face,
		const std::vector<vertex_t*>	&vertices,
		carve::geom3d::Vector			normal, 
		bool							poly_a, 
		carve::csg::FaceClass			face_class,
		carve::csg::CSG::Hooks			&hooks
	)  {
		// Invalidate current flag.
		flag_are_up_to_date = false;
		
		// case 0: face belongs to domain #0 and is forwarded to output #0
		if ( poly_a )
			FWD0( orig_face, vertices, normal, poly_a, face_class, hooks );
		
		// case 1: face belongs to domain #1 and is forwarded to output #1
		else
			FWD1( orig_face, vertices, normal, poly_a, face_class, hooks );
	}
	
	// ---------------------------------------------------------------------------- //
	void CMultiDomainCollector::collect(
		carve::csg::FaceLoopGroup* grp, 
		carve::csg::CSG::Hooks& hooks
	) {

		std::list< carve::csg::ClassificationInfo>& cinfo = ( grp->classification );

		if ( cinfo.empty() ) {
			std::cerr << "WARNING! group " << grp << " has no classification info!"
				<< std::endl;
			return;
		}

		carve::csg::FaceClass fc = carve::csg::FACE_UNCLASSIFIED;

		unsigned fc_closed_bits = 0;
		unsigned fc_open_bits = 0;
		unsigned fc_bits = 0;

		for ( 
			std::list< carve::csg::ClassificationInfo >::const_iterator i = grp->classification.begin(), e = grp->classification.end();
			i != e; 
			++i
		) {
			if ( (*i).intersected_mesh == nullptr ) 
			{
				fc_closed_bits = class_to_class_bit((*i).classification);
				break;
			}

			if ( (*i).classification == carve::csg::FACE_UNCLASSIFIED )
				continue;
			if ( (*i).intersectedMeshIsClosed() )
				fc_closed_bits |= carve::csg::class_to_class_bit( (*i).classification );
			else
				fc_open_bits |= carve::csg::class_to_class_bit( (*i).classification );
		}

		if ( fc_closed_bits )
			fc_bits = fc_closed_bits;
		else
			fc_bits = fc_open_bits;

		fc = carve::csg::class_bit_to_class(fc_bits);

		// handle the complex cases where a group is classified differently with
		// respect to two or more closed manifolds.
		if ( fc == carve::csg::FACE_UNCLASSIFIED )
		{
			unsigned inout_bits = fc_bits & carve::csg::FACE_NOT_ON_BIT;
			unsigned on_bits = fc_bits & carve::csg::FACE_ON_BIT;

			// both in and out. indicates an invalid manifold embedding.
			if ( inout_bits == ( carve::csg::FACE_IN_BIT | carve::csg::FACE_OUT_BIT ) )
				goto out;

			// on, both orientations. could be caused by two manifolds touching at a
			// face.
			if (on_bits == ( carve::csg::FACE_ON_ORIENT_IN_BIT | carve::csg::FACE_ON_ORIENT_OUT_BIT ) )
				goto out;

			// in or out, but also on (with orientation). the on classification takes
			// precedence.
			fc = carve::csg::class_bit_to_class(on_bits);
		}

	out:

		if ( fc == carve::csg::FACE_UNCLASSIFIED )
		{
			std::cerr << "group " << grp << " is unclassified!" << std::endl;

			#if defined(CARVE_DEBUG_WRITE_PLY_DATA)
				static int uc_count = 0;

				std::vector<carve::mesh::MeshSet<3>::face_t*> faces;

				for (FaceLoop* f = grp->face_loops.head; f; f = f->next) {
					carve::mesh::MeshSet<3>::face_t* temp =
						f->orig_face->create(f->vertices.begin(), f->vertices.end(), false);
					faces.push_back(temp);
				}

				carve::mesh::MeshSet<3>* p = new carve::mesh::MeshSet<3>(faces);

				std::ostringstream filename;
				filename << "classifier_fail_" << ++uc_count << ".ply";
				std::string out(filename.str().c_str());
				::writePLY(out, p, false);

				delete p;
			#endif

			return;
		}

		bool is_poly_a = grp->src == src0;

		for ( carve::csg::FaceLoop* f = grp->face_loops.head; f; f = f->next )
		{
			collect(
				f->orig_face, 
				f->vertices, 
				f->orig_face->plane.N, 
				is_poly_a, 
				fc,
				hooks
			);
		}
	}

	// ---------------------------------------------------------------------------- //
	CMultiDomainPreprocessor::MeshPtrType CMultiDomainCollector::done(
		carve::csg::CSG::Hooks& hooks
	) {
		return nullptr;
	}

	// ---------------------------------------------------------------------------- //
	void CMultiDomainCollector::UpdateFlags()
	{
		// faces0
		for ( auto f = faces0.begin(); f != faces0.end(); ++f )
		{
			if ( flag0->getAttribute( f->face ) == ERASE_FACE::PENDING )
				flag0->setAttribute( f->face, erase_new_marked_face0 );
		}
		// faces1
		for ( auto f = faces1.begin(); f != faces1.end(); ++f )
		{
			if ( flag1->getAttribute( f->face ) == ERASE_FACE::PENDING )
				flag1->setAttribute( f->face, erase_new_marked_face1 );
		}
		flag_are_up_to_date = true;
	}

	// Setter(s) ================================================================== //
	
	// ---------------------------------------------------------------------------- //
	void CMultiDomainCollector::SetFlags0( FlagSysType *_flag0 )
	{
		flag0 = _flag0;
	}
	
	// ---------------------------------------------------------------------------- //
	void CMultiDomainCollector::SetFlags1(FlagSysType *_flag1)
	{
		flag1 = _flag1;
	}
	
	// ---------------------------------------------------------------------------- //
	void CMultiDomainCollector::SetTags0(TagSysType *_tags0)
	{
		tags0 = _tags0;
	}

	// ---------------------------------------------------------------------------- //
	void CMultiDomainCollector::SetTags1( TagSysType *_tags1 )
	{
		tags1 = _tags1;
	}
			
	// ---------------------------------------------------------------------------- //
	void CMultiDomainCollector::SetID0( DomainIDType _ID0)
	{
		ID0 = _ID0;
	}
	
	// ---------------------------------------------------------------------------- //
	void CMultiDomainCollector::SetID1( DomainIDType _ID1 )
	{
		ID1 = _ID1;
	}
	
	// ---------------------------------------------------------------------------- //
	void CMultiDomainCollector::SetDomains( const DomainDescriptorList *_domains )
	{
		domains = _domains;
	}

	// ---------------------------------------------------------------------------- //
	void CMultiDomainCollector::SetMesh0( const MeshType *_src0 )
	{
		src0 = _src0;
	}

	// ---------------------------------------------------------------------------- //
	void CMultiDomainCollector::SetMesh1( const MeshType *_src1 )
	{
		src1 = _src1;
	}

	// ---------------------------------------------------------------------------- //
	void CMultiDomainCollector::EnforceUserPriorities( bool enforce )
	{
		enforce_priorities = enforce;
	}
	
	// Getter(s) ================================================================== //
	
	// ---------------------------------------------------------------------------- //
	CMultiDomainPreprocessor::MeshType* CMultiDomainCollector::GetOutput0( carve::csg::CSG::Hooks hooks /*= carve::csg::CSG::Hooks()*/ )
	{
		std::vector< face_t* > f;
		f.reserve( faces0.size() );
		for ( std::list<face_data_t>::iterator i = faces0.begin(); i != faces0.end(); ++i )
			f.push_back( (*i).face );

		MeshType* p = new MeshType(f);

		if ( hooks.hasHook( carve::csg::CSG::Hooks::RESULT_FACE_HOOK ) )
		{
			for ( std::list< face_data_t >::iterator i = faces0.begin(); i != faces0.end(); ++i )
				hooks.resultFace((*i).face, (*i).orig_face, (*i).flipped);
		}

		return p;
	}
	
	// ---------------------------------------------------------------------------- //
	CMultiDomainPreprocessor::MeshType* CMultiDomainCollector::GetCleanOutput0( carve::csg::CSG::Hooks hooks /*= carve::csg::CSG::Hooks()*/ )
	{
		// If user's priorities are enforced then there are no faces in a PENDING state.
		// Therefore, calling UpdateFlags is not necessary and we can spare a loop on mesh faces.
		if ( ( !flag_are_up_to_date ) && !enforce_priorities )
			UpdateFlags();

		std::vector< face_t* > f;
		f.reserve( faces0.size() );
		for ( std::list<face_data_t>::iterator i = faces0.begin(); i != faces0.end(); ++i )
		{
			if ( flag0->getAttribute( i->face ) == ERASE_FACE::KEEP )
				f.push_back( (*i).face );
		}

		MeshType* p = new MeshType(f);

		if ( hooks.hasHook( carve::csg::CSG::Hooks::RESULT_FACE_HOOK ) )
		{
			for ( std::list< face_data_t >::iterator i = faces0.begin(); i != faces0.end(); ++i )
				hooks.resultFace((*i).face, (*i).orig_face, (*i).flipped);
		}

		return p;
	}
	
	// ---------------------------------------------------------------------------- //
	CMultiDomainPreprocessor::MeshType* CMultiDomainCollector::GetOutput1( carve::csg::CSG::Hooks hooks /*= carve::csg::CSG::Hooks()*/ )
	{
		// Scope variables
		std::vector< face_t* > f;

		f.reserve( faces1.size() );
		for ( std::list< face_data_t >::iterator i = faces1.begin(); i != faces1.end(); ++i )
			f.push_back( (*i).face );

		MeshType* p = new MeshType( f );

		if ( hooks.hasHook( carve::csg::CSG::Hooks::RESULT_FACE_HOOK ) )
		{
			for ( std::list< face_data_t >::iterator i = faces1.begin(); i != faces1.end(); ++i )
				hooks.resultFace((*i).face, (*i).orig_face, (*i).flipped);
		}

		return p;
	}
	
	// ---------------------------------------------------------------------------- //
	CMultiDomainPreprocessor::MeshType* CMultiDomainCollector::GetCleanOutput1( carve::csg::CSG::Hooks hooks /*= carve::csg::CSG::Hooks()*/ )
	{
		// If user's priorities are enforced then there are no faces in a PENDING state.
		// Therefore, calling UpdateFlags is not necessary and we can spare 2 loops on the mesh faces.
		if ( ( !flag_are_up_to_date ) && !enforce_priorities )
			UpdateFlags();

		// Scope variables
		std::vector< face_t* > f;

		f.reserve( faces1.size() );
		for ( std::list< face_data_t >::iterator i = faces1.begin(); i != faces1.end(); ++i )
		{
			if ( flag1->getAttribute( i->face ) == ERASE_FACE::KEEP )
				f.push_back( (*i).face );
		}
		MeshType* p = new MeshType( f );

		if ( hooks.hasHook( carve::csg::CSG::Hooks::RESULT_FACE_HOOK ) )
		{
			for ( std::list< face_data_t >::iterator i = faces1.begin(); i != faces1.end(); ++i )
				hooks.resultFace( (*i).face, (*i).orig_face, (*i).flipped );
		}

		return p;
	}
	
	// ---------------------------------------------------------------------------- //
	void CMultiDomainCollector::Finalize()
	{
		if ( !enforce_priorities )
			UpdateFlags();
	}
		
	// ============================================================================ //
	// IMPLEMENTATION OF METHODS FOR CLASS CMultiDomainPreprocessor					//
	// ============================================================================ //

	// Static constant(s) ========================================================= //

	// ---------------------------------------------------------------------------- //
	const std::string CMultiDomainPreprocessor::TAG_NAME = "domain_tag";
	const CMultiDomainPreprocessor::DomainIDType CMultiDomainPreprocessor::NULL_DOMAIN_ID = -1;
	
	// Constructor(s) ============================================================= //

	// ---------------------------------------------------------------------------- //
	CMultiDomainPreprocessor::CMultiDomainPreprocessor() 
		
	= default;

	
	// Destructor(s) ============================================================== //

	// ---------------------------------------------------------------------------- //
	CMultiDomainPreprocessor::~CMultiDomainPreprocessor()
	{
		// Clear the list of input/output domains and release memory
		for (size_t i = 0, N = m_domains.size(); i < N; ++i)
		{
			delete std::get<0>( m_domains[i] );
			delete m_outputs.at(i);
		}
		
		// Clear tags
		delete [] m_tags;
	}
	
	// Configuration ============================================================== //

	// ---------------------------------------------------------------------------- //
	void CMultiDomainPreprocessor::SetPriority( DomainIDType domainID, PriorityScore priority )
	{
		// Check that the input ID is a valid domain ID.
		if ( ( domainID < 0 ) || ( domainID >= m_domains.size() ) )
			return;

		// Since priorities have changed, the output needs to be updated
		m_is_up_to_date = false;

		// Assign the new priority to the specified domain
		std::get<1>( m_domains.at( domainID ) ) = priority;
	}

	// ---------------------------------------------------------------------------- //
	CMultiDomainPreprocessor::DomainIDType CMultiDomainPreprocessor::AddDomain( const vtkSmartPointer< vtkPolyData > &domain, PriorityScore priority )
	{
		// Add the input domain to the internal list of domains to be processed.
		m_domains.push_back( std::make_tuple( nullptr, priority ) );
		
		// If the input vtkPolyData could not be converted into a carve mesh object, 
		// a message error is generated and stored internally.
		if ( !XCore::Modeling::CopyVtkPolyDataToCarveMesh( domain, std::get<0>( m_domains.back() ) ) )
		{
			// Generate and store the error message
			std::stringstream msg;
			msg << "CMultiDomainPreprocessor::AddDomain: ** ERROR ** Failed to convert input vtkPolyData to carve mesh container\n";
			m_errors += msg.str();
		}
		return (int) m_domains.size() - 1;
	}
	
	// ---------------------------------------------------------------------------- //
	CMultiDomainPreprocessor::DomainIDType CMultiDomainPreprocessor::AddDomain( const QGenTriangleMesh &domain, PriorityScore priority )
	{
		// Add the input domain to the internal list of domains to be processed.
		m_domains.push_back( std::make_tuple( nullptr, priority ) );
		
		// If the input QGenTriangleMesh could not be converted into a carve mesh object,
		// a message error is generated and stored internally.
		if ( !XCore::Modeling::CopyQGenTriangleMeshToCarveMesh( domain, std::get<0>( m_domains.back() ) ) )
		{
			// Generate the error message
			std::stringstream msg;
			msg << "CMultiDomainPreprocessor::AddDomain: ** ERROR ** Failed to convert input domain to carve mesh format\n";
			m_errors += msg.str();
		}
		return (int) m_domains.size() - 1;
	}
	
	// Getter(s) ================================================================== //

	// ---------------------------------------------------------------------------- //
	bool CMultiDomainPreprocessor::Warning() const
	{
		return ( m_warnings.length() != 0 );
	}
	
	// ---------------------------------------------------------------------------- //
	std::string CMultiDomainPreprocessor::GetWarningMsg() const
	{
		return m_warnings;
	}
	
	// ---------------------------------------------------------------------------- //
	bool CMultiDomainPreprocessor::Error() const
	{
		return ( m_errors.length() != 0 );
	}
	
	// ---------------------------------------------------------------------------- //
	std::string CMultiDomainPreprocessor::GetErrorMsg() const
	{
		return m_errors;
	}

	// ---------------------------------------------------------------------------- //
	bool CMultiDomainPreprocessor::Empty() const
	{
		return ( m_domains.empty() );
	}

	// ---------------------------------------------------------------------------- //
	const CMultiDomainPreprocessor::MeshType* CMultiDomainPreprocessor::GetOutputDomain( DomainIDType domainID ) const
	{
		// Check for empty preprocessor
		if ( Empty() )
			return nullptr;

		// Trigger an update if the output is not up-to-date
		if ( !m_is_up_to_date )
		{
			const_cast< CMultiDomainPreprocessor* >( this )->Update();
			if ( Error() )
				return nullptr;
		}

		return m_outputs.at( domainID );
	}

	// ---------------------------------------------------------------------------- //
	const CMultiDomainPreprocessor::TagSysType* CMultiDomainPreprocessor::GetOutputDomainTags( DomainIDType domainID ) const
	{
		// Check for empty preprocessor
		if ( Empty() )
			return nullptr;

		// Trigger an update if the output is not up-to-date
		if ( !m_is_up_to_date )
		{
			const_cast< CMultiDomainPreprocessor* >( this )->Update();
			if ( Error() )
				return nullptr;
		}

		return &( m_tags[ domainID ] );
	}
	
	// ---------------------------------------------------------------------------- //
	void CMultiDomainPreprocessor::GetMergedDomains( carve::mesh::MeshSet<3> *&merged, carve::interpolate::FaceAttr< TagType > *&tags ) const
	{
		// Check for empty preprocessor
		if ( this->Empty() )
			return;

		// Trigger an update if the output is not up-to-date
		if ( !m_is_up_to_date )
		{
			const_cast< CMultiDomainPreprocessor* >(this)->Update();
			if ( Error() )
				return;
		}

		// Clear the content of input mesh.
		if ( merged )
		{
			delete merged;
			merged = nullptr;
		}
		if ( tags )
		{
			delete tags;
			tags = nullptr;
		}

		// Typedef(s)
		typedef carve::mesh::MeshSet<3>::vertex_t::vector_t	vec3;
		using VertexList = std::vector<vec3>;
		using FaceList = std::vector<int>;
		using VertexIDMap = std::unordered_map<int, int>;

		// Scope variables
		VertexList		verts;
		FaceList		faces;
		size_t			ncells = 0;

		// Initializations
		{
			// Estimate the nr. of vertexes and faces
			size_t	nv = 0, nf = 0;
			for ( const auto &mesh : m_outputs )
			{
				nv += mesh->vertex_storage.size();
				for ( const auto &m : mesh->meshes )
					nf += m->faces.size();
			}

			// Reserve memory for vertex/faces insertion
			verts.reserve( nv );
			faces.reserve( 4*nf );
		}

		// Loop over all meshes and extract faces
		{
			// Scope variables
			size_t			domain_counter = 0;
			size_t			vertex_counter = 0;
			size_t			nv;

			// Loop over processed domains
			for ( const auto &mesh : m_outputs )
			{
				// Scope variables
				auto			VBase = mesh->vertex_storage.data();
				VertexIDMap		vmap;

				// Loop over faces in the processed domain
				auto f = mesh->faceBegin(), fe = mesh->faceEnd();
				for ( ; f != fe; ++f )
				{
					// Scope variables
					const auto	face = (*f);
					auto		edge = face->edge;

					// Face details
					nv = face->nVertices();

					// Collect face connectivity data
					faces.push_back( (int) nv );
					for ( int j = 0; j < nv; ++j, edge = edge->next )
					{
						// Scope variables
						auto vertId = edge->v1() - VBase;

						// Case 0: this vertex has not been imported yet
						if ( vmap.count( vertId ) == 0 )
						{
							// Add vertex to the vertex list
							verts.push_back( mesh->vertex_storage[vertId].v );

							// Update id mapper
							vmap[vertId] = (int) ( verts.size()-1 );

							// Update face connectivity
							faces.push_back( (int) ( verts.size()-1 ) );
						}

						// Case 1: this vertex has already been imported
						else
							faces.push_back( vmap.at( vertId ) );

					} //next j

					// Update cell counter
					++ncells;

				} //next f

			} //next mesh
		}

		// Instantiate a new mesh object
		merged	= new carve::mesh::MeshSet<3>( verts, ncells, faces );
		tags	= new carve::interpolate::FaceAttr< TagType >();

		// Collect tags
		{
			// Scope variables
			size_t	domain_counter = 0;
			auto	g = merged->faceBegin();

			// Loop over output meshes
			for ( const auto &mesh : m_outputs )
			{
				// Scope variables
				auto &tags_in = m_tags[ domain_counter++ ];
				auto f = mesh->faceBegin(), fe = mesh->faceEnd();

				// Loop over faces and collect tags
				for ( ; f != fe; ++f )
					tags->setAttribute( *(g++), tags_in.getAttribute( *f ) );

			}
		}

		return;
	}

	// ---------------------------------------------------------------------------- //
	void CMultiDomainPreprocessor::GetMergedDomains( vtkSmartPointer< vtkPolyData > &merged ) const
	{
		// Check for empty preprocessor
		if ( Empty() )
			return;

		// Initialize input
		merged = vtkSmartPointer< vtkPolyData >::New();

		// Scope variables
		carve::mesh::MeshSet<3>					*mesh 	= nullptr;
		carve::interpolate::FaceAttr< TagType > *tags 	= nullptr;
		std::stringstream						msg;

		// Extract merged domains
		this->GetMergedDomains( mesh, tags );
		if ( Error() )
		{
			msg << "CMultiDomainPreprocessor::GetMergedDomains: ** ERROR ** Failed to merge output domains\n";
			m_errors += msg.str();
			msg.str("");
			return;
		}

		// Converts to vtkPolyData
		if ( !XCore::Modeling::CopyCarveMeshToVtkPolyData( *mesh, merged, *tags ) )
		{
			msg << "CMultiDomainPreprocessor::GetMergedDomains: ** ERROR ** Failed to convert carve mesh to vtkPolyData\n";
			m_errors += msg.str();
			msg.str("");
		}
		
		// Clean up
		delete mesh;
		delete tags;
		return;
	}

	// ---------------------------------------------------------------------------- //
	void CMultiDomainPreprocessor::GetMergedDomains( QGenTriangleMesh &merged, std::vector< std::array< int, 2 > > &tags_out ) const
	{
		// Check for empty preprocessor
		if ( Empty() )
			return;

		// Scope variables
		carve::mesh::MeshSet<3>					*mesh = nullptr;
		carve::interpolate::FaceAttr< TagType >	*tags = nullptr;
		std::stringstream						msg;

		// Extract merged domains
		this->GetMergedDomains( mesh, tags );
		if ( Error() )
		{
			msg << "CMultiDomainPreprocessor::GetMergedDomains: ** ERROR ** Failed to extract merged domains\n";
			m_errors += msg.str();
			msg.str("");
			return;
		}

		// Copy mesh and tags
		if ( !XCore::Modeling::CopyCarveMeshToQGenTriangleMesh( *mesh, merged, *tags, tags_out ) )
		{
			msg << "CMultiDomainPreprocessor::GetMergedDomains: ** ERROR ** Failed to convert carve mesh to QGenTriangleMesh\n";
			m_errors += msg.str();
			msg.str("");
		}

		// Clean up
		delete mesh;
		delete tags;		
		return;
	}

	// Modifier(s) ================================================================ //

	// ---------------------------------------------------------------------------- //
	void CMultiDomainPreprocessor::Clear()
	{
		// Clear the list input domain
		for ( size_t i = 0, N= m_domains.size(); i < N; ++i )
			delete std::get<0>( m_domains[i] );
		
		for ( size_t i = 0, N= m_outputs.size(); i < N; ++i )
			delete m_outputs[i];
		
		// Clear tags
		delete [] m_tags;
		
		// Clear warnings/errors
		m_errors.clear();
		m_warnings.clear();
		
		m_is_up_to_date = true;
	}

	// ---------------------------------------------------------------------------- //
	bool CMultiDomainPreprocessor::Update()
	{
		// Clear previous error/warnings
		m_errors.clear();
		m_warnings.clear();

		// Update flags
		m_is_up_to_date = true;

		// Check for empty domain list
		if ( this->Empty() )
		{
			m_warnings += "CMultiDomainPreprocessor:Update ** WARNING ** the list of input domain is empty\n";
			return false;
		}

		// Check for errors in the inputs
		CheckInput();
		if ( Error() )
			return false;

		// Scope variables
		std::vector< DomainIDType >	domainIDs;
		std::vector< FlagSysType >	meshFlags;
		size_t						nDomains = m_domains.size();
		std::stringstream			msg;

		// Assign priorities
		{
			// Scope variables
			DomainIDType			counter = 0;

			// Adjust priorities based on topological properties of the ensemble
			if ( !m_enforce_user_priorities )
				AssignPrioritiesBaseOnTopology();

			// Fix priority conflicts
			domainIDs.reserve( nDomains );
			for ( const auto &domain : m_domains )
				domainIDs.push_back( counter++ );
			std::sort( 
				domainIDs.begin(), 
				domainIDs.end(), 
				[&]( const DomainIDType &id0, const DomainIDType &id1 )
				{
					return ( std::get<1>( m_domains.at(id0) ) < std::get<1>( m_domains.at(id1) ) );
				}
			);
			for ( size_t i = 1; i < nDomains; ++i )
			{
				DomainIDType	id0 = domainIDs.at(i-1),
								id1 = domainIDs.at(i);
				PriorityScore	&p0 = std::get<1>( m_domains.at( id0 ) ),
								&p1 = std::get<1>( m_domains.at( id1 ) );
				if ( p1 <= p0 )
				{
					msg << "CMultiDomainPreprocessor::Update ** WARNING Priority of domain " << id1 << " has been changed from " << p1 << " to " << p0 + 1 << "\n";
					m_warnings += msg.str();
					msg.str("");
					p1 = p0 + 1;
				}
			}
			
			// Sort domain based on priorities (from highest priority to the lowest)
			std::sort(
				domainIDs.begin(),
				domainIDs.end(),
				[&]( const DomainIDType &id0, const DomainIDType &id1 )
				{
					return ( std::get<1>( m_domains.at(id0) ) > std::get<1>( m_domains.at(id1) ) );
				}
			);

		}

		// Initialize meshes and tags
		{
			// Clean output from previous computations
			for ( auto &output : m_outputs )
			{
				if ( output )
					delete output;
			}
			std::vector< MeshPtrType >( nDomains, nullptr ).swap( m_outputs );
			
			// Clean tags from previous computations
			if ( m_tags )
				delete [] m_tags;
			m_tags = new TagSysType[nDomains];
			
			// Resize flags.
			meshFlags.resize( nDomains );

			// Loop on each input domain and initialize each output.
			for ( size_t i = 0; i < nDomains; ++i )
			{
				// Scope variables
				MeshType *&mesh = m_outputs.at(i);

				// Extract all the faces from the input mesh to create a copy of the i-th input domain
				XCore::Modeling::ExtractFaceSubSet( *( std::get<0>( m_domains[i] ) ), m_outputs[i], []( const DomainType::face_t * ){ return true; } );

				// Initialize tags and flags for the i-th domain
				TagSysType 	&tag = m_tags[i];
				FlagSysType &flag = meshFlags.at(i);
				for ( auto f = mesh->faceBegin(); f != mesh->faceEnd(); ++f )
				{
					tag.setAttribute( *f, std::make_pair( (int) i, NULL_DOMAIN_ID ) );
					flag.setAttribute( *f, CMultiDomainCollector::ERASE_FACE::KEEP );
				}
			} //next i

			// Quit if initialization failed
			if ( Error() )
				return false;
		}

		// Process each pair of domains (ii, jj) following the priority order
		{
			for ( size_t ii = 0; ii < nDomains; ++ii )
			{
				// Scope variables
				auto						i = domainIDs[ii];
				const DomainDescriptor		&domain0 = m_domains.at(i);
				MeshPtrType					&mesh0	 = m_outputs.at(i);
				TagSysType					&tag0	 = m_tags[i];
				FlagSysType					&flag0	 = meshFlags.at(i);

				for ( size_t jj = ii+1; jj < nDomains; ++jj )
				{

					// Scope variables
					auto						j = domainIDs[jj];
					const DomainDescriptor		&domain1 = m_domains.at(j);
					MeshPtrType					&mesh1	= m_outputs.at(j);
					TagSysType					&tag1	= m_tags[j];
					FlagSysType					&flag1  = meshFlags.at(j);

					// Check if the input domains actually intersect
					// (this check is not required, carve will already check and quit prematurely
					// if the input domains do not intersect)
					/*if ( intersect )*/
					{
						// Scope variables
						MeshType				*out0 = nullptr, 
												*out1 = nullptr;

						// Process domain pair (ii,jj)
						{
							// Scope variables
							CMultiDomainCollector	collector( mesh0, mesh1 );
							carve::csg::CSG			csg;

							// Initialize collector
							collector.SetID0( i );
							collector.SetTags0( &tag0 );
							collector.SetFlags0( &flag0 );
							collector.SetID1( j );
							collector.SetTags1( &tag1 );
							collector.SetFlags1( &flag1 );
							collector.SetDomains( &m_domains );
							collector.EnforceUserPriorities( m_enforce_user_priorities );

							// Cut each domain with the other
							csg.compute( mesh0, mesh1, collector, nullptr, carve::csg::CSG::CLASSIFY_NORMAL);
							
							// Finalize computation
							collector.Finalize();

							// Collect outputs
							out0 = collector.GetOutput0();
							out1 = collector.GetOutput1();
							
							// Check for valid output
							if ( !out0 || !out1 )
							{
								msg << "CMultiDomainPreprocessor::Update "" ERROR ** failed to process domains (" << i << ", " << j << ")\n";
								m_errors += msg.str();
								msg.str("");
								return false;
							}
						}

						// Overwrite the output domains from previous iterations with the new version
						delete mesh0;
						delete mesh1;
						mesh0 = out0;
						mesh1 = out1;
					}
				} //next j
			} //next i
		}

		// Post process output meshes
		// For each output domain, only the faces marked as CMultiDomainCollector::ERASE_FACE::KEEP are retained.
		for ( size_t i = 0; i < nDomains; ++i )
		{
			// Scope variables
			FlagSysType		&flag = meshFlags.at(i);
			TagSysType		&tags = m_tags[i];
			MeshPtrType		&mesh = m_outputs.at(i);
			MeshPtrType		out = nullptr;
			TagSysType		*tags_out_ptr = nullptr;

			// Functor(s)
			auto rule = [ &flag ]( const MeshType::face_t *face_ptr )
			{
				return ( flag.getAttribute(face_ptr) == CMultiDomainCollector::ERASE_FACE::KEEP );
			};
			
			// If extraction fails move to the next domain and leave the i-th domain untouched.
			if ( !XCore::Modeling::ExtractFaceSubSet( *mesh, out, rule )
			  || !XCore::Modeling::ExtractTagSubSet( *mesh, *out, tags, tags_out_ptr, rule ) )
			{
				msg << "CMultiDomainPreprocessor::Update ** ERROR ** Failed to retrieve domain " << i << "\n";
				m_errors += msg.str();
				msg.str("");
			}
			else
			{
				// Clear the current i-th domain
				delete mesh;
			
				// Overwrite the current i-th domain with the "cleaned" version
				mesh = out;
				
				// Copy the tags
				tags = *tags_out_ptr;

				// Clear the current i-th tags
				delete tags_out_ptr;
			}
		} //next i

		// Check Output
		CheckOutput();
		if ( Error() )
			return false;

		return true;
	}

	// ---------------------------------------------------------------------------- //
	void CMultiDomainPreprocessor::EnforceUserPriority( bool enforce )
	{
		// The priority schema has been changed so the current output might not be valid anymore.
		m_is_up_to_date = false;
		
		// Set flag value
		m_enforce_user_priorities = enforce;
	}

	// Private method(s) ========================================================== //

	// ---------------------------------------------------------------------------- //
	bool CMultiDomainPreprocessor::CheckInput() const
	{
		// Scope variables
		std::stringstream	msg;
		bool				no_errors = true;
		
		// Check input domains
		{
			for ( const auto &domain : m_domains )
			{
				DomainPtrType domain_ptr = std::get<0>(domain);

				// Pointer is not a valid memory location
				if ( !domain_ptr )
				{
					msg << "CMultiDomainPreprocessor::CheckInput ** ERROR ** domain: " << std::get<1>(domain) << " has not been initialized\n";
					m_errors += msg.str();
					msg.str("");
					no_errors = false;
				}

				// Pointer points to a empty domain
				else
				{
					if ( ( domain_ptr->meshes.empty() ) || ( domain_ptr->meshes.front()->faces.empty() ) )
					{
						msg << "CMultiDomainPreprocessor::CheckInput ** ERROR ** domain: " << std::get<1>(domain) << " is empty\n";
						m_errors += msg.str();
						msg.str("");
						no_errors = false;
					}
				}
			} //next domain
		}

		return no_errors;
	}

	// ---------------------------------------------------------------------------- //
	bool CMultiDomainPreprocessor::CheckOutput() const
	{
		// Trigger an update if outputs are not up-to-date
		if ( !m_is_up_to_date )
		{
			const_cast< CMultiDomainPreprocessor* >( this )->Update();
			if ( Error() )
				return false;
		}

		// Check domains
		{
			DomainIDType		i = 0;
			std::stringstream	msg;
			bool				no_errors = true;
			for ( const auto &output : m_outputs )
			{
				// Output domain points to a invalid memory location
				if ( !output )
				{
					msg << "CMultiDomainPreprocessor::CheckOutput ** ERROR ** invalid output for domain " << i << "\n";
					m_errors += msg.str();
					msg.str("");
					no_errors = false;
				}
					
				// Check for empty output
				else
				{
					if ( ( output->meshes.empty() ) || ( output->meshes.front()->faces.empty() ) )
					{
						msg << "CMultiDomainPreprocessor::CheckOutput ** WARNING ** Some domains were lost due to the priority scheme. Output domain " << i << " is empty\n";
						m_warnings += msg.str();
						msg.str("");
						no_errors = true;
					}
				}

				++i;
			}
			return no_errors;
		}
	}

	// ---------------------------------------------------------------------------- //
	bool CMultiDomainPreprocessor::AssignPrioritiesBaseOnTopology()
	{
		// Typedef(s)
		using AABB = carve::geom::aabb<3>;
		using AABBList = std::vector<AABB>;

		// Scope variables
		size_t				nDomains = m_domains.size();
		std::stringstream	msg;
		AABBList			aabbs( nDomains );

		// Functor(s)
		auto volumeAABB = []( const AABB &x )
		{
			return x.volume();
		};
		auto containedAABBs = [&volumeAABB]( const AABB &x, const AABB &y )
		{
			// Scope variables
			double			zmin, zmax;
			AABB::vector_t	center, extent;

			// Check for intersection between input AABBs
			if ( !x.intersects( y ) )
				return -1;

			// Compute the AABB of the intersection
			const auto	xmin = x.min(), xmax = x.max(),
						ymin = y.min(), ymax = y.max();
			for ( int i = 0; i < 3; ++i )
			{
				if ( xmin[i] > ymin[i] )
				{
					zmin = xmin[i];
					zmax = std::min( xmax[i], ymax[i] );
				}
				else
				{
					zmin = ymin[i];
					zmax = std::min( xmax[i], ymax[i] );
				}
				center[i] = 0.5 * ( zmin + zmax );
				extent[i] = center[i] - zmin;
			}
			AABB z( center, extent );

			// Detect whether x is partially included in y or the opposite.
			auto volx = x.volume();
			auto voly = y.volume();
			auto voli = z.volume();
			if ( voli / volx > 0.75 )
				return 1;
			if ( voli / voly > 0.75 )
				return 0;

			return -1;
		};

		// Compute and store domains AABBs
		for ( size_t i = 0; i < nDomains; ++i )
		{
			auto &domain = std::get<0>( m_domains[i] );
			aabbs[i] = domain->getAABB();
		}

		// Assign priority based on partial inclusion of domain's AABBs
		for ( size_t i = 0; i < nDomains; ++i )
		{
			auto &p0 = std::get<1>( m_domains[i] );

			for ( size_t j = i+1; j < nDomains; ++j )
			{
				auto &p1 = std::get<1>( m_domains[j] );
				auto flag = containedAABBs( aabbs[i], aabbs[j] );

				// aabbs[i] contains aabbs[j]
				if ( flag == 0 )
				{
					if ( p0 > p1 )
					{
						msg << "CMultiDomainPreprocessor::AssignPrioritiesBaseOnTopology ** WARNING ** Priority of domain " << j << " has been changed from " << p1 << " to " << p0 + 1 << "\n";
						p1 = p0 + 1;
						m_warnings += msg.str();
						msg.str("");
					}
				}

				// aabbs[j] contains aabbs[i]
				else if ( flag == 1 )
				{
					if ( p1 > p0 )
					{
						msg << "CMultiDomainPreprocessor::AssignPrioritiesBaseOnTopology ** WARNING ** Priority of domain " << i << " has been changed from " << p0 << " to " << p1 + 1 << "\n";
						p0 = p1 + 1;
						m_warnings += msg.str();
						msg.str("");
					}
				}

			} //next j
		} //next i

		return true;
	}
	
	// ============================================================================ //
	// CONVERSION OEPRATOR(S)                      									//
	// ============================================================================ //

	// vtk <-> carve ============================================================== //

	// ---------------------------------------------------------------------------- //
	bool	CopyVtkPolyDataToCarveMesh( 
		const vtkSmartPointer<vtkPolyData>	&mesh_in, 
		carve::mesh::MeshSet<3>				*&mesh_out
	) {
		// Trivial case: input vtkPolyData is empty
		if ( mesh_in->GetNumberOfCells() == 0 )
			return false;

		// Clean carve mesh
		if ( mesh_out )
		{
			delete mesh_out;
			mesh_out = nullptr;
		}

		// Typedef(s)
		using vec3 = carve::mesh::MeshSet<3>::vertex_t::vector_t;

		// Scope variables
		std::vector<vec3>	verts;
		std::vector<int>	face_ids;

		// Populate carve mesh with polydata vertexes
		{
			// Scope variables
			vec3				vert;
			int					nverts = mesh_in->GetNumberOfPoints();
			double				coords[3];

			// Create vertex list
			verts.reserve( nverts );
			for ( auto i = 0; i < nverts; ++i )
			{
				mesh_in->GetPoint( i, coords );
				vert[0] = coords[0];
				vert[1] = coords[1];
				vert[2] = coords[2];
				verts.push_back( vert );
			}
		}

		// Populate carve mesh with polydata faces
		{
			// Scope variables
			int					npolys = mesh_in->GetNumberOfCells();
			vtkIdType			nverts;
			vtkIdList			*ids = vtkIdList::New();

			// Create faces
			face_ids.reserve( 4 * npolys );
			for ( int i = 0; i < npolys; ++i )
			{
				mesh_in->GetCellPoints( i, ids );
				nverts = ids->GetNumberOfIds();
				face_ids.push_back( nverts );
				for ( auto j = 0; j < nverts; ++j )
					face_ids.push_back( ids->GetId( j ) );
			} //next i
			ids->Delete();

		}

		// Create a new carve mesh
		mesh_out = new carve::mesh::MeshSet<3>( verts, static_cast<size_t>( mesh_in->GetNumberOfCells() ), face_ids );
		return true;
	}
	
	// ---------------------------------------------------------------------------- //
	bool	CopyCarveMeshToVtkPolyData(
		const carve::mesh::MeshSet<3>		&mesh_in, 
		vtkSmartPointer<vtkPolyData>		&mesh_out,
		const carve::interpolate::FaceAttr< CMultiDomainPreprocessor::TagType > &tags_in
	) {
		// Scope variables
		size_t nV, nT;
		
		// Check inputs
		{
			// Nr of mesh verts
			nV = mesh_in.vertex_storage.size();

			// Nr of mesh polygons
			nT = 0;
			for ( const auto &m : mesh_in.meshes )
				nT += m->faces.size();
			if ( nV == 0 || nT == 0 )
				return false;
		}
		
		// Clean polys
		mesh_out = vtkSmartPointer< vtkPolyData >::New();
		vtkSmartPointer< vtkIntArray > tags_out = vtkSmartPointer< vtkIntArray >::New();
		
		// Configure tags
		tags_out->SetName( CMultiDomainPreprocessor::TAG_NAME.c_str() );
		tags_out->SetNumberOfComponents(2);
		
		// Create points
		{
			// Scope variables
			auto points = vtkSmartPointer< vtkPoints >::New();
			double	coords[3];

			// Populate points
			points->SetNumberOfPoints( nV );
			for ( auto i = 0; i < nV; ++i )
			{
				const auto &v = mesh_in.vertex_storage[i].v;
				coords[0] = v[0];
				coords[1] = v[1];
				coords[2] = v[2];
				points->SetPoint( i, coords );
			} //next i
			mesh_out->SetPoints( points );
		}

		// Create triangles
		{
			// Scope variables
			auto		f = mesh_in.faceBegin(), fe = mesh_in.faceEnd();
			auto		i = 0;
			auto		ids = vtkSmartPointer<vtkIdList>::New();
			const		carve::mesh::MeshSet<3>::vertex_t* VBase = mesh_in.vertex_storage.data();

			// Allocate memory for cells
			mesh_out->Allocate( nT );

			// Create polygonal cells
			ids->SetNumberOfIds(3);
			for ( ; f != fe; ++f, ++i )
			{
				// Scope variables
				const auto	face = *f;
				auto		edge = face->edge;
				CMultiDomainPreprocessor::TagType tag;

				// Retrieve tag associated to the current face
				tag = const_cast< CMultiDomainPreprocessor::TagSysType& >( tags_in ).getAttribute( *f );

				// Case 0: face is triangular
				if ( face->nVertices() == 3 )
				{
					auto *e = face->edge;
					for ( int j = 0; j < 3; j++, e = e->next )
						ids->SetId( j, static_cast<vtkIdType>( e->v1() - VBase ) );
					mesh_out->InsertNextCell( VTK_TRIANGLE, ids );
					tags_out->InsertNextTuple2( tag.first, tag.second );
				}

				// Case 1: face is polygonal
				else
				{
					// <<<<<<<<<<<< use carve triangulate >>>>>>>>>>> //
					{
						// // Scope variables
						// std::vector< carve::mesh::MeshSet<3>::vertex_t * > vloop;
						// std::vector< carve::triangulate::tri_idx > result;

						// // Triangulate the current face
						// face->getVertices(vloop);
						// carve::triangulate::triangulate(
							// carve::mesh::MeshSet<3>::face_t::projection_mapping(face->project),
							// vloop,
							// result
						// );
						// carve::triangulate::improve(
							// carve::mesh::MeshSet<3>::face_t::projection_mapping(face->project),
							// vloop,
							// carve::mesh::vertex_distance(),
							// result
						// );

						// // Add new triangles
						// for (auto i = 0; i < result.size(); ++i)
						// {
							// ids->SetId( 0, static_cast<bit32>(vloop[result[i].a] - VBase) );
							// ids->SetId( 1, static_cast<bit32>(vloop[result[i].b] - VBase) );
							// ids->SetId( 2, static_cast<bit32>(vloop[result[i].c] - VBase) );
							// mesh_out->InsertNextCell( VTK_TRIANGLE, ids );
							// tags_out->InsertNextTuple2( tag.first, tag.second );
						// } //next i

					}

					// <<<<<<<<<<<< use XCore CVertexLoopMesher >>>>>> //
					{
						// Typedefs
						using vert_t = std::array<double, 3>;

						// Scope variables
						std::vector< carve::mesh::MeshSet<3>::vertex_t * > vloop;
						face->getVertices(vloop);
						std::vector< vert_t > verts( vloop.size() );
						std::vector< const vert_t* > verts_ptr( vloop.size() );
						Meshing::CVertexLoopMesher< vert_t >::TriaList trias;

						// Conversion
						size_t v_id = 0;
						for ( auto &v : vloop )
						{
							verts[v_id][0] = v->v[0];
							verts[v_id][1] = v->v[1];
							verts[v_id][2] = v->v[2];
							verts_ptr[v_id] = verts.data() + v_id;
							++v_id;
						}

						// Mesh vertex loop
						Meshing::CVertexLoopMesher< vert_t > mesher( verts_ptr );
						if ( mesher.IsGood() )
							trias = mesher.Mesh();
						else
							return false;

						// Add triangles
						for ( const auto &tris : trias )
						{
							vtkIdType id;
							id = static_cast<vtkIdType>( vloop.at( tris[0] - verts.data() ) - VBase );
							ids->SetId( 0, id );
							id = static_cast<vtkIdType>( vloop.at( tris[1] - verts.data() ) - VBase );
							ids->SetId( 1, id );
							id = static_cast<vtkIdType>( vloop.at( tris[2] - verts.data() ) - VBase );
							ids->SetId( 2, id );
							mesh_out->InsertNextCell( VTK_TRIANGLE, ids );
							tags_out->InsertNextTuple2( tag.first, tag.second );
						}

					}
				}

			} //next f
			mesh_out->GetCellData()->AddArray( tags_out );
		}
		
		return true;
		
		// Create polygons
		{
			// Scope variables
			auto		f = mesh_in.faceBegin(), fe = mesh_in.faceEnd();
			auto		i = 0;
			vtkIdList	*ids = vtkIdList::New();
			size_t		nverts, ncells = 0;
			const		carve::mesh::MeshSet<3>::vertex_t* VBase = mesh_in.vertex_storage.data();

			// Initializations
			for ( const auto &m : mesh_in.meshes )
				ncells += m->faces.size();
			mesh_out->Allocate( ncells );
			tags_out->SetNumberOfValues( 2 * ncells );

			// Create polygonal cells
			for ( ; f != fe; ++f, ++i )
			{
				// Scope variables
				const auto	face = *f;
				auto		edge = face->edge;
				const auto 	&tag = const_cast< carve::interpolate::FaceAttr< CMultiDomainPreprocessor::TagType >& >( tags_in ).getAttribute( *f );
				nverts = face->nVertices();
				ids->SetNumberOfIds( nverts );
				for ( auto j = 0; j < nverts; ++j )
				{
					ids->SetId( j, static_cast<vtkIdType>( edge->v1() - VBase ) );
					edge = edge->next;
				}
				mesh_out->InsertNextCell( VTK_POLYGON, ids );
				tags_out->SetTuple2( i, tag.first, tag.second );
			} //next f
		}
		
		// Pack output
		mesh_out->GetCellData()->AddArray( tags_out );

		return true;
	}

	// QTech <-> carve ============================================================ //

	// ---------------------------------------------------------------------------- //
	bool CopyQGenTriangleMeshToCarveMesh( 
		const QGenTriangleMesh &mesh_in, 
		carve::mesh::MeshSet<3> *&mesh_out 
	) {
		// Check input
		if ( mesh_out )
		{
			delete mesh_out;
			mesh_out = nullptr;
		}
		
		// Typedef(s)
		using vec3 = carve::mesh::MeshSet<3>::vertex_t::vector_t;
		
		// Scope variables
		std::vector< vec3 > verts;
		std::vector< int >	faces;
		size_t				ncells = 0;

		// Collect vertexes
		{
			verts.reserve( mesh_in.verts_count() );
			vec3 v;
			for ( const auto &p : mesh_in.get_verts_buffer() )
			{
				v = p;
				verts.push_back(v);
			} //next point
		}

		// Collect faces
		{
			faces.reserve( mesh_in.tris_count() * 4 );
			for ( auto& t : mesh_in.get_tris_buffer() )
			{
				faces.push_back( 3 );
				++ncells;
				faces.push_back( t.inds[0] );
				faces.push_back( t.inds[1] );
				faces.push_back( t.inds[2] );
			} //next t
		}

		// Instantiate a new carve mesh object
		mesh_out = new carve::mesh::MeshSet<3>( verts, ncells, faces );

		return true;
	}

	// ---------------------------------------------------------------------------- //
	bool CopyCarveMeshToQGenTriangleMesh(
		const carve::mesh::MeshSet<3>			&mesh_in, 
		QGenTriangleMesh						&mesh_out
	) {

		// Scope variables
		QGenTriangleMesh::vec3_buffer	verts;
		QGenTriangleMesh::tri_buffer	tris;

		// Add vertexes
		{
			verts.reserve( mesh_in.vertex_storage.size() );
			for ( auto& p : mesh_in.vertex_storage )
				verts.emplace_back( p.v.x, p.v.y, p.v.z );
		}

		// Add faces
		{
			// Scope variables
			const carve::mesh::MeshSet<3>::vertex_t* VBase = mesh_in.vertex_storage.data();
			auto f = mesh_in.faceBegin(), fe = mesh_in.faceEnd();
			
			// Reserve memory for face/tag insertion
			size_t ncells = 0;
			for ( auto &m : mesh_in.meshes )
				ncells += m->faces.size();
			tris.reserve( ncells );
			
			// Loop over faces and extract triangular faces/triangulate polygonal faces
			for ( ; f != fe; ++f )
			{
				// Scope variables
				const	carve::mesh::MeshSet<3>::face_t *face = *f;


				// Case 0: face is triangular
				if ( face->nVertices() == 3 ) 
				{
					QGenTriangleMesh::tri t;
					auto *e = face->edge;
					for ( int i = 0; i < 3; i++, e = e->next)
						t.inds[i] = static_cast<bit32>( e->v1() - VBase );
					tris.push_back(t);
				}

				// Case 1: face is polygonal
				else
				{
					// <<<<<<<<<<<<<<< use carve::triangulate >>>>>>>> //
					{
						// // Scope variables
						// std::vector< carve::triangulate::tri_idx > result;
						// std::vector< carve::mesh::MeshSet<3>::vertex_t * > vloop;
						
						// // Get the vertex loop corresponding to the current face.
						// face->getVertices( vloop );

						// // Triangulate the current face
						// carve::triangulate::triangulate(
							// carve::mesh::MeshSet<3>::face_t::projection_mapping( face->project ),
							// vloop,
							// result
						// );

						// // Add new faces to tris buffer.
						// QGenTriangleMesh::tri t;
						// for ( auto i = 0; i < result.size(); ++i )
						// {
							// t.inds[0] = static_cast<bit32>( vloop[result[i].a] - VBase );
							// t.inds[1] = static_cast<bit32>( vloop[result[i].b] - VBase );
							// t.inds[2] = static_cast<bit32>( vloop[result[i].c] - VBase );
							// tris.push_back(t);
							// tags_out.push_back( {{ tag.first, tag.second }} );
						// } //next i
					}
					
					// <<<<<<<<<<<< use XCore CVertexLoopMesher >>>>>> //
					{
						// Typedefs
						using vert_t = std::array<double, 3>;

						// Scope variables
						std::vector< carve::mesh::MeshSet<3>::vertex_t * > vloop;
						face->getVertices(vloop);
						std::vector< vert_t > verts( vloop.size() );
						std::vector< const vert_t* > verts_ptr( vloop.size() );
						Meshing::CVertexLoopMesher< vert_t >::TriaList trias;

						// Conversion
						size_t v_id = 0;
						for ( auto &v : vloop )
						{
							verts[v_id][0] = v->v[0];
							verts[v_id][1] = v->v[1];
							verts[v_id][2] = v->v[2];
							verts_ptr[v_id] = verts.data() + v_id;
							++v_id;
						}

						// Mesh vertex loop
						Meshing::CVertexLoopMesher< vert_t > mesher( verts_ptr );
						if ( mesher.IsGood() )
							trias = mesher.Mesh();
						else
							return false;

						// Add triangles
						QGenTriangleMesh::tri t;
						for ( auto i = 0; i < trias.size(); ++i )
						{
							t.inds[0] = static_cast<bit32>( vloop.at( trias[i][0] - verts.data() ) - VBase );
							t.inds[1] = static_cast<bit32>( vloop.at( trias[i][1] - verts.data() ) - VBase );
							t.inds[2] = static_cast<bit32>( vloop.at( trias[i][2] - verts.data() ) - VBase );
							tris.push_back(t);
						} //next i

					}
				
				}
			} //next f
		}

		// Set vertexes and triangles in the QGenTriangleMesh object
		mesh_out.set_verts_buffer(verts);
		mesh_out.set_tris_buffer(tris);

		return true;
	}

	// ---------------------------------------------------------------------------- //
	bool CopyCarveMeshToQGenTriangleMesh( 
		const carve::mesh::MeshSet<3>			&mesh_in, 
		QGenTriangleMesh						&mesh_out,
		const carve::interpolate::FaceAttr< CMultiDomainPreprocessor::TagType > &tags_in,
		std::vector< std::array< int, 2 > >		&tags_out
	) {

		// Scope variables
		QGenTriangleMesh::vec3_buffer	verts;
		QGenTriangleMesh::tri_buffer	tris;

		// Add vertexes
		{
			verts.reserve( mesh_in.vertex_storage.size() );
			for ( auto& p : mesh_in.vertex_storage )
				verts.emplace_back( p.v.x, p.v.y, p.v.z );
		}

		// Add faces
		{
			// Scope variables
			const carve::mesh::MeshSet<3>::vertex_t* VBase = mesh_in.vertex_storage.data();
			auto f = mesh_in.faceBegin(), fe = mesh_in.faceEnd();
			
			// Reserve memory for face/tag insertion
			size_t ncells = 0;
			for ( auto &m : mesh_in.meshes )
				ncells += m->faces.size();
			tags_out.clear();
			tags_out.reserve( ncells );
			tris.reserve( ncells );
			
			// Loop over faces and extract triangular faces/triangulate polygonal faces
			for ( ; f != fe; ++f )
			{
				// Scope variables
				const	carve::mesh::MeshSet<3>::face_t *face = *f;
				CMultiDomainPreprocessor::TagType *faceTag = nullptr;
				const auto &tag = const_cast< carve::interpolate::FaceAttr< CMultiDomainPreprocessor::TagType >& >( tags_in ).getAttribute( *f );

				// Case 0: face is triangular
				if ( face->nVertices() == 3 ) 
				{
					QGenTriangleMesh::tri t;
					auto *e = face->edge;
					for ( int i = 0; i < 3; i++, e = e->next)
						t.inds[i] = static_cast<bit32>( e->v1() - VBase );
					tris.push_back(t);
					tags_out.push_back( {{ tag.first, tag.second }} );
				}

				// Case 1: face is polygonal
				else
				{
					// <<<<<<<<<<<<<<< use carve::triangulate >>>>>>>> //
					{
						// // Scope variables
						// std::vector< carve::triangulate::tri_idx > result;
						// std::vector< carve::mesh::MeshSet<3>::vertex_t * > vloop;
						
						// // Get the vertex loop corresponding to the current face.
						// face->getVertices( vloop );

						// // Triangulate the current face
						// carve::triangulate::triangulate(
							// carve::mesh::MeshSet<3>::face_t::projection_mapping( face->project ),
							// vloop,
							// result
						// );

						// // Add new faces to tris buffer.
						// QGenTriangleMesh::tri t;
						// for ( auto i = 0; i < result.size(); ++i )
						// {
							// t.inds[0] = static_cast<bit32>( vloop[result[i].a] - VBase );
							// t.inds[1] = static_cast<bit32>( vloop[result[i].b] - VBase );
							// t.inds[2] = static_cast<bit32>( vloop[result[i].c] - VBase );
							// tris.push_back(t);
							// tags_out.push_back( {{ tag.first, tag.second }} );
						// } //next i
					}
					
					// <<<<<<<<<<<< use XCore CVertexLoopMesher >>>>>> //
					{
						// Typedefs
						using vert_t = std::array<double, 3>;

						// Scope variables
						std::vector< carve::mesh::MeshSet<3>::vertex_t * > vloop;
						face->getVertices(vloop);
						std::vector< vert_t > verts( vloop.size() );
						std::vector< const vert_t* > verts_ptr( vloop.size() );
						Meshing::CVertexLoopMesher< vert_t >::TriaList trias;

						// Conversion
						size_t v_id = 0;
						for ( auto &v : vloop )
						{
							verts[v_id][0] = v->v[0];
							verts[v_id][1] = v->v[1];
							verts[v_id][2] = v->v[2];
							verts_ptr[v_id] = verts.data() + v_id;
							++v_id;
						}

						// Mesh vertex loop
						Meshing::CVertexLoopMesher< vert_t > mesher( verts_ptr );
						if ( mesher.IsGood() )
							trias = mesher.Mesh();
						else
							return false;

						// Add triangles
						QGenTriangleMesh::tri t;
						for ( auto i = 0; i < trias.size(); ++i )
						{
							t.inds[0] = static_cast<bit32>( vloop.at( trias[i][0] - verts.data() ) - VBase );
							t.inds[1] = static_cast<bit32>( vloop.at( trias[i][1] - verts.data() ) - VBase );
							t.inds[2] = static_cast<bit32>( vloop.at( trias[i][2] - verts.data() ) - VBase );
							tris.push_back(t);
							tags_out.push_back( {{ tag.first, tag.second }} );
						} //next i

					}
				
				}
			} //next f
		}

		// Set vertexes and triangles in the QGenTriangleMesh object
		mesh_out.set_verts_buffer(verts);
		mesh_out.set_tris_buffer(tris);

		return true;
	}

	// ============================================================================ //
	// UTILITIES                                   									//
	// ============================================================================ //
	
	// ---------------------------------------------------------------------------- //
	bool ImprintMeshes( const QGenTriangleMesh &mesh0, const QGenTriangleMesh &mesh1, QGenTriangleMesh &out0, QGenTriangleMesh &out1, bool symmetric)
	{
		// Check input(s)
		{
			auto nc0 = mesh0.tris_count(),
				 nc1 = mesh1.tris_count();
			if ( ( nc0 == 0 ) || ( nc1 == 0 ) )
				return false;
		}

		// Scope variables
		carve::mesh::MeshSet<3>		*m0 = nullptr,
									*m1 = nullptr,
									*o0 = nullptr,
									*o1 = nullptr;

		// Conversions
		if ( !CopyQGenTriangleMeshToCarveMesh( mesh0, m0 ) )
			return false;
		if ( !CopyQGenTriangleMeshToCarveMesh( mesh1, m1 ) )
			return false;

		// Symmetric imprint
		{
			// Typedef(s)
			using DomainDescriptorList = CMultiDomainPreprocessor::DomainDescriptorList;
			using Collector = CMultiDomainCollector;
			using TagSysType = CMultiDomainPreprocessor::TagSysType;
			using FlagSysType = CMultiDomainPreprocessor::FlagSysType;

			// Scope variables
			DomainDescriptorList	descriptors(2);
			TagSysType				tags0, tags1;
			FlagSysType				flags0, flags1;
			std::vector< std::array< int, 2 > > tags_out;

			// Initialize data
			std::get<0>( descriptors[0] ) = m0;
			std::get<1>( descriptors[0] ) = 0;
			std::get<0>( descriptors[1] ) = m1;
			std::get<1>( descriptors[1] ) = 1;

			// Scope variables
			Collector				collector( m0, m1 );
			carve::csg::CSG			csg;

			// Initialize collector
			collector.SetID0( 0 );
			collector.SetTags0( &tags0 );
			collector.SetFlags0( &flags0 );
			collector.SetID1( 1 );
			collector.SetTags1( &tags1 );
			collector.SetFlags1( &flags1 );
			collector.SetDomains( &descriptors );

			// Cut domain i with domain j and viceversa
			csg.compute( m0, m1, collector, nullptr, carve::csg::CSG::CLASSIFY_NORMAL);
			collector.Finalize();

			// Collect outputs
			if (symmetric)
			{
				o0 = collector.GetOutput0();
				o1 = collector.GetOutput1();
			}
			else
			{
				o0 = collector.GetCleanOutput0();
				o1 = collector.GetCleanOutput1();
			}
		}		
		
		// Converts back to QGenTriangleMesh
		if ( !CopyCarveMeshToQGenTriangleMesh( *o0, out0 ) )
			return false;
		if ( !CopyCarveMeshToQGenTriangleMesh( *o1, out1 ) )
			return false;

		return true;
	}

	
} }
