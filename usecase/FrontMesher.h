// ================================================================================ //
//								** Front mesher **									//
//																					//
// Utilities for meshing non-convex vertex loops.									//
// ================================================================================ //
/*!
 *	@file	VkiMultiDomainPreprocessing.h
 *	@author	Alessandro Alaia (mail: alaia@itis.swiss.ch)
 *	@brief	Utilities for meshing non-convex vertex loops.
*/
# pragma once

// ================================================================================ //
// INCLUDES																			//
// ================================================================================ //

// Standard Template Library
# include <exception>
# include <limits>
# include <vector>
# include <array>


namespace  XCore { namespace Modeling
{
	namespace Meshing
	{
		// ======================================================================== //
		// STATIC CONSTANT(S)														//
		// ======================================================================== //
		static const double pi = 3.14159265358979;

		// ======================================================================== //
		// DEFINITION OF CLASS CVertexLoop											//
		// ======================================================================== //
		/*!
		 *	@brief Data structure holding a vertex loop.
		*/
		template< class T >
		class CVertexLoop
		{
			// Typedef(s) --------------------------------------------------------- //
		public:
			using VertexObj = T;													///< vertex object
			using VertexObjPtr = const VertexObj*;									///< vertex pointer
			using VertexList = std::vector< VertexObj >;							///< vertex list
			using VertexPtrList = std::vector< VertexObjPtr >;						///< vertex pointer list
		private:
			using SelfType = CVertexLoop<T>;										///< self type

			// Constructor(s) ----------------------------------------------------- //
		public:
			/*!
			 *	@brief Initialize a empty vertex loop
			*/
			CVertexLoop()
			= default;
			/*!
			 *	@brief Initialize a loop from the list of vertex pointers.
			*/
			CVertexLoop( const VertexPtrList &verts ) :
				m_loop( verts )
			{
				ComputeLoopNormal();
				ComputeLoopAngles();
			}

			// Static method(s) --------------------------------------------------- //
		public:
			/*!
			 *	@brief Compute the oriented angle between two vectors on the plane with specified normal.
			 *
			 *	@note The order of argument is relevant, since the oriented angle between y and x is the complementary of
			 *	the oriented angle between x and y.
			 *
			 *	@param [in]		x 		1st vector
			 *	@param [in]		y 		2nd vector
			 *	@param [in]		n 		plane normal
			*/
			static double		ComputeOrientedAngle( VertexObj x, VertexObj y, const VertexObj &n )
			{
				// Normalize input vectors
				SelfType::Normalize( x );
				SelfType::Normalize( y );
				
				// If the two input vectors are parallel returns pi
				double dp1 = SelfType::ScalarProduct( x, y );
				if ( std::abs( dp1 + 1. ) < 1.e-12 )
					return pi;
				
				// Compute the oriented angle
				// Compute the oriented angle
				auto t = SelfType::CrossProduct( x, y );
				SelfType::Normalize( t );
				double theta = std::acos( std::min( 1., std::max( -1., dp1 ) ) );
				double dp2 = SelfType::ScalarProduct( t, n );
				
				// Correct angle for non-convex cases
				if ( dp2 < 0. )
					theta = 2. * pi - theta;
				return theta;
			}
			/*!
			 *	@brief Compute the Euclidean norm of the input vector.
			*/
			static double		Norm( const VertexObj &p )
			{
				return sqrt( p[0] * p[0]  + p[1] * p[1] + p[2] * p[2] );
			}
			/*!
			 *	@brief Normalize the input vector.
			*/
			static void			Normalize( VertexObj &p )
			{
				auto n = SelfType::Norm(p);
				p[0] /= n;
				p[1] /= n;
				p[2] /= n;
			}
			/*!
			 *	@brief Divide the input vector by the input scalar.
			*/
			static VertexObj	Divide( const VertexObj &p, double s )
			{
				VertexObj o(p);
				o[0] /= s;
				o[1] /= s;
				o[2] /= s;
				return o;
			}
			/*!
			 *	@brief Multiply the input vector by the input scalar.
			*/
			static VertexObj	Mult( const VertexObj &p, double s )
			{
				VertexObj o(p);
				o[0] *= s;
				o[1] *= s;
				o[2] *= s;
				return o;
			}
			/*!
			 *	@brief Compute the vector sum of the input vectors.
			*/
			static VertexObj	Add( const VertexObj &p, const VertexObj &q )
			{
				VertexObj o;
				o[0] = p[0] + q[0];
				o[1] = p[1] + q[1];
				o[2] = p[2] + q[2];
				return o;
			}
			/*!
			 *	@brief Compute the vector difference of the input vectors.
			*/
			static VertexObj	Diff( const VertexObj &p, const VertexObj &q )
			{
				VertexObj o;
				o[0] = p[0] - q[0];
				o[1] = p[1] - q[1];
				o[2] = p[2] - q[2];
				return o;
			}
			/*!
			 *	@brief Compute the scalar product of the input vectors.
			*/
			static double		ScalarProduct( const VertexObj &p, const VertexObj &q )
			{
				return p[0] * q[0] + p[1] * q[1] + p[2] * q[2];
			}
			/*!
			 *	@brief Compute the cross-product of the input vectors.
			*/
			static VertexObj	CrossProduct( const VertexObj &p, const VertexObj &q )
			{
				VertexObj o;
				o[0] = p[1] * q[2] - p[2] * q[1];
				o[1] = p[2] * q[0] - p[0] * q[2];
				o[2] = p[0] * q[1] - p[1] * q[0];
				return o;
			}

			// Getter(s) ---------------------------------------------------------- //
		public:
			/*!
			 *	@brief Returns the nr. of vertexes in the loop.
			*/
			size_t				GetNumberOfVertexes() const
			{
				return m_loop.size();
			}
			/*!
			 *	@brief Returns the normal unit vector to the plane supporting the vertex loop.
			*/
			const VertexObj&	GetLoopNormal() const
			{
				return m_normal;
			}
			/*!
			 *	@brief Returns the angle at the i-th vertex of the loop.
			*/
			double				GetVertexAngle( size_t i ) const
			{
				return m_angles.at(i);
			}
			/*!
			 *	@brief Returns the coordinates of the i-th vertex in the loop.
			*/
			const VertexObj&	GetVertex( size_t i ) const
			{
				return *(m_loop.at(i));
			}
			/*!
			 *	@brief Returns const pointer to the i-th vertex in the loop.
			*/
			VertexObjPtr		GetVertexPtr( size_t i ) const
			{
				return m_loop.at(i);
			}
			/*!
			 *	@brief Returns the subloop obtained by moving from i to j (following the loop direction)
			*/
			SelfType			ExtractSubLoop( size_t i, size_t j )
			{
				VertexPtrList verts;
				for ( ; i != j; i = Right(i) )
					verts.push_back( m_loop[i] );
				verts.push_back( m_loop[j] );
				return SelfType( verts );
			}

			// Modifier(s) -------------------------------------------------------- //
		public:
			/*!
			 *	@brief Set the vertex loop
			*/
			void				SetLoop( const VertexPtrList &verts )
			{
				VertexPtrList( verts ).swap( m_loop );
				ComputeLoopNormal();
			}
			/*!
			 *	@brief Remove a vertex from the loop
			*/
			void				EraseVertex( size_t i )
			{
				// Remove vertex from the loop
				m_loop.erase( m_loop.begin() + i );
				m_angles.erase( m_angles.begin() + i );

				// Update angles
				if ( m_loop.size() >= 3 )
				{
					i = i % m_loop.size();
					auto l = Left(i);
					m_angles[i] = ComputeAngle(i);
					m_angles[l] = ComputeAngle(l);
				}
			}
			/*!
			 *	@brief Insert vertex in the loop before the specified position
			*/
			void				InsertVertex( const VertexObjPtr *vert, size_t i )
			{
				m_loop.insert( m_loop.begin() + i, vert );
				m_angles.insert( m_angles.begin() + i, 0. );
				if ( m_loop.size() >= 3 )
				{
					auto l = Left(i);
					auto r = Right(i);
					m_angles[i] = ComputeAngle(i);
					m_angles[l] = ComputeAngle(l);
					m_angles[r] = ComputeAngle(r);
				}
			}
			/*!
			 *	@brief Clear the content of this vertex loop
			*/
			void				Clear()
			{
				m_normal[0] = m_normal[1] = m_normal[2] = 0.;
				m_loop.clear();
				m_angles.clear();
			}

			// Info --------------------------------------------------------------- //
		public:
			/*!
			 *	@brief Returns the left neighbor of the i-th vertex in the loop.
			*/
			size_t				Left( size_t i ) const
			{
				return ( m_loop.size() + i - 1 ) % m_loop.size();
			}
			/*!
			 *	@brief Returns the right neighbor of the i-th vertex in the loop.
			*/
			size_t				Right( size_t i ) const
			{
				return ( i + 1 ) % m_loop.size();
			}
			/*!
			 *	@brief Returns (true) if the vertex loop is in a valid state
			*/
			bool				IsGood() const
			{
				if ( std::abs( CVertexLoop::Norm( m_normal ) - 1. ) > 1.e-12 )
					return false;
				return ( m_loop.size() >= 3 );
			}
			/*!
			 *	@brief Returns (true) if the vertex loop is convex (i.e. all the vertex angles
			 *	are below pi - tol).
			 *
			 *	@param [in]		tol		(default = .0001) tolerance for the convexity check
			*/
			bool				IsConvex( double tol = 0.0001 ) const
			{
				if ( Empty() )
					return false;
				for ( const auto &angle : m_angles )
				{
					if ( angle > pi - tol )
						return false;
				}
				return true;
			}
			/*!
			 *	@brief Returns (true) if the loop is empty
			*/
			bool				Empty() const
			{
				return m_loop.empty();
			}

			// Member(s) ---------------------------------------------------------- //
		protected:
			/*!
			 *	@brief Compute the normal unit vector of the plane supporting the vertex loop
			 *	using the Newell's method.
			 *
			 *	The normal vector is computed using least square.
			*/
			void				ComputeLoopNormal( )
			{
				// <<<<<<<<<<<<<< new implementation >>>>>>>>>>>>>>>>>>>> //
				// Note: this algorithm uses the Newell's which is robust
				// w.r.t. off-plane deviations, parallel consecutive edges
				// Requirements: all the vertexes should be distinct.
				// Reset normal
				m_normal[0] = m_normal[1] = m_normal[2] = 0.;
				for ( size_t i = 0; i < m_loop.size(); ++ i )
				{
					const auto &v = *( m_loop.at(i) );
					const auto &r = *( m_loop.at( Right(i) ) );
					m_normal[0] += ( v[1] - r[1] ) * ( v[2] + r[2] );
					m_normal[1] += ( v[2] - r[2] ) * ( v[0] + r[0] );
					m_normal[2] += ( v[0] - r[0] ) * ( v[1] + r[1] );
				}
				SelfType::Normalize( m_normal );
				return;

				// <<<<<<<<<<<<<< old implementation >>>>>>>>>>>>>>>>>>>> //
				// Note: This algorithms estimate the plane's normal using
				// the least square method. Consequently the plane normal 
				// might be inverted.
				// Requirements: the vertexes should not lie on a straight line.
				// Check input
				if ( m_loop.size() < 3 )
					return;

				// Scope variables
				size_t		n = m_loop.size();
				VertexObj	centroid;
				double		xx = 0., yy = 0., zz = 0., xy = 0., xz = 0., yz = 0.;

				// Compute centroid coordinates
				centroid[0] = centroid[1] = centroid[2] = 0.;
				for ( const auto &v : m_loop )
				{
					centroid[0] += (*v)[0];
					centroid[1] += (*v)[1];
					centroid[2] += (*v)[2];
				}
				centroid[0] /= (double) n;
				centroid[1] /= (double) n;
				centroid[2] /= (double) n;

				// Compute covariance matrix
				for ( const auto &v : m_loop )
				{
					xx += ( (*v)[0] - centroid[0] ) * ( (*v)[0] - centroid[0] );
					yy += ( (*v)[1] - centroid[1] ) * ( (*v)[1] - centroid[1] );
					zz += ( (*v)[2] - centroid[2] ) * ( (*v)[2] - centroid[2] );
					xy += ( (*v)[0] - centroid[0] ) * ( (*v)[1] - centroid[1] );
					xz += ( (*v)[0] - centroid[0] ) * ( (*v)[2] - centroid[2] );
					yz += ( (*v)[1] - centroid[1] ) * ( (*v)[2] - centroid[2] );
				}
				double det_x = yy*zz - yz*yz;
				double det_y = xx*zz - xz*xz;
				double det_z = xx*yy - xy*xy;
				double det_max = std::max( std::max( det_x, det_y ), det_z );
				if ( det_max <= 0.0 )
				{
					m_normal[0] = m_normal[1] = m_normal[2] = 0.;
					return;
				}

				// Pick path with best conditioning:
				if ( det_max == det_x )
				{
					m_normal[0] = det_x;
					m_normal[1] = xz*yz - xy*zz;
					m_normal[2] = xy*yz - xz*yy;
				}
				else if ( det_max == det_y )
				{
					m_normal[0] = xz*yz - xy*zz;
					m_normal[1] = det_y;
					m_normal[2] = xy*xz - yz*xx;
				}
				else
				{
					m_normal[0] = xy*yz - xz*yy;
					m_normal[1] = xy*xz - yz*xx;
					m_normal[2] = det_z;
				}
				SelfType::Normalize( m_normal );
			}
			/*!
			 *	@brief Compute the angle at the i-th vertex in the loop.
			*/
			double				ComputeAngle( size_t i ) const
			{
				// Scope variables
				const auto &v = *( m_loop.at(i) );
				const auto &l = *( m_loop.at( Left(i) ) );
				const auto &r = *( m_loop.at( Right(i) ) );

				// Compute the angle 
				auto ul = SelfType::Diff( l, v );
				auto ur = SelfType::Diff( r, v );
				return SelfType::ComputeOrientedAngle( ur, ul, m_normal );			
			}
			/*!
			 *	@brief Compute angles for all vertexes in the front
			*/
			void				ComputeLoopAngles()
			{
				m_angles.resize( m_loop.size() );
				for ( size_t i = 0; i < m_loop.size(); ++i )
					m_angles[i] = ComputeAngle(i);
			}

			// Member(s) ---------------------------------------------------------- //
		protected:
			VertexPtrList		m_loop;												///< vertex loop
			std::vector<double>	m_angles;											///< vertex angles
			VertexObj			m_normal;											///< normal unit vector to the plane supporting the loop

		};
		
		// ======================================================================== //
		// DEFINITION OF CLASS CConvexLoopMesher										//
		// ======================================================================== //
		/*!
		 *	@brief Generate a triangulation for a convex planar vertex loop.
		 *
		 *	@tparam			T		vertex object type
		*/
		template< class T >
		class CConvexLoopMesher
		{
			// Typedef(s) --------------------------------------------------------- //
		public:
			using VertexObj = T;													///< vertex object
			using VertexObjPtr = const VertexObj*;									///< vertex pointer
			using TriaObj = std::array< VertexObjPtr, 3 >;							///< tria object
			using TriaList = std::vector< TriaObj >;								///< tria list

			// Constructor(s) ----------------------------------------------------- //
		public:
			/*!
			 *	@brief Default constructor.
			 *
			 *	Initialize a empty mesher for a vertex loop.
			*/
			CConvexLoopMesher( )
			= default;
			/*!
			 *	@brief Initialize a mesher for the input vertex loop.
			*/
			CConvexLoopMesher( const CVertexLoop<VertexObj> &loop )
			{
				if ( loop.IsConvex() )
					m_loop = loop;
			}

			// Method(s) ---------------------------------------------------------- //
		protected:
			/*!
			 *	@brief Returns the next vertex in the loop to be processed
			*/
			size_t				Next() const
			{
				// Scope variables
				double min_angle = std::numeric_limits< double >::max();
				size_t j = -1;

				for ( auto i = 0; i < m_loop.GetNumberOfVertexes(); ++i )
				{
					double angle = m_loop.GetVertexAngle(i);
					if ( angle < min_angle )
					{
						min_angle = angle;
						j = i;
					}
				} //next i
				return j;
			}
			public:
			/*!
			 *	@brief Returns (true) if the settings are valid.
			*/
			bool				IsGood() const
			{
				if ( m_loop.Empty() )
					return false;
				return m_loop.IsGood();
			}
			/*!
			 *	@brief Set the loop
			*/
			void				SetLoop( const CVertexLoop<VertexObj> &loop )
			{
				if ( loop.IsConvex() )
					m_loop = loop;
			}
			/*!
			 *	@brief Returns a triangulation of the loop
			*/
			TriaList			Mesh()
			{
				if ( !IsGood() )
					return TriaList();

				// Scope variables
				TriaList	trias;

				// Generate the triangulation
				while ( m_loop.GetNumberOfVertexes() > 3 )
				{
					auto i = Next();
					auto l = m_loop.Left(i);
					auto r = m_loop.Right(i);
					trias.push_back( {{ m_loop.GetVertexPtr(l), m_loop.GetVertexPtr(i), m_loop.GetVertexPtr(r) }} );
					m_loop.EraseVertex(i);
				}

				// Closing triangles
				trias.push_back( {{ m_loop.GetVertexPtr(0), m_loop.GetVertexPtr(1), m_loop.GetVertexPtr(2) }} );

				return trias;
			}
			
			// Member(s) ---------------------------------------------------------- //
		protected:
			CVertexLoop<VertexObj>	m_loop;											///< vertex loop
		};

		// ======================================================================== //
		// DEFINITION OF CLASS CVertexLoopMesher										//
		// ======================================================================== //
		/*!
		 *	@brief Mesher for non-convex vertex loops
		*/
		template< class T >
		class CVertexLoopMesher
		{
			// Typedef(s) --------------------------------------------------------- //
		public:
			using VertexObj = T;													///< vertex object
			using VertexObjPtr = VertexObj*;										///< vertex pointer
			using VertexPtrList = std::vector< VertexObjPtr >;						///< vertex pointer list
			using TriaList = typename CConvexLoopMesher<T>::TriaList;				///< triangle list
			using LoopList = std::vector< CVertexLoop<T> >;							///< list of vertex loops
		private:
			using SelfType = CVertexLoopMesher<T>;									///< self type

			// Constructor(s) ----------------------------------------------------- //
		public:
			/*!
			 *	@brief Initialize a empty mesher.
			*/
			CVertexLoopMesher()
			= default;
			/*!
			 *	@brief Initialize a mesher for the input loop
			*/
			CVertexLoopMesher( const CVertexLoop<VertexObj> &loop )
			{
				m_loops.push_back( loop );
				if ( m_loops.front().IsGood() )
					BreakIntoConvexLoops();
			}

			// Method(s) ---------------------------------------------------------- //
		public:
			/*!
			 *	@brief Returns (true) if the mesher is empty.
			*/
			bool		Empty() const
			{
				return m_loops.empty();
			}
			/*!
			 *	@brief Returns (true) if the mesher is in good state.
			*/
			bool		IsGood() const
			{
				if ( Empty() )
					return false;
				for ( const auto &loop : m_loops )
				{
					if ( !loop.IsGood() )
						return false;
				}
				return true;
			}
			/*!
			 *	@brief Set the vertex loop
			*/
			void		SetLoop( const CVertexLoop<VertexObj> &verts )
			{
				LoopList( 1, verts ).swap( m_loops );
				if ( m_loops.front().IsGood() )
					BreakIntoConvexLoops( );
			}
			/*!
			 *	@brief Generate the mesh for the input loop
			*/
			TriaList	Mesh()
			{
				// Check if empty
				if ( Empty() )
					return TriaList();

				// Check if mesher is in good state
				if ( !IsGood() )
					return TriaList();

				// Generate mesh
				TriaList	 trias;
				for ( const auto &loop : m_loops )
				{
					CConvexLoopMesher<VertexObj> cmesher( loop );
					if ( cmesher.IsGood() )
					{
						TriaList tmp = cmesher.Mesh();
						trias.insert( trias.end(), tmp.cbegin(), tmp.cend() );
					}
				}

				return trias;
			}
			/*!
			 *	@brief Returns the nr. of convex loops
			*/
			size_t		GetNumberOfConvexLoops() const
			{
				return m_loops.size();
			}
		protected:
			/*!
			 *	@brief Break the input loop into convex loops.
			*/
			void		BreakIntoConvexLoops( )
			{
				// Constant(s)
				const double tol = 1.e-12;

				// Check for empty mesher
				if ( Empty() )
					return;

				// Check if loop is in good state
				if ( !m_loops.front().IsGood() )
					return;

				// If input loop is convex there is nothing to do
				size_t i = 0;
				while ( i < m_loops.size() )
				{
					// Scope variables
					auto		&loop = m_loops[i];

					// Loop is not convex
					if ( ( loop.GetNumberOfVertexes() > 3 ) && !loop.IsConvex() )
					{
						// Scope variables
						size_t		j;
						size_t		m = -1;
						size_t		k = -1;
						const auto	&loop_normal = loop.GetLoopNormal();

						// Locate the max concave angle
						{
							double angle;
							double max_angle = -std::numeric_limits< double >::max();
							for ( j = 0; j < loop.GetNumberOfVertexes(); ++j )
							{
								angle = loop.GetVertexAngle(j);
								if ( angle > max_angle )
								{
									max_angle = angle;
									k = j;
								}
							} //next j
						}

						// Locate the opposite vertex
						{
							// Scope variables
							double max_score = -std::numeric_limits<double>::max();
							double min_dist = std::numeric_limits<double>::max();

							// Compute the ideal bisection
							size_t		l = loop.Left(k),
										r = loop.Right(k);
							const auto	&vl = loop.GetVertex(l),
										&vr = loop.GetVertex(r),
										&vk = loop.GetVertex(k);
							VertexObj	rk = CVertexLoop<T>::Diff( vr, vk ),
										kl = CVertexLoop<T>::Diff( vk, vl ),
										bl, br, b;
							CVertexLoop<T>::Normalize(rk);
							CVertexLoop<T>::Normalize(kl);
							bl = CVertexLoop<T>::CrossProduct( loop_normal, kl );
							br = CVertexLoop<T>::CrossProduct( loop_normal, rk );
							b = CVertexLoop<T>::Add( bl, br );
							CVertexLoop<T>::Normalize(b);

							// Loop over front vertexes
							j = loop.Right( loop.Right(k) );
							l = loop.Left(k);
							while ( j != l )
							{
								// Scope variables
								const auto	&vj = loop.GetVertex(j);
								VertexObj	jk = CVertexLoop<T>::Diff( vj, vk );
								double		dist = CVertexLoop<T>::Norm(jk);

								// Current splitting direction
								jk = CVertexLoop<T>::Divide( jk, dist );
								double score = CVertexLoop<T>::ScalarProduct( b, jk );

								// Update max score
								if ( score > max_score + tol )
								{
									min_dist = dist;
									max_score = score;
									m = j;
								}
								else if ( std::abs( score - max_score ) <= tol )
								{
									if ( dist < min_dist )
									{
										min_dist = dist;
										max_score = score;
										m = j;
									}
								}

								// Move to next vertex
								j = loop.Right(j);
							}
						}

						// Split loop at location m
						{
							auto subloop1 = loop.ExtractSubLoop( k, m );
							auto subloop2 = loop.ExtractSubLoop( m, k );
							if ( subloop1.Empty() )
								std::cout << "subloop1 is empty" << std::endl;
							if ( subloop2.Empty() )
								std::cout << "subloop2 is empty" << std::endl;
							m_loops[i] = subloop1;
							m_loops.push_back( subloop2 );
							--i;
						}
					}

					// Update counter
					++i;
				} //next i

			}

			// Member(s) ---------------------------------------------------------- //
		protected:
			std::vector< CVertexLoop<VertexObj> >	m_loops;

		};
		
	}
} }