include_directories("${carve_SOURCE_DIR}/include")

file(GLOB carve_HEADERS *.hpp)

add_library(carve ${carve_HEADERS}
            aabb.cpp
            carve.cpp
            convex_hull.cpp
            csg.cpp
            csg_collector.cpp
            edge.cpp
            face.cpp
            geom.cpp
            geom2d.cpp
            geom3d.cpp
            intersect.cpp
            intersect_classify_edge.cpp
            intersect_classify_group.cpp
            intersect_debug.cpp
            intersect_face_division.cpp
            intersect_group.cpp
            intersect_half_classify_group.cpp
            intersection.cpp
            math.cpp
            mesh.cpp
            octree.cpp
            pointset.cpp
            polyhedron.cpp
            polyline.cpp
            tag.cpp
            timing.cpp
            triangulator.cpp
            triangle_intersection.cpp
			shewchuk_predicates.cpp
)

# set compile properties for predicates, to avoid bad behavior
if(MSVC)
	set_source_files_properties(shewchuk_predicates.cpp PROPERTIES COMPILE_FLAGS "/Od /fp:strict")
else(MSVC)
	set_source_files_properties(shewchuk_predicates.cpp PROPERTIES COMPILE_FLAGS -O0)
endif(MSVC)

set_target_properties(carve PROPERTIES
                      VERSION   "${carve_VERSION_MAJOR}.${carve_VERSION_MINOR}.${carve_VERSION_PATCH}"
                      SOVERSION "${carve_VERSION_MAJOR}.${carve_VERSION_MINOR}")

install(TARGETS carve
        LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/lib"
        ARCHIVE DESTINATION "${CMAKE_INSTALL_PREFIX}/lib")
