#pragma once

#include <vector>

using bit32 = unsigned int;

class QGenTriangleMesh
{
public:
	struct tri
	{
		bit32 inds[3];
	};
	struct vec3
	{
		union {
			float v[3];
			struct {float x, y, z; };
		};

		vec3(float v0, float v1, float v2) : x(v0), y(v1), z(v2) {}

		float& operator[](std::size_t idx) { return v[idx]; }
		const float& operator[](std::size_t idx) const { return v[idx]; }
	};
	using vec3_buffer = std::vector<vec3>;
	using tri_buffer = std::vector<tri>;

	bit32 verts_count() const { return static_cast<bit32>(_verts.size()); }
	const std::vector<vec3>& get_verts_buffer() const { return _verts; }
	void set_verts_buffer(const std::vector<vec3>& v) { _verts = v; }

	bit32 tris_count() const { return static_cast<bit32>(_tris.size()); }
	const std::vector<tri>& get_tris_buffer() const { return _tris; }
	void set_tris_buffer(const std::vector<tri>& t) { _tris = t; }

private:
	vec3_buffer _verts;
	tri_buffer _tris;
};
