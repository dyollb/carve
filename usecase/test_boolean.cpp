#include "Carve.h"

#include "QGenTriangleMesh.h"

#include <vtkCellArray.h>
#include <vtkDataSetReader.h>
#include <vtkPoints.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

QGenTriangleMesh load_mesh(const std::string& fname)
{
	auto r = vtkSmartPointer<vtkDataSetReader>::New();
	r->SetFileName(fname.c_str());
	r->Update();
	auto pd = r->GetPolyDataOutput();

	// todo: copy
	QGenTriangleMesh::vec3_buffer verts;
	for (vtkIdType i = 0; i < pd->GetNumberOfPoints(); ++i)
	{
		double x[3];
		pd->GetPoint(i, x);
		verts.push_back(QGenTriangleMesh::vec3(x[0], x[1], x[2]));
	}

	QGenTriangleMesh::tri_buffer tris;
	vtkIdType npts, *pts;
	auto polys = pd->GetPolys();
	for (polys->InitTraversal(); polys->GetNextCell(npts, pts); )
	{
		QGenTriangleMesh::tri t;
		t.inds[0] = pts[0];
		t.inds[1] = pts[1];
		t.inds[2] = pts[2];
		tris.push_back(t);
	}

	QGenTriangleMesh mesh;
	mesh.set_verts_buffer(verts);
	mesh.set_tris_buffer(tris);
	return mesh;
};

int main(int argc, char** argv)
{
	// load meshes t1, t2
	std::string fname1 = "/Users/lloyd/Code/carve/usecase/_block.vtk";
	std::string fname2 = "/Users/lloyd/Code/carve/usecase/_helix.vtk";

	auto t1 = load_mesh(fname1);
	auto t2 = load_mesh(fname2);

	auto result = XCore::Boolean(t1, t2, XCore::eBoolean::kA_Minus_B);

	if (result)
	{
		std::cerr << result->verts_count() << std::endl;
	}

	return 0;
}
