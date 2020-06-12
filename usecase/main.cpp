#include <carve/mesh.hpp>

#include "VkiMultiDomainPreprocessing.h"

#include <vtkPLYReader.h>
#include <vtkPLYWriter.h>

int main(int argc, char** argv)
{
	using namespace XCore::Modeling;

	{
	std::string dir = "/Users/lloyd/Models/BooleanEval/";
	//dir = "/home/lloyd/Data/";
	//dir = "Z:/Models/BooleanEval/";
	std::string f1 = dir + "A.ply";
	std::string f2 = dir + "B.ply";
	std::string of = dir + "A_B_prepro.ply";

	auto r1 = vtkSmartPointer<vtkPLYReader>::New();
	r1->SetFileName(f1.c_str());
	r1->Update();
	vtkSmartPointer<vtkPolyData> mesh1 = r1->GetOutput();

	auto r2 = vtkSmartPointer<vtkPLYReader>::New();
	r2->SetFileName(f2.c_str());
	r2->Update();
	vtkSmartPointer<vtkPolyData> mesh2 = r2->GetOutput();

	CMultiDomainPreprocessor prepro;
	prepro.EnforceUserPriority(true);
	auto id1 = prepro.AddDomain(mesh1, 1);
	auto id2 = prepro.AddDomain(mesh2, 2);
	bool ok = prepro.Update();
	assert(ok == true);

	assert(prepro.Error() == false);
	assert(prepro.Warning() == false);

	vtkSmartPointer<vtkPolyData> out;
	prepro.GetMergedDomains(out);

	auto w = vtkSmartPointer<vtkPLYWriter>::New();
	w->SetInputData(out);
	w->SetFileName(of.c_str());
	w->SetFileTypeToBinary();
	w->Write();
	
	}
	
	return 0;
}
