// Copyright 2006-2015 Tobias Sargeant (tobias.sargeant@gmail.com).
//
// This file is part of the Carve CSG Library (http://carve-csg.com/)
//
// Permission is hereby granted, free of charge, to any person
// obtaining a copy of this software and associated documentation
// files (the "Software"), to deal in the Software without
// restriction, including without limitation the rights to use, copy,
// modify, merge, publish, distribute, sublicense, and/or sell copies
// of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be
// included in all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
// EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
// MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
// NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
// BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
// ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
// CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.


#include <carve/math.hpp>
#include <carve/matrix.hpp>

#include <iostream>

#define D(x) strtod(x, NULL)

int main(int argc, char** argv)
{
	carve::math::Matrix3 m;
	m._11 = D(argv[1]);
	m._12 = D(argv[2]);
	m._13 = D(argv[3]);
	m._21 = D(argv[2]);
	m._22 = D(argv[4]);
	m._23 = D(argv[5]);
	m._31 = D(argv[3]);
	m._32 = D(argv[5]);
	m._33 = D(argv[6]);

	double l1, l2, l3;
	carve::geom::vector<3> e1, e2, e3;

	carve::math::eigSolveSymmetric(m, l1, e1, l2, e2, l3, e3);
	std::cout << l1 << " " << e1 << std::endl;
	std::cout << l2 << " " << e2 << std::endl;
	std::cout << l3 << " " << e3 << std::endl;

	std::cout << m * e1 - l1 * e1 << "  " << (m * e1 - l1 * e1).isZero()
						<< std::endl;
	std::cout << m * e2 - l2 * e2 << "  " << (m * e2 - l2 * e2).isZero()
						<< std::endl;
	std::cout << m * e3 - l3 * e3 << "  " << (m * e3 - l3 * e3).isZero()
						<< std::endl;

	eigSolve(m, l1, l2, l3);

	std::cout << l1 << " " << e1 << std::endl;
	std::cout << l2 << " " << e2 << std::endl;
	std::cout << l3 << " " << e3 << std::endl;
}
