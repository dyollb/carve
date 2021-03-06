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

#include <gtest/gtest.h>


#include <carve/carve.hpp>
#include <carve/exact.hpp>

using namespace carve::exact;

inline std::ostream& operator<<(std::ostream& out,
		const std::vector<double>& p)
{
	out << '{';
	out << p[0];
	for (size_t i = 1; i < p.size(); ++i)
	{
		out << ';' << p[i];
	}
	out << '}';
	return out;
}

std::vector<double> dblvec(double a)
{
	std::vector<double> r;
	r.push_back(a);
	return r;
}

std::vector<double> dblvec(double a, double b)
{
	std::vector<double> r;
	r.push_back(a);
	r.push_back(b);
	return r;
}

void code(double* a, double* b, double* r)
{
	detail::op<1, 1>::add_fast(a, b, r);
}

TEST(ExactTest, ExactTest)
{
	double a, b, c;

	a = 4;
	b = 3;
	c = 1e60;

	EXPECT_EQ((detail::add<1, 1>(&a, &b)), exact_t(0.0, 7.0));
	EXPECT_EQ((detail::add<1, 1>(&a, &c)), exact_t(4.0, 1.0e60));

	EXPECT_EQ((detail::sub<1, 1>(&a, &b)), exact_t(0.0, 1.0));
	EXPECT_EQ((detail::sub<1, 1>(&c, &a)), exact_t(-4.0, 1e60));

	EXPECT_EQ((detail::sub<2, 1>(detail::sub<1, 1>(&c, &a), &a)), exact_t(0.0, -8.0, 1e60));

	double s1 = 4;
	double s2 = 3 * exp2(60);
	double s3 = 2 * exp2(120);
	double s4 = 1 * exp2(180);

	exact_t s;

	s = detail::add<2, 2>(detail::add<1, 1>(&s1, &s2), detail::add<1, 1>(&s3, &s4));
	EXPECT_EQ(s, exact_t(s1, s2, s3, s4));

	s = detail::add<2, 2>(detail::add<1, 1>(&s4, &s3), detail::add<1, 1>(&s2, &s1));
	EXPECT_EQ(s, exact_t(s1, s2, s3, s4));

	double add1 = 3;
	EXPECT_EQ((detail::add<4, 1>(s, &add1)).compressed(), exact_t(0, 7, s2, s3, s4).compressed());

	// (c-a) - (b-c) == (c+c) - (a+b)
	EXPECT_EQ(
			(detail::sub<2, 2>(detail::sub<1, 1>(&c, &a), detail::sub<1, 1>(&b, &c))).compressed(),
			(detail::sub<2, 2>(detail::add<1, 1>(&c, &c), detail::add<1, 1>(&a, &b))).compressed());

	//   EXPECT_EQ(detail::prod(a, b), detail::mkdbl(12.0, 0.0));
	//   EXPECT_EQ(detail::square(a),  detail::mkdbl(16.0, 0.0));

	//   EXPECT_EQ(detail::prod(a, c), detail::mkdbl(4e60, 0.0));

	//   EXPECT_EQ(detail::sum(a, c),  detail::sum(c, a));
	//   EXPECT_EQ(detail::prod(a, c), detail::prod(c, a));

	//   a = 805306457.0;
	//   b = 1610612741.0;
	//   c = 3221225473.0;
	//   d = 4294967291.0;

	//   EXPECT_EQ(detail::prod(detail::sum(a,b),  detail::sum(c,d)),
	//   detail::mkdbl(18158514394416283648.0,
	//                                                                               -376.0,
	//                                                                               0,0,0,0,0,0));
	//   EXPECT_EQ(detail::prod(detail::prod(a,b), detail::sum(c,d)),
	//   detail::mkdbl(9748778911853561743662907392.0,
	//                                                                               -543313364724.0,
	//                                                                               0,0,0,0,0,0));
	//   EXPECT_EQ(detail::prod(detail::prod(a,b), detail::prod(c,d)),
	//   detail::mkdbl(17944579966664105723432610040410800128.0,
	//                                                                               -312874079305504653312.0,
	//                                                                               -2048.0,
	//                                                                               -177.0,
	//                                                                               0, 0, 0, 0));

	//   EXPECT_EQ(detail::square(a), detail::mkdbl(648518489685892864.0, -15.0));

	//   EXPECT_EQ(exact_t(detail::square(detail::square(a))).compressed(),
	//             exact_t(detail::prod(detail::prod(a,a),
	//             detail::prod(a,a))).compressed());

	//   EXPECT_EQ(exact_t(detail::square(detail::square(a))).compressed(),
	//             exact_t(detail::prod(a, detail::prod(a, detail::prod(a,
	//             a)))).compressed());

	//   EXPECT_EQ(exact_t(detail::prod(detail::prod(a,b),
	//   detail::prod(c,d))).compressed(),
	//             exact_t(detail::mkdbl(17944579966664105723432610040410800128.0,
	//                                   -312874079305504653312.0,
	//                                   -2225.0)));

	//   EXPECT_EQ(exact_t(detail::prod(detail::prod(a,b),
	//   detail::prod(c,d))).compressed(),
	//             exact_t(detail::prod(detail::prod(d,b),
	//             detail::prod(c,a))).compressed());
}

TEST(ExactTest, SumZeroelim)
{
	//   exact_t result;
	//   sum_zeroelim(detail::mkdbl(18158514394416283648.0, -376.0),
	//                detail::mkdbl(9748778911853561743662907392.0,
	//                -543313364724.0),
	//                result);
	//   std::cerr << result << std::endl;
}
