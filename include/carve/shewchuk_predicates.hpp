
#pragma once

#include "carve.hpp"

namespace shewchuk {
CARVE_API double orient2dfast(const double* pa, const double* pb, const double* pc);
CARVE_API double orient2dexact(const double* pa, const double* pb, const double* pc);
CARVE_API double orient2dslow(const double* pa, const double* pb, const double* pc);
CARVE_API double orient2dadapt(const double* pa, const double* pb, const double* pc, double detsum);
CARVE_API double orient2d(const double* pa, const double* pb, const double* pc);

CARVE_API double orient3dfast(const double* pa, const double* pb, const double* pc, const double* pd);
CARVE_API double orient3dexact(const double* pa, const double* pb, const double* pc, const double* pd);
CARVE_API double orient3dslow(const double* pa, const double* pb, const double* pc, const double* pd);
CARVE_API double orient3dadapt(const double* pa, const double* pb, const double* pc, const double* pd, double permanent);
CARVE_API double orient3d(const double* pa, const double* pb, const double* pc, const double* pd);

CARVE_API double incirclefast(const double* pa, const double* pb, const double* pc, const double* pd);
CARVE_API double incircleexact(const double* pa, const double* pb, const double* pc, const double* pd);
CARVE_API double incircleslow(const double* pa, const double* pb, const double* pc, const double* pd);
CARVE_API double incircleadapt(const double* pa, const double* pb, const double* pc, const double* pd, double permanent);
CARVE_API double incircle(const double* pa, const double* pb, const double* pc, const double* pd);

CARVE_API double inspherefast(const double* pa, const double* pb, const double* pc, const double* pd, const double* pe);
CARVE_API double insphereexact(const double* pa, const double* pb, const double* pc, const double* pd, const double* pe);
CARVE_API double insphereslow(const double* pa, const double* pb, const double* pc, const double* pd, const double* pe);
CARVE_API double insphereadapt(const double* pa, const double* pb, const double* pc, const double* pd, const double* pe, double permanent);
CARVE_API double insphere(const double* pa, const double* pb, const double* pc, const double* pd, const double* pe);
} // namespace shewchuk
