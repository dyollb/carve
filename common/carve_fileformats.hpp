#pragma once

#if defined(WIN32) && !defined(CARVE_STATIC)
#  if defined(carve_fileformats_EXPORTS)
#    define CARVE_IO_API __declspec(dllexport)
#  else
#    define CARVE_IO_API __declspec(dllimport)
#  endif
#else
#    define CARVE_IO_API
#endif
