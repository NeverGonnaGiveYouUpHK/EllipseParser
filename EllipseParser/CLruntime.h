#pragma once

#include <CL/cl.hpp>
#include <iostream>
#include <fstream>
#include <sstream>

inline cl_int initCL();
inline cl_int getPlatforms();
inline cl_int getDevices();

inline cl_int compileKernel(cl::Context context, const char* path);