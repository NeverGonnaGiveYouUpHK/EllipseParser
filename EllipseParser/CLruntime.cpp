#include "CLruntime.h"

inline cl_int initCL() {
	cl_int err;

	cl_context_properties clContextProperties = {};

	cl::Context clContext(
		CL_DEVICE_GPI
	);
}

inline cl_int compileKernel(cl::Context context, const char* path) {

	std::ifstream kernelFileStream;

	kernelFileStream.exceptions(std::ifstream::badbit | std::ifstream::failbit);

	std::string kernelOutput;

	try {
		kernelFileStream.open(path);

		std::stringstream kernelStringStream;
		kernelStringStream << kernelFileStream.rdbuf();

		// close file handlers 
		kernelFileStream.close();

		// Convert stream into GLchar array 
		kernelOutput = kernelStringStream.str();
	}
	catch (std::ifstream::failure e) {
		return CL_LINK_PROGRAM_FAILURE;
	}

	cl::Program::Sources kernelSource(
		1,
		std::make_pair(kernelOutput.c_str(), kernelOutput.size() + 1)
	);


}