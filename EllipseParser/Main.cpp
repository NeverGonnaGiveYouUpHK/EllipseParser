#define _USE_MATH_DEFINES
#include <libtiff/tiffio.h>

#include <iostream>
#include <fstream>
#include <chrono>
#include <math.h>

#include <filesystem>
namespace fs = std::filesystem;
#include <Eigen/Dense>
using namespace Eigen;


struct Point {
	uint32 x;
	uint32 y;
};

struct ellipseABCDEF {
	double A;
	double B;
	double C;
	double D;
	double E;
	double F;
};

struct ellipseABCDEF fitEllipse(struct Point p1, struct Point p2, struct Point p3, struct Point p4, struct Point p5) {
	Matrix<double, 5, 6> mat56(5, 6);
	mat56 << 
		(double)p1.x * (double)p1.x, (double)p1.x * -(double)p1.y, (double)p1.y * (double)p1.y, (double)p1.x, -(double)p1.y, 1,
		(double)p2.x * (double)p2.x, (double)p2.x * -(double)p2.y, (double)p2.y * (double)p2.y, (double)p2.x, -(double)p2.y, 1,
		(double)p3.x * (double)p3.x, (double)p3.x * -(double)p3.y, (double)p3.y * (double)p3.y, (double)p3.x, -(double)p3.y, 1,
		(double)p4.x * (double)p4.x, (double)p4.x * -(double)p4.y, (double)p4.y * (double)p4.y, (double)p4.x, -(double)p4.y, 1,
		(double)p5.x * (double)p5.x, (double)p5.x * -(double)p5.y, (double)p5.y * (double)p5.y, (double)p5.x, -(double)p5.y, 1;
	
	CompleteOrthogonalDecomposition<Matrix<double, Dynamic, Dynamic> > cod;
	cod.compute(mat56);

	// Find URV^T
	MatrixXd V = cod.matrixZ().transpose();
	MatrixXd Null_space = V.block(0, cod.rank(), V.rows(), V.cols() - cod.rank());
	MatrixXd P = cod.colsPermutation();
	Null_space = P * Null_space; // Unpermute the columns


	return {
		Null_space(0) / Null_space(5),
		Null_space(1) / Null_space(5),
		Null_space(2) / Null_space(5),
		Null_space(3) / Null_space(5),
		Null_space(4) / Null_space(5),
		1
	};
}


static uint16 finalAverage;
static std::ofstream outStream("out.csv");

typedef unsigned long long uint64;



struct followReturn {
	bool success;
	struct Point point;
};

bool validateEdge(bool* considered, struct Point point, uint32 width, uint32 height) {
	struct Point nearby[4] = {
		{point.x + 0, point.y - 1},
		{point.x - 1, point.y + 0},
		{point.x + 1, point.y + 0},
		{point.x + 0, point.y + 1}
	};

	for (uint32 i = 0; i < 4; i++) {
		struct Point tested = nearby[i];

		if (tested.x != 0xFFFFFFFF && tested.y != 0xFFFFFFFF && tested.x != width && tested.y != height) {
			uint32 index = tested.y * width + tested.x;

			if (!considered[index]) {
				return true;
			}
		}
	}

	return false;
}

struct EllipseParams {
	double fCenterX, fCenterY;
	double fLongLength, fShortLength;
	double fAngle;
};

EllipseParams toEllipseParams(ellipseABCDEF a) {

	EllipseParams eOutParams;

	eOutParams.fCenterX = ((2 * a.C * a.D) - (a.B * a.E)) / ((a.B * a.B) - (4 * a.A * a.C));
	eOutParams.fCenterY = ((2 * a.A * a.E) - (a.B * a.D)) / ((a.B * a.B) - (4 * a.A * a.C));

	eOutParams.fLongLength = -sqrt(2 * (a.A * pow(a.E, 2) + a.C * pow(a.D, 2) - a.B * a.D * a.E + (pow(a.B, 2) - 4 * a.A * a.C) * a.F) * ((a.A + a.C) + sqrt(pow((a.A - a.C), 2) + pow(a.B, 2)))) /
							 ((a.B * a.B) - (4 * a.A * a.C));

	eOutParams.fShortLength = -sqrt(2 * (a.A * pow(a.E, 2) + a.C * pow(a.D, 2) - a.B * a.D * a.E + (pow(a.B, 2) - 4 * a.A * a.C) * a.F) * ((a.A + a.C) - sqrt(pow((a.A - a.C), 2) + pow(a.B, 2)))) /
							 ((a.B * a.B) - (4 * a.A * a.C));

	if (a.B != 0)
		eOutParams.fAngle = atan((a.C - a.A - sqrt(pow((a.A - a.C), 2) + pow(a.B, 2))) / a.B);
	else if
		(a.A <= a.C) eOutParams.fAngle = 0;
	else
		eOutParams.fAngle = M_PI / 2;

	return eOutParams;
}

double rad2Deg(double rad) {
	return (rad * 180) / M_PI;
}

//filename, ellipse_center_x, ellipse_center_y, ellipse_majoraxis, ellipse_minoraxis, ellipse_angle, elapsed_time

void saveEllipse(bool empty, EllipseParams params, float elapsedTime, char* fileName){
	
	if (empty) {

		outStream << fileName << ",";
		outStream << ",";
		outStream << ",";
		outStream << ",";
		outStream << ",";
		outStream << ",";
		outStream << (int)(elapsedTime * 1000) << std::endl;
		return;
	}

	outStream << fileName << ",";
	outStream << params.fCenterX << ",";
	outStream << -params.fCenterY << ",";
	outStream << params.fLongLength << ",";
	outStream << params.fShortLength << ",";
	outStream << fmod(180 - rad2Deg(params.fAngle), 360) << ",";
	outStream << (int)(elapsedTime * 1000) << std::endl;

}


double ratePoint(struct Point point, struct EllipseParams ellipse) {
	double xi = (double)point.x;
	double yi = (double)point.y;

	xi = xi - ellipse.fCenterX;
	yi = yi + ellipse.fCenterY;

	double x = xi * cos(ellipse.fAngle) + yi * sin(ellipse.fAngle);
	double y = xi * -sin(ellipse.fAngle) + yi * cos(ellipse.fAngle);

	double angle = atan(y / x);

	double distanceFromCenter = hypot(x, y);

	double distanceFromCenterToEdge = ellipse.fLongLength * ellipse.fShortLength / (sqrt(pow(ellipse.fShortLength * cos(angle), 2) + pow(ellipse.fLongLength * sin(angle), 2)));

	double finalDistance = abs(distanceFromCenter - distanceFromCenterToEdge);

	//some gaussian trickery right here
	double rating = exp(-(finalDistance * finalDistance) / 8);
	return rating;
}


void ParseEllipse(std::string path) { 
	TIFF* tif = TIFFOpen(path.c_str(), "r");

	if (tif == nullptr) {
		std::cout << path << " is not a valid TIFF image!" << std::endl << std::endl;
		return;
	}

	uint32 width, height;

	TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);           // uint32 width;
	TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);        // uint32 height;

	std::cout << std::endl << "Image: " << path << std::endl;
	printf("Width: %d\nHeight: %d\n", width, height);

	uint16* raster = (uint16*)_TIFFmalloc(width * height * sizeof(uint16));

	for (uint32 row = 0; row < height; row++) {
		TIFFReadScanline(tif, raster + width * row, row, 0);
	}

	//Time Taken
	auto timeStart = std::chrono::system_clock::now();

	uint64 average = 0;
	for (uint32 i = 0; i < width * height; i++) {
		average += uint64(raster[i]);
	}

	finalAverage = average / (uint64(width) * uint64(height)) * 1;

	/*uint32 maxLuminance = 0;

	for (uint32 i = 0; i < width * height; i++) {
		if (raster[i] > maxLuminance) {
			maxLuminance = raster[i];
		}
	}

	double fmaxLuminance = (double)maxLuminance;

	double average = 0;

	for (uint32 i = 0; i < width * height; i++) {
		average += (1 - pow(exp(-(pow((double)raster[i] / fmaxLuminance, 2) / 0.5)), 64)) * fmaxLuminance;
	}

	average /= width * height;


	finalAverage = (uint16)average;*/

	bool* considered = (bool*)malloc(width * height * sizeof(bool));
	bool* alreadyTaken = (bool*)calloc(width * height, sizeof(bool));

	for (uint32 y = 0; y < height; y++) {
		for (uint32 x = 0; x < width; x++) {
			const uint32 index = y * width + x;

			//considered[index] = bool(raster[index] < fmaxLuminance * 0.45);
			considered[index] = bool(raster[index] < finalAverage);
		}
	}


	std::vector<std::vector<struct Point>> sequences;

	std::vector<struct Point> paths[2] = {
		std::vector<struct Point>(),
		std::vector<struct Point>()
	};



	for (uint32 y = 0; y < height; y++) {
		for (uint32 x = 0; x < width; x++) {
			uint32 indexInitial = y * width + x;

			if (alreadyTaken[indexInitial] || !considered[indexInitial]) {
				continue;
			}

			if (!validateEdge(considered, { x, y }, width, height)) {
				continue;
			}

			uint32 newX = x, newY = y;

			struct Point nearbyInitial[8] = {
				{newX - 1, newY - 1},
				{newX + 0, newY - 1},
				{newX + 1, newY - 1},
				{newX - 1, newY + 0},
				{newX + 1, newY + 0},
				{newX - 1, newY + 1},
				{newX + 0, newY + 1},
				{newX + 1, newY + 1}
			};

			uint32 part = 0;
			uint32 firstPointIndex = 0;



			while (part < 2 && firstPointIndex < 8) {
				newX = nearbyInitial[firstPointIndex].x;
				newY = nearbyInitial[firstPointIndex].y;
				firstPointIndex++;

				if (newX == 0xFFFFFFFF || newY == 0xFFFFFFFF || newX == width || newY == height) {
					continue;
				}
				if (!validateEdge(considered, { newX, newY }, width, height)) {
					continue;
				}

				bool endPath = false;
				do {
					uint32 indexThis = newY * width + newX;


					struct Point nearby[8] = {
						{newX - 1, newY - 1},
						{newX + 0, newY - 1},
						{newX + 1, newY - 1},
						{newX - 1, newY + 0},
						{newX + 1, newY + 0},
						{newX - 1, newY + 1},
						{newX + 0, newY + 1},
						{newX + 1, newY + 1}
					};

					endPath = true;

					for (uint32 i = 0; i < 8; i++) {
						struct Point tested = nearby[i];

						if (tested.x != 0xFFFFFFFF && tested.y != 0xFFFFFFFF && tested.x != width && tested.y != height) {
							if (!validateEdge(considered, tested, width, height)) {
								continue;
							}

							uint32 testedIndex = tested.y * width + tested.x;

							if (alreadyTaken[testedIndex] || !considered[testedIndex]) {
								continue;
							}


							paths[part].push_back({ newX, newY });


							newX = tested.x;
							newY = tested.y;

							endPath = false;

							alreadyTaken[testedIndex] = true;
							break;
						}
					}
					

				} while (!endPath);

				part++;
			}

			uint32 sizePart1 = paths[0].size();
			uint32 sizePart2 = paths[1].size();

			if (sizePart1 + sizePart2 >= 10) {
				std::vector<struct Point> finalPath = std::vector<struct Point>();

				finalPath.reserve(sizePart1 + sizePart2);

				for (int i = sizePart1 - 1; i >= 0; i--) {
					finalPath.push_back(paths[0][i]);

				}

				int a = 0;

				for (int i = 0; i < sizePart2; i++) {
					finalPath.push_back(paths[1][i]);
				}

				paths[0] = std::vector<struct Point>();
				paths[1] = std::vector<struct Point>();

				sequences.push_back(finalPath);
			}
			else {
				paths[0].clear();
				paths[1].clear();
			}
		}
	}

	if (sequences.size() == 0) {

		auto timeEnd = std::chrono::system_clock::now();

		std::chrono::duration<float> elapsedTime = timeEnd - timeStart;
		std::cout << "No ellipse found!" << std::endl;

		saveEllipse(true, {}, elapsedTime.count(), (char*)fs::path(path).filename().string().c_str());
		return;
	}

	uint32 maxLength = 0;
	uint32 longestIndex = 0;
	for (uint32 i = 0; i < sequences.size(); i++) {
		if (sequences[i].size() > maxLength) {
			maxLength = sequences[i].size();
			longestIndex = i;
		}
	}

	std::vector<struct Point> longest = sequences[longestIndex];



	//Sampling
	Point pSamples[5];
	const int nMaxSampleItterations = 4; //Linear addition

	struct EllipseParams bestFit;
	double bestFitRating = 0;

	//Iterate through
	for (int nCurrItteration = 0; nCurrItteration < nMaxSampleItterations; nCurrItteration++) {

		//Populate samples 
		int nPointDistance = longest.size() / (5 + 2 * nCurrItteration);

		if (nPointDistance == 0) {
			break;
		}

		for (int i = 0; i < 5; i++) {
			pSamples[i] = longest[i * nPointDistance];
		}


		int nCurrOffset = 0;
		//Go through the path

		while ((4 * nPointDistance + nCurrOffset) < longest.size()) {

			for (int i = 0; i < 5; i++) {
				pSamples[i] = longest[i * nPointDistance + nCurrOffset];
			}

			EllipseParams eFitEllipse = toEllipseParams(fitEllipse(pSamples[0], pSamples[1], pSamples[2], pSamples[3], pSamples[4]));
			


			double rating = 0;
			for (int i = 0; i < longest.size(); i++) {
				rating += ratePoint(longest[i], eFitEllipse);
			}

			if (rating > bestFitRating) {
				bestFitRating = rating;
				bestFit = eFitEllipse;
			}

			//Update offset for next iteration
			nCurrOffset++;
		}
	}

	std::cout << "Best fit got rated: " << bestFitRating << " (out of " << longest.size() << ")" << std::endl;

	auto timeEnd = std::chrono::system_clock::now();

	std::chrono::duration<float> elapsedTime = timeEnd - timeStart;
	std::cout << "Finding the best ellipse took: " << elapsedTime.count() << " s" << std::endl;

	saveEllipse(false, bestFit, elapsedTime.count(), (char*)fs::path(path).filename().string().c_str());

	TIFFClose(tif);
	//free mallocs

	_TIFFfree(raster);
	free(considered);
	free(alreadyTaken);
}

int main(int argc, char** argv) {

	outStream << "filename,ellipse_center_x,ellipse_center_y,ellipse_majoraxis,ellipse_minoraxis,ellipse_angle,elapsed_time" << std::endl; //write the header

	if (argc == 1) {
		std::cout << "No input directory provided! Exiting..." << std::endl;
		return 0;
	}

	std::string path(argv[1]);
	try {
		for (const auto& entry : fs::recursive_directory_iterator(path))
			ParseEllipse(entry.path().string());
	}
	catch (std::exception e) {
		std::cout << argv[1] << " is not a valid folder!" << std::endl;
	}

	outStream.close();

	return 0;
}