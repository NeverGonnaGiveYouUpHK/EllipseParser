#include <libtiff/tiffio.h>

#include <iostream>

#include <math.h>

#include <Eigen/Dense>
using namespace Eigen;

#define OLC_PGE_APPLICATION
#include "olcPixelGameEngine.h"

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
		p1.x * p1.x, p1.y * p1.y, p1.x * p1.y, p1.x, p1.y, 1,
		p2.x * p2.x, p2.y * p2.y, p2.x * p2.y, p2.x, p2.y, 1,
		p3.x * p3.x, p3.y * p3.y, p3.x * p3.y, p3.x, p3.y, 1,
		p4.x * p4.x, p4.y * p4.y, p4.x * p4.y, p4.x, p4.y, 1,
		p5.x * p5.x, p5.y * p5.y, p5.x * p5.y, p5.x, p5.y, 1;
	
	CompleteOrthogonalDecomposition<Matrix<double, 5, 6> > cod;
	cod.compute(mat56);
	std::cout << "rank : " << cod.rank() << "\n";
	// Find URV^T
	MatrixXd V = cod.matrixZ().transpose();
	MatrixXd Null_space = V.block(0, cod.rank(), V.rows(), V.cols() - cod.rank());
	MatrixXd P = cod.colsPermutation();
	Null_space = P * Null_space; // Unpermute the columns
	// The Null space:
	std::cout << "The null space: \n" << Null_space << "\n";
	// Check that it is the null-space:
	std::cout << "mat56 * Null_space = \n" << mat56 * Null_space << '\n';

	return {
		Null_space(0),
		Null_space(1),
		Null_space(2),
		Null_space(3),
		Null_space(4),
		Null_space(5)
	};
}


static uint16 finalAverage;

class sampletext : public olc::PixelGameEngine {
public:
	uint16* content;

	sampletext(uint16* raster) {
		content = raster;
	}
	
	virtual bool OnUserUpdate(float fElapsedTime) override {
		return true;
	}

	virtual bool OnUserCreate() override {
		for (uint32 y = 0; y < ScreenHeight(); y++) {
			for (uint32 x = 0; x < ScreenWidth(); x++) {
				Draw(x, y, olc::Pixel(content[y * ScreenWidth() + x] > finalAverage ? 0 : 255, 0, 0));
			}
		}
		
		return true;
	}
};

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


int main(int argc, char** argv) {
	TIFF* tif = TIFFOpen("./test/2018-02-15 19.22.13.141000.tiff", "r");

	uint32 width, height;

	TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &width);           // uint32 width;
	TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &height);        // uint32 height;

	printf("width: %d\nheight: %d\n", width, height);

	uint16* raster = (uint16*)_TIFFmalloc(width * height * sizeof(uint16));

	

	for (uint32 row = 0; row < height; row++) {
		TIFFReadScanline(tif, raster + width * row, row, 0);
	}
	
	uint64 average = 0;
	for (uint32 i = 0; i < width * height; i++) {
		average += uint64(raster[i]);
	}

	std::cout << average / (uint64(width) * uint64(height)) << std::endl;

	finalAverage = average / (uint64(width) * uint64(height)) * 1;

	std::cout << finalAverage << std::endl;



	bool* considered = (bool*)malloc(width * height * sizeof(bool));
	bool* alreadyTaken = (bool*)calloc(width * height, sizeof(bool));

	for (uint32 y = 0; y < height; y++) {
		for (uint32 x = 0; x < width; x++) {
			const uint32 index = y * width + x;
			considered[index] = bool(raster[index] < finalAverage);
		}
	}

	std::vector<std::vector<struct Point>> sequences;

	std::vector<struct Point> paths[2] = {
		std::vector<struct Point>(),
		std::vector<struct Point>()
	};

	std::cout << validateEdge(considered, { 393, 341 }, width, height);


	for (uint32 y = 0; y < height; y++) {
		for (uint32 x = 0; x < width; x++) {
			uint32 indexInitial = y * width + x;

			//std::cout << (alreadyTaken[indexInitial] || !considered[indexInitial]) << std::endl;

			if (alreadyTaken[indexInitial] || !considered[indexInitial]) {
				continue;
			}
			else {
				//printf("sample\n");
			}

			if (!validateEdge(considered, { x, y }, width, height)) {
				continue;
			}

			uint32 newX = x, newY = y;

			printf("newX %lu newY %lu\n", newX, newY);

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


							printf("adding %lu %lu\n", newX, newY);

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
					printf("%lu\n", i);
					finalPath.push_back(paths[0][i]);
					
				}

				int a = 0;

				for (int i = 0; i < sizePart2; i++) {
					finalPath.push_back(paths[1][i]);
				}

				paths[0] = std::vector<struct Point>();
				paths[1] = std::vector<struct Point>();
		
				sequences.push_back(finalPath);
			} else {
				paths[0].clear();
				paths[1].clear();
			}
		}
	}

	printf("%d paths\n", sequences.size());

	if (sequences.size() == 0) {
		//NO ELLIPSE AT ALL
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







	sampletext game = sampletext(raster);

	if (game.Construct(width, height, 1, 1, false, false))
		game.Start();

	TIFFClose(tif);
	return 0;
}

