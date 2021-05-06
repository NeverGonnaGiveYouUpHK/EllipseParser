#include <libtiff/tiffio.h>

#include <iostream>

#include <math.h>

#define OLC_PGE_APPLICATION
#include "olcPixelGameEngine.h"

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

struct Point {
	uint32 x;
	uint32 y;
};

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
	TIFF* tif = TIFFOpen("./test/2018-02-15 17.55.41.631000.tiff", "r");

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

	finalAverage = average / (uint64(width) * uint64(height)) * 28;

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

	sequences.push_back(std::vector<struct Point>());

	int pathCount = 0;

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
			bool isEmptyFlag = true;
			uint32 part = 0;
			uint32 firstPointIndex = 0;
			while (part < 2 && firstPointIndex < 8) {
				newX = nearbyInitial[firstPointIndex].x;
				newY = nearbyInitial[firstPointIndex].y;
				firstPointIndex++;

				if (newX == 0xFFFFFFFF || newY == 0xFFFFFFFF || newX == width || newY == height) {
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

							newX = tested.x;
							newY = tested.y;

							printf("%lu %lu\n", newX, newY);
							isEmptyFlag = false;

							endPath = false;

							alreadyTaken[testedIndex] = true;
							break;
						}
					}


				} while (!endPath);

				part++;
			}
			//printf("%d paths\n", pathCount);
			if (!isEmptyFlag) {
				printf("%d paths\n", pathCount);
				pathCount++;
			}
			
		}
	}

	printf("%d paths\n", pathCount);






	sampletext game = sampletext(raster);

	if (game.Construct(width, height, 1, 1, false, false))
		game.Start();

	TIFFClose(tif);
	return 0;
}

