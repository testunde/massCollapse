/*
 * main.cpp
 *
 *  Created on: Dec 16, 2018
 *      Author: root
 */

#include <stdio.h>
#include <opencv4/opencv2/opencv.hpp>
#include <list>
#include <stdlib.h> // rand
#include <time.h> // seed for rand
#include <random> // normal_distribution

#include "Global.h"
#include "Droplet.h"

using namespace std;

void update() {
	// grow condensation
	Droplet *growDrop = Droplet::bigger_h;
	while (growDrop != nullptr) {
		growDrop->growCondensation();
		growDrop = growDrop->smaller;
	}

	// (plain) falling
	Droplet *fallDrop = Droplet::below_h;
	while (fallDrop != nullptr) {
		Droplet *nextFallDrop = fallDrop->above;
		fallDrop->fallBy(fallDrop->getVelocity());
		fallDrop = nextFallDrop;
	}
//	Droplet::sortListHeight();
}

void clearEnvironment() {
	Droplet::left_h->recursiveDeleteWidthList();
}

double getRandom() {
	return ((double) random()) / ((double) RAND_MAX);
}

int main(int, char**) {
	srand(time(nullptr));
	normal_distribution<double> distribution(ENVIRONMENT_SPAWN_DROP_SIZE, ENVIRONMENT_SPAWN_DROP_SIZE_STD_2 * .5);
	default_random_engine rnd_gen;

	// initialize environment
	printf("Initialize environment + generate droplets...\n");

	// generate droplets
	int dCount = 0;
	double tSmax = DBL_MIN, tSmin = DBL_MAX;
	for (int c = 0; c < ENVIRONMENT_SPAWN_DROPS_TOTAL; c++) {
		double tempCoord[] = {getRandom() * ENVIRONMENT_WIDTH,
				getRandom() * ENVIRONMENT_SPAWN_THICKNESS};
		double tempSize = 0.;
		while (abs(tempSize - ENVIRONMENT_SPAWN_DROP_SIZE) > ENVIRONMENT_SPAWN_DROP_SIZE_STD_2) {
			tempSize = distribution(rnd_gen);
		}

		if (tempSize > tSmax)
			tSmax = tempSize;
		if (tempSize < tSmin)
			tSmin = tempSize;

		Droplet * tempDrop = new Droplet(tempSize / 2., tempCoord);
		Droplet::dropList->push_back(tempDrop);

		dCount++;
		if (dCount % 100 == 0) {
			printf("%.4f%%\r", 100. * dCount / ENVIRONMENT_SPAWN_DROPS_TOTAL);
			fflush(stdout);
		}
	}
	printf("min: %f | max: %f [µm]\n", tSmin * 1.E6, tSmax * 1.E6);
	printf("Sorting...\n");
	Droplet::sortListWidth();
	Droplet::sortListHeight();
	Droplet::sortListSize();

	// init openCV
	cv::namedWindow("env_simu", cv::WINDOW_AUTOSIZE);
	int visu_width = ENVIRONMENT_WIDTH * VISU_WIDTH_PX_PER_METER;
	int visu_height = ENVIRONMENT_HEIGHT * VISU_HEIGHT_PX_PER_METER;
	cv::Mat envVisu(visu_height, visu_width, CV_8UC3);

	// simulation start
	Droplet *statDrop = Droplet::bigger_h;
	double meanTotalMassPerPx = DENSITY_WATER * pow(ENVIRONMENT_SPAWN_DROP_SIZE * .5, 3.) * (M_PI * 4. / 3.) * ENVIRONMENT_SPAWN_DROPS_TOTAL / (visu_width * visu_height);
	printf("starting simulation... (meanTotalMassPerPx/biggestMass = %f)\n", meanTotalMassPerPx/Droplet::bigger_h->getMass());
	for (int t = 0; t <= SIMULATION_TIME_MAX; t++) {
		update();

		if (t % SIMULATION_TIME_STATS_UPDATE == 0) {
			// count mass per pixel
			double total_mass[visu_width][visu_height] = {0.};
			Droplet *currentDrop = Droplet::left_h;
			while (currentDrop != nullptr) {
				int idx[2] = {(int) (currentDrop->getCoord(0) * VISU_WIDTH_PX_PER_METER), (int) (currentDrop->getCoord(1) * VISU_HEIGHT_PX_PER_METER)};
				total_mass[idx[0]][(idx[1] > visu_height) ? visu_height : idx[1]] += currentDrop->getMass();

				currentDrop = currentDrop->right;
			}

			// visualization
			double avgColor = 0.;
			for (int w = 0; w < visu_width; w++) {
				for (int h = 0; h < visu_height; h++) {
					// TODO: insert color variations depending on avgRadius distribution
					double massRatio = total_mass[w][h] / Droplet::bigger_h->getMass();
					double scaledRatio = (log10(massRatio + 0.1) + 1) / (log10(1.1) + 1);
					int colorMass = (int) (.5 * 255. * scaledRatio);//meanTotalMassPerPx);
					int colorFinal = (colorMass > 255) ? 255 : colorMass;
					envVisu.at<cv::Vec3b>(h, w) = cv::Vec3b(colorFinal, colorFinal, colorFinal);
					avgColor += colorFinal;
				}
			}
			avgColor /= visu_width * visu_height;
			cv::imshow("env_simu", envVisu);
			if (cv::waitKey(30) >= 0)
				break;

			// TODO: insert "remaining drop-count"
			printf("time: %d [s] | drop size: %f [mm] | velocity: %f [cm/s] | height: %f [m] | avgColor: %f [0-256[\n", t, statDrop->getRadius() * 2. * 1000., statDrop->getVelocity() * 100., statDrop->getCoord(1), avgColor);
			fflush(stdout);
		}
	}

	// clear simulation environment
	clearEnvironment();
}
