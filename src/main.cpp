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
}

void clearEnvironment() {
	Droplet::left_h->recursiveDeleteWidthList();
}

double getRandom() {
	return ((double) random()) / ((double) RAND_MAX);
}

void insertInListLinks(Droplet *drp) {
	// left/right list links & tail/heads
	Droplet *preRight = nullptr, *currentRight = Droplet::left_h;
	while (currentRight != nullptr && currentRight->getCoord(0) < drp->getCoord(0)) {
		preRight = currentRight;
		currentRight = currentRight->right;
	}
	if (preRight == nullptr) {
		Droplet::left_h = drp;
	} else {
		preRight->right = drp;
	}
	if (currentRight == nullptr) {
		Droplet::right_h = drp;
	} else {
		currentRight->left = drp;
	}
	drp->right = currentRight;
	drp->left = preRight;

	// above/below list links & tail/heads
	Droplet *preBelow = nullptr, *currentBelow = Droplet::above_h;
	while (currentBelow != nullptr && currentBelow->getCoord(1) < drp->getCoord(1)) {
		preBelow = currentBelow;
		currentBelow = currentBelow->below;
	}
	if (preBelow == nullptr) {
		Droplet::above_h = drp;
	} else {
		preBelow->below = drp;
	}
	if (currentBelow == nullptr) {
		Droplet::below_h = drp;
	} else {
		currentBelow->above = drp;
	}
	drp->below = currentBelow;
	drp->above = preBelow;

	// bigger/smaller list links & tail/heads
	Droplet *preSmaller = nullptr, *currentSmaller = Droplet::bigger_h;
	while (currentSmaller != nullptr && currentSmaller->getRadius() > drp->getRadius()) {
		preSmaller = currentSmaller;
		currentSmaller = currentSmaller->smaller;
	}
	if (preSmaller == nullptr) {
		Droplet::bigger_h = drp;
	} else {
		preSmaller->smaller = drp;
	}
	if (currentSmaller == nullptr) {
		Droplet::smaller_h = drp;
	} else {
		currentSmaller->bigger = drp;
	}
	drp->smaller = currentSmaller;
	drp->bigger = preSmaller;
}

int main(int, char**) {
	srand(time(nullptr));
	normal_distribution<double> distribution(ENVIRONMENT_SPAWN_DROP_SIZE, ENVIRONMENT_SPAWN_DROP_SIZE_STD_2 * .5);
	default_random_engine rnd_gen;

	// initialize environment
	printf("Initialize environment + generate droplets...\n");

	// generate droplets
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
		insertInListLinks(tempDrop);
	}
	printf("min: %f | max: %f [µm]\n", tSmin * 1.E6, tSmax * 1.E6);

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

			printf("time: %d [s] | drop size: %f [mm] | velocity: %f [cm/s] | height: %f [m] | avgColor: %f [0-256[\n", t, statDrop->getRadius() * 2. * 1000., statDrop->getVelocity() * 100., statDrop->getCoord(1), avgColor);
			fflush(stdout);
		}
	}

	// clear simulation environment
	clearEnvironment();
}
