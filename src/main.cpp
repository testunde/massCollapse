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
	Droplet *fallDrop = Droplet::below_h; // it does not matter which dimension
	while (fallDrop != nullptr) {
		fallDrop->fallBy(fallDrop->getVelocity());
		fallDrop = fallDrop->above;;
	}
	Droplet::sortListHeight();

	// merging
	Droplet *mergeDrop = Droplet::bigger_h;
	while (mergeDrop != nullptr) {
		list<Droplet*> potentialDrps;

		// to the left
		Droplet *tempDrop = mergeDrop->left;
		double bound = mergeDrop->getCoord(0) - 2. * mergeDrop->getRadius();
		while (tempDrop != nullptr && (tempDrop->getCoord(0) >= bound)) {
			potentialDrps.push_back(tempDrop);
			tempDrop = tempDrop->left;
		}

		// to the right
		tempDrop = mergeDrop->right;
		bound = mergeDrop->getCoord(0) + 2. * mergeDrop->getRadius();
		while (tempDrop != nullptr && (tempDrop->getCoord(0) <= bound)) {
			potentialDrps.push_back(tempDrop);
			tempDrop = tempDrop->right;
		}

		// calculate closest distance within the last time step between mergeDrop and potential drops
		list<Droplet*> toMerge;
		double mD_velo = mergeDrop->getCoord(1) - mergeDrop->getCoordPre(1); // [m/2] assuming only height change and always falling down (+)
		for (Droplet *d : potentialDrps) {
			double d_velo = d->getCoord(1) - d->getCoordPre(1); // [m/2]
			double timeClosest = (mergeDrop->getCoordPre(1) - d->getCoordPre(1)) /
					(d_velo - mD_velo); // [s]
			// limit time point to [0s, 1s]
			timeClosest = (timeClosest < 0.) ? 0 : ((timeClosest > 1.) ? 1. : timeClosest);

			double distance = sqrt(pow((mergeDrop->getCoordPre(1) + timeClosest * mD_velo) - (d->getCoordPre(1) + timeClosest * d_velo), 2.)
					+ pow(mergeDrop->getCoordPre(0) - d->getCoordPre(0), 2.)); // [m]

			if (distance <= mergeDrop->getRadius() + d->getRadius())
				toMerge.push_back(d);
		}

		mergeDrop->merge(&toMerge);
		mergeDrop = mergeDrop->smaller;
	}
	Droplet::sortListSize();
	Droplet::sortListWidth();
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
	double meanTotalMassPerPx = DENSITY_WATER * pow(ENVIRONMENT_SPAWN_DROP_SIZE * .5, 3.) * (M_PI * 4. / 3.) * ENVIRONMENT_SPAWN_DROPS_TOTAL / (visu_width * visu_height);
	printf("starting simulation... (meanTotalMassPerPx/biggestDropMass = %f)\n", meanTotalMassPerPx/Droplet::bigger_h->getMass());
	for (int t = 0; t <= SIMULATION_TIME_MAX; t++) {
		update();

		if (t % SIMULATION_TIME_STATS_UPDATE == 0) {
			// count mass per pixel
			double total_mass[visu_width][visu_height] = {0.};
			double max_mass[visu_width][visu_height] = {0.};
			bool contain_drop[visu_width][visu_height] = {false};
			double avgSize = 0.;
			Droplet *currentDrop = Droplet::left_h;
			while (currentDrop != nullptr) {
				int idx[2] = {(int) (currentDrop->getCoord(0) * VISU_WIDTH_PX_PER_METER), (int) (currentDrop->getCoord(1) * VISU_HEIGHT_PX_PER_METER)};
				idx[1] = (idx[1] > visu_height) ? visu_height : idx[1];
				total_mass[idx[0]][idx[1]] += currentDrop->getMass();
				max_mass[idx[0]][idx[1]] = max(max_mass[idx[0]][idx[1]], currentDrop->getMass());
				contain_drop[idx[0]][idx[1]] = true;

				avgSize += currentDrop->getRadius();

				currentDrop = currentDrop->right;
			}
			avgSize *= 2. / (double) (ENVIRONMENT_SPAWN_DROPS_TOTAL - Droplet::mergeCount);

			// visualization
			double avgColor = 0.;
			for (int w = 0; w < visu_width; w++) {
				for (int h = 0; h < visu_height; h++) {
					double totalMassRatio = total_mass[w][h] / meanTotalMassPerPx;//Droplet::bigger_h->getMass();
					totalMassRatio = (log10(totalMassRatio + 0.1) + 1) / (log10(1.1) + 1);
					int colorTotalMass = (int) (.68 * 255. * totalMassRatio);
					colorTotalMass = (colorTotalMass > 255) ? 255 : colorTotalMass;

					double maxMassRatio = max_mass[w][h] / Droplet::bigger_h->getMass();
					maxMassRatio = (log10(maxMassRatio + 0.1) + 1) / (log10(1.1) + 1);
					int colorMaxMass = (int) (255. * maxMassRatio);
					colorMaxMass = (colorMaxMass > 255) ? 255 : colorMaxMass;

					int colorContainDrop = contain_drop[w][h] ? 255 : 0;

					envVisu.at<cv::Vec3b>(h, w) = cv::Vec3b(colorMaxMass, colorMaxMass, colorMaxMass);
					avgColor += sqrt(pow(envVisu.at<cv::Vec3b>(h, w)[0], 2.) +
							pow(envVisu.at<cv::Vec3b>(h, w)[1], 2.) +
							pow(envVisu.at<cv::Vec3b>(h, w)[2], 2.));
				}
			}
			avgColor /= visu_width * visu_height;

			cv::imshow("env_simu", envVisu);
			if (cv::waitKey(30) >= 0)
				break;

			printf("time: %d [s] | avg drop size: %f [mm] | avgColor: %f [0-256[ | remaining drops: %d | lowest height: %f [m]\n",
					t, avgSize, avgColor, ((int) ENVIRONMENT_SPAWN_DROPS_TOTAL) - Droplet::mergeCount, ENVIRONMENT_HEIGHT - Droplet::below_h->getCoord(1));
			fflush(stdout);
		}
	}

	// clear simulation environment
	clearEnvironment();
}
