/*
 * main.cpp
 *
 *  Created on: Dec 16, 2018
 */

#include <stdio.h>
#include <list>
#include <stdlib.h> // rand
#include <time.h> // seed for rand
#include <random> // normal_distribution
#include <float.h> // DLB_MIN + DBL_MAX

#include "Global.h"
#include "Droplet.h"

#if USE_OPENCV
#include <opencv4/opencv2/opencv.hpp>
#endif

using namespace std;

long currentMicroSec() {
    struct timespec t;
    clock_gettime(CLOCK_MONOTONIC, &t);
    return t.tv_sec * 1E6L + t.tv_nsec / 1E3L;
}

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
		Droplet *nextDrop = fallDrop->above;
		if (fallDrop->getCoordPre(1) >= ENVIRONMENT_HEIGHT)
			fallDrop->deleteInstance();

		fallDrop = nextDrop;
	}
	Droplet::sortListHeight();

	// merging
	Droplet *mergeDrop = Droplet::bigger_h;
	while (mergeDrop != nullptr) {
		list<Droplet*> potentialDrps;

		// above
		Droplet *tempDrop = mergeDrop->above;
		double bound = mergeDrop->getCoordPre(1) - 2. * mergeDrop->getRadius();
		while (tempDrop != nullptr && (tempDrop->getCoord(1) >= bound)) {
			potentialDrps.push_back(tempDrop);
			tempDrop = tempDrop->above;
		}

		// below
		tempDrop = mergeDrop->below;
		bound = mergeDrop->getCoord(1) + 2. * mergeDrop->getRadius();
		while (tempDrop != nullptr && (tempDrop->getCoordPre(1) <= bound)) {
			potentialDrps.push_back(tempDrop);
			tempDrop = tempDrop->below;
		}

		// calculate closest distance within the last time step between mergeDrop and potential drops
		list<Droplet*> toMerge;
		double mD_velo = mergeDrop->getCoord(1) - mergeDrop->getCoordPre(1); // [m/2] assuming only height change and always falling down (+)
		double widthBound[2] = {mergeDrop->getCoord(0) - 2. * mergeDrop->getRadius(), mergeDrop->getCoord(0) + 2. * mergeDrop->getRadius()};
		for (Droplet *d : potentialDrps) {
			// check if drop even in horizontal range
			if ((d->getCoord(0) < widthBound[0]) || (d->getCoord(0) > widthBound[1]))
				continue;

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
	printf("Clearing environment...\n");
	fflush(stdout);
	for (Droplet * d: *Droplet::dropList) {
		delete d;
	}
	Droplet::dropList->clear();

	Droplet::bigger_h = nullptr; // head
	Droplet::smaller_h = nullptr; // tail
	Droplet::above_h = nullptr; // head
	Droplet::below_h = nullptr; // tail
	Droplet::left_h = nullptr; // head
	Droplet::right_h = nullptr; // tail
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
		if (dCount % ((int) ENVIRONMENT_SPAWN_DROPS_TOTAL / 10) == 0) {
			printf("%.0f%%\r", 100. * dCount / ENVIRONMENT_SPAWN_DROPS_TOTAL);
			fflush(stdout);
		}
	}
	printf("min: %f | max: %f [µm]\n", tSmin * 1.E6, tSmax * 1.E6);
	printf("Sorting...\n");
	Droplet::sortListWidth();
	Droplet::sortListHeight();
	Droplet::sortListSize();

	// init openCV
	int visu_width = ENVIRONMENT_WIDTH * VISU_WIDTH_PX_PER_METER;
	int visu_height = ENVIRONMENT_HEIGHT * VISU_HEIGHT_PX_PER_METER;
#if USE_OPENCV
	cv::namedWindow("env_simu", cv::WINDOW_AUTOSIZE);
	cv::Mat envVisu(visu_height, visu_width, CV_8UC3);
	cv::VideoWriter video(OPENCV_VIDEO_FILENAME, cv::VideoWriter::fourcc(OPENCV_VIDEO_FORMAT), OPENCV_VIDEO_SECONDS_PER_SECOND, cv::Size(visu_width, visu_height));
#endif

	// simulation start
	double meanTotalMassPerPx = DENSITY_WATER * pow(ENVIRONMENT_SPAWN_DROP_SIZE * .5, 3.) * (M_PI * 4. / 3.) * ENVIRONMENT_SPAWN_DROPS_TOTAL / (visu_width * visu_height);
	printf("starting simulation... (meanTotalMassPerPx/biggestDropMass = %f)\n", meanTotalMassPerPx/Droplet::bigger_h->getMass());
	long startTimeStamp = currentMicroSec();
	for (int t = 0; t <= SIMULATION_TIME_MAX; t++) { // [s]
		update();
		if (Droplet::disposedDrops >= ENVIRONMENT_SPAWN_DROPS_TOTAL) {
			printf("No drops remaining. Exiting...\n");
			break;
		}

		if (t % SIMULATION_TIME_STATS_UPDATE == 0) {
			// count mass per pixel
			double total_mass[visu_width][visu_height] = {0.};
			double max_mass[visu_width][visu_height] = {0.};
			bool contain_drop[visu_width][visu_height] = {false};
			// reset first column (weird bug?! because it contains the values of the prior visualization step even with init=0.)
			for (int h = 0; h < visu_height; h++) {
				total_mass[0][h] = 0.;
				max_mass[0][h] = 0.;
				contain_drop[0][h] = false;
			}
			double avgSize = 0.;
			Droplet *currentDrop = Droplet::left_h;
			while (currentDrop != nullptr) {
				int idx[2] = {(int) (currentDrop->getCoord(0) * VISU_WIDTH_PX_PER_METER), (int) (currentDrop->getCoord(1) * VISU_HEIGHT_PX_PER_METER)};
				idx[1] = (idx[1] >= visu_height) ? (visu_height - 1) : idx[1];
				total_mass[idx[0]][idx[1]] += currentDrop->getMass();
				max_mass[idx[0]][idx[1]] = max(max_mass[idx[0]][idx[1]], currentDrop->getMass());
				contain_drop[idx[0]][idx[1]] = true;

				avgSize += currentDrop->getRadius();

				currentDrop = currentDrop->right;
			}
			avgSize *= 2. * 1.E3 / (double) Droplet::remainingDrops();

			// visualization
			double avgColor = 0.;
#if USE_OPENCV
			envVisu.setTo(cv::Scalar(0, 0, 0));
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

			video.write(envVisu);
			cv::imshow("env_simu", envVisu);
			if (cv::waitKey(30) >= 0)
				break;
#endif

			long currentTimeStamp = currentMicroSec();
			long eta_seconds = (long) ((SIMULATION_TIME_MAX - (long) t) * (currentTimeStamp - startTimeStamp) / (1E6L * (long) t));
			printf("time: %d [s] | drop size (avg/max): %f/%f [mm] | avgColor: %f [0-256[ | remaining drops: %d | lowest height: %f [m] (ETA: %ldm %lds)\n",
					t, avgSize, Droplet::bigger_h->getRadius() * 2. * 1.E3, avgColor, Droplet::remainingDrops(), (Droplet::below_h == nullptr) ? -1 : (ENVIRONMENT_HEIGHT - Droplet::below_h->getCoord(1)),
					eta_seconds / 60, eta_seconds % 60);
			fflush(stdout);
		}
	}

	// clear simulation environment
#if USE_OPENCV
	video.release();
#endif
	clearEnvironment();
}
