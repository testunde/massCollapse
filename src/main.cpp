/*
 * main.cpp
 *
 *  Created on: Jan 12, 2019
 */

#include <float.h> // DLB_MIN + DBL_MAX
#include <random>  // normal_distribution
#include <stdio.h>
#include <stdlib.h> // rand
#include <time.h>   // seed for rand
#include <vector>

#include "Global.h"
#include "Particle.h"

#ifdef USE_OPENCV
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
    Particle *growDrop = Particle::bigger_h;
    while (growDrop != nullptr) {
        growDrop->growCondensation();
        growDrop = growDrop->smaller;
    }

    // (plain) falling
    vector<Particle *> to_delete{};
    for (Particle *fallDrop : *Particle::dropList) {
        fallDrop->fallBy(fallDrop->getVelocity());
        if (fallDrop->getCoordPre(1) >= ENVIRONMENT_HEIGHT) {
            to_delete.push_back(fallDrop);
        }
    }
    for (Particle *fallDrop : to_delete) {
        fallDrop->deleteInstance();
    }
    Particle::sortListHeight();

    // merging
    Particle *mergeDrop = Particle::bigger_h;
    while (mergeDrop != nullptr) {
        vector<Particle *> potentialDrps;

        // above
        Particle *tempDrop = mergeDrop->above;
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

        // calculate closest distance within the last time step between
        // mergeDrop and potential drops
        vector<Particle *> toMerge;
        // [m/2] assuming only height change and always falling down (+)
        double mD_velo = mergeDrop->getCoord(1) - mergeDrop->getCoordPre(1);
        double widthBound[2] = {
            mergeDrop->getCoord(0) - 2. * mergeDrop->getRadius(),
            mergeDrop->getCoord(0) + 2. * mergeDrop->getRadius()};
        for (Particle *d : potentialDrps) {
            // check if drop even in horizontal range
            if ((d->getCoord(0) < widthBound[0]) ||
                (d->getCoord(0) > widthBound[1]))
                continue;

            double d_velo = d->getCoord(1) - d->getCoordPre(1); // [m/2]
            double timeClosest =
                (mergeDrop->getCoordPre(1) - d->getCoordPre(1)) /
                (d_velo - mD_velo); // [s]
            // limit time point to [0s, 1s]
            timeClosest = (timeClosest < 0.)
                              ? 0
                              : ((timeClosest > 1.) ? 1. : timeClosest);

            double distHeight =
                (mergeDrop->getCoordPre(1) + timeClosest * mD_velo) -
                (d->getCoordPre(1) + timeClosest * d_velo);
            double distWidth = mergeDrop->getCoordPre(0) - d->getCoordPre(0);
            double distance =
                sqrt(distHeight * distHeight + distWidth * distWidth); // [m]

            if (distance <= mergeDrop->getRadius() + d->getRadius())
                toMerge.push_back(d);
        }

        mergeDrop->merge(&toMerge);
        mergeDrop = mergeDrop->smaller;
    }
    Particle::sortListSize();
    Particle::sortListWidth();
}

void clearEnvironment() {
    printf("Clearing environment...\n");
    fflush(stdout);
    for (Particle *d : *Particle::dropList) {
        delete d;
    }
    Particle::dropList->clear();

    Particle::bigger_h = nullptr;  // head
    Particle::smaller_h = nullptr; // tail
    Particle::above_h = nullptr;   // head
    Particle::below_h = nullptr;   // tail
    Particle::left_h = nullptr;    // head
    Particle::right_h = nullptr;   // tail
}

double getRandom() { return ((double)random()) / ((double)RAND_MAX); }

template <class T> vector<vector<T>> init_matrix(int a, int b, T v) {
    vector<vector<T>> matrix;
    for (int i = 0; i < a; i++) {
        matrix.emplace_back(b, v);
    }
    return matrix;
}

int main(int, char **) {
    srand(time(nullptr));
    normal_distribution<double> distribution(
        ENVIRONMENT_SPAWN_DROP_SIZE, ENVIRONMENT_SPAWN_DROP_SIZE_STD_2 * .5);
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
        while (abs(tempSize - ENVIRONMENT_SPAWN_DROP_SIZE) >
               ENVIRONMENT_SPAWN_DROP_SIZE_STD_2) {
            tempSize = distribution(rnd_gen);
        }

        if (tempSize > tSmax)
            tSmax = tempSize;
        if (tempSize < tSmin)
            tSmin = tempSize;

        Particle *tempDrop = new Particle(tempSize / 2., tempCoord);
        Particle::dropList->push_back(tempDrop);

        dCount++;
        if (dCount % ((int)ENVIRONMENT_SPAWN_DROPS_TOTAL / 10) == 0) {
            printf("%.0f%%\r", 100. * dCount / ENVIRONMENT_SPAWN_DROPS_TOTAL);
            fflush(stdout);
        }
    }
    printf("min: %f | max: %f [µm]\n", tSmin * 1.E6, tSmax * 1.E6);
    printf("Sorting...\n");
    Particle::sortListWidth();
    Particle::sortListHeight();
    Particle::sortListSize();

    // init openCV
    int visu_width = ENVIRONMENT_WIDTH * VISU_WIDTH_PX_PER_METER;
    int visu_height = ENVIRONMENT_HEIGHT * VISU_HEIGHT_PX_PER_METER;
#ifdef USE_OPENCV
    cv::namedWindow("env_simu", cv::WINDOW_AUTOSIZE);
    cv::Mat envVisu(visu_height, visu_width, CV_8UC3);
    cv::VideoWriter video(
        OPENCV_VIDEO_FILENAME, cv::VideoWriter::fourcc(OPENCV_VIDEO_FORMAT),
        OPENCV_VIDEO_SECONDS_PER_SECOND, cv::Size(visu_width, visu_height));
    visu_width /= OPENCV_VIDEO_SCALE;
    visu_height /= OPENCV_VIDEO_SCALE;
#endif

    // simulation start
    double meanTotalMassPerPx =
        DENSITY_WATER *
        (ENVIRONMENT_SPAWN_DROP_SIZE * ENVIRONMENT_SPAWN_DROP_SIZE *
         ENVIRONMENT_SPAWN_DROP_SIZE * .125) *
        (M_PI * 4. / 3.) * ENVIRONMENT_SPAWN_DROPS_TOTAL /
        (visu_width * visu_height);
    printf("starting simulation... (meanTotalMassPerPx/biggestDropMass = %f)\n",
           meanTotalMassPerPx / Particle::bigger_h->getMass());
    long startTimeStamp = currentMicroSec();
    for (int t = 0; t <= SIMULATION_TIME_MAX; t++) { // [s]
        update();
        if (Particle::disposedDrops >= ENVIRONMENT_SPAWN_DROPS_TOTAL) {
            printf("No drops remaining. Exiting...\n");
            break;
        }

        if (t % SIMULATION_TIME_STATS_UPDATE == 0) {
            // count mass per pixel
            auto total_mass = init_matrix<double>(visu_width, visu_height, 0.);
            auto max_mass = init_matrix<double>(visu_width, visu_height, 0.);
            auto contain_drop =
                init_matrix<bool>(visu_width, visu_height, false);
            auto count_drop = init_matrix<int>(visu_width, visu_height, 0);

            double avgSize = 0.;
            Particle *currentDrop = Particle::left_h;
            int pixels_with_droplets = 0;
            int count_drop_max = 0;
            while (currentDrop != nullptr) {
                int idx[2] = {
                    (int)(currentDrop->getCoord(0) * VISU_WIDTH_PX_PER_METER /
                          OPENCV_VIDEO_SCALE),
                    (int)(currentDrop->getCoord(1) * VISU_HEIGHT_PX_PER_METER /
                          OPENCV_VIDEO_SCALE)};
                idx[1] = (idx[1] >= visu_height) ? (visu_height - 1) : idx[1];
                total_mass[idx[0]][idx[1]] += currentDrop->getMass();
                max_mass[idx[0]][idx[1]] =
                    max(max_mass[idx[0]][idx[1]], currentDrop->getMass());
                if (contain_drop[idx[0]][idx[1]] == false)
                    pixels_with_droplets++;
                contain_drop[idx[0]][idx[1]] = true;
                count_drop[idx[0]][idx[1]]++;
                count_drop_max =
                    max(count_drop_max, count_drop[idx[0]][idx[1]]);

                avgSize += currentDrop->getRadius();

                currentDrop = currentDrop->right;
            }
            avgSize *= 2. * 1.E3 / (double)Particle::remainingDrops();

            // visualization
            double avgColor = 0.;
#ifdef USE_OPENCV
            envVisu.setTo(cv::Scalar(0, 0, 0));
            for (int w = 0; w < visu_width; w++) {
                for (int h = 0; h < visu_height; h++) {
                    double totalMassRatio =
                        total_mass[w][h] /
                        meanTotalMassPerPx; // Droplet::bigger_h->getMass();
                    // totalMassRatio =
                    //    (log10(totalMassRatio + 0.1) + 1) / (log10(1.1) + 1);
                    int colorTotalMass = (int)(.68 * 255. * totalMassRatio);
                    colorTotalMass =
                        (colorTotalMass > 255) ? 255 : colorTotalMass;

                    double maxMassRatio =
                        max_mass[w][h] / Particle::bigger_h->getMass();
                    // maxMassRatio =
                    //    (log10(maxMassRatio + 0.1) + 1) / (log10(1.1) + 1);
                    int colorMaxMass = (int)(255. * maxMassRatio);
                    colorMaxMass = (colorMaxMass > 255) ? 255 : colorMaxMass;

                    int colorContainDrop = contain_drop[w][h] ? 255 : 0;

                    int colorCountDrop =
                        (int)(255. * ((double)count_drop[w][h]) /
                              (double)count_drop_max);

                    cv::Vec3b finalColor(colorCountDrop, colorMaxMass, 0);
                    for (int ww = 0; ww < OPENCV_VIDEO_SCALE; ww++)
                        for (int hh = 0; hh < OPENCV_VIDEO_SCALE; hh++)
                            envVisu.at<cv::Vec3b>(h * OPENCV_VIDEO_SCALE + hh,
                                                  w * OPENCV_VIDEO_SCALE + ww) =
                                finalColor;
                    avgColor +=
                        (finalColor[0] + finalColor[1] + finalColor[2]) / 3.;
                }
            }
            avgColor /= visu_width * visu_height;

            video.write(envVisu);
            cv::imshow("env_simu", envVisu);
            if (cv::waitKey(30) >= 0)
                break;
#endif

            long currentTimeStamp = currentMicroSec();
            long eta_seconds =
                (long)((SIMULATION_TIME_MAX - (long)t) *
                       (currentTimeStamp - startTimeStamp) / (1E6L * (long)t));
            printf(
                "time: %d [s] | drop size (avg/max): %f/%f [mm] | avgColor: %f "
                "[0-256[ | remaining drops: %d"
                " | lowest height: %f [m] (ETA: %ldm %lds)\n",
                t, avgSize, Particle::bigger_h->getRadius() * 2. * 1.E3,
                avgColor, Particle::remainingDrops(),
                (Particle::below_h == nullptr)
                    ? -1
                    : (ENVIRONMENT_HEIGHT - Particle::below_h->getCoord(1)),
                eta_seconds / 60, eta_seconds % 60);
            fflush(stdout);
        }
    }

    // clear simulation environment
#ifdef USE_OPENCV
    video.release();
#endif
    clearEnvironment();
}
