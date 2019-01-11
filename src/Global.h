/*
 * Global.h
 *
 *  Created on: Dec 16, 2018
 */

#define DENSITY_WATER 997.                             // [kg/m^3]
#define DROP_GROWTH_CONDENSATION_S 1.001               // [1]
#define DROP_GROWTH_CONDENSATION_F 1.E+10              // [s/m^2]
#define DROP_GROWTH_CONDENSATION_LIMIT_SIZE 0.31623E-6 // [m]

#define ENVIRONMENT_WIDTH .05           // [m]
#define ENVIRONMENT_HEIGHT 2000         // [m]
#define ENVIRONMENT_SPAWN_THICKNESS 100 // [m]
#define ENVIRONMENT_SPAWN_DROPS_TOTAL                                          \
    (ENVIRONMENT_WIDTH * ENVIRONMENT_SPAWN_THICKNESS *                         \
     2.1E4) // hint: accretion (merging) occur beginning at 1E4 drop/m^2
#define ENVIRONMENT_SPAWN_DROP_SIZE 2.E-6        // [m] (diameter)
#define ENVIRONMENT_SPAWN_DROP_SIZE_STD_2 1.8E-6 // [m] (2*std; diameter)

#define SIMULATION_TIME_MAX 3600       //(5 * 3600) // [s]
#define SIMULATION_TIME_STATS_UPDATE 1 // [s]

#define VISU_WIDTH_PX_PER_METER 5000         // [px/m]
#define VISU_HEIGHT_PX_PER_METER (2.5 * .18) // [px/m]

#ifdef USE_OPENCV
#define OPENCV_VIDEO_FILENAME "video.avi"
#define OPENCV_VIDEO_FORMAT 'H', 'F', 'Y', 'U' // raw; as for .avi format
#define OPENCV_VIDEO_SECONDS_PER_SECOND                                        \
    60                       // [s/s] simulation seconds per video second (FPS)
#define OPENCV_VIDEO_SCALE 1 // only positive int!
#else
#define OPENCV_VIDEO_SCALE 1 // don't change!
#endif
