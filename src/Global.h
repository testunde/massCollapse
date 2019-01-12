/*
 * Global.h
 *
 *  Created on: Jan 12, 2019
 */

#define GRAVITAIONAL_CONSTANT 6.67408E-11 // [m^3 / (kg * s^2)]

#define ENVIRONMENT_WIDTH 1                         // [m]
#define ENVIRONMENT_HEIGHT 1                        // [m]
#define ENVIRONMENT_SPAWN_PARTICLES_TOTAL 2         //.5E4
#define ENVIRONMENT_SPAWN_PARTICLE_MASS 5.E2        // [kg]
#define ENVIRONMENT_SPAWN_PARTICLE_MASS_STD_2 1.8E1 // [kg] (2*std)
#define ENVIRONMENT_SPAWN_START_VELO                                           \
    -2.E-3 // [m/2] (orthogonal towards the center)

#define SIMULATION_TIME_MAX 3600       // [s]
#define SIMULATION_TIME_STATS_UPDATE 1 // [s]

#define VISU_WIDTH_PX_PER_METER 800  // [px/m]
#define VISU_HEIGHT_PX_PER_METER 800 // [px/m]

#ifdef USE_OPENCV
#define OPENCV_VIDEO_FILENAME "video.avi"
#define OPENCV_VIDEO_FORMAT 'H', 'F', 'Y', 'U' // raw; as for .avi format
#define OPENCV_VIDEO_SECONDS_PER_SECOND                                        \
    60                       // [s/s] simulation seconds per video second (FPS)
#define OPENCV_VIDEO_SCALE 1 // only positive int!
#else
#define OPENCV_VIDEO_SCALE 1 // don't change!
#endif
