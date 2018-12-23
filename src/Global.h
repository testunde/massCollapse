/*
 * Global.h
 *
 *  Created on: Dec 16, 2018
 *      Author: root
 */

#define DENSITY_WATER 997. // [kg/m^3]
#define DROP_GROWTH_CONDENSATION_S 1.001 // [1]
#define DROP_GROWTH_CONDENSATION_F 1.E+10 // [s/m^2]
#define DROP_GROWTH_CONDENSATION_LIMIT_SIZE 0.31623E-6 // [m]

#define ENVIRONMENT_WIDTH 0.1 // [m]
#define ENVIRONMENT_HEIGHT 3. // [m]
#define ENVIRONMENT_SPAWN_THICKNESS 0.3 // [m]
#define ENVIRONMENT_SPAWN_DROPS_TOTAL 1E4
#define ENVIRONMENT_SPAWN_DROP_SIZE 2.E-6 // [m] (diameter)
#define ENVIRONMENT_SPAWN_DROP_SIZE_STD_2 1.8E-6 // [m] (2*std; diameter)

#define SIMULATION_TIME_MAX 1200 // [s]
#define SIMULATION_TIME_STATS_UPDATE 1 // [s]

#define VISU_WIDTH_PX_PER_METER 1000 // [px/m]
#define VISU_HEIGHT_PX_PER_METER 300 // [px/m]
