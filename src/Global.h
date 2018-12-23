/*
 * Global.h
 *
 *  Created on: Dec 16, 2018
 *      Author: root
 */

#define DENSITY_WATER 997. // [kg/m^3]
#define DROP_GROWTH_CONDENSATION_S 1.001 // [1]
#define DROP_GROWTH_CONDENSATION_F 1.E+10 // [s/m^2]

#define ENVIRONMENT_CELL_SIZE 1.E-3 // [m]
#define ENVIRONMENT_CELLS_WIDTH 100 // 1 cell = 1 [mm]
#define ENVIRONMENT_CELLS_HEIGHT 1000 // 1 cell = 1 [mm]
#define ENVIRONMENT_CELLS_SPAWN_THICKNESS 300 // 1 cell = 1 [mm]
#define ENVIRONMENT_SPAWN_DROPS_PER_CELL 10
#define ENVIRONMENT_SPAWN_DROP_SIZE 2.E-6 // [m] (diameter)
#define ENVIRONMENT_SPAWN_DROP_SIZE_STD_2 1.8E-6 // [m] (2*std; diameter)

#define SIMULATION_TIME_MAX 1200 // [s]
#define SIMULATION_TIME_STATS_UPDATE 1 // [s]
