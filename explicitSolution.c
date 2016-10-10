//
//  explicitSolution.c
//  HeatEquation
//
//  Created by Mallory Lai on 10/7/16.
//  Copyright © 2016 Mallory Lai. All rights reserved.
//

#include <stdio.h>
#include <math.h>

int main(int argc, const char * argv[]) {
    
    // Physical parameters
    double Length = 38.1; //length of muffin radius in mm
    int TotalTime = 100; //duration of time in minutes
    
    // Numerical parameters
    double DeltaT = 0.45; //change in time
    double DeltaX = 1; //change in x
    int rows = (TotalTime/DeltaT) + 1;
    int cols = (Length/DeltaX) + 1;
    
    // Counter variables
    int i;
    int j;
    
    // Create grid
    double TemperatureGrid[rows][cols];
    // Initialize boundary conditions of grid
    
    // Left boundary
    
    // Right boundary
    
    // Top boundary
    
    
    // Calculation
    for(i = 2; i <= rows; i++){
        for(j = 2; j < cols; j++){
    TemperatureGrid[i,j] <- TemperatureGrid[i-1, j] + DeltaT*((TemperatureGrid[i-1, j-1] - 2*TemperatureGrid[i-1,j] +
         TemperatureGrid[i-1, j+1])/DeltaX*DeltaX)
        }
    }
    
    printf("Rows = %d\nCols = %d\n", rows, cols);
    
    //free(rows)
    
    return 0;
}
