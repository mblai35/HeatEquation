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
    
    // Create grid
    double Temperature[rows][cols];
    // Initialize boundary conditions of grid
    // Left boundary
    // Bottom boundary
    
    
    
    printf("Rows = %d\nCols = %d\n", rows, cols);
    
    //free(rows)
    
    return 0;
}
