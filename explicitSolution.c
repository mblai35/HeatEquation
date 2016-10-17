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
    //firstCol = 81.09*exp(-.09*Lincrements) + 92.93*exp(-.00217*Lincrements)
    // Right boundary
    //lastCol <- 80.35*exp(-.1156*Lincrements) + 93.69*exp(-.002442*Lincrements)
    // Top boundary
    // ???
    
    /*
    // Calculation
    for(i = 2; i <= rows; i++){
        for(j = 2; j < cols; j++){
    TemperatureGrid[i,j] <- TemperatureGrid[i-1, j] + DeltaT*((TemperatureGrid[i-1, j-1] - 2*TemperatureGrid[i-1,j] +
         TemperatureGrid[i-1, j+1])/DeltaX*DeltaX)
        }
    }
    */
    
    /*
    for(i = 1; i <= rows; i++){
        for(j = 1; j <= cols; j++){
            if(j = rows){
                printf(" %d\n", TemperatureGrid[i,j]);
            }
            else{
                printf(" %d ", TemperatureGrid[i,j]);
            }
        }
     }
     */
    
    printf("Rows = %d\nCols = %d\n", rows, cols);
    
    //free(rows)
    
    return 0;
}
