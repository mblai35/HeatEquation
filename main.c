//
//  main.c
//  ExplicitMain.c
//
//  Created by Mallory Lai on 11/13/16.
//  Copyright © 2016 Mallory Lai. All rights reserved.
//

#include <stdio.h>
#include <math.h>

int main(int argc, const char * argv[]) {
    
    // Physical parameters
    double Length = 38.1; //length of muffin radius in mm
    int TotalTime = 100; //duration of time in minutes
    
    // Numerical parameters
    double DeltaT = 0.45; //change in time in min
    double DeltaX = 1; //change in x
    int rows = (TotalTime/DeltaT) + 1; //number of rows in grid
    int cols = (Length/DeltaX) + 1; //number of columns in grid
    
    // Counter variables
    int i;
    int j;
    
    // Create temperature grid with dimensions matching the number of columns and rows
    double TemperatureGrid[rows][cols];
    
    // Find the appropriate row increments to plug into exponential function for boundary conditions.
    double rowIncrements[rows];
    
    // Row increments
    for(i = 0; i < rows; i++){
        rowIncrements[i] = i * DeltaT;
    }
    
    // Initialize boundary conditions of grid
    // Left boundary; Populate first column of TemperatureGrid using exponential line of best fit.
    for(i = 0; i < rows; i++){
        TemperatureGrid[i][0] = 81.09*exp(-.09*rowIncrements[i]) + 92.93*exp(-.00217*rowIncrements[i]);
    }
    
    // Right boundary; Populate last column of TemperatureGrid using exponential line of best fit.
    for(i = 0; i < rows; i++){
        TemperatureGrid[i][cols-1] = 80.35*exp(-.1156*rowIncrements[i]) + 93.69*exp(-.002442*rowIncrements[i]);
    }
    
    // Top boundary; Linear interpolation of first row.
    for(i = 0; i < cols; i++){
        TemperatureGrid[0][i+1] = TemperatureGrid[0][i] + (TemperatureGrid[0][cols-1]-TemperatureGrid[0][0])/(cols-1);
    }
    
    
    // Calculation; explicit solution.
    for(i = 1; i < rows; i++){
        for(j = 1; j < (cols - 1); j++){
            TemperatureGrid[i][j] = TemperatureGrid[i-1][j] + DeltaT*((TemperatureGrid[i-1][j-1] - 2*TemperatureGrid[i-1][j] +
                                                                       TemperatureGrid[i-1][j+1])/DeltaX*DeltaX);
        }
    }
    
    
    // Write TemperatureGrid to csv file.
    FILE * grid=fopen("Grid.csv", "w+");
    
    for(i = 0; i < rows; i++){
        for(j = 0; j < cols; j++){
            if(j == (cols-1)){
                fprintf(grid, " %f\n", TemperatureGrid[i][j]);
            }
            else{
                fprintf(grid, " %f,", TemperatureGrid[i][j]);
            }
        }
    }
    

    
    fclose(grid);
    
    return 0;
}