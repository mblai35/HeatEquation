# explicitSoln.R
# R version 3.2.2 (2015-08-14)
# October 8, 2016. Mallory B. Lai.
# Reviewed by: TODO (Mallory B. Lai) : Find reviewer to proofread
# Explicit Solution to 1D heat equation for English muffin.

# Collaborators: Xiukun Hu, Geeta Monpara 

#------------------------------------------------------------------------------

# Define duration of temperature recordings in minutes:
L <- 100
# Define English muffin radius in mm: 
x <- 38.1
# Define delta t: 
dt <- .45
# Define delta x:
dx <- 1
# Define initial external boundary temperature.
initExt <- 156

# Initialize first row with initial boundary conditions. 
# Find increments of x.
Xincrements <- seq(0, x, dx)
# Create even sequence from left to right boundary with the appropriate length.
row1 <- seq(col1[1], initExt, length.out = length(Xincrements))

# Initialize first column with initial boundary conditions. 
# Find increments of L. 
Lincrements <- seq(0, L, dt)
# Use the exponential function found by Geeta to populate field with correct 
# boundary temperature. 
col1 <- 81.09*exp(-.09*Lincrements) + 92.93*exp(-.00217*Lincrements)

# Initialize last column with initial boundary conditions. 
# Use the exponential function found by Geeta to populate field with correct 
# boundary temperature. 
lastCol <- 80.35*exp(-.1156*Lincrements) + 93.69*exp(-.002442*Lincrements)

# Create a datatable with zeros.
TemperatureGrid <- matrix(0, nrow = length(Lincrements), 
                          ncol = length(Xincrements))
# Update matrix to include initial boundary conditions
TemperatureGrid[1, ] <- row1
TemperatureGrid[, 1] <- col1
TemperatureGrid[, length(Xincrements)] <- lastCol

# Update values of grid using the explicit solution.
for (i in 2:(dim(TemperatureGrid)[1])){
  for(j in 2:(dim(TemperatureGrid)[2]-1)){
    TemperatureGrid[i,j] <- TemperatureGrid[i-1, j] + 
      dt*((TemperatureGrid[i-1, j-1] - 2*TemperatureGrid[i-1,j] +
        TemperatureGrid[i-1, j+1])/dx^2)
  }
}
  


