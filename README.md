# ESN-Code
ESN = "Echo State Network" neural network

cdw.cpp -
  "main" code, reads in data then trains an ESN on that data if desired and/or makes predictions based on that data

dyanmics.cpp -
  functions to help model a lattice of sites (grid of points) in a particular dynamical system

esn.cpp -
  defines the ESN as a class along with its initialization, and functions used for training and testing


simulation comparsion.png - lattice snapshots to compare target with predictions from two different models over three time steps

correlation function.png - plots of calculated "correlation function" values vs. distance between lattice sites, each plot represents a different time step, used for analysis based on theory
