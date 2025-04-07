# CFD

This project is a CFD model that solves the conservation equations of mass, momentum, and energy. It is currently a 1D model. The model is set up to handle axially symmetric geometries such as converging, diverging, and converging-diverging nozzles. The user can specify a sonic input (sonic conditions at inlet) for some geometries. The model compares the CFD result to the isentropic solution for the same geometry.

# Features
- Solves the conservation equations (mass, momentum, and energy) for compressible flow in a 1D, axially symmetric geometry
- Capable of modeling different nozzle shapes:
    - Converging geometry
    - Diverging geometry 
    - Converging-diverging geometry
- Isentropic solution and comparison to CFD
- Numerical solution uses finite difference approach

# Installation
1. Clone the respository:  
git clone https://github.com/tylerbarrows/CFD.git  
cd CFD  
3. Install required dependencies
    - Numpy
    - Math
    - Matplotlib

# Usage
Input parameters can be modified in CFD_Main.py  
Main input parameters are:
- inlet conditions (pressure, temperature)
- outlet conditions (pressure)
- initial conditions
- geometry shape
- inlet area, throat area, outlet area
- dx (incremental distance used in solution)
- length of nozzle (L)

# Running the model
To run the model, run the following from the command line in the correct directory:  
python3 CFD_Main.py
