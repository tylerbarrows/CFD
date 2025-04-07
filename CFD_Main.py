## Main CFD Script
# Notes
# maybe put all of this below in an inputs script??

## ---------------------------------------------------------- ##
## Import
import numpy as np
import math 
from CFD_MeshGeneration import meshGenerationFunction
from CFD_SS_Solution import SSSolutionFunction
from CFD_Plot import plotSSResultFunction
from CFD_IsentropicSolution import isentropicSolutionFunc
from CFD_PreExamine import PreExamineFunction
## --  --  --  --  --  --  --  --  --  --  --  --  --  --  -- ##

## ---------------------------------------------------------- ##
## Inputs
# Geometry, fluid properties, intial and boundary conditions, fluid inputs, solver information, plotting requests
class geometryInputs:
    # select geometry shape: straight, converging, diverging, converging-diverging
    # diverging note: have to be careful with intial conditions
    shape = 'converging'
    shape = 'diverging'
    # shape = 'converging-diverging'
    # select 1D or 2D solution
    dimension = '1D'
    # geometries are specified in meters (m)
    # length of segment
    length = 0.15
    # inlet width
    heightInlet = 0.1
    # heightInlet = 0.1 # use for diverging
    # throat width
    heightThroat = 0.1
    # outlet width
    heightOutlet = 0.15
    # dx
    dx = 0.001
    # dx = 0.01
    # dy 
    dy = 0.0005 

# Fluid properties
class fluidProperties:
    # Gas constant, units J/kg/K
    R = 287
    # Specific heat ratio
    k = 1.4
    # Specific heat (constant pressure and constant volume, J/kg/K)
    cp = k/(k-1)*R
    cv = cp/k
    # thermal conductivity (W/m/K)
    kcond = 0.026
    # kcond = 0

# Initial Conditions
class initialConditions:
    class inlet:
        # P - bar, T - deg C, u - m/s
        P = 1.16
        T = 96.5
        u = 100
    class outlet:
        # P - bar, T - deg C, u - m/s
        P = 1.19
        # P = 1.22 # use for diverging
        T = 100
        u = 100
    class core:
        # P - bar, T - deg C, u - m/s
        P = 1.16
        T = 100
        u = 100
    # set sonicInputPin to 1 for pinning input to be sonic (will be checked again in pre-examine)
    sonicInputPin = 0
    fixedInputTestPin = 0 # only use if want to do something that is not really valie

# Boundary conditions
class boundaryConditions:
    # eventually put something here for heat transfer, etc.
    # heat transfer boundary condition
    heatBoundary = 'isothermal'
    Twall = 500
    # no slip (u=0 at boundary) (or perpendicular velocity = 0)
    noSlip = 'on'
    # friction
    includeFriction = 'no'
    # friction related variables
    fd = 0
    # heat transfer
    class q:
        # heat transfer shape: constant, linear, exponential, etc.
        shape = 'constant'
        qMax = 0#10
        qMin = 0#10
        # y = q*e^(ax)
        a = 0
    class intletConditions: 
        # Need to specify two inlet conditions for subsonic (3 for sonic)
        one = 'pressure'
        two = 'temperature'

# Solver Information
# steady state, step, ramp, etc.
class solverInformation:
    # run type - steady state (SS), step, ramp, etc.
    runType = 'SS'
    # solver type - finite element (FE), finite volume (FV)
    solverType = 'FE'
    # input type - mass flow step, pressure step, density step, etc.
    inputType = 'mass flow step'
    # acceptable error (%) - not really error but used in covergence
    error_step = 1/100000000
    error_SS = 1/10000000
    error_Isen = 1/10000
    error_pre = 1/1000
    # weights for semi-implicit solver
    X,Y = 1,99 # Y > X implies more implicit solver # problem with X and Y weights?
    # factor for finding dt (C should be greater than 1)
    C = 1.2
    C = 1.8
    # number of digits in dt
    numDigits = 10
    # max iterations
    class iterationInformation:
        # max number of iterations for finding next step
        nextStepMax = 30
        # number of time steps for steady state
        SSMax = 10000#10000
        # SSMax = 50000
        # number of steps for isentropic convergence
        IsenMax = 100
        # max number of iterations for pre_exam
        preExamMax = 100

class checkAgainstIsentropicSolution:
    # either yes or no
    check = 'yes'
    # include comparison in plots (yes/no)
    includeInPlots = 'yes'

# Plotting requests
class plottingRequests:
    # number of plots
    plotNumber = 1
    # plot type
    plotOne = 'P_vs_x'
## --  --  --  --  --  --  --  --  --  --  --  --  --  --  -- ##

## ---------------------------------------------------------- ##
## Mesh setup
meshGeneration = meshGenerationFunction(geometryInputs)
# print(meshGeneration.h_vector)

## --  --  --  --  --  --  --  --  --  --  --  --  --  --  -- ##

## ---------------------------------------------------------- ##
## Pre-solution exmaination
PreExamResults = PreExamineFunction(geometryInputs,fluidProperties,initialConditions,boundaryConditions,solverInformation)

## --  --  --  --  --  --  --  --  --  --  --  --  --  --  -- ##

## ---------------------------------------------------------- ##
## Valid inputs
flagValid = PreExamResults.flagValid

## --  --  --  --  --  --  --  --  --  --  --  --  --  --  -- ##

match flagValid:
    case 'valid':

        ## ---------------------------------------------------------- ##
        ## Solver setup
        print('SonicInput: ',PreExamResults.sonicInput,'0-off,1-on')
        SS_Solution = SSSolutionFunction(meshGeneration,geometryInputs,fluidProperties,initialConditions,boundaryConditions,solverInformation,PreExamResults)
        # print(SS_Solution.int_track[-1][0])
        # print(SS_Solution.P_track)
        print('Numerical iterations: ',SS_Solution.SSStep)

        ## --  --  --  --  --  --  --  --  --  --  --  --  --  --  -- ##

        ## ---------------------------------------------------------- ##
        ## Isentropic solution
        P_isen, T_isen, rho_isen, u_isen, M_isen, flag_isen = isentropicSolutionFunc(meshGeneration,geometryInputs,fluidProperties,initialConditions,boundaryConditions,solverInformation,SS_Solution,PreExamResults)
        IsentropicSolution = np.vstack((P_isen,T_isen))
        IsentropicSolution = np.vstack((IsentropicSolution,rho_isen))
        IsentropicSolution = np.vstack((IsentropicSolution,u_isen))
        IsentropicSolution = np.vstack((IsentropicSolution,M_isen))

        ## --  --  --  --  --  --  --  --  --  --  --  --  --  --  -- ##

        ## ---------------------------------------------------------- ##
        ## Plotting SS results
        plotSSResultFunction(SS_Solution,plottingRequests,IsentropicSolution,boundaryConditions)

        ## --  --  --  --  --  --  --  --  --  --  --  --  --  --  -- ##
    case 'not valid':
        print('Inputs are not valid for finding solution')

    case _:
        print('Error')



# use this template below for more functions
## ---------------------------------------------------------- ##
## 

## --  --  --  --  --  --  --  --  --  --  --  --  --  --  -- ##