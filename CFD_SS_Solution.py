def SSSolutionFunction(meshGeneration,geometryInputs,fluidProperties,initialConditions,boundaryConditions,solverInformation,PreExamResults):
    # first check if 1D or 2D solution (only including 1D for now)
    # first build P,rho,T vectors that will turn into matrices
    # need time vector as well
    # find dt based on dx and c and v (estimate v)
    # right now, P is specified at end (unless supersonic), everything else floats
    # if supersonic, P needs to float as well
    # P and rho specified at inlet, let u float

    # Import
    import math
    import numpy as np
    from CFD_NextStepResult import nextStepFunction_1D

    class SS_Solution_result:
        ## ---------------------------------------------------------------- ##
        ## Import 
        # Geometry
        x_vector = meshGeneration.x_vector
        h_vector = meshGeneration.h_vector
        dx = meshGeneration.dx
        dimension = geometryInputs.dimension

        # Fluid Properties
        k = fluidProperties.k
        R = fluidProperties.R
        cp = fluidProperties.cp
        cv = fluidProperties.cv

        # Solver related inputs
        error_SS = solverInformation.error_SS
        C = solverInformation.C
        numDigits = solverInformation.numDigits
        maxSSSteps = solverInformation.iterationInformation.SSMax
        # PreExamResults = PreExamResults

        ## ---------------------------------------------------------------- ##
        # P,T,u vectors
        P_vector = np.linspace(initialConditions.core.P,initialConditions.core.P,len(x_vector))
        P_vector[0] = initialConditions.inlet.P
        P_vector[-1] = initialConditions.outlet.P

        T_vector = np.linspace(initialConditions.core.T,initialConditions.core.T,len(x_vector))
        T_vector[0] = initialConditions.inlet.T
        T_vector[-1] = initialConditions.outlet.T

        u_vector = np.linspace(initialConditions.core.u,initialConditions.core.u,len(x_vector))
        u_vector[0] = initialConditions.inlet.u
        u_vector[-1] = initialConditions.outlet.u

        # Converting units
        P_vector = P_vector * 100000
        T_vector = T_vector + 273.15

        # rho vector
        rho_vector = P_vector/(R*T_vector)

        # tracking vectors
        P_track = P_vector
        T_track = T_vector
        u_track = u_vector
        rho_track = rho_vector
        int_track = np.array([0])
        SSStep_track = np.array([0])
        time_track = np.array([0])

        # finding dt step
        c = math.sqrt(k*R*max(T_vector))
        M = max(np.absolute(u_vector)/(np.sqrt(k*R*T_vector)))
        dt = round(dx/(C*(M+1)*c),numDigits) # need this to update if solution becomes supersonic

        # friction, heat transfer boundary condition, velocity boundary condition
        # nothing here right now

        ## ---------------------------------------------------------------- ##
        ## Solving
        match dimension:
            case '1D':
                # First step
                rho_vector, P_vector, T_vector, u_vector, iterations, dt = nextStepFunction_1D(x_vector,h_vector,P_vector,T_vector,rho_vector,u_vector,dx,dt,fluidProperties,solverInformation,initialConditions,boundaryConditions,PreExamResults)
                # stacking arrays
                rho_track = np.vstack((rho_track,rho_vector))
                P_track = np.vstack((P_track,P_vector))
                T_track = np.vstack((T_track,T_vector))
                u_track = np.vstack((u_track,u_vector))
                int_track = np.vstack((int_track,iterations))
                SSStep_track = np.vstack((SSStep_track,np.array([0])))
                time_track = np.vstack((time_track,time_track[-1]+dt))
                # While loop
                SSStep = 2
                # while (max(abs((P_track[-1,:]-P_track[-2,:])/P_track[-2,:])) > error_SS or max(abs((T_track[-1,:]-T_track[-2,:])/T_track[-2,:])) > error_SS or max(abs((rho_track[-1,:]-rho_track[-2,:])/rho_track[-2,:])) > error_SS) and SSStep < maxSSSteps:
                for SSStep in range(3,maxSSSteps):
                    # SSStep = SSStep + 1
                    c = math.sqrt(k*R*max(T_vector))
                    M = max(np.absolute(u_vector)/(np.sqrt(k*R*T_vector)))
                    dt = round(dx/(C*(M+1)*c),numDigits) # need this to update if solution becomes supersonic
                    dt = dx/(C*(M+1)*c)
                    rho_vector, P_vector, T_vector, u_vector, iterations, dt = nextStepFunction_1D(x_vector,h_vector,P_vector,T_vector,rho_vector,u_vector,dx,dt,fluidProperties,solverInformation,initialConditions,boundaryConditions,PreExamResults)

                    rho_track = np.vstack((rho_track,rho_vector))
                    P_track = np.vstack((P_track,P_vector))
                    T_track = np.vstack((T_track,T_vector))
                    u_track = np.vstack((u_track,u_vector))
                    int_track = np.vstack((int_track,iterations))
                    SSStep_track = np.vstack((SSStep_track,SSStep))
                    time_track = np.vstack((time_track,time_track[-1][0]+dt))

                    if min(rho_vector) < 0 or min(P_vector) < 0 or min(T_vector) < 0:
                        print('rho=',min(rho_vector),'at x=',x_vector[np.where(rho_vector==min(rho_vector))[0][0]])
                        print('P=',min(P_vector),'at x=',x_vector[np.where(P_vector==min(P_vector))[0][0]])
                        print('T=',min(T_vector),'at x=',x_vector[np.where(T_vector==min(T_vector))[0][0]])
                        print('u=',min(u_vector),'at x=',x_vector[np.where(u_vector==min(u_vector))[0][0]])
                        print("breaking...")
                        break
                    elif (max(abs((P_track[-1,:]-P_track[-2,:])/P_track[-2,:])) < error_SS and max(abs((T_track[-1,:]-T_track[-2,:])/T_track[-2,:])) < error_SS and max(abs((rho_track[-1,:]-rho_track[-2,:])/rho_track[-2,:])) < error_SS):
                        break
            case '2D':
                flag = 1
            case _:
                flag = 1

    return SS_Solution_result