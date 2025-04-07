def nextStepFunction_1D(x_vector,h_vector,P_vector,T_vector,rho_vector,u_vector,dx,dt,fluidProperties,solverInformation,initialConditions,boundaryConditions,PreExamResults):

    # Import
    import math 
    import numpy as np
    from numpy import inf

    # Initial conditions
    # only need P and u for now
    P_inlet = initialConditions.inlet.P * 100000
    P_outlet = initialConditions.outlet.P * 100000
    u_inlet = initialConditions.inlet.u
    T_inlet = initialConditions.inlet.T + 273.15

    # Fluid Properties
    k = fluidProperties.k
    R = fluidProperties.R
    cp = fluidProperties.cp
    cv = fluidProperties.cv
    kcond = fluidProperties.kcond
    fd = boundaryConditions.fd

    # Heat input
    qMin = boundaryConditions.q.qMin
    qMax = boundaryConditions.q.qMax
    if qMin > 0 and qMax > 0:
        a = 1/x_vector[-1]*math.log(qMax/qMin)
    else:
        a = 0

    q_array = np.linspace(1,1,len(x_vector))
    temp_counter = 0
    for i in q_array:
        q_i = qMin * math.exp(a*x_vector[temp_counter])
        q_array[temp_counter] = q_i
        temp_counter = temp_counter + 1

    # Max Iterations
    maxIterationsStep = solverInformation.iterationInformation.nextStepMax

    # Step Convergence
    convergenceError = solverInformation.error_step

    # Weights
    X = solverInformation.X
    Y = solverInformation.Y

    # Other solver information
    sonicInput = PreExamResults.sonicInput
    c = math.sqrt(k*R*T_inlet)
    fixedInputTestPin = initialConditions.fixedInputTestPin
    if fixedInputTestPin == 1:
        c = u_inlet
        sonicInput = 1
    C = 0.9 # used to help make sure variables don't become negative
    Cx = 0.1 # only using Cx in 1D case
    # Cx = 0.3
    Cy = 0.1

    # First Step Estimate

    # First finding core vectors
    P_core = P_vector[1:-1]
    T_core = T_vector[1:-1]
    rho_core = rho_vector[1:-1]
    u_core = u_vector[1:-1]
    M_inlet = u_vector[1]/math.sqrt(k*R*T_vector[1])
    M_outlet = u_vector[-1]/math.sqrt(k*R*T_vector[-1])

    # Used in calculations
    lenV = len(P_vector)

    # First time derivatives
    drho_dt = -u_core * (rho_vector[2:lenV]-rho_vector[0:-2]) / (2*dx) \
        - rho_core * (u_vector[2:lenV]-u_vector[0:-2]) / (2*dx) \
        - rho_core * u_core * (np.log(h_vector[2:lenV])-np.log(h_vector[0:-2])) / (2*dx)
    du_dt = -u_core * (u_vector[2:lenV]-u_vector[0:-2]) / (2*dx) \
        - 1/rho_core * (P_vector[2:lenV]-P_vector[0:-2]) / (2*dx) \
        - u_core/rho_core * fd * 1 / h_vector[1:-1]
    dT_dt = -u_core * (T_vector[2:lenV]-T_vector[0:-2]) / (2*dx) \
        - 1/cv*P_core/rho_core * (u_vector[2:lenV]-u_vector[0:-2]) / (2*dx) \
        - 1/cv*P_core*u_core/rho_core * (np.log(h_vector[2:lenV])-np.log(h_vector[0:-2])) / (2*dx) \
        + 1/cv/rho_core * kcond * (T_vector[2:lenV]-2*T_core+T_vector[0:-2])/(math.pow(dx,2)) \
        + 1/cv/rho_core*kcond*(T_vector[2:lenV]-T_vector[0:-2]) / (2*dx) * (np.log(h_vector[2:lenV])-np.log(h_vector[0:-2])) / (2*dx) \
        + 1/cv/rho_core*q_array[1:-1]

    # Artificial viscosity
    S_rho = np.abs(P_vector[2:lenV]-2*P_core+P_vector[0:-2])/(P_vector[2:lenV]+2*P_core+P_vector[0:-2]) \
        * (rho_vector[2:lenV]-2*rho_core+rho_vector[0:-2])
    S_u = np.abs(P_vector[2:lenV]-2*P_core+P_vector[0:-2])/(P_vector[2:lenV]+2*P_core+P_vector[0:-2]) \
        * (u_vector[2:lenV]-2*u_core+u_vector[0:-2])
    S_T = np.abs(P_vector[2:lenV]-2*P_core+P_vector[0:-2])/(P_vector[2:lenV]+2*P_core+P_vector[0:-2]) \
        * (T_vector[2:lenV]-2*T_core+T_vector[0:-2])
    # S_rho = np.nan_to_num(S_rho,nan=0.0)
    # S_rho[S_rho==inf]=0
    # S_u = np.nan_to_num(S_u,nan=0.0)
    # S_u[S_u==inf]=0
    # S_T = np.nan_to_num(S_T,nan=0.0)
    # S_T[S_T==inf]=0
    
    # First estimate at new step
    rho_core_next = rho_core + drho_dt * dt + S_rho
    u_core_next = u_core + du_dt * dt + S_u
    T_core_next = T_core + dT_dt * dt + S_T
    P_core_next = rho_core_next * R * T_core_next

    # P at inlet is specified, need to check about outlet
    # Specify P at inlet, let density float, find T using ideal gas law
    # Check if P floats at outlet or not, let density float, find T using ideal gas law
    if M_outlet < 1.0:
        P_vector_next = np.append(np.append(np.array([P_inlet]),P_core_next),np.array([P_outlet]))
    else:
        P_vector_next = np.append(np.append(np.array([P_inlet]),P_core_next),np.array([2*P_core_next[-1]-P_core_next[-2]]))
    if sonicInput == 1:
        u_vector_next = np.append(np.append(np.array([c]),u_core_next),np.array([2*u_core_next[-1]-u_core_next[-2]]))
    else:
        u_vector_next = np.append(np.append(np.array([2*u_core_next[0]-u_core_next[1]]),u_core_next),np.array([2*u_core_next[-1]-u_core_next[-2]]))
    T_vector_next = np.append(np.append(np.array([T_inlet]),T_core_next),np.array([2*T_core_next[-1]-T_core_next[-2]]))
    rho_vector_next = np.append(np.append(np.array([P_vector_next[0]/R/T_vector_next[0]]),rho_core_next),np.array([P_vector_next[-1]/R/T_vector_next[-1]]))
    
    # new time derivative with weights
    # X weights
    # drho_dt_X = -u_core * (rho_vector[2:lenV]-rho_vector[0:-2]) / (2*dx) \
    #     - rho_core * (u_vector[2:lenV]-u_vector[0:-2]) / (2*dx) \
    #     - rho_core * u_core * (np.log(h_vector[2:lenV])-np.log(h_vector[0:-2])) / (2*dx)
    # du_dt_X = -u_core * (u_vector[2:lenV]-u_vector[0:-2]) / (2*dx) \
    #     - 1/rho_core * (P_vector[2:lenV]-P_vector[0:-2]) / (2*dx) \
    #     - u_core/rho_core * fd * 1 / h_vector[1:-1]
    # dT_dt_X = -u_core * (T_vector[2:lenV]-T_vector[0:-2]) / (2*dx) \
    #     - 1/cv*P_core/rho_core * (u_vector[2:lenV]-u_vector[0:-2]) / (2*dx) \
    #     - 1/cv*P_core*u_core/rho_core * (np.log(h_vector[2:lenV])-np.log(h_vector[0:-2])) / (2*dx) \
    #     + 1/cv/rho_core * kcond * (T_vector[2:lenV]-2*T_core+T_vector[0:-2])/(math.pow(dx,2)) \
    #     + 1/cv/rho_core*kcond*(T_vector[2:lenV]-T_vector[0:-2]) / (2*dx) * (np.log(h_vector[2:lenV])-np.log(h_vector[0:-2])) / (2*dx) 
    drho_dt_X = drho_dt
    du_dt_X = du_dt
    dT_dt_X = dT_dt
    # Artificial viscosity
    # S_rho_X = np.abs(P_vector[2:lenV]-2*P_core+P_vector[0:-2])/(P_vector[2:lenV]+2*P_core+P_vector[0:-2]) \
    #     * (rho_vector[2:lenV]-2*rho_core+rho_vector[0:-2])
    # S_u_X = np.abs(P_vector[2:lenV]-2*P_core+P_vector[0:-2])/(P_vector[2:lenV]+2*P_core+P_vector[0:-2]) \
    #     * (u_vector[2:lenV]-2*u_core+u_vector[0:-2])
    # S_T_X = np.abs(P_vector[2:lenV]-2*P_core+P_vector[0:-2])/(P_vector[2:lenV]+2*P_core+P_vector[0:-2]) \
    #     * (T_vector[2:lenV]-2*T_core+T_vector[0:-2])
    S_rho_X = S_rho
    S_u_X = S_u
    S_T_X = S_T
    # S_rho_X = np.nan_to_num(S_rho_X,nan=0.0)
    # S_rho_X[S_rho_X==inf]=0
    # S_u_X = np.nan_to_num(S_u_X,nan=0.0)
    # S_u_X[S_u_X==inf]=0
    # S_T_X = np.nan_to_num(S_T_X,nan=0.0)
    # S_T_X[S_T_X==inf]=0
    # Y weights
    drho_dt_Y = -u_core_next * (rho_vector_next[2:lenV]-rho_vector_next[0:-2]) / (2*dx) \
        - rho_core_next * (u_vector_next[2:lenV]-u_vector_next[0:-2]) / (2*dx) \
        - rho_core_next * u_core_next * (np.log(h_vector[2:lenV])-np.log(h_vector[0:-2])) / (2*dx)
    du_dt_Y = -u_core_next * (u_vector_next[2:lenV]-u_vector_next[0:-2]) / (2*dx) \
        - 1/rho_core_next * (P_vector_next[2:lenV]-P_vector_next[0:-2]) / (2*dx) \
        - u_core_next/rho_core_next * fd * 1 / h_vector[1:-1]
    dT_dt_Y = -u_core_next * (T_vector_next[2:lenV]-T_vector_next[0:-2]) / (2*dx) \
        - 1/cv*P_core_next/rho_core_next * (u_vector_next[2:lenV]-u_vector_next[0:-2]) / (2*dx) \
        - 1/cv*P_core_next*u_core_next/rho_core_next * (np.log(h_vector[2:lenV])-np.log(h_vector[0:-2])) / (2*dx) \
        + 1/cv/rho_core_next * kcond * (T_vector_next[2:lenV]-2*T_core_next+T_vector_next[0:-2])/(math.pow(dx,2)) \
        + 1/cv/rho_core_next*kcond*(T_vector_next[2:lenV]-T_vector_next[0:-2]) / (2*dx) * (np.log(h_vector[2:lenV])-np.log(h_vector[0:-2])) / (2*dx) \
        + 1/cv/rho_core_next*q_array[1:-1]
    # Artificial viscosity
    S_rho_Y = np.abs(P_vector_next[2:lenV]-2*P_core_next+P_vector_next[0:-2])/(P_vector_next[2:lenV]+2*P_core_next+P_vector_next[0:-2]) \
        * (rho_vector_next[2:lenV]-2*rho_core_next+rho_vector_next[0:-2])
    S_u_Y = np.abs(P_vector_next[2:lenV]-2*P_core_next+P_vector_next[0:-2])/(P_vector_next[2:lenV]+2*P_core_next+P_vector_next[0:-2]) \
        * (u_vector_next[2:lenV]-2*u_core_next+u_vector_next[0:-2])
    S_T_Y = np.abs(P_vector_next[2:lenV]-2*P_core_next+P_vector_next[0:-2])/(P_vector_next[2:lenV]+2*P_core_next+P_vector_next[0:-2]) \
        * (T_vector_next[2:lenV]-2*T_core_next+T_vector_next[0:-2])
    # S_rho_Y = np.nan_to_num(S_rho_Y,nan=0.0)
    # S_rho_Y[S_rho_Y==inf]=0
    # S_u_Y = np.nan_to_num(S_u_Y,nan=0.0)
    # S_u_Y[S_u_Y==inf]=0
    # S_T_Y = np.nan_to_num(S_T_Y,nan=0.0)
    # S_T_Y[S_T_Y==inf]=0
    # Weighted derivatives
    drho_dt = (X*drho_dt_X+Y*drho_dt_Y)/(X+Y)
    du_dt = (X*du_dt_X+Y*du_dt_Y)/(X+Y)
    dT_dt = (X*dT_dt_X+Y*dT_dt_Y)/(X+Y)

    # Estimate at next step
    rho_core_next_next = rho_core + drho_dt * dt + (X*S_rho_X+Y*S_rho_Y)/(X+Y)
    u_core_next_next = u_core + du_dt * dt + (X*S_u_X+Y*S_u_Y)/(X+Y)
    T_core_next_next = T_core + dT_dt * dt + (X*S_T_X+Y*S_T_Y)/(X+Y)
    P_core_next_next = rho_core_next_next * R * T_core_next_next

    # Appending arrays 
    # P at inlet is specified, need to check about outlet
    # Specify P at inlet, let density float, find T using ideal gas law
    # Check if P floats at outlet or not, let density float, find T using ideal gas law
    if M_outlet < 1.0:
        P_vector_next_next = np.append(np.append(np.array([P_inlet]),P_core_next_next),np.array([P_outlet]))
    else:
        P_vector_next_next = np.append(np.append(np.array([P_inlet]),P_core_next_next),np.array([2*P_core_next_next[-1]-P_core_next_next[-2]]))
    if sonicInput == 1:
        u_vector_next_next = np.append(np.append(np.array([c]),u_core_next_next),np.array([2*u_core_next_next[-1]-u_core_next_next[-2]]))
    else:
        u_vector_next_next = np.append(np.append(np.array([2*u_core_next_next[0]-u_core_next_next[1]]),u_core_next_next),np.array([2*u_core_next_next[-1]-u_core_next_next[-2]]))
    T_vector_next_next = np.append(np.append(np.array([T_inlet]),T_core_next_next),np.array([2*T_core_next_next[-1]-T_core_next_next[-2]]))
    rho_vector_next_next = np.append(np.append(np.array([P_vector_next_next[0]/R/T_vector_next_next[0]]),rho_core_next_next),np.array([P_vector_next_next[-1]/R/T_vector_next_next[-1]]))

    # while loop for converging (comparing first iteration next to next next until solution stops changing)
    current_iteration = 0
    while (max(abs((P_vector_next_next-P_vector_next)/P_vector_next)) > convergenceError or max(abs((T_vector_next_next-T_vector_next)/T_vector_next)) > convergenceError or max(abs((rho_vector_next_next-rho_vector_next)/rho_vector_next))) > convergenceError and current_iteration < maxIterationsStep:
        if min(P_vector_next_next) < 0 or min(rho_vector_next_next) < 0 or min(T_vector_next_next) < 0:
            dt = math.pow(dt,C)
        
        # print(min(rho_core_next))
        current_iteration = current_iteration + 1
        # core vector update
        rho_core_next = rho_core_next_next
        P_core_next = P_core_next_next
        T_core_next = T_core_next_next
        u_core_next = u_core_next_next
        # vector update
        rho_vector_next = rho_vector_next_next
        P_vector_next = P_vector_next_next
        T_vector_next = T_vector_next_next
        u_vector_next = u_vector_next_next
        # new time derivative with weights
        # # X weights
        # drho_dt_X = -u_core * (rho_vector[2:lenV]-rho_vector[0:-2]) / (2*dx) \
        #     - rho_core * (u_vector[2:lenV]-u_vector[0:-2]) / (2*dx) \
        #     - rho_core * u_core * (np.log(h_vector[2:lenV])-np.log(h_vector[0:-2])) / (2*dx)
        # du_dt_X = -u_core * (u_vector[2:lenV]-u_vector[0:-2]) / (2*dx) \
        #     - 1/rho_core * (P_vector[2:lenV]-P_vector[0:-2]) / (2*dx) \
        #     - u_core/rho_core * fd * 1 / h_vector[1:-1]
        # dT_dt_X = -u_core * (T_vector[2:lenV]-T_vector[0:-2]) / (2*dx) \
        #     - 1/cv*P_core/rho_core * (u_vector[2:lenV]-u_vector[0:-2]) / (2*dx) \
        #     - 1/cv*P_core*u_core/rho_core * (np.log(h_vector[2:lenV])-np.log(h_vector[0:-2])) / (2*dx) \
        #     + 1/cv/rho_core * kcond * (T_vector[2:lenV]-2*T_core+T_vector[0:-2])/(math.pow(dx,2)) \
        #     + 1/cv/rho_core*kcond*(T_vector[2:lenV]-T_vector[0:-2]) / (2*dx) * (np.log(h_vector[2:lenV])-np.log(h_vector[0:-2])) / (2*dx) 
        # Y weights
        drho_dt_Y = -u_core_next * (rho_vector_next[2:lenV]-rho_vector_next[0:-2]) / (2*dx) \
            - rho_core_next * (u_vector_next[2:lenV]-u_vector_next[0:-2]) / (2*dx) \
            - rho_core_next * u_core_next * (np.log(h_vector[2:lenV])-np.log(h_vector[0:-2])) / (2*dx)
        du_dt_Y = -u_core_next * (u_vector_next[2:lenV]-u_vector_next[0:-2]) / (2*dx) \
            - 1/rho_core_next * (P_vector_next[2:lenV]-P_vector_next[0:-2]) / (2*dx) \
            - u_core_next/rho_core_next * fd * 1 / h_vector[1:-1]
        dT_dt_Y = -u_core_next * (T_vector_next[2:lenV]-T_vector_next[0:-2]) / (2*dx) \
            - 1/cv*P_core_next/rho_core_next * (u_vector_next[2:lenV]-u_vector_next[0:-2]) / (2*dx) \
            - 1/cv*P_core_next*u_core_next/rho_core_next * (np.log(h_vector[2:lenV])-np.log(h_vector[0:-2])) / (2*dx) \
            + 1/cv/rho_core_next * kcond * (T_vector_next[2:lenV]-2*T_core_next+T_vector_next[0:-2])/(math.pow(dx,2)) \
            + 1/cv/rho_core_next*kcond*(T_vector_next[2:lenV]-T_vector_next[0:-2]) / (2*dx) * (np.log(h_vector[2:lenV])-np.log(h_vector[0:-2])) / (2*dx) \
            + 1/cv/rho_core_next*q_array[1:-1]
        # Artificial viscosity
        S_rho_Y = np.abs(P_vector_next[2:lenV]-2*P_core_next+P_vector_next[0:-2])/(P_vector_next[2:lenV]+2*P_core_next+P_vector_next[0:-2]) \
            * (rho_vector_next[2:lenV]-2*rho_core_next+rho_vector_next[0:-2])
        S_u_Y = np.abs(P_vector_next[2:lenV]-2*P_core_next+P_vector_next[0:-2])/(P_vector_next[2:lenV]+2*P_core_next+P_vector_next[0:-2]) \
            * (u_vector_next[2:lenV]-2*u_core_next+u_vector_next[0:-2])
        S_T_Y = np.abs(P_vector_next[2:lenV]-2*P_core_next+P_vector_next[0:-2])/(P_vector_next[2:lenV]+2*P_core_next+P_vector_next[0:-2]) \
            * (T_vector_next[2:lenV]-2*T_core_next+T_vector_next[0:-2])
        # S_rho_Y = np.nan_to_num(S_rho_Y,nan=0.0)
        # S_rho_Y[S_rho_Y==inf]=0
        # S_u_Y = np.nan_to_num(S_u_Y,nan=0.0)
        # S_u_Y[S_u_Y==inf]=0
        # S_T_Y = np.nan_to_num(S_T_Y,nan=0.0)
        # S_T_Y[S_T_Y==inf]=0
        # Weighted derivatives
        drho_dt = (X*drho_dt_X+Y*drho_dt_Y)/(X+Y)
        du_dt = (X*du_dt_X+Y*du_dt_Y)/(X+Y)
        dT_dt = (X*dT_dt_X+Y*dT_dt_Y)/(X+Y)

        # Estimate at next step
        rho_core_next_next = rho_core + drho_dt * dt + (X*S_rho_X+Y*S_rho_Y)/(X+Y)
        u_core_next_next = u_core + du_dt * dt + (X*S_u_X+Y*S_u_Y)/(X+Y)
        T_core_next_next = T_core + dT_dt * dt + (X*S_T_X+Y*S_T_Y)/(X+Y)
        P_core_next_next = rho_core_next_next * R * T_core_next_next

        # Appending arrays 
        # P at inlet is specified, need to check about outlet
        # Specify P at inlet, let density float, find T using ideal gas law
        # Check if P floats at outlet or not, let density float, find T using ideal gas law
        if M_outlet < 1.0:
            P_vector_next_next = np.append(np.append(np.array([P_inlet]),P_core_next_next),np.array([P_outlet]))
        else:
            P_vector_next_next = np.append(np.append(np.array([P_inlet]),P_core_next_next),np.array([2*P_core_next_next[-1]-P_core_next_next[-2]]))
        if sonicInput == 1:
            u_vector_next_next = np.append(np.append(np.array([c]),u_core_next_next),np.array([2*u_core_next_next[-1]-u_core_next_next[-2]]))
        else:
            u_vector_next_next = np.append(np.append(np.array([2*u_core_next_next[0]-u_core_next_next[1]]),u_core_next_next),np.array([2*u_core_next_next[-1]-u_core_next_next[-2]]))
        T_vector_next_next = np.append(np.append(np.array([T_inlet]),T_core_next_next),np.array([2*T_core_next_next[-1]-T_core_next_next[-2]]))
        rho_vector_next_next = np.append(np.append(np.array([P_vector_next_next[0]/R/T_vector_next_next[0]]),rho_core_next_next),np.array([P_vector_next_next[-1]/R/T_vector_next_next[-1]]))

    rho_vector = rho_vector_next_next
    P_vector = P_vector_next_next
    T_vector = T_vector_next_next
    u_vector = u_vector_next_next

    return rho_vector, P_vector, T_vector, u_vector, current_iteration, dt