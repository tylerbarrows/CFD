def nextStepFunction_1D(x_vector,h_vector,P_vector,T_vector,rho_vector,u_vector,dx,dt,fluidProperties,solverInformation,initialConditions,boundaryConditions,PreExamResults):

    # Import
    import math 
    import numpy as np

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
        + 1/cv/rho_core*kcond*(T_vector[2:lenV]-T_vector[0:-2]) / (2*dx) * (np.log(h_vector[2:lenV])-np.log(h_vector[0:-2])) / (2*dx) 
    
    # First estimate at new step
    C = 0.9
    rho_core_next = rho_core + drho_dt * dt
    u_core_next = u_core + du_dt * dt
    T_core_next = T_core + dT_dt * dt
    # ints = 1
    # while min(rho_core_next) < 0:
    #     ints = ints + 1
    #     rho_core_next = rho_core + math.pow(C,ints) * drho_dt * dt
    # ints = 1
    # while min(u_core_next) < 0:
    #     ints = ints + 1
    #     u_core_next = u_core + math.pow(C,ints) * du_dt * dt
    # ints = 1
    # while min(T_core_next) < 0:
    #     ints = ints + 1
    #     T_core_next = T_core + math.pow(C,ints) * dT_dt * dt

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
    drho_dt = -(X*u_core+Y*u_core_next)/(X+Y) * ((X*rho_vector[2:lenV]+Y*rho_vector_next[2:lenV])/(X+Y)-(X*rho_vector[0:-2]+Y*rho_vector_next[0:-2])/(X+Y)) / (2*dx) \
        - (X*rho_core+Y*rho_core_next)/(X+Y) * ((X*u_vector[2:lenV]+Y*u_vector_next[2:lenV])/(X+Y)-(X*u_vector[0:-2]+Y*u_vector_next[0:-2])/(X+Y)) / (2*dx) \
        - (X*rho_core+Y*rho_core_next)/(X+Y) * (X*u_core+Y*u_core_next)/(X+Y) * (np.log(h_vector[2:lenV])-np.log(h_vector[0:-2])) / (2*dx)
    du_dt = -(X*u_core+Y*u_core_next)/(X+Y) * ((X*u_vector[2:lenV]+Y*u_vector_next[2:lenV])/(X+Y)-(X*u_vector[0:-2]+Y*u_vector_next[0:-2])/(X+Y)) / (2*dx) \
        - 1/(X*rho_core+Y*rho_core_next)*(X+Y) * ((X*P_vector[2:lenV]+Y*P_vector_next[2:lenV])/(X+Y)-(X*P_vector[0:-2]+Y*P_vector_next[0:-2])/(X+Y)) / (2*dx) \
        - 1/(X*rho_core+Y*rho_core_next)*(X*u_core+Y*u_core_next)*fd*1/h_vector[1:-1]
    dT_dt = -(X*u_core+Y*u_core_next)/(X+Y) * ((X*T_vector[2:lenV]+Y*T_vector_next[2:lenV])/(X+Y)-(X*T_vector[0:-2]+Y*T_vector_next[0:-2])/(X+Y)) / (2*dx) \
        - 1/cv*(X*P_core+Y*P_core_next)/(X*rho_core+Y*rho_core_next) * ((X*u_vector[2:lenV]+Y*u_vector_next[2:lenV])/(X+Y)-(X*u_vector[0:-2]+Y*u_vector_next[0:-2])/(X+Y)) / (2*dx) \
        - 1/cv*(X*P_core+Y*P_core_next)/(X*rho_core+Y*rho_core_next) * (X*u_core+Y*u_core_next)/(X+Y) * (np.log(h_vector[2:lenV])-np.log(h_vector[0:-2])) / (2*dx) \
        + 1/cv/(X*rho_core+Y*rho_core_next)*(X+Y)*kcond*((X*T_vector[2:lenV]+Y*T_vector_next[2:lenV])/(X+Y)-2*(X*T_core+Y*T_core_next)/(X+Y)+(X*T_vector[0:-2]+Y*T_vector_next[0:-2])/(X+Y))/(math.pow(dx,2)) \
        + 1/cv/(X*rho_core+Y*rho_core_next)*(X+Y)*kcond*((X*T_vector[2:lenV]+Y*T_vector_next[2:lenV])/(X+Y)-(X*T_vector[0:-2]+Y*T_vector_next[0:-2])/(X+Y))/(2*dx) * (np.log(h_vector[2:lenV])-np.log(h_vector[0:-2])) / (2*dx) 

    # Estimate at next step
    rho_core_next_next = rho_core + drho_dt * dt
    u_core_next_next = u_core + du_dt * dt
    T_core_next_next = T_core + dT_dt * dt
    # ints = 1
    # while min(rho_core_next_next) < 0:
    #     ints = ints + 1
    #     rho_core_next_next = rho_core + math.pow(C,ints) * drho_dt * dt
    # ints = 1
    # while min(u_core_next_next) < 0:
    #     ints = ints + 1
    #     u_core_next_next = u_core + math.pow(C,ints) * du_dt * dt
    # ints = 1
    # while min(T_core_next_next) < 0:
    #     ints = ints + 1
    #     T_core_next_next = T_core + math.pow(C,ints) * dT_dt * dt
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
        drho_dt = -(X*u_core+Y*u_core_next)/(X+Y) * ((X*rho_vector[2:lenV]+Y*rho_vector_next[2:lenV])/(X+Y)-(X*rho_vector[0:-2]+Y*rho_vector_next[0:-2])/(X+Y)) / (2*dx) \
            - (X*rho_core+Y*rho_core_next)/(X+Y) * ((X*u_vector[2:lenV]+Y*u_vector_next[2:lenV])/(X+Y)-(X*u_vector[0:-2]+Y*u_vector_next[0:-2])/(X+Y)) / (2*dx) \
            - (X*rho_core+Y*rho_core_next)/(X+Y) * (X*u_core+Y*u_core_next)/(X+Y) * (np.log(h_vector[2:lenV])-np.log(h_vector[0:-2])) / (2*dx)
        du_dt = -(X*u_core+Y*u_core_next)/(X+Y) * ((X*u_vector[2:lenV]+Y*u_vector_next[2:lenV])/(X+Y)-(X*u_vector[0:-2]+Y*u_vector_next[0:-2])/(X+Y)) / (2*dx) \
            - 1/(X*rho_core+Y*rho_core_next)*(X+Y) * ((X*P_vector[2:lenV]+Y*P_vector_next[2:lenV])/(X+Y)-(X*P_vector[0:-2]+Y*P_vector_next[0:-2])/(X+Y)) / (2*dx) \
            - 1/(X*rho_core+Y*rho_core_next)*(X*u_core+Y*u_core_next)*fd*1/h_vector[1:-1]
        dT_dt = -(X*u_core+Y*u_core_next)/(X+Y) * ((X*T_vector[2:lenV]+Y*T_vector_next[2:lenV])/(X+Y)-(X*T_vector[0:-2]+Y*T_vector_next[0:-2])/(X+Y)) / (2*dx) \
            - 1/cv*(X*P_core+Y*P_core_next)/(X*rho_core+Y*rho_core_next) * ((X*u_vector[2:lenV]+Y*u_vector_next[2:lenV])/(X+Y)-(X*u_vector[0:-2]+Y*u_vector_next[0:-2])/(X+Y)) / (2*dx) \
            - 1/cv*(X*P_core+Y*P_core_next)/(X*rho_core+Y*rho_core_next)*(X*u_core+Y*u_core_next)/(X+Y) * (np.log(h_vector[2:lenV])-np.log(h_vector[0:-2])) / (2*dx) \
            + 1/cv/(X*rho_core+Y*rho_core_next)*(X+Y)*kcond*((X*T_vector[2:lenV]+Y*T_vector_next[2:lenV])/(X+Y)-2*(X*T_core+Y*T_core_next)/(X+Y)+(X*T_vector[0:-2]+Y*T_vector_next[0:-2])/(X+Y))/(math.pow(dx,2)) \
            + 1/cv/(X*rho_core+Y*rho_core_next)*(X+Y)*kcond*((X*T_vector[2:lenV]+Y*T_vector_next[2:lenV])/(X+Y)-(X*T_vector[0:-2]+Y*T_vector_next[0:-2])/(X+Y))/(2*dx) * (np.log(h_vector[2:lenV])-np.log(h_vector[0:-2])) / (2*dx) \
            + 1/cv/(X*rho_core+Y*rho_core_next)*np.power((X*u_core+Y*u_core_next),2)/(X+Y)

        # Estimate at next step
        rho_core_next_next = rho_core + drho_dt * dt
        u_core_next_next = u_core + du_dt * dt
        T_core_next_next = T_core + dT_dt * dt
        # ints = 1
        # while min(rho_core_next_next) < 0:
        #     ints = ints + 1
        #     rho_core_next_next = rho_core + math.pow(C,ints) * drho_dt * dt
        # ints = 1
        # while min(u_core_next_next) < 0:
        #     ints = ints + 1
        #     u_core_next_next = u_core + math.pow(C,ints) * du_dt * dt
        # ints = 1
        # while min(T_core_next_next) < 0:
        #     ints = ints + 1
        #     T_core_next_next = T_core + math.pow(C,ints) * dT_dt * dt
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

    return rho_vector, P_vector, T_vector, u_vector, current_iteration