def isentropicSolutionFunc(meshGeneration,geometryInputs,fluidProperties,initialConditions,boundaryConditions,solverInformation,SS_Solution,PreExamResults):
    # First need to check if converging or diverging vs converging-diverging
    # if converging or diverging, can calculate things directly
    # if converging-diverging, need to check At/Ae and Po/Pe to see if throat is choked 
    # bisection search to find shock location (if applicable)
    # isentropic upstream and downstream of shock

    # import
    import math 
    import numpy as np

    # Fluid Properties
    k = fluidProperties.k
    R = fluidProperties.R
    cp = fluidProperties.cp
    cv = fluidProperties.cv
    kcond = fluidProperties.kcond

    # Inlet conditions and outlet
    P_inlet = initialConditions.inlet.P * 100000
    P_outlet = initialConditions.outlet.P * 100000
    u_inlet = initialConditions.inlet.u
    T_inlet = initialConditions.inlet.T + 273.15

    # Geometry vectors
    x_vector = meshGeneration.x_vector
    h_vector = meshGeneration.h_vector

    # Error
    error = solverInformation.error_Isen
    maxInt = solverInformation.iterationInformation.IsenMax

    # Geometry - areas
    Ai = geometryInputs.heightInlet
    At = geometryInputs.heightThroat
    Ae = geometryInputs.heightOutlet
    shape = geometryInputs.shape

    # Solver error
    errorFindingThroatChoked = solverInformation.error_SS

    # CFD Data
    P_track = SS_Solution.P_track
    rho_track = SS_Solution.rho_track
    T_track = SS_Solution.T_track
    u_track = SS_Solution.u_track
    M_track = u_track[-1,:] / np.sqrt(k*R*T_track[-1,:])

    # flag
    flag_isen = 0 # not sure what I'm doing with this
    if PreExamResults.sonicInput == 1:
        inputFlag = 'sonic'
    else:
        inputFlag = 'subsonic'

    match shape:

        case 'converging':
            Pi = P_inlet
            Ti = T_inlet
            Pe = P_outlet
            # Astar bounds, iterate on Astar until finding Mi and Po that give Pe
            Astar_l = 0
            Astar_u = Ae
            for i in range(0,maxInt):
                Astar = 1/2*(Astar_l+Astar_u)
                Ai_Astar = Ai/Astar
                Mi_l = 0
                Mi_u = 1
                Mi = 1/2*(Mi_l+Mi_u)
                Ai_Astar_fM = 1/Mi*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Mi,2)),(k+1)/2/(k-1))
                for x in range(0,maxInt):
                    if Ai_Astar_fM > Ai_Astar:
                        Mi_u = Mi_u
                        Mi_l = Mi
                    else:
                        Mi_u = Mi
                        Mi_l = Mi_l
                    Mi = 1/2*(Mi_u+Mi_l)
                    Ai_Astar_fM = 1/Mi*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Mi,2)),(k+1)/2/(k-1))
                Po = Pi * math.pow(1+(k-1)/2*math.pow(Mi,2),k/(k-1))
                To = Ti * math.pow(Po/Pi,(k-1)/k)
                Ae_Astar = Ae/Astar
                Me_l = 0
                Me_u = 1
                Me = 1/2*(Me_l+Me_u)
                Ae_Astar_fM = 1/Me*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Me,2)),(k+1)/2/(k-1))
                for x in range(0,maxInt):
                    if Ae_Astar_fM > Ae_Astar:
                        Me_u = Me_u
                        Me_l = Me
                    else:
                        Me_u = Me
                        Me_l = Me_l
                    Me = 1/2*(Me_u+Me_l)
                    Ae_Astar_fM = 1/Me*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Me,2)),(k+1)/2/(k-1))
                Pe_calc = Po / math.pow(1+(k-1)/2*math.pow(Me,2),k/(k-1))
                if Pe_calc > Pe:
                    Astar_u = Astar_u
                    Astar_l = Astar
                else:
                    Astar_u = Astar
                    Astar_l = Astar_l
            # now that we have Astar, find all the values we want
            A_array = np.linspace(Ai,Ae,len(x_vector))
            M_array = np.linspace(1,1,len(x_vector))
            T_array = np.linspace(1,1,len(x_vector))
            P_array = np.linspace(1,1,len(x_vector))
            u_array = np.linspace(1,1,len(x_vector))
            A_Astar_array = A_array/Astar
            for i in range(0,len(M_array)+1):
                A_Astar = A_Astar_array[i-1]
                M_l = 0
                M_u = 1
                M = 1/2*(M_l+M_u)
                A_Astar_fM = 1/M*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(M,2)),(k+1)/2/(k-1))
                for x in range(0,maxInt):
                    if A_Astar_fM > A_Astar:
                        M_u = M_u
                        M_l = M
                    else:
                        M_u = M
                        M_l = M_l
                    M = 1/2*(M_u+M_l)
                    A_Astar_fM = 1/M*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(M,2)),(k+1)/2/(k-1))
                P = Po / math.pow(1+(k-1)/2*math.pow(M,2),k/(k-1))
                T = To / (1+(k-1)/2*math.pow(M,2))
                u = M * math.sqrt(k*R*T)
                P_array[i-1] = P
                T_array[i-1] = T
                u_array[i-1] = u
                M_array[i-1] = M

        case 'diverging':
            # two cases to consider, sonic input and not sonic input
            # do case thing with sonic and subsonic
            # sonic 
            # As bounds, inlet and outlet
            # iterate until find Pe that is input
            # loop over A_array treating A < As without shock and A > As with shock
            # subsonic
            # Astar bounds
            # iterate until find Pe that is input
            # loop over A_array 

            match inputFlag:
                case 'subsonic':
                    Pi = P_inlet
                    Ti = T_inlet
                    Pe = P_outlet
                    Astar_l = 0
                    Astar_u = Ai
                    for i in range(0,maxInt):
                        Astar = 1/2*(Astar_l+Astar_u)
                        Ai_Astar = Ai/Astar
                        Mi_l = 0
                        Mi_u = 1
                        Mi = 1/2*(Mi_l+Mi_u)
                        Ai_Astar_fM = 1/Mi*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Mi,2)),(k+1)/2/(k-1))
                        for x in range(0,maxInt):
                            if Ai_Astar_fM > Ai_Astar:
                                Mi_u = Mi_u
                                Mi_l = Mi
                            else:
                                Mi_u = Mi
                                Mi_l = Mi_l
                            Mi = 1/2*(Mi_u+Mi_l)
                            Ai_Astar_fM = 1/Mi*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Mi,2)),(k+1)/2/(k-1))
                        Po = Pi * math.pow(1+(k-1)/2*math.pow(Mi,2),k/(k-1))
                        To = Ti * math.pow(Po/Pi,(k-1)/k)
                        Ae_Astar = Ae/Astar
                        Me_l = 0
                        Me_u = 1
                        Me = 1/2*(Me_l+Me_u)
                        Ae_Astar_fM = 1/Me*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Me,2)),(k+1)/2/(k-1))
                        for x in range(0,maxInt):
                            if Ae_Astar_fM > Ae_Astar:
                                Me_u = Me_u
                                Me_l = Me
                            else:
                                Me_u = Me
                                Me_l = Me_l
                            Me = 1/2*(Me_u+Me_l)
                            Ae_Astar_fM = 1/Me*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Me,2)),(k+1)/2/(k-1))
                        Pe_calc = Po / math.pow(1+(k-1)/2*math.pow(Me,2),k/(k-1))
                        if Pe_calc > Pe:
                            Astar_u = Astar
                            Astar_l = Astar_l
                        else:
                            Astar_u = Astar_u
                            Astar_l = Astar
                    
                    A_array = np.linspace(Ai,Ae,len(x_vector))
                    M_array = np.linspace(1,1,len(x_vector))
                    T_array = np.linspace(1,1,len(x_vector))
                    P_array = np.linspace(1,1,len(x_vector))
                    u_array = np.linspace(1,1,len(x_vector))
                    A_Astar_array = A_array/Astar

                    for i in range(0,len(M_array)+1):
                        A_Astar = A_Astar_array[i-1]
                        M_l = 0
                        M_u = 1
                        M = 1/2*(M_l+M_u)
                        A_Astar_fM = 1/M*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(M,2)),(k+1)/2/(k-1))
                        for x in range(0,maxInt):
                            if A_Astar_fM > A_Astar:
                                M_u = M_u
                                M_l = M
                            else:
                                M_u = M
                                M_l = M_l
                            M = 1/2*(M_u+M_l)
                            A_Astar_fM = 1/M*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(M,2)),(k+1)/2/(k-1))
                        P = Po / math.pow(1+(k-1)/2*math.pow(M,2),k/(k-1))
                        T = To / (1+(k-1)/2*math.pow(M,2))
                        u = M * math.sqrt(k*R*T)
                        P_array[i-1] = P
                        T_array[i-1] = T
                        u_array[i-1] = u
                        M_array[i-1] = M
                
                case 'sonic':
                    Pi = P_inlet
                    Ti = T_inlet
                    Pe = P_outlet
                    Mi = 1
                    Astar = Ai
                    Po = Pi * math.pow(1+(k-1)/2*math.pow(Mi,2),k/(k-1))
                    To = Ti * math.pow(Po/Pi,(k-1)/k)
                    As_u = Ae
                    As_l = Ai
                    for i in range(0,maxInt):
                        As = 1/2*(As_l+As_u)
                        As_Astar = As/Astar
                        Ms_l = 1
                        Ms_u = 100
                        Ms = 1/2*(Ms_l+Ms_u)
                        As_Astar_fM = 1/Ms*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Ms,2)),(k+1)/2/(k-1))
                        for x in range(0,maxInt): 
                            if As_Astar_fM > As_Astar:
                                Ms_u = Ms
                                Ms_l = Ms_l
                            else:
                                Ms_u = Ms_u
                                Ms_l = Ms
                            Ms = 1/2*(Ms_u+Ms_l)
                            As_Astar_fM = 1/Ms*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Ms,2)),(k+1)/2/(k-1))
                        Ps = Po / math.pow(1+(k-1)/2*math.pow(Ms,2),k/(k-1))
                        Ts = To * math.pow(Ps/Po,(k-1)/k)
                        rhos = Ps*100000/Ts/R
                        us = Ms * math.sqrt(k*R*Ts)
                        c1 = rhos*us
                        c2 = Ps*100000+rhos*math.pow(us,2)
                        c3 = k/(k-1)*Ps*100000/rhos+1/2*math.pow(us,2)
                        if math.pow(k/(k-1)*c2,2)-4*c3*1/2*(k+1)/(k-1)*math.pow(c1,2) < 0:
                            rho =  (k/(k-1)*c2+0)/2/c3
                        else:
                            rho = (k/(k-1)*c2+math.sqrt(math.pow(k/(k-1)*c2,2)-4*c3*1/2*(k+1)/(k-1)*math.pow(c1,2)))/2/c3
                        u = c1/rho
                        P = (c2-rho*math.pow(u,2))/100000
                        T = P*100000/R/rho
                        M = u/math.sqrt(k*R*T)
                        Pe_o = P * math.pow(1+(k-1)/2*math.pow(M,2),k/(k-1))
                        Te_o = T * math.pow(Pe_o/P,(k-1)/k)
                        As_Astar = 1/M*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(M,2)),(k+1)/2/(k-1))
                        Astar_o = As/As_Astar
                        Ae_Astar = Ae/Astar_o
                        Me_u = 1
                        Me_l = 0
                        Me = 1/2*(Me_u+Me_l)
                        Ae_Astar_fM = 1/Me*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Me,2)),(k+1)/2/(k-1))
                        for x in range(0,maxInt):
                            if Ae_Astar_fM > Ae_Astar:
                                Me_u = Me_u
                                Me_l = Me
                            else:
                                Me_u = Me
                                Me_l = Me_l
                            Me = 1/2*(Me_u+Me_l)
                            Ae_Astar_fM = 1/Me*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Me,2)),(k+1)/2/(k-1))
                        Pe_calc = Pe_o / math.pow(1+(k-1)/2*math.pow(Me,2),k/(k-1))
                        if Pe_calc > Pe:
                            As_u = As_u
                            As_l = As
                        else:
                            As_u = As
                            As_l = As_l

                    A_array = np.linspace(Ai,Ae,len(x_vector))
                    M_array = np.linspace(1,1,len(x_vector))
                    T_array = np.linspace(1,1,len(x_vector))
                    P_array = np.linspace(1,1,len(x_vector))
                    u_array = np.linspace(1,1,len(x_vector))
                    
                    index = 0
                    for A in A_array:
                        if A > As:
                            # here put shock things
                            A_Astar = A/Astar_o
                            M_l = 0
                            M_u = 1
                            for i in range(0,maxInt):
                                M = 1/2*(M_l+M_u)
                                A_Astar_fM = 1/M*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(M,2)),(k+1)/2/(k-1))
                                if A_Astar_fM > A_Astar:
                                    M_u = M_u
                                    M_l = M
                                else:
                                    M_u = M
                                    M_l = M_l
                            P = Pe_o / math.pow(1+(k-1)/2*math.pow(M,2),k/(k-1))
                            T = Te_o * math.pow(P/Pe_o,(k-1)/k)
                            u = M*math.sqrt(k*R*T)

                            M_array[index] = M
                            T_array[index] = T
                            P_array[index] = P
                            u_array[index] = u
                        else:
                            # only isentropic expansion
                            A_Astar = A/Astar
                            M_l = 1
                            M_u = 100
                            for i in range(0,maxInt):
                                M = 1/2*(M_l+M_u)
                                A_Astar_fM = 1/M*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(M,2)),(k+1)/2/(k-1))
                                if A_Astar_fM > A_Astar:
                                    M_u = M
                                    M_l = M_l
                                else:
                                    M_u = M_u
                                    M_l = M
                            P = Po / math.pow(1+(k-1)/2*math.pow(M,2),k/(k-1))
                            T = To * math.pow(P/Po,(k-1)/k)
                            u = M*math.sqrt(k*R*T)

                            M_array[index] = M
                            T_array[index] = T
                            P_array[index] = P
                            u_array[index] = u

                        index = index + 1
                
                case _:
                    print('Should never see this message')

        case 'converging-diverging':
            Aint = np.where(h_vector==At)
            # Mt = M_track[-1,Aint[0][0]]
            Mt = M_track[Aint[0][0]]
            Mt = max(M_track)
            # print(Mt)
            # print(max(M_track))
            if Mt * (1+errorFindingThroatChoked) > 1:
                inputFlag = 'sonic'
            match inputFlag:
                case 'subsonic':
                    Pi = P_inlet
                    Ti = T_inlet
                    Pe = P_outlet
                    Astar_l = 0
                    Astar_u = At
                    for i in range(0,maxInt):
                        Astar = 1/2*(Astar_l+Astar_u)
                        Ai_Astar = Ai/Astar
                        Mi_l = 0
                        Mi_u = 1
                        Mi = 1/2*(Mi_l+Mi_u)
                        Ai_Astar_fM = 1/Mi*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Mi,2)),(k+1)/2/(k-1))
                        for x in range(0,maxInt):
                            if Ai_Astar_fM > Ai_Astar:
                                Mi_u = Mi_u
                                Mi_l = Mi
                            else:
                                Mi_u = Mi
                                Mi_l = Mi_l
                            Mi = 1/2*(Mi_u+Mi_l)
                            Ai_Astar_fM = 1/Mi*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Mi,2)),(k+1)/2/(k-1))
                        Po = Pi * math.pow(1+(k-1)/2*math.pow(Mi,2),k/(k-1))
                        To = Ti * math.pow(Po/Pi,(k-1)/k)
                        Ae_Astar = Ae/Astar
                        Me_l = 0
                        Me_u = 1
                        Me = 1/2*(Me_l+Me_u)
                        Ae_Astar_fM = 1/Me*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Me,2)),(k+1)/2/(k-1))
                        for x in range(0,maxInt):
                            if Ae_Astar_fM > Ae_Astar:
                                Me_u = Me_u
                                Me_l = Me
                            else:
                                Me_u = Me
                                Me_l = Me_l
                            Me = 1/2*(Me_u+Me_l)
                            Ae_Astar_fM = 1/Me*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Me,2)),(k+1)/2/(k-1))
                        Pe_calc = Po / math.pow(1+(k-1)/2*math.pow(Me,2),k/(k-1))
                        if Pe_calc > Pe:
                            # Astar_u = Astar
                            # Astar_l = Astar_l
                            Astar_u = Astar_u
                            Astar_l = Astar
                        else:
                            # Astar_u = Astar_u
                            # Astar_l = Astar
                            Astar_u = Astar
                            Astar_l = Astar_l
                    # print(Pe_calc)
                    # A_array = np.linspace(Ai,Ae,len(x_vector))
                    A_array = h_vector
                    M_array = np.linspace(1,1,len(x_vector))
                    T_array = np.linspace(1,1,len(x_vector))
                    P_array = np.linspace(1,1,len(x_vector))
                    u_array = np.linspace(1,1,len(x_vector))
                    A_Astar_array = A_array/Astar

                    for i in range(0,len(M_array)+1):
                        A_Astar = A_Astar_array[i-1]
                        M_l = 0
                        M_u = 1
                        M = 1/2*(M_l+M_u)
                        A_Astar_fM = 1/M*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(M,2)),(k+1)/2/(k-1))
                        for x in range(0,maxInt):
                            if A_Astar_fM > A_Astar:
                                M_u = M_u
                                M_l = M
                            else:
                                M_u = M
                                M_l = M_l
                            M = 1/2*(M_u+M_l)
                            A_Astar_fM = 1/M*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(M,2)),(k+1)/2/(k-1))
                        P = Po / math.pow(1+(k-1)/2*math.pow(M,2),k/(k-1))
                        T = To / (1+(k-1)/2*math.pow(M,2))
                        u = M * math.sqrt(k*R*T)
                        P_array[i-1] = P
                        T_array[i-1] = T
                        u_array[i-1] = u
                        M_array[i-1] = M

                case 'sonic':
                    Pi = P_inlet
                    Ti = T_inlet
                    Pe = P_outlet
                    Mt = 1
                    Astar = At
                    Ai_Astar = Ai/Astar
                    Mi_l = 0
                    Mi_u = 1
                    for i in range(0,maxInt):
                        Mi = 1/2*(Mi_l+Mi_u)
                        Ai_Astar_fM = 1/Mi*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Mi,2)),(k+1)/2/(k-1))
                        if Ai_Astar_fM > Ai_Astar:
                            Mi_u = Mi_u
                            Mi_l = Mi
                        else:
                            Mi_u = Mi
                            Mi_l = Mi_l
                    Po = Pi * math.pow(1+(k-1)/2*math.pow(Mi,2),k/(k-1))
                    To = Ti * math.pow(Po/Pi,(k-1)/k)
                    As_u = Ae
                    As_l = At
                    for i in range(0,maxInt):
                        As = 1/2*(As_l+As_u)
                        As_Astar = As/Astar
                        Ms_l = 1
                        Ms_u = 100
                        Ms = 1/2*(Ms_l+Ms_u)
                        As_Astar_fM = 1/Ms_l*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Ms,2)),(k+1)/2/(k-1))
                        for x in range(0,maxInt): 
                            if As_Astar_fM > As_Astar:
                                Ms_u = Ms
                                Ms_l = Ms_l
                            else:
                                Ms_u = Ms_u
                                Ms_l = Ms
                            Ms = 1/2*(Ms_u+Ms_l)
                            As_Astar_fM = 1/Ms*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Ms,2)),(k+1)/2/(k-1))
                        Ps = Po / math.pow(1+(k-1)/2*math.pow(Ms,2),k/(k-1))
                        Ts = To * math.pow(Ps/Po,(k-1)/k)
                        rhos = Ps*100000/Ts/R
                        us = Ms * math.sqrt(k*R*Ts)
                        c1 = rhos*us
                        c2 = Ps*100000+rhos*math.pow(us,2)
                        c3 = k/(k-1)*Ps*100000/rhos+1/2*math.pow(us,2)
                        if math.pow(k/(k-1)*c2,2)-4*c3*1/2*(k+1)/(k-1)*math.pow(c1,2) < 0:
                            rho =  (k/(k-1)*c2+0)/2/c3
                        else:
                            rho = (k/(k-1)*c2+math.sqrt(math.pow(k/(k-1)*c2,2)-4*c3*1/2*(k+1)/(k-1)*math.pow(c1,2)))/2/c3
                        u = c1/rho
                        P = (c2-rho*math.pow(u,2))/100000
                        T = P*100000/R/rho
                        M = u/math.sqrt(k*R*T)
                        Pe_o = P * math.pow(1+(k-1)/2*math.pow(M,2),k/(k-1))
                        Te_o = T * math.pow(Pe_o/P,(k-1)/k)
                        As_Astar = 1/M*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(M,2)),(k+1)/2/(k-1))
                        Astar_o = As/As_Astar
                        Ae_Astar = Ae/Astar_o
                        Me_u = 1
                        Me_l = 0
                        Me = 1/2*(Me_u+Me_l)
                        Ae_Astar_fM = 1/Me*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Me,2)),(k+1)/2/(k-1))
                        for x in range(0,maxInt):
                            if Ae_Astar_fM > Ae_Astar:
                                Me_u = Me_u
                                Me_l = Me
                            else:
                                Me_u = Me
                                Me_l = Me_l
                            Me = 1/2*(Me_u+Me_l)
                            Ae_Astar_fM = 1/Me*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Me,2)),(k+1)/2/(k-1))
                        Pe_calc = Pe_o / math.pow(1+(k-1)/2*math.pow(Me,2),k/(k-1))
                        if Pe_calc > Pe:
                            As_u = As_u
                            As_l = As
                        else:
                            As_u = As
                            As_l = As_l

                    # A_array = np.linspace(Ai,Ae,len(x_vector))
                    A_array = h_vector
                    M_array = np.linspace(1,1,len(x_vector))
                    T_array = np.linspace(1,1,len(x_vector))
                    P_array = np.linspace(1,1,len(x_vector))
                    u_array = np.linspace(1,1,len(x_vector))
                    
                    index = 0
                    for A in A_array:
                        if A > As and index > Aint[0][0]:
                            # here put shock things
                            A_Astar = A/Astar_o
                            M_l = 0
                            M_u = 1
                            for i in range(0,maxInt):
                                M = 1/2*(M_l+M_u)
                                A_Astar_fM = 1/M*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(M,2)),(k+1)/2/(k-1))
                                if A_Astar_fM > A_Astar:
                                    M_u = M_u
                                    M_l = M
                                else:
                                    M_u = M
                                    M_l = M_l
                            P = Pe_o / math.pow(1+(k-1)/2*math.pow(M,2),k/(k-1))
                            T = Te_o * math.pow(P/Pe_o,(k-1)/k)
                            u = M*math.sqrt(k*R*T)

                            M_array[index] = M
                            T_array[index] = T
                            P_array[index] = P
                            u_array[index] = u
                        else:
                            if index > Aint[0][0]: #A > As:
                                # only isentropic expansion
                                A_Astar = A/Astar
                                M_l = 1
                                M_u = 100
                                for i in range(0,maxInt):
                                    M = 1/2*(M_l+M_u)
                                    A_Astar_fM = 1/M*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(M,2)),(k+1)/2/(k-1))
                                    if A_Astar_fM > A_Astar:
                                        M_u = M
                                        M_l = M_l
                                    else:
                                        M_u = M_u
                                        M_l = M
                                P = Po / math.pow(1+(k-1)/2*math.pow(M,2),k/(k-1))
                                T = To * math.pow(P/Po,(k-1)/k)
                                u = M*math.sqrt(k*R*T)

                                M_array[index] = M
                                T_array[index] = T
                                P_array[index] = P
                                u_array[index] = u
                            else:
                                # only isentropic expansion
                                A_Astar = A/Astar
                                M_l = 0
                                M_u = 1
                                for i in range(0,maxInt):
                                    M = 1/2*(M_l+M_u)
                                    A_Astar_fM = 1/M*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(M,2)),(k+1)/2/(k-1))
                                    if A_Astar_fM > A_Astar:
                                        M_u = M_u
                                        M_l = M
                                    else:
                                        M_u = M
                                        M_l = M_l
                                P = Po / math.pow(1+(k-1)/2*math.pow(M,2),k/(k-1))
                                T = To * math.pow(P/Po,(k-1)/k)
                                u = M*math.sqrt(k*R*T)

                                M_array[index] = M
                                T_array[index] = T
                                P_array[index] = P
                                u_array[index] = u

                        index = index + 1
        case _:
            flag = 1

    return P_array, T_array, P_array/R/T_array, u_array, M_array, flag_isen