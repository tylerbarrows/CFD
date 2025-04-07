def PreExamineFunction(geometryInputs,fluidProperties,initialConditions,boundaryConditions,solverInformation):
    # Function that checks inputs and geometry that could influence how solution is found
    
    # Basic imports
    import math

    # Importing
    shape = geometryInputs.shape
    dimension = geometryInputs.dimension
    P_inlet = initialConditions.inlet.P
    P_outlet = initialConditions.outlet.P
    T_inlet = initialConditions.inlet.T
    k = fluidProperties.k
    R = fluidProperties.R
    Ai = geometryInputs.heightInlet
    At = geometryInputs.heightThroat
    Ae = geometryInputs.heightOutlet
    error = solverInformation.error_pre
    maxInt = solverInformation.iterationInformation.preExamMax
    sonicInputPin = initialConditions.sonicInputPin
    
    class PreExamResults:
        # return PreExamResults for modifying how solution works
        match dimension:
            case '1D':
                match shape:
                    case 'converging':
                        # maybe want to add check to make sure that Pe isn't too low that 
                        # what does this look like
                        # Astar = Ae -> find Pe that corresponds to this 
                        flagValid = 'valid'
                        if sonicInputPin == 1:
                            sonicInput = 1
                        else:
                            sonicInput = 0
                        if Ae > Ai or P_outlet > P_inlet:
                            flagValid = 'not valid'
                        
                        soincInput = 0 # turning this off
                        
                        Astar = Ae
                        Ai_Astar = Ai/Astar
                        Pi = P_inlet
                        Pe = P_outlet
                        Mi_u = 1
                        Mi_l = 0
                        for i in range(0,maxInt):
                            Mi = 1/2*(Mi_l+Mi_u)
                            Ai_Astar_fM = 1/Mi*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Mi,2)),(k+1)/2/(k-1))
                            if Ai_Astar_fM > Ai_Astar:
                                Mi_u = Mi_u
                                Mi_l = Mi
                            else:
                                Mi_u = Mi
                                Mi_l = Mi_l
                        Po_calc = Pi * math.pow(1+(k-1)/2*math.pow(Mi,2),k/(k-1))
                        Me = 1
                        Pe_calc = Po_calc / math.pow(1+(k-1)/2*math.pow(Me,2),k/(k-1))
                        if Pe < Pe_calc:
                            flagValid = 'not valid'

                    case 'diverging':
                        subsonicSolutionFoundFlag = 0
                        sonicSolutionFoundFlag = 0
                        flagValid = 'valid'

                        Pi = P_inlet
                        Pe = P_outlet
                        Ti = T_inlet
                        
                        # looking for subsonic solution
                        Astar_u = Ai
                        Astar_l = 0.01
                        # first finding upper bound where Astar = Ai
                        Ai_Astar = Ai/Astar_u
                        Mi_u = 1
                        Mi_l = 0
                        Mi = 1/2*(Mi_u+Mi_l)
                        Ai_Astar_fM = 1/Mi*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Mi,2)),(k+1)/2/(k-1))
                        for i in range(0,maxInt):
                            if Ai_Astar_fM > Ai_Astar:
                                Mi_u = Mi_u
                                Mi_l = Mi
                            else:
                                Mi_u = Mi
                                Mi_l = Mi_l
                            Mi = 1/2*(Mi_u+Mi_l)
                            Ai_Astar_fM = 1/Mi*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Mi,2)),(k+1)/2/(k-1))
                        Po_calc = Pi * math.pow(1+(k-1)/2*math.pow(Mi,2),k/(k-1))
                        Ae_Astar = Ae/Astar_u
                        Me_u = 1
                        Me_l = 0
                        Me = 1/2*(Me_l+Me_u)
                        Ae_Astar_fM = 1/Me*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Me,2)),(k+1)/2/(k-1))
                        for i in range(0,maxInt):
                            if Ae_Astar_fM > Ae_Astar:
                                Me_u = Me_u
                                Me_l = Me
                            else:
                                Me_u = Me
                                Me_l = Me_l
                            Me = 1/2*(Me_u+Me_l)
                            Ae_Astar_fM = 1/Me*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Me,2)),(k+1)/2/(k-1))
                        Pe_subsonic_u = Po_calc / math.pow(1+(k-1)/2*math.pow(Me,2),k/(k-1))
                        # next finding lower bound where Astar = 0.01
                        Ai_Astar = Ai/Astar_l
                        Mi_u = 1
                        Mi_l = 0
                        Mi = 1/2*(Mi_u+Mi_l)
                        Ai_Astar_fM = 1/Mi*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Mi,2)),(k+1)/2/(k-1))
                        for i in range(0,maxInt):
                            if Ai_Astar_fM > Ai_Astar:
                                Mi_u = Mi_u
                                Mi_l = Mi
                            else:
                                Mi_u = Mi
                                Mi_l = Mi_l
                            Mi = 1/2*(Mi_u+Mi_l)
                            Ai_Astar_fM = 1/Mi*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Mi,2)),(k+1)/2/(k-1))
                        Po_calc = Pi * math.pow(1+(k-1)/2*math.pow(Mi,2),k/(k-1))
                        Ae_Astar = Ae/Astar_l
                        Me_u = 1
                        Me_l = 0
                        Me = 1/2*(Me_l+Me_u)
                        Ae_Astar_fM = 1/Me*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Me,2)),(k+1)/2/(k-1))
                        for i in range(0,maxInt):
                            if Ae_Astar_fM > Ae_Astar:
                                Me_u = Me_u
                                Me_l = Me
                            else:
                                Me_u = Me
                                Me_l = Me_l
                            Me = 1/2*(Me_u+Me_l)
                            Ae_Astar_fM = 1/Me*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Me,2)),(k+1)/2/(k-1))
                        Pe_subsonic_l = Po_calc / math.pow(1+(k-1)/2*math.pow(Me,2),k/(k-1))
                        if Pe > Pe_subsonic_l and Pe < Pe_subsonic_u:
                            subsonicSolutionFoundFlag = 1
                        
                        # next check for sonic solution
                        # looking for sonic solution
                        Ai_Astar = 1
                        Mi = 1
                        Po = Pi * math.pow(1+(k-1)/2*math.pow(Mi,2),k/(k-1))
                        To = Ti * math.pow(Po/Pi,(k-1)/k)
                        # bounds for shock
                        As_u = Ae
                        As_l = Ai
                        # first finding upper bound on Pe
                        As = As_l # higher Pe when shock is closer to inlet
                        As_Astar = As/Ai
                        Ms_u = 100
                        Ms_l = 1
                        Ms = 1/2*(Ms_u+Ms_l)
                        As_Astar_fM = 1/Ms*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Ms,2)),(k+1)/2/(k-1))
                        for i in range(0,maxInt): 
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
                        As_Astar = 1/M*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(M,2)),(k+1)/2/(k-1))
                        Astar = As/As_Astar
                        Ae_Astar = Ae/Astar
                        Me_u = 1
                        Me_l = 0
                        Me = 1/2*(Me_u+Me_l)
                        Ae_Astar_fM = 1/Me*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Me,2)),(k+1)/2/(k-1))
                        for i in range(0,maxInt):
                            if Ae_Astar_fM > Ae_Astar:
                                Me_u = Me_u
                                Me_l = Me
                            else:
                                Me_u = Me
                                Me_l = Me_l
                            Me = 1/2*(Me_u+Me_l)
                            Ae_Astar_fM = 1/Me*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Me,2)),(k+1)/2/(k-1))
                        Pe_sonic_u = Pe_o / math.pow(1+(k-1)/2*math.pow(Me,2),k/(k-1))
                        # next finding lower bound on 
                        As = As_u # higher Pe when shock is closer to inlet
                        As_Astar = As/Ai
                        Ms_u = 100
                        Ms_l = 1
                        Ms = 1/2*(Ms_u+Ms_l)
                        As_Astar_fM = 1/Ms*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Ms,2)),(k+1)/2/(k-1))
                        for i in range(0,maxInt): 
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
                        As_Astar = 1/M*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(M,2)),(k+1)/2/(k-1))
                        Astar = As/As_Astar
                        Ae_Astar = Ae/Astar
                        Me_u = 1
                        Me_l = 0
                        Me = 1/2*(Me_u+Me_l)
                        Ae_Astar_fM = 1/Me*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Me,2)),(k+1)/2/(k-1))
                        for i in range(0,maxInt):
                            if Ae_Astar_fM > Ae_Astar:
                                Me_u = Me_u
                                Me_l = Me
                            else:
                                Me_u = Me
                                Me_l = Me_l
                            Me = 1/2*(Me_u+Me_l)
                            Ae_Astar_fM = 1/Me*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Me,2)),(k+1)/2/(k-1))
                        Pe_sonic_l = Pe_o / math.pow(1+(k-1)/2*math.pow(Me,2),k/(k-1))
                        if Pe > Pe_sonic_l and Pe < Pe_sonic_u :
                            sonicSolutionFoundFlag = 1
                        
                        if subsonicSolutionFoundFlag == 1 and sonicInputPin == 0:
                            sonicInput = 0
                            flagValid = 'valid'
                        elif sonicInputPin == 1 and sonicSolutionFoundFlag == 1:
                            sonicInput = 1
                            flagValid = 'valid'
                        elif subsonicSolutionFoundFlag == 1 and sonicInputPin == 0 and sonicSolutionFoundFlag == 1:
                            sonicInput = 0
                            flagValid = 'valid'
                        elif subsonicSolutionFoundFlag == 0 and sonicInputPin == 0 and sonicSolutionFoundFlag == 1:
                            sonicInput = 1
                            flagValid = 'valid'
                        else :
                            sonicInput = 0
                            flagValid = 'not valid'

                        if Ae < Ai :
                            flagValid = 'not valid'

                        if sonicInputPin == 0 and sonicInput == 1:
                            print('Finding solution but with sonic input assumption')
                        
                    case 'converging-diverging':
                        subsonicSolutionFoundFlag = 0
                        sonicSolutionFoundFlag = 0
                        flagValid = 'valid'

                        Pi = P_inlet
                        Pe = P_outlet
                        Ti = T_inlet
                        
                        # looking for subsonic solution
                        Astar_u = At
                        Astar_l = 0.01
                        # first finding upper bound where Astar = Ai
                        Ai_Astar = Ai/Astar_u
                        Mi_u = 1
                        Mi_l = 0
                        Mi = 1/2*(Mi_u+Mi_l)
                        Ai_Astar_fM = 1/Mi*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Mi,2)),(k+1)/2/(k-1))
                        for i in range(0,maxInt):
                            if Ai_Astar_fM > Ai_Astar:
                                Mi_u = Mi_u
                                Mi_l = Mi
                            else:
                                Mi_u = Mi
                                Mi_l = Mi_l
                            Mi = 1/2*(Mi_u+Mi_l)
                            Ai_Astar_fM = 1/Mi*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Mi,2)),(k+1)/2/(k-1))
                        Po_calc = Pi * math.pow(1+(k-1)/2*math.pow(Mi,2),k/(k-1))
                        Ae_Astar = Ae/Astar_u
                        Me_u = 1
                        Me_l = 0
                        Me = 1/2*(Me_l+Me_u)
                        Ae_Astar_fM = 1/Me*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Me,2)),(k+1)/2/(k-1))
                        for i in range(0,maxInt):
                            if Ae_Astar_fM > Ae_Astar:
                                Me_u = Me_u
                                Me_l = Me
                            else:
                                Me_u = Me
                                Me_l = Me_l
                            Me = 1/2*(Me_u+Me_l)
                            Ae_Astar_fM = 1/Me*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Me,2)),(k+1)/2/(k-1))
                        Pe_subsonic_l = Po_calc / math.pow(1+(k-1)/2*math.pow(Me,2),k/(k-1))
                        # next finding lower bound where Astar = 0.01
                        Ai_Astar = Ai/Astar_l
                        Mi_u = 1
                        Mi_l = 0
                        Mi = 1/2*(Mi_u+Mi_l)
                        Ai_Astar_fM = 1/Mi*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Mi,2)),(k+1)/2/(k-1))
                        for i in range(0,maxInt):
                            if Ai_Astar_fM > Ai_Astar:
                                Mi_u = Mi_u
                                Mi_l = Mi
                            else:
                                Mi_u = Mi
                                Mi_l = Mi_l
                            Mi = 1/2*(Mi_u+Mi_l)
                            Ai_Astar_fM = 1/Mi*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Mi,2)),(k+1)/2/(k-1))
                        Po_calc = Pi * math.pow(1+(k-1)/2*math.pow(Mi,2),k/(k-1))
                        Ae_Astar = Ae/Astar_l
                        Me_u = 1
                        Me_l = 0
                        Me = 1/2*(Me_l+Me_u)
                        Ae_Astar_fM = 1/Me*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Me,2)),(k+1)/2/(k-1))
                        for i in range(0,maxInt):
                            if Ae_Astar_fM > Ae_Astar:
                                Me_u = Me_u
                                Me_l = Me
                            else:
                                Me_u = Me
                                Me_l = Me_l
                            Me = 1/2*(Me_u+Me_l)
                            Ae_Astar_fM = 1/Me*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Me,2)),(k+1)/2/(k-1))
                        Pe_subsonic_u = Po_calc / math.pow(1+(k-1)/2*math.pow(Me,2),k/(k-1))
                        if Pe > Pe_subsonic_l and Pe < Pe_subsonic_u:
                            subsonicSolutionFoundFlag = 1
                        # print(Pe)
                        # print(Pe_subsonic_l)
                        # print(Pe_subsonic_u)
                        # print(subsonicSolutionFoundFlag)
                        # next check for sonic solution
                        # looking for sonic solution
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
                        # bounds for shock
                        As_u = Ae
                        As_l = At
                        # first finding upper bound on Pe
                        As = As_l # higher Pe when shock is closer to inlet
                        As_Astar = As/At
                        Ms_u = 100
                        Ms_l = 1
                        Ms = 1/2*(Ms_u+Ms_l)
                        As_Astar_fM = 1/Ms*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Ms,2)),(k+1)/2/(k-1))
                        for i in range(0,maxInt): 
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
                        As_Astar = 1/M*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(M,2)),(k+1)/2/(k-1))
                        Astar = As/As_Astar
                        Ae_Astar = Ae/Astar
                        Me_u = 1
                        Me_l = 0
                        Me = 1/2*(Me_u+Me_l)
                        Ae_Astar_fM = 1/Me*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Me,2)),(k+1)/2/(k-1))
                        for i in range(0,maxInt):
                            if Ae_Astar_fM > Ae_Astar:
                                Me_u = Me_u
                                Me_l = Me
                            else:
                                Me_u = Me
                                Me_l = Me_l
                            Me = 1/2*(Me_u+Me_l)
                            Ae_Astar_fM = 1/Me*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Me,2)),(k+1)/2/(k-1))
                        Pe_sonic_u = Pe_o / math.pow(1+(k-1)/2*math.pow(Me,2),k/(k-1))
                        # next finding lower bound on 
                        As = As_u # higher Pe when shock is closer to inlet
                        As_Astar = As/At
                        Ms_u = 100
                        Ms_l = 1
                        Ms = 1/2*(Ms_u+Ms_l)
                        As_Astar_fM = 1/Ms*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Ms,2)),(k+1)/2/(k-1))
                        for i in range(0,maxInt): 
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
                        As_Astar = 1/M*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(M,2)),(k+1)/2/(k-1))
                        Astar = As/As_Astar
                        Ae_Astar = Ae/Astar
                        Me_u = 1
                        Me_l = 0
                        Me = 1/2*(Me_u+Me_l)
                        Ae_Astar_fM = 1/Me*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Me,2)),(k+1)/2/(k-1))
                        for i in range(0,maxInt):
                            if Ae_Astar_fM > Ae_Astar:
                                Me_u = Me_u
                                Me_l = Me
                            else:
                                Me_u = Me
                                Me_l = Me_l
                            Me = 1/2*(Me_u+Me_l)
                            Ae_Astar_fM = 1/Me*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Me,2)),(k+1)/2/(k-1))
                        Pe_sonic_l = Pe_o / math.pow(1+(k-1)/2*math.pow(Me,2),k/(k-1))
                        # print(Pe_sonic_l)
                        # print(Pe_sonic_u)
                        if Pe > Pe_sonic_l and Pe < Pe_sonic_u :
                            sonicSolutionFoundFlag = 1
                        # print(sonicSolutionFoundFlag)
                        if subsonicSolutionFoundFlag == 1 and sonicInputPin == 0:
                            sonicInput = 0
                            flagValid = 'valid'
                        elif sonicInputPin == 1 and sonicSolutionFoundFlag == 1:
                            sonicInput = 1
                            flagValid = 'valid'
                        elif subsonicSolutionFoundFlag == 1 and sonicInputPin == 0 and sonicSolutionFoundFlag == 1:
                            sonicInput = 0
                            flagValid = 'valid'
                        elif subsonicSolutionFoundFlag == 0 and sonicInputPin == 0 and sonicSolutionFoundFlag == 1:
                            sonicInput = 1
                            flagValid = 'valid'
                        else :
                            sonicInput = 0
                            flagValid = 'not valid'

                        if Ae < At or Ai < At :
                            flagValid = 'not valid'

                        sonicInput = 0
                        if sonicInputPin == 0 and sonicInput == 1:
                            print('Finding solution but with sonic input assumption')

                    case _:
                        flag = 1
            
            case '2D':
                flag = 0
            
            case _:
                flag = 1

    return PreExamResults
