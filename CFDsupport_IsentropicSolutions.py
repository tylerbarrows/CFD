# import
import math
import numpy as np
import matplotlib.pyplot as plt

# fluid properties
k = 1.4
R = 287

# geometry
shape = 'converging-diverging'
# shape = 'diverging'
# shape = 'converging'
# Ai = 0.25
Ai = 0.2
At = 0.1
Ae = 0.25

# inlet conditions
Pi = 1.2 # bar
Ti = 300 # K

# convergence error
error = 0.1/100
maxInt = 50

match shape:

    case 'converging':
        ########## Diverging ##########
        # really only want to check what happens when input is subsonic vs sonic
        # 2025/02/09 - don't know how to handle M=1 at converging inlet
        # possible to handle supersonic case but would have minimum supersonic speed I think
        # will have to expand to handle supersonic input eventually 

        # to fix CFD_IsentropicSolution.py - don't love the way we find converging isentropic solution
        N = 500
        Astar_range = np.linspace(0.01,Ae,N)
        Mi_range = np.linspace(-1,-1,N)
        Po_range = np.linspace(-1,-1,N)
        Pe_range = np.linspace(-1,-1,N)
        Me_range = np.linspace(-1,-1,N)
        index = 0
        for Astar in Astar_range:
            Ai_Astar = Ai/Astar
            Mi_u = 1
            Mi_l = 0
            Mi = 1/2*(Mi_l+Mi_u)
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
            Po = Pi / math.pow(1+(k-1)/2*math.pow(Mi,2),k/(k-1))
            # finding upper Pe (Pe needs to be this low for shock)
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
            Pe = Po / math.pow(1+(k-1)/2*math.pow(Me,2),k/(k-1))

            Mi_range[index] = Mi
            Po_range[index] = Po
            Pe_range[index] = Pe
            Me_range[index] = Me

            index = index + 1
        
        fig1, axs1 = plt.subplots()
        axs1.scatter(Mi_range,Pe_range,s=10)
        axs1.set_xlabel("Inlet mach #")
        axs1.set_ylabel("Pe (bar)")

        fig2, axs2 = plt.subplots()
        axs2.scatter(Astar_range,Pe_range,s=10)
        axs2.set_xlabel("Astar")
        axs2.set_ylabel("Pe (bar)")

        fig3, axs3 = plt.subplots()
        axs3.scatter(Mi_range,Me_range,s=10)
        axs3.set_xlabel("Inlet mach #")
        axs3.set_ylabel("Outlet mach #")

        plt.show()

    case 'diverging':
        ########## Diverging ##########
        # assuming choked throat, what are critical Pe
        u_choked = math.sqrt(k*R*Ti)
        M_choked = 1
        Po = Pi * math.pow(1+(k-1)/2*math.pow(M_choked,2),k/(k-1))
        To = Ti * math.pow(Po/Pi,(k-1)/k)
        # finding upper Pe (Pe needs to be this low for shock)
        Ae_Astar = Ae/Ai
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
        Me_critical_u = Me
        Pe_critical_u = Po / math.pow(1+(k-1)/2*math.pow(Me,2),(k/(k-1)))
        Te_critical_u = To * math.pow(Pe_critical_u/Po,(k-1)/k)
        # finding upper Pe (Pe needs to be this low for shock)
        Ae_Astar = Ae/Ai
        Me_u = 100
        Me_l = 1
        Me = 1/2*(Me_u+Me_l)
        Ae_Astar_fM = 1/Me*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Me,2)),(k+1)/2/(k-1))
        for i in range(0,maxInt):
            if Ae_Astar_fM > Ae_Astar:
                Me_u = Me
                Me_l = Me_l
            else:
                Me_u = Me_u
                Me_l = Me
            Me = 1/2*(Me_u+Me_l)
            Ae_Astar_fM = 1/Me*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Me,2)),(k+1)/2/(k-1))
        Me_critical_l = Me
        Pe_critical_l = Po / math.pow(1+(k-1)/2*math.pow(Me,2),(k/(k-1)))
        Te_critical_l = To * math.pow(Pe_critical_l/Po,(k-1)/k)
        rhoe_critical_l = Pe_critical_l * 100000 / Te_critical_l / R
        ue_critical_l = Me_critical_l * math.sqrt(k*R*Te_critical_l)
        c1 = rhoe_critical_l*ue_critical_l
        c2 = Pe_critical_l*100000+rhoe_critical_l*math.pow(ue_critical_l,2)
        c3 = k/(k-1)*Pe_critical_l*100000/rhoe_critical_l+1/2*math.pow(ue_critical_l,2)
        if math.pow(k/(k-1)*c2,2)-4*c3*1/2*(k+1)/(k-1)*math.pow(c1,2) < 0:
            rhoe_critical_l_new =  (k/(k-1)*c2+0)/2/c3
        else:
            rhoe_critical_l_new = (k/(k-1)*c2+math.sqrt(math.pow(k/(k-1)*c2,2)-4*c3*1/2*(k+1)/(k-1)*math.pow(c1,2)))/2/c3
        ue_critical_l_new = c1/rhoe_critical_l_new
        Pe_critical_l_new = (c2-rhoe_critical_l_new*math.pow(ue_critical_l_new,2))/100000
        Te_critical_l_new = Pe_critical_l_new*100000/R/rhoe_critical_l_new
        Me_critical_l_new = ue_critical_l_new/math.sqrt(k*R*Te_critical_l_new)

        # let's find relationship between Astar and Pe
        Astar_range = np.linspace(0.01,Ai,100)
        Pe_Astar_range = np.linspace(1,1,100)
        Po_Astar_range = np.linspace(1,1,100)
        index = 0
        for Astar in Astar_range:
            Ai_Astar = Ai/Astar
            Mi_u = 1
            Mi_l = 0
            Mi = 1/2*(Mi_l+Mi_u)
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
            Ae_Astar = Ae/Astar
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
            Pe_calc = Po_calc / math.pow(1+(k-1)/2*math.pow(Me,2),k/(k-1))
            Po_Astar_range[index] = Po_calc
            Pe_Astar_range[index] = Pe_calc
            index = index + 1
        
        # let's check what the As to Pe relationship is
        As_range_sonic_sweep = np.linspace(Ai,Ae,100)
        Pe_range_sonic_sweep = np.linspace(1,1,100)
        index = 0
        for As in As_range_sonic_sweep:
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
            Pe_calc = Pe_o / math.pow(1+(k-1)/2*math.pow(Me,2),k/(k-1))
            Pe_range_sonic_sweep[index] = Pe_calc
            index = index + 1
        
        # finding number of solutions
        N = 500
        Pe_range = np.linspace(Pe_critical_l,Pe_critical_u,N)
        # Pe_range = np.append(Pe_range,np.linspace(Pe_critical_u+(Pe_critical_u-Pe_critical_l)/100,Po,100))
        Pe_range = np.append(Pe_range,np.linspace(Pe_critical_u,Po,2*N))
        Pi_range = np.linspace(Pi,Pi,len(Pe_range))
        As_range_sonic = np.linspace(0,0,len(Pe_range))
        Po_track_sonic = np.linspace(0,0,len(Pe_range))
        Mi_track_sonic = np.linspace(0,0,len(Pe_range))
        error_track = np.linspace(0,0,len(Pe_range))
        error_track_subsonic = np.linspace(0,0,len(Pe_range))
        Po_track_subsonic = np.linspace(0,0,len(Pe_range))
        Mi_track_subsonic = np.linspace(0,0,len(Pe_range))
        Me_track_subsonic = np.linspace(0,0,len(Pe_range))
        index = 0
        for Pe in Pe_range:
            # first look for subsonic solution
            Astar_u = Ai
            Astar_l = 0
            Astar = 1/2*(Astar_u+Astar_l)
            Ai_Astar = Ai/Astar
            Mi_u = 1
            Mi_l = 0
            Mi = 1/2*(Mi_l+Mi_u)
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
            Ae_Astar = Ae/Astar
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
            Pe_calc = Po_calc / math.pow(1+(k-1)/2*math.pow(Me,2),k/(k-1))
            for i in range(0,maxInt):
                if Pe_calc > Pe:
                    Astar_u = Astar
                    Astar_l = Astar_l
                else:
                    Astar_u = Astar_u
                    Astar_l = Astar
                Astar = 1/2*(Astar_u+Astar_l)
                Ai_Astar = Ai/Astar
                Mi_u = 1
                Mi_l = 0
                Mi = 1/2*(Mi_l+Mi_u)
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
                Ae_Astar = Ae/Astar
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
                Pe_calc = Po_calc / math.pow(1+(k-1)/2*math.pow(Me,2),k/(k-1))
            if abs(Pe_calc-Pe)/Pe > error:
                Po_track_subsonic[index] = -1
                Mi_track_subsonic[index] = -1
                Me_track_subsonic[index] = -1
            else:
                Po_track_subsonic[index] = Po_calc
                Mi_track_subsonic[index] = Mi
                Me_track_subsonic[index] = Me
            
            error_track_subsonic[index] = abs(Pe_calc-Pe)/Pe

            # next look for sonic solution
            Astar = Ai
            # find relationship again between As and Pe
            # upper and lower As
            # find solution upstream and downstream of shock
            # expand to exit and find outlet condition
            # update As_u and As_l and iterate 
            # check error if solution is found or not
            Astar = Ai # need choked flow at throat
            As_u = Ae
            As_l = Ai
            for l in range(0,maxInt):
                Astar = Ai
                As = 1/2*(As_u+As_l)
                As_Astar = As/Astar
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
                # print(abs(As_Astar_fM-As_Astar)/As_Astar)
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
                Pe_calc = Pe_o / math.pow(1+(k-1)/2*math.pow(Me,2),k/(k-1))
                if Pe_calc > Pe:
                    As_u = As_u
                    As_l = As
                    # As_u = As
                    # As_l = As_l
                else:
                    As_u = As
                    As_l = As_l
                    # As_u = As_u
                    # As_l = As
            error_track[index] = abs(Pe_calc-Pe)/Pe
            if abs(Pe_calc-Pe)/Pe > error:
                Mi_track_sonic[index] = -1
                As_range_sonic[index] = -1
                Po_track_sonic[index] = -1
                Mi_track_sonic[index] = -1
            else:
                Mi_track_sonic[index] = 1
                As_range_sonic[index] = As
                Po_track_sonic[index] = Po
                Mi_track_sonic[index] = math.pow(2/(k-1)*(math.pow(Po/Pi,(k-1)/k)-1),1/2)
            
            index = index + 1

        # print(Po)
        # print(Pe_critical_l)
        # print(Pe_critical_l_new)
        # print(Pe_critical_u)

        ## Plotting
        fig1, axs1 = plt.subplots()
        axs1.scatter(Pe_range,Mi_track_sonic,s=10,label="Sonic")
        axs1.scatter(Pe_range,Mi_track_subsonic,s=10,label="Subsonic")
        axs1.set_xlabel("Pe (bar)")
        axs1.set_ylabel("Mach #")
        fig1.legend()
        plt.ylim(0,1)
        plt.title('Mt vs Pe (sonic and subsonic)')

        fig2, axs2 = plt.subplots()
        axs2.scatter(Pe_range,As_range_sonic,s=10,label="Shock area")
        axs2.set_xlabel("Pe (bar)")
        axs2.set_ylabel("Shock area")
        plt.ylim(Ai,Ae)
        plt.title('Pe vs shock area (sonic)')

        fig3, axs3 = plt.subplots()
        axs3.scatter(Pe_range,error_track,s=10,label="Sonic error")
        axs3.set_xlabel("Pe (bar)")
        axs3.set_ylabel("Error")
        plt.title('Pe vs error (sonic)')

        fig4, axs4 = plt.subplots()
        axs4.scatter(As_range_sonic_sweep,Pe_range_sonic_sweep,s=10,label="Pe")
        axs4.set_xlabel("As")
        axs4.set_ylabel("Pe")
        plt.title('As vs Pe (sonic sweep)')

        fig5, axs5 = plt.subplots()
        axs5.scatter(Astar_range,Pe_Astar_range,s=10,label="Pe")
        axs5.set_xlabel("Astar")
        axs5.set_ylabel("Pe")
        plt.title('Pe vs Astar (subsonic)')

        fig6, axs6 = plt.subplots()
        axs6.scatter(Pe_range,error_track_subsonic,s=10,label="Pe")
        axs6.set_xlabel("Pe")
        axs6.set_ylabel("Error")
        plt.title('Pe vs error (subsonic)')

        Pi_range_check = Po_track_subsonic / np.power(1+(k-1)/2*np.power(Mi_track_subsonic,2),k/(k-1))
        Pe_range_check = Po_track_subsonic / np.power(1+(k-1)/2*np.power(Me_track_subsonic,2),k/(k-1))
        for i in range(0,len(Pi_range_check)):
            if Mi_track_subsonic[i] == -1:
                Pi_range_check[i] = -1
                Pe_range_check[i] = -1
        fig7, axs7 = plt.subplots()
        axs7.scatter(Pe_range,Pi_range_check,s=10,label="Pi")
        axs7.scatter(Pe_range,Pe_range_check,s=10,label="Pe")
        axs7.set_xlabel("Pe")
        axs7.set_ylabel("Pi and Pe")
        plt.title('Pi and Pe check (sonic)')
        fig7.legend()
        plt.ylim(0,2*Pi)
        # Nbuffer = 20
        # plt.xlim(Pe_range[min(np.argwhere(Pe_range_check>0))[0]-Nbuffer],Pe_range[max(np.argwhere(Pe_range_check>0))[0]+Nbuffer])
        # available vectors
        # Po_track_subsonic
        # Mi_track_subonic
        # Me_track_subsonic
        # Pe_range

        fig8, axs8 = plt.subplots()
        axs8.scatter(Pi_range/Pe_range,Mi_track_sonic,s=10,label="Sonic")
        axs8.scatter(Pi_range/Pe_range,Mi_track_subsonic,s=10,label="Subsonic")
        axs8.set_xlabel("Pi/Pe")
        axs8.set_ylabel("Valve inlet Mach #")
        fig8.legend()
        plt.ylim(0,1)
        plt.xlim(0,2)
        plt.title('Inlet Mach # vs valve pressure ratio (sonic and subsonic)')

        fig9, axs9 = plt.subplots()
        axs9.scatter(Po_track_sonic/Pe_range,Mi_track_sonic,s=10,label="Sonic")
        axs9.scatter(Po_track_subsonic/Pe_range,Mi_track_subsonic,s=10,label="Subsonic")
        axs9.set_xlabel("Po/Pe")
        axs9.set_ylabel("Valve inlet Mach #")
        fig9.legend()
        plt.ylim(0,1)
        plt.xlim(1,2)
        plt.title('Inlet Mach # vs valve pressure ratio (sonic and subsonic)')

        # print(Pi/Pe_range)
        fig10, axs10 = plt.subplots()
        axs10.scatter(Pi/Pe_range,Po_track_sonic/Pe_range,s=10,label="Sonic")
        axs10.scatter(Pi/Pe_range,Po_track_subsonic/Pe_range,s=10,label="Subsonic")
        axs10.set_xlabel("Pi/Pe")
        axs10.set_ylabel("Po/Pe")
        fig10.legend()
        plt.ylim(1,2)
        plt.xlim(0,2)
        plt.title('Inlet Mach # vs valve pressure ratio (sonic and subsonic)')

        plt.show()

    case 'converging-diverging':
        ########## Converging-Diverging ##########
        # assuming choked throat, what are critical Pe
        # first find Po assuming choked throat
        Ai_Astar = Ai/At
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
        Po = Pi * math.pow(1+(k-1)/2*math.pow(Mi,2),(k/(k-1)))
        To = Ti * math.pow(Po/Pi,(k-1)/k)
        # finding upper Pe (Pe needs to be this low for shock)
        Ae_Astar = Ae/At
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
        Me_critical_u = Me
        Pe_critical_u = Po / math.pow(1+(k-1)/2*math.pow(Me,2),(k/(k-1)))
        Te_critical_u = To * math.pow(Pe_critical_u/Po,(k-1)/k)
        # finding upper Pe (Pe needs to be this low for shock)
        Ae_Astar = Ae/At
        Me_u = 100
        Me_l = 1
        Me = 1/2*(Me_u+Me_l)
        Ae_Astar_fM = 1/Me*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Me,2)),(k+1)/2/(k-1))
        for i in range(0,maxInt):
            if Ae_Astar_fM > Ae_Astar:
                Me_u = Me
                Me_l = Me_l
            else:
                Me_u = Me_u
                Me_l = Me
            Me = 1/2*(Me_u+Me_l)
            Ae_Astar_fM = 1/Me*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Me,2)),(k+1)/2/(k-1))
        Me_critical_l = Me
        Pe_critical_l = Po / math.pow(1+(k-1)/2*math.pow(Me,2),(k/(k-1)))
        Te_critical_l = To * math.pow(Pe_critical_l/Po,(k-1)/k)
        rhoe_critical_l = Pe_critical_l * 100000 / Te_critical_l / R
        ue_critical_l = Me_critical_l * math.sqrt(k*R*Te_critical_l)
        c1 = rhoe_critical_l*ue_critical_l
        c2 = Pe_critical_l*100000+rhoe_critical_l*math.pow(ue_critical_l,2)
        c3 = k/(k-1)*Pe_critical_l*100000/rhoe_critical_l+1/2*math.pow(ue_critical_l,2)
        if math.pow(k/(k-1)*c2,2)-4*c3*1/2*(k+1)/(k-1)*math.pow(c1,2) < 0:
            rhoe_critical_l_new =  (k/(k-1)*c2+0)/2/c3
        else:
            rhoe_critical_l_new = (k/(k-1)*c2+math.sqrt(math.pow(k/(k-1)*c2,2)-4*c3*1/2*(k+1)/(k-1)*math.pow(c1,2)))/2/c3
        ue_critical_l_new = c1/rhoe_critical_l_new
        Pe_critical_l_new = (c2-rhoe_critical_l_new*math.pow(ue_critical_l_new,2))/100000
        Te_critical_l_new = Pe_critical_l_new*100000/R/rhoe_critical_l_new
        Me_critical_l_new = ue_critical_l_new/math.sqrt(k*R*Te_critical_l_new)

        # let's check what the As to Pe relationship is
        As_range_sonic_sweep = np.linspace(At,Ae,100)
        Pe_range_sonic_sweep = np.linspace(1,1,100)
        index = 0
        for As in As_range_sonic_sweep:
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
            # print(P/Ps)
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
            Pe_calc = Pe_o / math.pow(1+(k-1)/2*math.pow(Me,2),k/(k-1))
            Pe_range_sonic_sweep[index] = Pe_calc
            index = index + 1

        # let's find relation between At_Astar and Pe
        Astar_range = np.linspace(0.01,At,100)
        Pe_Astar_range = np.linspace(1,1,100)
        Po_Astar_range = np.linspace(1,1,100)
        index = 0
        for Astar in Astar_range:
            Ai_Astar = Ai/Astar
            Mi_u = 1
            Mi_l = 0
            Mi = 1/2*(Mi_l+Mi_u)
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
            # print(Po_calc)
            Ae_Astar = Ae/Astar
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
            Pe_calc = Po_calc / math.pow(1+(k-1)/2*math.pow(Me,2),k/(k-1))
            Po_Astar_range[index] = Po_calc
            Pe_Astar_range[index] = Pe_calc
            index = index + 1

        # finding number of solutions
        N = 500
        Pe_range = np.linspace(Pe_critical_l,Pe_critical_u,N)
        # Pe_range = np.append(Pe_range,np.linspace(Pe_critical_u+(Pe_critical_u-Pe_critical_l)/100,Po,100))
        Pe_range = np.append(Pe_range,np.linspace(Pe_critical_u,Po,2*N))
        Pi_range = np.linspace(Pi,Pi,len(Pe_range))
        Mt_range_subsonic = np.linspace(0,0,len(Pe_range))
        Mt_range_sonic = np.linspace(0,0,len(Pe_range))
        As_range_sonic = np.linspace(0,0,len(Pe_range))
        Po_track_sonic = np.linspace(0,0,len(Pe_range))
        Mi_track_sonic = np.linspace(0,0,len(Pe_range))
        error_track = np.linspace(0,0,len(Pe_range))
        error_track_subsonic = np.linspace(0,0,len(Pe_range))
        Po_track_subsonic = np.linspace(0,0,len(Pe_range))
        Mi_track_subsonic = np.linspace(0,0,len(Pe_range))
        Me_track_subsonic = np.linspace(0,0,len(Pe_range))
        index = 0
        for Pe in Pe_range:
            # first look for subsonic solution
            # Astar bounds
            # Ai/Astar gives Mi
            # with Mi and Pi, get Po
            # Ae/Astar gives Me
            # Po and Me give Pe
            # Adjust Astar based on Pe
            Astar_u = At
            Astar_l = 0
            Astar = 1/2*(Astar_u+Astar_l)
            Ai_Astar = Ai/Astar
            Mi_u = 1
            Mi_l = 0
            Mi = 1/2*(Mi_l+Mi_u)
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
            Ae_Astar = Ae/Astar
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
            Pe_calc = Po_calc / math.pow(1+(k-1)/2*math.pow(Me,2),k/(k-1))
            for i in range(0,maxInt):
                if Pe_calc > Pe:
                    Astar_u = Astar
                    Astar_l = Astar_l
                else:
                    Astar_u = Astar_u
                    Astar_l = Astar
                Astar = 1/2*(Astar_u+Astar_l)
                Ai_Astar = Ai/Astar
                Mi_u = 1
                Mi_l = 0
                Mi = 1/2*(Mi_l+Mi_u)
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
                Ae_Astar = Ae/Astar
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
                Pe_calc = Po_calc / math.pow(1+(k-1)/2*math.pow(Me,2),k/(k-1))
            if abs(Pe_calc-Pe)/Pe > error:
                Mt_range_subsonic[index] = -1
                Po_track_subsonic[index] = -1
                Mi_track_subsonic[index] = -1
                Me_track_subsonic[index] = -1
            else:
                At_Astar = At/Astar
                Mt_u = 1
                Mt_l = 0
                Mt = 1/2*(Mt_l+Mt_u)
                At_Astar_fM = 1/Mt*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Mt,2)),(k+1)/2/(k-1))
                for i in range(0,maxInt):
                    if At_Astar_fM > At_Astar:
                        Mt_u = Mt_u
                        Mt_l = Mt
                    else:
                        Mt_u = Mt
                        Mt_l = Mt_l
                    Mt = 1/2*(Mt_u+Mt_l)
                    At_Astar_fM = 1/Mt*math.pow(2/(k+1)*(1+(k-1)/2*math.pow(Mt,2)),(k+1)/2/(k-1))
                Mt_range_subsonic[index] = Mt
                Po_track_subsonic[index] = Po_calc
                Mi_track_subsonic[index] = Mi
                Me_track_subsonic[index] = Me
            
            error_track_subsonic[index] = abs(Pe_calc-Pe)/Pe

            # now look for sonic solution
            Astar = At # need choked flow at throat
            As_u = Ae
            As_l = At
            for l in range(0,maxInt):
                Astar = At
                As = 1/2*(As_u+As_l)
                As_Astar = As/Astar
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
                # print(abs(As_Astar_fM-As_Astar)/As_Astar)
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
                Pe_calc = Pe_o / math.pow(1+(k-1)/2*math.pow(Me,2),k/(k-1))
                if Pe_calc > Pe:
                    As_u = As_u
                    As_l = As
                    # As_u = As
                    # As_l = As_l
                else:
                    As_u = As
                    As_l = As_l
                    # As_u = As_u
                    # As_l = As
            error_track[index] = abs(Pe_calc-Pe)/Pe
            if abs(Pe_calc-Pe)/Pe > error:
                Mt_range_sonic[index] = -1
                As_range_sonic[index] = -1
                Po_track_sonic[index] = -1
                Mi_track_sonic[index] = -1
            else:
                Mt_range_sonic[index] = 1
                As_range_sonic[index] = As
                Po_track_sonic[index] = Po
                Mi_track_sonic[index] = math.pow(2/(k-1)*(math.pow(Po/Pi,(k-1)/k)-1),1/2)
            index = index + 1

        # print(Po)
        # print(Pe_critical_l)
        # print(Pe_critical_l_new)
        # print(Pe_critical_u)

        ## Plotting
        fig1, axs1 = plt.subplots()
        axs1.scatter(Pe_range,Mt_range_sonic,s=10,label="Sonic")
        axs1.scatter(Pe_range,Mt_range_subsonic,s=10,label="Subsonic")
        axs1.set_xlabel("Pe (bar)")
        axs1.set_ylabel("Mach #")
        fig1.legend()
        plt.ylim(0,1)
        plt.title('Mt vs Pe (sonic and subsonic)')

        fig2, axs2 = plt.subplots()
        axs2.scatter(Pe_range,As_range_sonic,s=10,label="Shock area")
        axs2.set_xlabel("Pe (bar)")
        axs2.set_ylabel("Shock area")
        plt.ylim(At,Ae)
        plt.title('Pe vs shock area (sonic)')

        fig3, axs3 = plt.subplots()
        axs3.scatter(Pe_range,error_track,s=10,label="Sonic error")
        axs3.set_xlabel("Pe (bar)")
        axs3.set_ylabel("Error")
        plt.title('Pe vs error (sonic)')

        fig4, axs4 = plt.subplots()
        axs4.scatter(As_range_sonic_sweep,Pe_range_sonic_sweep,s=10,label="Pe")
        axs4.set_xlabel("As")
        axs4.set_ylabel("Pe")
        plt.title('As vs Pe (sonic sweep)')

        fig5, axs5 = plt.subplots()
        axs5.scatter(Astar_range,Pe_Astar_range,s=10,label="Pe")
        axs5.set_xlabel("Astar")
        axs5.set_ylabel("Pe")
        plt.title('Pe vs Astar (subsonic)')

        fig6, axs6 = plt.subplots()
        axs6.scatter(Pe_range,error_track_subsonic,s=10,label="Pe")
        axs6.set_xlabel("Pe")
        axs6.set_ylabel("Error")
        plt.title('Pe vs error (subsonic)')

        Pi_range_check = Po_track_subsonic / np.power(1+(k-1)/2*np.power(Mi_track_subsonic,2),k/(k-1))
        Pe_range_check = Po_track_subsonic / np.power(1+(k-1)/2*np.power(Me_track_subsonic,2),k/(k-1))
        for i in range(0,len(Pi_range_check)):
            if Mt_range_subsonic[i] == -1:
                Pi_range_check[i] = -1
                Pe_range_check[i] = -1
        fig7, axs7 = plt.subplots()
        axs7.scatter(Pe_range,Pi_range_check,s=10,label="Pi")
        axs7.scatter(Pe_range,Pe_range_check,s=10,label="Pe")
        axs7.set_xlabel("Pe")
        axs7.set_ylabel("Pi and Pe")
        plt.title('Pi and Pe check (sonic)')
        fig7.legend()
        plt.ylim(0,2*Pi)
        # Nbuffer = 20
        # plt.xlim(Pe_range[min(np.argwhere(Pe_range_check>0))[0]-Nbuffer],Pe_range[max(np.argwhere(Pe_range_check>0))[0]+Nbuffer])
        # available vectors
        # Po_track_subsonic
        # Mi_track_subonic
        # Me_track_subsonic
        # Pe_range

        fig8, axs8 = plt.subplots()
        axs8.scatter(Pi_range/Pe_range,Mi_track_sonic,s=10,label="Sonic")
        axs8.scatter(Pi_range/Pe_range,Mi_track_subsonic,s=10,label="Subsonic")
        axs8.set_xlabel("Pi/Pe")
        axs8.set_ylabel("Valve inlet Mach #")
        fig8.legend()
        plt.ylim(0,1)
        plt.xlim(0,2)
        plt.title('Inlet Mach # vs valve pressure ratio (sonic and subsonic)')

        fig9, axs9 = plt.subplots()
        axs9.scatter(Po_track_sonic/Pe_range,Mi_track_sonic,s=10,label="Sonic")
        axs9.scatter(Po_track_subsonic/Pe_range,Mi_track_subsonic,s=10,label="Subsonic")
        axs9.set_xlabel("Po/Pe")
        axs9.set_ylabel("Valve inlet Mach #")
        fig9.legend()
        plt.ylim(0,1)
        plt.xlim(1,2)
        plt.title('Inlet Mach # vs valve pressure ratio (sonic and subsonic)')

        # print(Pi/Pe_range)
        fig10, axs10 = plt.subplots()
        axs10.scatter(Pi/Pe_range,Po_track_sonic/Pe_range,s=10,label="Sonic")
        axs10.scatter(Pi/Pe_range,Po_track_subsonic/Pe_range,s=10,label="Subsonic")
        axs10.set_xlabel("Pi/Pe")
        axs10.set_ylabel("Po/Pe")
        fig10.legend()
        plt.ylim(1,2)
        plt.xlim(0,2)
        plt.title('Inlet Mach # vs valve pressure ratio (sonic and subsonic)')

        plt.show()

    case _:
        print('No acceptable shape entered')