def plotSSResultFunction(SS_Solution,plottingRequests,IsentropicSolution,boundaryConditions):
    # Import
    import math
    import numpy as np
    import matplotlib.pyplot as plt

    P_track = SS_Solution.P_track
    rho_track = SS_Solution.rho_track
    T_track = SS_Solution.T_track
    u_track = SS_Solution.u_track
    x_vector = SS_Solution.x_vector
    h_vector = SS_Solution.h_vector
    int_track = SS_Solution.int_track
    SSStep_track = SS_Solution.SSStep_track
    time_track = SS_Solution.time_track
    cp = SS_Solution.cp
    # heat input
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

    # Isentropic solution
    P_isen = IsentropicSolution[0,:]
    # print(P_isen)
    T_isen = IsentropicSolution[1,:]
    rho_isen = IsentropicSolution[2,:]
    u_isen = IsentropicSolution[3,:]
    M_isen = IsentropicSolution[4,:]
    
    # Steady state result
    fig1, axs1 = plt.subplots(2,2)
    axs1[0,0].plot(x_vector,P_track[-1,:]/100000,label="CFD")
    axs1[0,0].plot(x_vector,P_isen/100000,label="Isentropic")
    axs1[0,0].set_xlabel("x (m)")
    axs1[0,0].set_ylabel("Pressure (bar)")
    fig1.legend()
    
    axs1[0,1].plot(x_vector,rho_track[-1,:]) 
    axs1[0,1].plot(x_vector,rho_isen) 
    axs1[0,1].set_xlabel("x (m)")
    axs1[0,1].set_ylabel("Density (kg/m3)")

    axs1[1,0].plot(x_vector,T_track[-1,:]-273.15) 
    axs1[1,0].plot(x_vector,T_isen-273.15) 
    axs1[1,0].set_xlabel("x (m)")
    axs1[1,0].set_ylabel("Temperature (deg C)")

    axs1[1,1].plot(x_vector,u_track[-1,:]) 
    axs1[1,1].plot(x_vector,u_isen) 
    axs1[1,1].set_xlabel("x (m)")
    axs1[1,1].set_ylabel("Velocity (m/s)")

    # Mach number plot
    fig15, axs15 = plt.subplots()
    axs15.plot(x_vector,u_track[-1,:]/np.sqrt(SS_Solution.k*SS_Solution.R*T_track[-1,:]))
    axs15.set_xlabel("x (m)")
    axs15.set_ylabel("Mach Number")

    # Iterations within step
    fig2, axs2 = plt.subplots()
    axs2.plot(SSStep_track,int_track)
    axs2.set_xlabel("Steady step number")
    axs2.set_ylabel("Iterations at this step")

    # Time based results
    fig3, axs3 = plt.subplots(3,2)
    # fig3, axs3 = plt.subplots(3,2,figsize = (4,2))
    axs3[0,0].plot(time_track,P_track[:,0]/100000,label="x0")
    axs3[0,0].plot(time_track,P_track[:,1]/100000,label="x1")
    axs3[0,0].plot(time_track,P_track[:,-2]/100000,label="xN-1")
    axs3[0,0].plot(time_track,P_track[:,-1]/100000,label="xN")
    axs3[0,0].set_xlabel("Time (s)")
    axs3[0,0].set_ylabel("Pressure (bar)")
    # axs3[0,0].set_xticks([0,0.005,0.010,0.015])
    axs3[0,0].set_xticks([0,0.005,0.010])
    # axs3[0,0].set_xticks([0,0.01,0.02,0.03,0.04,0.05,0.06])
    fig3.legend()

    axs3[0,1].plot(time_track,rho_track[:,0])
    axs3[0,1].plot(time_track,rho_track[:,1])
    axs3[0,1].plot(time_track,rho_track[:,-2])
    axs3[0,1].plot(time_track,rho_track[:,-1])
    axs3[0,1].set_xlabel("Time (s)")
    axs3[0,1].set_ylabel("Density (kg/m3)")
    # axs3[0,1].set_xticks([0,0.005,0.010,0.015])
    axs3[0,1].set_xticks([0,0.005,0.010])
    # axs3[0,1].set_xticks([0,0.01,0.02,0.03,0.04,0.05,0.06])

    axs3[1,0].plot(time_track,T_track[:,0]-273.15)
    axs3[1,0].plot(time_track,T_track[:,1]-273.15)
    axs3[1,0].plot(time_track,T_track[:,-2]-273.15)
    axs3[1,0].plot(time_track,T_track[:,-1]-273.15)
    axs3[1,0].set_xlabel("Time (s)")
    axs3[1,0].set_ylabel("Temperature (K)")
    # axs3[1,0].set_xticks([0,0.005,0.010,0.015])
    axs3[1,0].set_xticks([0,0.005,0.010])
    # axs3[1,0].set_xticks([0,0.01,0.02,0.03,0.04,0.05,0.06])

    axs3[1,1].plot(time_track,u_track[:,0])
    axs3[1,1].plot(time_track,u_track[:,1])
    axs3[1,1].plot(time_track,u_track[:,-2])
    axs3[1,1].plot(time_track,u_track[:,-1])
    axs3[1,1].set_xlabel("Time (s)")
    axs3[1,1].set_ylabel("Velocity (m/s)")
    # axs3[1,1].set_xticks([0,0.005,0.010,0.015])
    axs3[1,1].set_xticks([0,0.005,0.010])
    # axs3[1,1].set_xticks([0,0.01,0.02,0.03,0.04,0.05,0.06])

    axs3[2,0].plot(time_track,cp*T_track[:,0]+1/2*np.power(u_track[:,0],2))
    axs3[2,0].plot(time_track,cp*T_track[:,1]+1/2*np.power(u_track[:,1],2))
    axs3[2,0].plot(time_track,cp*T_track[:,-2]+1/2*np.power(u_track[:,-2],2))
    axs3[2,0].plot(time_track,cp*T_track[:,-1]+1/2*np.power(u_track[:,-1],2))
    axs3[2,0].set_xlabel("Time (s)")
    axs3[2,0].set_ylabel("Enthalpy (J/kg)")
    # axs3[2,0].set_xticks([0,0.005,0.010,0.015])
    axs3[2,0].set_xticks([0,0.005,0.010])
    # axs3[2,0].set_xticks([0,0.01,0.02,0.03,0.04,0.05,0.06])

    axs3[2,1].plot(time_track,u_track[:,0]/np.sqrt(SS_Solution.k*SS_Solution.R*T_track[:,0]))
    axs3[2,1].plot(time_track,u_track[:,1]/np.sqrt(SS_Solution.k*SS_Solution.R*T_track[:,1]))
    axs3[2,1].plot(time_track,u_track[:,-2]/np.sqrt(SS_Solution.k*SS_Solution.R*T_track[:,-2]))
    axs3[2,1].plot(time_track,u_track[:,-1]/np.sqrt(SS_Solution.k*SS_Solution.R*T_track[:,-1]))
    axs3[2,1].set_xlabel("Time (s)")
    axs3[2,1].set_ylabel("Mach Number")
    # axs3[2,1].set_xticks([0,0.005,0.010,0.015])
    axs3[2,1].set_xticks([0,0.005,0.010])
    # axs3[2,1].set_xticks([0,0.01,0.02,0.03,0.04,0.05,0.06])

    # troubleshooting plot
    # axs3[2,1].plot(time_track,2*P_track[:,-1]/100000-P_track[:,-2]/100000)

    # Geometry
    fig4, axs4 = plt.subplots()
    axs4.plot(SS_Solution.x_vector,SS_Solution.h_vector)
    axs4.set_xlabel("x (m)")
    axs4.set_ylabel("y (m)")

    # Heat input
    fig5,axs5 = plt.subplots()
    axs5.plot(x_vector,q_array)
    axs5.set_xlabel("x (m)")
    axs5.set_ylabel("q (W/m2)")

    # Geometry 2
    h_u_vector = np.linspace(1,1,len(x_vector))
    h_d_vector = np.linspace(1,1,len(x_vector))
    temp_counter = 0
    for i in h_u_vector:
        h_max = max(h_vector)
        h_u_vector[temp_counter] = h_max - (h_max-h_vector[temp_counter])/2
        h_d_vector[temp_counter] = (h_max-h_vector[temp_counter])/2
        temp_counter = temp_counter + 1

    fig6,axs6 = plt.subplots()
    axs6.plot(x_vector,h_u_vector,color='blue')
    axs6.plot(x_vector,h_d_vector,color='blue')
    axs6.set_xlabel("x (m)")
    axs6.set_ylabel("y (m)")

    plt.show()