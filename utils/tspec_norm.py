def calc_Tspec_norm(N_temps, T_low, T_high, norms, mode, tlscp):

    temps = np.linspace(T_low, T_high, N_temps)
    temps1 = temps[:int(N_temps/2)]
    temps2 = temps[int(N_temps/2):]

    #print(temps)
    #print(temps1)
    #print(temps2)

    temps_spec = {}
    temps_diff = {}

    dts = {}

    for normm in norms:

        temps_spec[str(normm)] = []
        temps_diff[str(normm)] = []

        tt = []

        dts[str(normm)] = []
        dd = []
        
        #temps_prom = np.zeros(N_temps)
   	#N_usr = 2
    	#for i in range(0, N_usr):
    	#    print(i, end=" ")

        for t1, t2 in zip(temps1, temps2):
            
            abundanc = 1.0
            if mode == 'lines':
                temp_test, flux_test, dt_left_test, dt_right_test = single_T_fit_lines(t2, t1, 2, abundanc, 
                                                                                       tlscp, normm, 
                                                                                       texp=40000, stpar=True, 
                                                                                       plot=False, Xplot=False)
            elif mode == 'continuum':
                temp_test, flux_test, dt_left_test, dt_right_test = single_T_fit_continuum(t2, t1, 2, 
                                                                                           tlscp, normm, 
                                                                                           texp=40000, stpar=True, 
                                                                                           plot=False, Xplot=False)
            elif mode == 'realistic':
                temp_test, flux_test, dt_left_test, dt_right_test = single_T_realistic(t2, t1, 2, abundanc, 
                                                                                       tlscp, normm, 
                                                                                       texp=40000, stpar=True, 
                                                                                       plot=False, Xplot=False)

            temps_spec[str(normm)].append(temp_test[0])
            tt.append(temp_test[1])

            if dt_left_test != []:
                dts[str(normm)].append((dt_left_test[0], dt_right_test[0]))
                dd.append((dt_left_test[1], dt_right_test[1]))

	#    temps_prom = [art+tra for art, tra in zip(temps_prom, [aaa/N_usr for aaa in temps_spec[str(normm)]])]
        #temps_spec[str(normm)] = temps_prom
        
        temps_spec[str(normm)] += tt
        temps_diff[str(normm)] = temps_spec[str(normm)]-temps

        dts[str(normm)] += dd

        print("norm = "+str(normm)+": done")
        
    return temps, temps_spec, temps_diff, dts
    
  
def draw_Tspec_norm(temps, temps_spec, temps_diff, dts, mode):

    plt.figure(figsize=(14.5, 7))
    plt.suptitle('Difference of $T_{original}$ and $T_{spec}$ for various normalization values: '+mode, fontsize = 15)

    for normm in (1, 0.1, 0.011):

        plt.subplot(1,2,1)
        plt.xlim(temps[0]-0.3, temps[-1]+0.3)
        plt.ylim(temps[0]-0.3, temps[-1]+0.3)
        if normm == 0.011:
            y_errors = [temps_spec[str(normm)], temps_spec[str(normm)]]-np.transpose(dts['0.011'])
            plt.errorbar(temps, temps_spec[str(normm)], yerr=np.abs(y_errors), label='norm = '+str(normm), capsize=2.5, barsabove=True)
        else:
            plt.plot(temps, temps_spec[str(normm)], label='norm = '+str(normm))
        plt.plot(temps, temps, linestyle = '--', color='black')
        plt.xlabel('$T_{orig}$, keV', fontsize = 15)
        plt.ylabel('$T_{spec}$, keV', fontsize = 15)
        plt.grid()
        plt.legend()

        plt.subplot(3,2,2)
        if normm == 0.011:
            plt.errorbar(temps, temps_diff[str(normm)], yerr=np.abs(y_errors), label='norm = '+str(normm), capsize=2.5)
        else:
            plt.plot(temps, temps_diff[str(normm)], label='norm = '+str(normm))
        plt.axhline(0, linestyle = '--', color='black')
        #plt.xlabel('$T_{original}$, keV', fontsize = 15)
        #plt.xticks([],[])
        plt.ylabel('$T_{spec}-T_{orig}$, keV', fontsize = 12)
        plt.grid()

        plt.subplot(3,2,4)
        if normm == 0.011:
            plt.errorbar(temps, [tttt*normm for tttt in temps_diff[str(normm)]], yerr=np.abs(y_errors)*normm, label='norm = '+str(normm), capsize=2)
        else:
            plt.plot(temps, [tttt*normm**(1/2) for tttt in temps_diff[str(normm)]], label='norm = '+str(normm))
        plt.axhline(0, linestyle = '--', color='black')
        #plt.xlabel('$T_{original}$, keV', fontsize = 15)
        #plt.xticks([], [])
        plt.ylabel('$(T_{spec}-T_{orig}) \cdot norm^{0.5}$, keV', fontsize = 10)
        plt.grid()

        plt.subplot(3,2,6)
        if normm == 0.011:
            plt.errorbar(temps, [a/b for a, b in zip(temps_diff[str(normm)], temps)], yerr=np.abs(y_errors)/temps, label='norm = '+str(normm), capsize=2)
        else:
            plt.plot(temps, [a/b for a, b in zip(temps_diff[str(normm)], temps)], label='norm = '+str(normm))
        plt.axhline(0, linestyle = '--', color='black')
        plt.xlabel('$T_{orig}$, keV', fontsize = 15)
        plt.ylabel('$(T_{spec}-T_{orig})/T_{orig}$', fontsize = 12)
        plt.grid()

        #plt.subplot(2,3,6)
        #plt.plot(temps, [tttt*normm for tttt in temps_diff[str(normm)]], label='norm = '+str(normm))
        #plt.axhline(0, linestyle = '--', color='black')
        #plt.xlabel('$T_{original}$, keV', fontsize = 15)
        #plt.ylabel('$(T_{spec}-T_{original}) \cdot$ norm, keV', fontsize = 10)
        #plt.grid()
        #plt.legend(title='APEC')

    plt.show()
