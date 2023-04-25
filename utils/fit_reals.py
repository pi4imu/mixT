def model_realistic(T_minnn, T_maxxx, f_minnn, f_maxxx, abund, nrm, cfluxxx=False):
    
    if not cfluxxx:
        
        mod = x.Model('phabs*(const*apec+const*apec)')
        mod.setPars(0.0, f_minnn, T_minnn, abund, 0, nrm, f_maxxx, T_maxxx, abund, 0, nrm)
        mod(9).link = "4"
        mod(10).link = "5"
        mod(11).link = "6"
    
        #mod = x.Model('phabs*(const*mekal+const*mekal)')
        #mod.setPars(0.0, f_minnn, T_minnn, 1e-6, abund, 0, 1, nrm, f_maxxx, T_maxxx, 1e-6, abund, 0, 1, nrm)
        #mod(12).link = "5"
        #mod(13).link = "6"
        #mod(15).link = "8"
    	
    x.AllModels.show()
    return mod

def single_T_realistic(T_minnn, T_maxxx, N_fmins, abund, telescope_name, nrm, texp, bnds, stpar=False, plot=False, Xplot=False):

    x.AllData.clear()
    x.AllData.removeDummyrsp()
    x.AllData.dummyrsp(lowE=0.1, highE=50.0, nBins=1024)
    x.Xset.addModelString("APEC_TRACE_ABUND", "0")
    
    if Xplot:
        x.Plot.device = "/xs"
    else:
        x.Plot.device = '/null'

    tspec_list = []
    flux_list=[]
    
    dt_lefts = []
    dt_rights = []

    for l in range(0, N_fmins):

        f_minnn = 0.0+1/(N_fmins-1)*(l)
        f_maxxx = 1-f_minnn
        
        # set model for fakeit
        mod = model_realistic(T_minnn, T_maxxx, f_minnn, f_maxxx, abund, nrm, cfluxxx=False)
        
        # calculating flux
        #x.AllModels.calcFlux('0.7 10.0')
        #fluxx = x.AllModels(1).flux[0]
        #flux_list.append(fluxx) # in units of ergs/cm2/s 
        #or use [3] in units of photons / s / cm^2
       
        # plot model
        if plot:
            plt.figure(figsize=(15, 17))
            plt.subplot(3,2,1)
            draw_model(nrm, linesandcont=False)
        
        # data from fake spectrum
        perform_fakeit(telescope_name, str(texp))
        x.AllData.ignore("**-0.7 10.0-**")             # IMPORTANT !
        x.AllData.show()
        
        bounds = bnds[l]
        
        # fitting
        x.AllModels.clear()
        x.AllData.removeDummyrsp()
        mod2fit = x.Model("phabs*apec")
        mod2fit.setPars(0.0, 1.0, abund, 0., nrm)
        mod2fit(1).frozen = True    # n_H 
        mod2fit(2).values = f"{(bounds[0]+bounds[1])/2}, 0.01, {bounds[0]}, {bounds[0]}, {2*bounds[1]}, {2*bounds[1]}" # temperature
        mod2fit(3).frozen = False   # abundance
        #mod2fit(4).frozen = False  # redshift   
        #mod2fit(5).frozen = True   # norm
        #mod2fit(5).values = f"{nrm}, -1, 0.0, 0.0, 1.1, 1.1"

        #x.AllModels.clear()
        #mod2fit = x.Model("phabs*mekal")
        #mod2fit.setPars(0.0, 1.0, 1e-6, abund, 0., 1, nrm)
        #mod2fit(1).frozen = True    # n_H 
        #mod2fit(4).frozen = False   # abundance
        ##mod2fit(5).frozen = False  # redshift   
        ##mod2fit(6).frozen = True   # norm
        ##mod2fit(6).values = f"{nrm}, -1, 0.0, 0.0, 1.1, 1.1"
        
        #x.AllData.ignore("bad")
        x.Fit.renorm('auto')
        x.Fit.nIterations = 100
        x.Fit.query = 'yes'
        x.Fit.weight = 'standard'
        x.Fit.statMethod = 'chi'
        x.Fit.perform()
        #x.Fit.goodness(10)
        #x.AllModels.show()
        x.Fit.show()

        # steppar
        if stpar:    
            N_steps = 20
            perform_steppar(mod2fit, 2, 0.05, 3, 0.5, N_steps)
            
        if nrm == 0.011:
            x.Fit.error('2')
            dleft, dright, _ = mod2fit(2).error
            #print(mod2fit(2).error)
            dt_lefts.append(dleft)
            dt_rights.append(dright)
            
        # return some parameters
        x.AllModels.show()
        best_kT = mod2fit(2).values[0]
        abund_from_fit = mod2fit(3).values[0]
        norm = mod2fit(5).values[0]
        tspec_list.append(best_kT)
        #print(best_kT)

        # calculating flux
        #fluxx = x.AllData(1).rate[0]
        #flux_list.append(fluxx) # in units of counts / s
        x.AllModels.calcFlux('0.7 10.0')
        fluxx = x.AllData(1).flux[0]
        flux_list.append(fluxx) # in units of ergs/cm2/s 
        # or use [3] in units of photons / s / cm^2      
        
        if plot:
            plt.suptitle(f'$T_{{min}}={T_minnn} \ keV, \ T_{{max}}={T_maxxx} \ keV, \ f_{{min}}={f_minnn:.2f}, \ f_{{max}}={f_maxxx:.2f}$ \n $t_{{exp}}={texp} \ s, \ Z_{{model}} ={abund}, \ Z_{{from \ fit}} = {abund_from_fit:.2f}, \ norm = {norm:.4f}$ \n $T_{{spec}}={best_kT:.4f} \ keV, \ $Flux$ = {fluxx*10**12:.3f}\cdot 10^{{-12}} \ ergs/cm^2/s$ \n ', fontsize = 20)
            plt.subplot(3,2,3)
            if stpar:
                plot_contours_from_steppar(N_steps, 2, 3, mod2fit, zoomin=True)
            else:
                draw_best_model(nrm, linesandcont=False)
                #draw_goodness()
            plt.subplot(3,2,2)
            draw_data_and_best_model(nrm, linesandcont=False)
            plt.subplot(3,2,1)
            draw_best_model(nrm, linesandcont=False)
            plt.show()
            
        #x.Plot.commands=()
        x.AllData.clear()
        x.AllModels.clear()

    return tspec_list, flux_list, dt_lefts, dt_rights
    
    
def f_line(telescope_name, temperature, mode, abundance):
    
    c__T = c_T(telescope_name, temperature, mode)
    l__T = l_T(telescope_name, temperature, mode)
    
    return l__T*abundance / (l__T*abundance + c__T)
