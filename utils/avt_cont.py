def get_Tspec_continuum_eq45(fminnn, Tmin, Tmax):
    
    # derive T_spec from given values 
    # of T_min, T_max and f_min
    
    Tspec = []
    
    for fmin in fminnn:
        
        temperatures = [Tmin, Tmax]
        weights = temperatures*0

        alpha = 0.75
        
        weights = [temperatures[i]**(-alpha) for i in range(0, len(temperatures))]

        weights = np.multiply(weights, [fmin, (1-fmin)])

        #print(weights)

        num = np.dot(weights, temperatures)
        denom = sum(weights)

        #print(num/denom)
        
        Tspec.append(num/denom)

    return Tspec

def calc_c_T(T, T_left, T_right, telescope_name, Xplot=False):
    
    # calculating photon count rate for continuum
    
    x.AllData.clear()
    x.AllData.removeDummyrsp()
    x.AllData.dummyrsp(lowE=0.1, highE=50.0, nBins=1024)
    x.Xset.addModelString("APEC_TRACE_ABUND", "0")
    
    if telescope_name == 'Chandra/ACIS-OLD':
        RMF_NAME = 'telescopes/chandra/djs50.ugc3957_v05.rmf' 
        ARF_NAME = 'telescopes/chandra/djs50.ugc3957_v05.arf' 
    elif telescope_name == 'SRG/eROSITA':
        RMF_NAME = 'telescopes/erosita/erosita_pirmf_v20210719.rmf'
        ARF_NAME = 'telescopes/erosita/tm1_arf_open_000101v02.fits'
    elif telescope_name == 'Chandra/ACIS-NEW':
        RMF_NAME = 'telescopes/acis/acisi.rmf'
        ARF_NAME = 'telescopes/acis/acisi_namp_qc.arf'
    elif telescope_name == 'XMM-Newton/MOS':
        RMF_NAME = 'telescopes/xmm-newton/m1_thin1v9q19t5r5_all_15.rsp'
        ARF_NAME = ''
    elif telescope_name == 'Chandra/ACIS-2002':
        RMF_NAME = 'telescopes/chandra-2002/acisf03243_000N022_r0087_rmf3.fits'
        ARF_NAME = 'telescopes/chandra-2002/acisf03243_000N022_r0087_arf3.fits'

    # set model for fakeit
    
    mod = x.Model('phabs*apec')
    mod.setPars(0.00, T, 0.0, 0, 1)
    x.AllModels.show()
    
    if Xplot:
        x.Plot.device = '/xs'
    else:
        x.Plot.device = '/null'

    # fake spectrum
    fs = x.FakeitSettings(response = RMF_NAME, 
                               arf = ARF_NAME, 
                        background = '', 
                          exposure = '40000', 
                        correction = '', 
                      backExposure = '', 
                          fileName = 'fakeit.pha')
    x.AllData.fakeit(nSpectra = 1, 
                     settings = fs, 
                   applyStats = True,
                   filePrefix = "",
                      noWrite = True)

    x.AllData.ignore(f"**-{T_left} {T_right}-**")             # IMPORTANT !
    
    x.Plot.xAxis = "keV"
    #x.AllData.show()
    x.Plot("ldata")
    xVals = x.Plot.x()
    yVals = x.Plot.y()
    
    s1 = x.AllData(1).rate[0]
    
    x.AllModels.calcFlux(f"{T_left} {T_right}")
    flx = x.AllData(1).flux[0]
    
    x.AllData.clear()
    x.AllModels.clear()
     
    return flx, s1, np.dot(xVals, yVals)/sum(yVals)
    
def c_T(telescope_name, temperature, mode):

    if telescope_name == 'Chandra/ACIS-OLD':
    	tt = 'Chandra_ACIS-OLD'
    if telescope_name == 'SRG/eROSITA':
    	tt = 'SRG_eROSITA'
    if telescope_name == 'XMM-Newton/MOS':
    	tt = 'XMM-Newton_MOS'
    if telescope_name == 'Chandra/ACIS-2002':
    	tt = 'Chandra_ACIS-2002'
    
    read_cT = pd.read_csv('c(T)/c(T)_'+str(tt)+'.csv', header=None, delimiter=' ')
    temps = read_cT[0].values
    flux_photons = read_cT[1].values
    count_rate = read_cT[2].values
    npdot = read_cT[3].values
    
    if mode == 'flux':
        return flux_photons[np.argmin(np.abs(temps - temperature))]
    elif mode == 'rate':
        return count_rate[np.argmin(np.abs(temps - temperature))]
    elif mode == 'npdot':
        return npdot[np.argmin(np.abs(temps - temperature))]


def get_Tspec_continuum_eq46(fminnn, Tmin, Tmax, alpha, telescope_name):
    
    # derive T_spec from given values 
    # of T_min, T_max and f_min
    
    Tspec = []
    temperatures = [Tmin, Tmax]
    
    for fmin in fminnn:
        
        weights = temperatures*0
        
        weights = [ts**(-alpha) for ts in temperatures]

        weights = np.multiply(weights, [fmin, (1-fmin)])
        
        #c_T_min = calc_c_T(Tmin, 0.7, 10.0, telescope_name, Xplot=False)
        #c_T_max = calc_c_T(Tmax, 0.7, 10.0, telescope_name, Xplot=False)
        
        md = 'rate'
        
        c_T_min = c_T(telescope_name, Tmin, md)
        c_T_max = c_T(telescope_name, Tmax, md)
        
        weights = np.multiply(weights, [c_T_min, c_T_max])

        #print(weights)

        num = np.dot(weights, temperatures)
        denom = sum(weights)

        #print(num/denom)
        
        Tspec.append(num/denom)

    return Tspec
    
def fancy_fig4():

	plt.ylim(0.1, 30)
	plt.yscale('log')
	plt.xticks(size=15)
	plt.yticks([0.1, 1, 10], [0.1, 1, 10], size=15)
	plt.xlabel('$f_{min}$', fontsize = 15)
	plt.ylabel('$T_{spec}$ (keV)', fontsize = 15)
	plt.title('Continuum-dominated spectra ('+telescope+') \n $\\alpha=$'+str(alpha_current), fontsize = 15)

	handles, labels = plt.gca().get_legend_handles_labels()
	#line_n = Line2D([], [], label='Naive weighting', color='black', linestyle='--', linewidth=2)
	line_e = Line2D([], [], label='$T_{spec}$ from eq. 4,5', color='black', linestyle=':', linewidth=3)
	dots_f = Line2D([], [], label='$T_{spec}$ from eq. 4,6', color='black', linestyle='-', linewidth=1)
	dots_T = Line2D([], [], label='Single-T fit', color='black', marker='.', linewidth=0, markersize=12)
	handles.extend([line_e, dots_f, dots_T])
	plt.legend(handles=handles, fontsize=15)
 
