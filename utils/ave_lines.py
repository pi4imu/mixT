def calc_l_T(T, telescope_name, Xplot=False):

    x.AllData.clear()
    x.AllData.removeDummyrsp()
    x.AllData.dummyrsp(lowE=0.1, highE=10.0, nBins=1024)
    x.Xset.addModelString("APEC_TRACE_ABUND", "0")

    if Xplot:
        x.Plot.device = '/xs'
    else:
        x.Plot.device = '/null'
    
    Ab = 1.0
    Norm = 1
    z = 0
    nH = 0.0
    
    mod = x.Model('phabs*(apec+const*apec)')
    mod.setPars(nH, T, Ab, z, Norm, 1., T, 0.0, z, Norm)
    mod(6).values =  "-1, 0.0001, -1, -1, 1, 1"
    mod(7).link = "2"
    mod(9).link = "4"
    mod(10).link = "5"
    
    x.AllModels.show()
    #x.Plot.show()
    #x.AllModels.setEnergies("0.1 10.0 10 log")
    #x.AllModels.setEnergies("reset")
    #x.Plot("model")
    
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
    
    fs = x.FakeitSettings(response = RMF_NAME, 
                               arf = ARF_NAME, 
                        background = '', 
                          exposure = 40000, 
                        correction = '', 
                      backExposure = '', 
                          fileName = 'fakeit.pha')
    x.AllData.fakeit(nSpectra = 1, 
                     settings = fs, 
                   applyStats = True,
                   filePrefix = "",
                      noWrite = True)

    x.AllData.ignore("**-0.7 10.0-**")             # IMPORTANT !
    #x.AllData.show()
    x.AllModels.setEnergies("reset")
    
    x.Plot("data")
    x.Plot.xAxis = "keV"
    xVals = x.Plot.x()
    yVals = x.Plot.y()
         
    x.AllModels.calcFlux('0.7 10.0')
    flx = x.AllData(1).flux[0]
    
    cr = x.AllData(1).rate[2]
    
    #channels = len(x.AllData(1).noticed)
    ens = x.AllData(1).energies
    #e_min = min(ens)[0]
    #e_max = max(ens)[1]
    
    E_i = np.zeros(len(ens))
    dE = np.zeros(len(ens))
    
    for i in ens:
               
        dE[ens.index(i)] = i[1]-i[0]
        E_i[ens.index(i)] = (i[0]+i[1])/2
    
    s_i = x.AllData(1).values
    
    return flx, cr, np.dot(xVals, yVals)/sum(yVals), np.dot(E_i, s_i)/cr


def l_T(telescope_name, temperature, mode):

    if telescope_name == 'Chandra/ACIS-OLD':
        tt = 'Chandra_ACIS-OLD'
    if telescope_name == 'SRG/eROSITA':
        tt = 'SRG_eROSITA'
    if telescope_name == 'XMM-Newton/MOS':
        tt = 'XMM-Newton_MOS'
    
    read_lT = pd.read_csv("l(T)/l(T)_"+str(tt)+".csv", header=None, delimiter=' ')
    temps1 = read_lT[0].values
    flux_photons1 = read_lT[1].values
    count_rate1 = read_lT[2].values
    npdot1 = read_lT[3].values
    av_en1 = read_lT[4].values
    
    if mode == 'flux':
        return flux_photons1[np.argmin(np.abs(temps1 - temperature))]
    elif mode == 'rate':
        return count_rate1[np.argmin(np.abs(temps1 - temperature))]
    elif mode == 'npdot':
        return npdot1[np.argmin(np.abs(temps1 - temperature))]
    elif mode == 'av_en':
        return av_en1[np.argmin(np.abs(temps1 - temperature))]
        
        
def e_from_t(TT, tlscp):
    
    return l_T(tlscp, TT, 'av_en')
    
    
def t_from_e(EE, telescope_name):
    
    if telescope_name == 'Chandra/ACIS-OLD':
        tt = 'Chandra_ACIS-OLD'
    if telescope_name == 'SRG/eROSITA':
        tt = 'SRG_eROSITA'
    if telescope_name == 'XMM-Newton/MOS':
        tt = 'XMM-Newton_MOS'
    
    table = pd.read_csv('l(T)/l(T)_'+tt+'.csv', sep = ' ', header=None, index_col=False)
    
    return table[0][np.argmin(np.abs(table[4]-EE))]
    
    
def calc_Tspec_from_avE(Tmin, Tmax, N_fmins, telescope_name):
    
    # derive T_spec from given values 
    # of T_min, T_max, f_min and for given <E>(T)
    
    if telescope_name == 'Chandra/ACIS-OLD':
        tt = 'Chandra_ACIS-OLD'
    if telescope_name == 'SRG/eROSITA':
        tt = 'SRG_eROSITA'
    if telescope_name == 'XMM-Newton/MOS':
        tt = 'XMM-Newton_MOS'
    
    table = pd.read_csv('l(T)/l(T)_'+tt+'.csv', sep = ' ', header=None, index_col=False)

    E_1 = e_from_t(Tmin, telescope_name)
    E_2 = e_from_t(Tmax, telescope_name)
    
    S_j_1 = table[2][np.argmin(np.abs(table[0] - Tmin))]
    S_j_2 = table[2][np.argmin(np.abs(table[0] - Tmax))]

    TTT = []

    for fmin in np.linspace(0, 1, N_fmins):
    
    	num =    fmin*S_j_1*E_1 + (1-fmin)*S_j_2*E_2
    	denom =  fmin*S_j_1     + (1-fmin)*S_j_2
       
    	E_tot = num/denom
    	TTT.append(t_from_e(E_tot, telescope_name))
    	#plt.scatter(fmin, TTT, color="blue")

    return TTT
    
    
def plot_Tspec_fmin_details():  
    
    plt.xticks(size=15)
    plt.yticks(size=15)
    plt.xlabel('$f_{min}$', fontsize = 15)
    plt.ylabel('$T_{spec}$ (keV)', fontsize = 15)
    
    handles, labels = plt.gca().get_legend_handles_labels()
    line_n = Line2D([], [], label='Naive weighting', color='black', linestyle='--', linewidth=2)
    line_e = Line2D([], [], label='$T_{spec}$ from eq. [1-3]', color='blue', linestyle='-', linewidth=2)
    dots_f = Line2D([], [], label='Single-T fit', color='black', marker='.', linestyle='None', markersize=12)
    handles.extend([line_n, line_e, dots_f])
    plt.legend(handles=handles, fontsize=15)
