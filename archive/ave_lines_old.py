# Defining functions to draw E(T) for lines:

def f(temperature, eMean, tList): 
    
    # <E> = f(T)
    #returns energy
    
    return eMean[np.argmin(np.abs(tList - temperature))]

    
def f_inv(energy, eMean, tList): 
    
    # T = f^(-1) (<E>)
    #returns temperature
    
    return tList[np.argmin(np.abs(eMean - energy))]
    
    
def get_data(dataName, show_table=False):
    
    # extracts data obtained from stats.sh
    # dataName is stats_NN.dat

    headers = [ 'Flux', 'Abund', 'T', 'z', 'n_H', 'Chnls', '$E_{min}$', '$E_{max}$', '$E_{sum}$', 'cs', 'ecs', 'rate'] 
    table = pd.read_csv(dataName, sep = ' ', names = headers, index_col=False)
    
    # adding column with E_mean
    
    table['$E_{mean}$'] = table['ecs']/table['cs']
    
    if show_table:
        display(table)
        
    # returns table as numpy array (?)
    
    return table.to_numpy().astype(float)
    

def plot_E_T(dataName, spectra_type, telescope):
    
    # plotting <E>(T) as in Fig.2 for given dataName
    # you should add spectra type and telescope name by yourself
    
    # possible values:
    # spectra_type = 'Line' or sectra_type = 'Continuum'
    # telescope = 'Chandra' or telescope = 'SRG/eROSITA'

    #print(dataName)
    #headers = [ 'Flux', 'Abund', 'T', 'z', 'n_H', 'Chnls', '$E_{min}$', '$E_{max}$', '$E_{sum}$', 'cs', 'ecs' ] 
    #table = pd.read_csv(dataName, sep = ' ', names = headers)
    #adding column with E_mean
    #table['$E_{mean}$'] = table['ecs']/table['cs'] 
    #data = table.to_numpy().astype(float)
    
    data = get_data(dataName, show_table=False)
    
    e_mean = data[:,12]
    temp = data[:,2]
    
    # taking first ever values of corresponding characteristics 
    
    #abundance = table['Abund'].to_numpy().astype(float)[0]
    #absorption = table['n_H'].to_numpy().astype(float)[0]
    #redshift = table['z'].to_numpy().astype(float)[0]
    
    abundance = data[0, 1]
    absorption = data[0, 4]
    redshift = data[0, 3]
    
    #from scipy.interpolate import make_interp_spline, BSpline
    #xnew = np.linspace(0.01, 10, 1001)
    #spl = make_interp_spline(temp, e_mean, k=5)  # type: BSpline
    #smooth = spl(xnew)
    #plt.plot(xnew, spl(xnew), linewidth = 3, label = 'spline')
    
    plt.plot(temp, e_mean, linewidth = 3, label = telescope) # label = spectra_type + ', ' + telescope
    plt.xlim(0.08, 15)
    plt.ylim(0, 3)
    plt.xlabel('Temperature (keV)', fontsize = 15)
    plt.ylabel('Average energy (keV)', fontsize = 15)
    plt.title(spectra_type + '-dominated spectra (' + telescope + 
              ') \n $n_H =' + str(absorption) + '\cdot 10^{22} \ cm^{-2}$; z = ' + 
              str(redshift), fontsize = 15) # +'; Z = '+str(abundance)+' Solar units')
    plt.xticks(size=15)
    plt.yticks(size=15)
    plt.xscale('log')
    #plt.yscale('log')
    #plt.grid()
    plt.legend(fontsize = 15, loc=2)
    
    plt.xticks([0.1, 1., 10.], [0.1, 1, 10], size=15)    

    #return data[:,11], data[:,2]
    #return e_mean, temp
    
def add_T(Tmin, Tmax, dataName):
    
    # adding vertical lines for given temperatures
    # and horizontal lines for corresponding energies
    
    data = get_data(dataName, show_table=False)
    
    eMean = data[:,12]
    tList = data[:,2]
    
    #E_min = e_mean[np.argmin(np.abs(temp - Tmin))]
    #E_max = e_mean[np.argmin(np.abs(temp - Tmax))]
    E_min = f(Tmin, eMean, tList)
    E_max = f(Tmax, eMean, tList)
    
    #print('E_max =', round(E_max,2), 'keV')
    #print('E_min =', round(E_min,2), 'keV')
    
    plt.axvline(Tmin, ymin=0, ymax=E_min/3.0, linewidth=2, linestyle="-.", color='black')
    plt.axvline(Tmax, ymin=0, ymax=E_max/3.0,  linewidth=2, linestyle="--", color='black')
    
    plt.axhline(E_min, xmax = 0.35, linewidth=2, linestyle="-.", color='black', label=f'$T_{{min}}  = {Tmin:.2f} \ keV, E_{{min}} = {E_min:.2f} \ keV$')
    plt.axhline(E_max, xmax = 0.72, linewidth=2, linestyle="--", color='black', label=f'$T_{{max}}  = {Tmax:.2f} \ keV, E_{{max}} = {E_max:.2f} \ keV$')
    
    plt.legend(fontsize=15, loc = 0)
    
# Алгоритм вычисления $T_{spec}$ по графикам $E \ (T)$ согласно формулам (1-3):

def get_Tspec_lines(fmin, Tmin, Tmax, Data):
    
    # derive T_spec from given values 
    # of T_min, T_max, f_min and for given <E>(T)

    E_1 = f(Tmin, Data[:,12], Data[:,2])
    E_2 = f(Tmax, Data[:,12], Data[:,2])

    S_j_1 = Data[np.argmin(np.abs(Data[:,2] - Tmin)), 9]
    S_j_2 = Data[np.argmin(np.abs(Data[:,2] - Tmax)), 9]
    
    #print(Data[np.argmin(np.abs(tList - Tmin)),9], Data[np.argmin(np.abs(tList - Tmax)),9])
    #print()

    num =    fmin*S_j_1*E_1 + (1-fmin)*S_j_2*E_2
    denom =  fmin*S_j_1     + (1-fmin)*S_j_2
    
    #print(num, denom)
    
    #Etot = num/denom
    #Tspec = f_inv(Etot, Data[:,11], Data[:,2])
    #plt.scatter(fmin, Tspec, color="blue")

    return f_inv(num/denom, Data[:,12], Data[:,2])


def plot_Tspec_fmin(Tmin, Tmax, fmin, Data, naive=False, lstyle = '-', cline = 'blue'):

    # draw plot like Fig.3

    #f_min = np.linspace(0, 1, N_fmins+1)
    
    #naive weighting
    if naive:
        plt.plot(fmin, fmin*T_min+(1-fmin)*T_max, linestyle = '--', linewidth=2, color='black')

    T_spec = fmin*0

    for i in range(0, len(fmin)):

        T_spec[i] = get_Tspec_lines(fmin[i], T_min, T_max, Data)

    #print('*************')

    plt.plot(fmin, T_spec, linewidth=2, linestyle = lstyle, color = cline)
    #, label='eq. [1-3]')#'from '+str(T_max)+' to '+str(T_min))

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
