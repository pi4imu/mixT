def draw_model(nrm, linesandcont):
    
    x.Plot("model")
    x.Plot.add = True
    #x.Plot.setID()
    xVals = x.Plot.x()
    modVals = x.Plot.model()
    yAdd1 = x.Plot.addComp(1)
    yAdd2 = x.Plot.addComp(2)
    plt.plot(xVals, modVals, linewidth = 3, label='Initial model', color='black')
    if linesandcont:
        yAdd3 = x.Plot.addComp(3)
        yAdd4 = x.Plot.addComp(4)
        yAdd12 = [a+b for a,b in zip(yAdd1, yAdd2)]
        yAdd34 = [c+d for c,d in zip(yAdd3, yAdd4)]
    else:
        yAdd12 = yAdd1
        yAdd34 = yAdd2
    plt.plot(xVals, yAdd12, linewidth = 2, linestyle = ":", 
             label=f'Low T', color='red')
    plt.plot(xVals, yAdd34, linewidth = 2, linestyle = "--", 
             label=f'High T', color='green')
    #plt.plot(xVals, yAdd12, label='sum 1 2')
    #plt.plot(xVals, yAdd34, label='sum 3 4')
    #plt.plot(xVals, yAdd3, label='3')
    #plt.plot(xVals, yAdd4, label='4')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(0.1, 14)
    plt.ylim(bottom=10**(-7), top=10**(1.5)*nrm)
    plt.legend(fontsize = 15, loc=1)
    add_plt_fancy()
    #plt.show()
    
    
def draw_best_model(nrm, linesandcont):
    
    x.Plot("model")
    modVals = x.Plot.model()
    xVals = x.Plot.x()
    if linesandcont:
        y1 = x.Plot.addComp(1)
        y2 = x.Plot.addComp(2)
        plt.plot(xVals, y1, label="APEC: Z=$Z_{{from \ fit}}$")
        #plt.plot(xVals, y2, label="APEC: Z=0")
        plt.plot(xVals, [-aa for aa in y2], label='APEC: Z=0', linestyle = '-.')
        #plt.plot(xVals, [c+d for c,d in zip(y1, y2)], label='sum')
    plt.plot(xVals, modVals, label=f"Best-fit", color='blue', alpha=0.3)
    plt.xscale('log')
    plt.yscale('log')
    plt.xlim(0.1, 14)
    plt.ylim(bottom=10**(-6), top=10**(1.5)*nrm)
    add_plt_fancy()
    #plt.title(f"Best-fit model (reduced $\\chi^2=$ {x.Fit.statistic/x.Fit.dof:.3f})", fontsize = 15)
    #plt.show()
    
    
def draw_goodness():
    
    x.Plot("goodness")
    xValsG = x.Plot.x()
    yValsG = x.Plot.y()
    #print(xValsG, yValsG)
    #plt.hist(yValsG, bins=xValsG)
    plt.bar(xValsG, height=yValsG, width = (np.min(xValsG)-np.max(xValsG))/len(xValsG))
    #plt.plot(xValsG, yValsG, label="goodness")
    plt.xlabel(x.Plot.labels()[0], fontsize = 15)
    plt.ylabel(x.Plot.labels()[1], fontsize = 15)
    plt.title(x.Plot.labels()[2], fontsize = 15)
    plt.xticks(size=15)
    plt.yticks(size=15)
    plt.grid()
    plt.axvline(x.Fit.statistic, linestyle = "--", color='red', linewidth = 3)
    
    
def draw_data_and_best_model(nrm, linesandcont):
    
    x.Plot("data")
    x.Plot.xAxis = "keV"
    x.Plot.add = True
    xVals = x.Plot.x()
    yVals = x.Plot.y()
    if linesandcont:
        y1 = x.Plot.addComp(1)
        y2 = x.Plot.addComp(2)
    modVals = x.Plot.model()
    #plt.plot(xVals, modVals, linewidth = 5, color = 'green')
    plt.yscale('log')
    #plt.plot(xVals, yVals, label='Data', color='black') 
    xErrs = x.Plot.xErr()
    yErrs = x.Plot.yErr()
    plt.errorbar(xVals, yVals, yErrs, xErrs, 
                 fmt = "none", ecolor = "black", label="Data with errors")
    if linesandcont:
        plt.plot(xVals, y1, label="APEC: Z=$Z_{{from \ fit}}$")
        plt.plot(xVals, [-aa for aa in y2], label="APEC: Z=0")
        summa = [c+d for c,d in zip(y1, y2)]
        plt.plot(xVals, summa, label='Best-fit', color='red') #modVals = summa
    else:
        plt.plot(xVals, modVals, label='Best-fit', color='red')
    add_plt_fancy()
    plt.xlabel("")
    plt.ylim(0.0001, 10**4*nrm)
    plt.title(f"Data and best-fit model (reduced $\\chi^2=$ {x.Fit.statistic/x.Fit.dof:.3f})", fontsize = 15)
    
    plt.subplot(6,2,6)
    x.Plot("resid")
    #plt.subplot(2,2,4)
    xValsR = x.Plot.x()
    yValsR = x.Plot.y()
    xErrsR = x.Plot.xErr()
    yErrsR = x.Plot.yErr()
    plt.errorbar(xValsR, yValsR, yErrsR, xErrsR, 
                 fmt = "none", ecolor = "black", label="errors")
    if not linesandcont:
        summa = modVals
    plt.scatter(xValsR, [a-b for a,b in zip(yVals, summa)], s=3, 
                color='blue', label='data minus \nbest fit')
    add_plt_fancy()
    plt.xlabel("")
    plt.ylabel("")
    plt.xticks([])

    plt.subplot(6,2,8)
    x.Plot("chi")
    #plt.subplot(2,2,4)
    xValsC = x.Plot.x()
    yValsC = x.Plot.y()
    #xErrsR = x.Plot.xErr()
    #yErrsR = x.Plot.yErr()
    plt.scatter(xValsC, yValsC, color = "black", label = "contribution to \n the fit statistic \n from each bin")
    add_plt_fancy()
    plt.title(x.Plot.labels()[1], fontsize = 14)
    plt.ylabel("")
    
            
def add_plt_fancy():
    
    plt.xlabel(x.Plot.labels()[0], fontsize = 14)
    plt.ylabel(x.Plot.labels()[1], fontsize = 14)
    plt.title(x.Plot.labels()[2], fontsize = 14)
    plt.xticks(size=15)
    plt.yticks(size=15)
    plt.grid()
    plt.xscale('log')
    plt.xticks([0.1, 1., 10.], [0.1, 1, 10])
    plt.legend(fontsize=15, loc=0)
    

def perform_fakeit(tname, expos):

    if tname == 'Chandra/ACIS-OLD':
        RMF_NAME = 'telescopes/chandra/djs50.ugc3957_v05.rmf' 
        ARF_NAME = 'telescopes/chandra/djs50.ugc3957_v05.arf' 
    elif tname == 'SRG/eROSITA':
        RMF_NAME = 'telescopes/erosita/erosita_pirmf_v20210719.rmf'
        ARF_NAME = 'telescopes/erosita/tm1_arf_open_000101v02.fits'
    elif tname == 'Chandra/ACIS-NEW':
        RMF_NAME = 'telescopes/acis/acisi.rmf'
        ARF_NAME = 'telescopes/acis/acisi_namp_qc.arf'
    elif tname == 'XMM-Newton/MOS':
        RMF_NAME = 'telescopes/xmm-newton/m1_thin1v9q19t5r5_all_15.rsp'
        ARF_NAME = ''    

    fs = x.FakeitSettings(response = RMF_NAME, 
                               arf = ARF_NAME, 
                        background = '', 
                          exposure = expos, 
                        correction = '', 
                      backExposure = '', 
                          fileName = 'fakeit.pha')
    x.AllData.fakeit(nSpectra = 1, 
                     settings = fs, 
                   applyStats = True,
                   filePrefix = "",
                      noWrite = True)
                      
    
def plot_contours_from_steppar(Nst, par_x_num, par_y_num, mmmodel, zoomin=True):

    #x.Plot.device = "/xs"
    x.Plot.device = "/null"
    x.Plot("contour")
    #x.Plot.device = "/null"
    chi2 = x.Plot.z()
    #print(x.Fit.statistic, x.Plot.contourLevels())

    #plt.subplot(1,2,1)
    
    par_x = mmmodel(par_x_num).values[0]
    par_y = mmmodel(par_y_num).values[0]

    #x11 = np.linspace(best_kT-0.005, best_kT+0.005, len(chi2))
    #y11 = np.linspace(abund_from_fit-0.01, abund_from_fit+0.01, len(chi2))
    #print(x11, y11)
    x11 = x.Fit.stepparResults(str(par_x_num))[0:Nst+1]
    y11 = x.Fit.stepparResults(str(par_y_num))[0::Nst+1]
    #print(x11, y11)
    X, Y = np.meshgrid(x11, y11)
    center = np.argwhere(chi2==np.min(chi2))[0]
    #plt.contour(X, Y, XCYC, levels=[x.Fit.statistic], colors='yellow')
    center_x = X[center[0]][center[1]]
    center_y = Y[center[0]][center[1]]

    #plt.contourf(X, Y, XCYC, 20, cmap='jet')

    contours = plt.contour(X, Y, chi2, levels=x.Plot.contourLevels(), colors='red')
    #plt.clabel(contours, inline=True, fontsize=8)
    
    exxxt = [x11[0], x11[-1], y11[0], y11[-1]]
    plt.xlim(x11[0], x11[-1])
    plt.ylim(y11[0], y11[-1])
    
    if zoomin:
        #for ii, seg in enumerate(contours.allsegs[2]):
        #    xleft, xright = np.min(seg), np.max(seg)
        #    print(np.min(seg), np.max(seg))
        #    plt.plot(seg[:,0], seg[:,1], '.-', label=ii)
        #plt.legend(fontsize=9, loc='best')

        #print(len(contours.allsegs), contours.allsegs)

        # finding borders of 3sigma level (or 2sigma) for fancier imshow
        seg = contours.allsegs[len(contours.allsegs)-1]
        x_left = np.min(seg[0][:,0])
        x_right = np.max(seg[0][:,0])
        y_left = np.min(seg[0][:,1])
        y_right = np.max(seg[0][:,1])

        scale_n = 2
        x_c = (center_x + par_x)/2
        x_l = x_c - scale_n*( x_c - np.min([par_x, center_x, x_left])) #(best_kT - x_left)
        x_r = x_c + scale_n*(-x_c + np.max([par_x, center_x, x_right])) #(x_right - best_kT)
        y_c = (center_y + par_y)/2
        y_l = y_c - scale_n*( y_c - np.min([par_y, center_y, y_left])) #(abund_from_fit - y_left)
        y_r = y_c + scale_n*(-y_c + np.max([par_y, center_y, y_right])) #(y_right - abund_from_fit)

        plt.xlim(x_l, x_r)
        plt.ylim(y_l, y_r)
        exxxt = x_l, x_r, y_l, y_r
        
    plt.imshow(chi2, extent=exxxt, aspect='auto', origin='lower', cmap='viridis')
    plt.colorbar(fraction=0.046, pad=0.04)

    plt.scatter(center_x, center_y, marker='+', c='cyan', label=f'$X_{{centered}}={center_x:.3f}, Y_{{centered}}={center_y:.3f}$')
    plt.axvline(center_x, color='cyan', alpha=0.3)
    plt.axhline(center_y, color='cyan', alpha=0.3)

    plt.axvline(par_x, color='yellow', alpha=0.3)
    plt.axhline(par_y, color='yellow', alpha=0.3)
    plt.scatter(par_x, par_y, marker='+', color='yellow', label=f'$X_{{min \ \\chi^2}}={par_x:.3f}, Y_{{min \ \\chi^2}}={par_y:.3f}$')

    #seg0 = contours.allsegs[0]
    #dx_left = np.min(seg[0][:,0]) - center_x
    #dx_right = np.max(seg[0][:,0]) - center_x
    #dy_left = np.min(seg[0][:,1]) - center_x
    #dy_right = np.max(seg[0][:,1]) - center_x
    
    plt.xlabel(x.Plot.labels()[0], fontsize = 12)
    plt.ylabel(x.Plot.labels()[1], fontsize = 12)
    #titlelabels = [float(f"{a:.3f}") for a in x.Plot.contourLevels()]
    #plt.title(x.Plot.labels()[2]+f"\n Cyan cross = {x.Fit.statistic:.3f}; yellow cross = {chi2[center[0]][center[1]]:.3f} \n levels = {*titlelabels,}", fontsize = 12)
    plt.legend()
    #plt.show()
    
    #return dx_left, dx_right, dy_left, dy_right 

    #plt.subplot(1,2,2)
    # 3d
    #ax = plt.axes(projection='3d')
    #ax.contour3D(X, Y, chi2, 100, cmap='viridis')#, rstride=1, cstride=1, edgecolor='none')
    #ax.set_xlabel(x.Plot.labels()[0])
    #ax.set_ylabel(x.Plot.labels()[1])
    #ax.set_zlabel('Chi-Squared')
    #ax.set_xlim(x_l, x_r)
    #ax.set_ylim(y_l, y_r)
    #plt.colorbar(fraction=0.046, pad=0.04)
    #ax.view_init(60, 35)
    #plt.show()
    

def perform_steppar(mmmodel, par_x_num, par_x_delta, par_y_num, par_y_delta, Nst):
            
    best_X = mmmodel(par_x_num).values[0]
    best_Y = mmmodel(par_y_num).values[0]
    x.Xset.parallel.steppar = 4
    N_steps = Nst
    par1_delta = par_x_delta
    par2_delta = par_y_delta
    x.Fit.steppar(f"{par_x_num} {best_X-par1_delta} {best_X+par1_delta} {N_steps} {par_y_num} {best_Y-par2_delta} {best_Y+par2_delta} {N_steps}")
    #print(x.Fit.stepparResults('2'))
    
   
#def add_background():
        #x.Plot("data")
        #bkg = x.Plot.backgroundVals()
        #plt.plot(xVals, yVals)
        
        
def print_parnames(MODEL):
    
    ncomp = len(MODEL.componentNames)
    for icomp in MODEL.componentNames:
        print (icomp, eval(f'MODEL.{icomp}.parameterNames'))
        

def universal_function(flin, delta1, delta2, beta):
     
    return np.exp( -(flin/delta1)**(2*beta) ) * np.exp( -(flin/delta2)**(8) )
