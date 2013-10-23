#!/usr/bin/env python
import sys
import matplotlib.pyplot as plt
from matplotlib import rc
from mesa import *
rc('font',**{'family':'sans-serif','sans-serif':['Helvetica'], 'size' : 20})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
rc('text', usetex=True)
import numpy

"""
Plots the evolution of a series of models with fixed initial masses and different periods.
Options are:
    mass_donor: string with the donor mass
    mass_accretor: string with the accretor mass
    periods: array of strings with periods, must be sorted in increasing order
    z_label: data to plot in z direction (colormap), options are:
        total_mass_lost
        omega_crit
        lg_mstar_dot_1
        lg_mstar_dot_2
        n14
    y_label:
        mass
        age
    folder: base folder where models are located

The folder with the data is assumed to be in:
folder+mass_donor+"_"+mass_accretor+"_"+period

Output is stored in:
mass_donor+"_"+mass_accretor+"_"+z_label+"_"+ylabel+".pdf"

"""
def plot_fixed_masses(mass_donor, mass_accretor, periods, \
      z_label = "lg_mstar_dot_1", y_label = "mass", folder ='.', num_levels=30):
    history = []
    binary_history = []
    y_axis = []
    z_axis = []
    min_s1_h1 = []
    min_s2_h1 = []
    setted_limits = False
    min_z=0
    max_z=0
    folder = folder + "/"

    #solve boundaries of columns in plot
    boundaries = zeros(len(periods)+1)
    for i in range(1,len(boundaries)-1):
        boundaries[i] = (float(periods[i-1]) + float(periods[i])) / 2
    boundaries[0] = 2*float(periods[0]) - boundaries[1]
    boundaries[len(boundaries)-1] = 2*float(periods[0]) - boundaries[1]
    
    #load data
    i = 0
    for period in periods:
        folder_name = folder+mass_donor+"_"+mass_accretor+"_"+period
        print folder_name
        history.append(history_data(folder_name + "/LOGS2", clean_starlog=False))
        binary_history.append(history_data(folder_name, slname="binary_history.data", clean_starlog=False))

        if y_label == "mass":
            y_axis.append(binary_history[i].get("star_1_mass"))
        elif y_label == "age":
            y_axis.append(binary_history[i].get("age"))
        else:
            print "Unknown y_label "+y_label
            return

        if z_label == "lg_mstar_dot_1":
            z_axis.append(binary_history[i].get("lg_mstar_dot_1"))
        elif z_label == "lg_mstar_dot_2":
            z_axis.append(binary_history[i].get("lg_mstar_dot_2"))
        elif z_label == "omega_crit":
            z_axis.append(history[i].get("surf_avg_omega_div_omega_crit"))
        elif z_label == "total_mass_lost":
            z_axis.append(binary_history[i].get("star_1_mass")+binary_history[i].get("star_2_mass"))
            z_axis[i] = z_axis[i][0]-z_axis[i]
        elif z_label == "n14":
            z_axis.append(history[i].get("surface_n14"))
        else:
            print "Unknown z_label "+z_label
            return
    
        min_s2_h1.append(min(history[i].get("center_h1")))
        min_s1_h1.append(min(history_data(folder_name + "/LOGS1", clean_starlog=False).get("center_h1")))
    
        #z_axis.append(12+numpy.log10(history[i].get("surface_n14")/history[i].get("surface_h1")))
        if setted_limits:
            min_z = min(min_z,min(z_axis[i]))
            max_z = max(max_z,max(z_axis[i]))
        else:
            min_z = min(z_axis[i])
            max_z = max(z_axis[i])
            setted_limits = True
        i = i + 1
    
    if z_label=="omega_crit":
        min_z = 0.0
        max_z = 1.0
    
    norm = matplotlib.colors.Normalize(vmin=min_z, vmax=max_z)
    levels = []
    for i in range(0,num_levels):
        levels.append(min_z+(max_z-min_z)/(num_levels-1)*i)
    
    #plot and save
    for i in range(0,len(history)):
        num_models = min(len(z_axis[i]),len(y_axis[i]))
        Xvalues = numpy.zeros((2,num_models))
        Yvalues = numpy.zeros((2,num_models))
        Zvalues = numpy.zeros((2,num_models))
        period = float(periods[i])
        Xvalues[0,:] = boundaries[i]
        Xvalues[1,:] = boundaries[i+1]
        for j in range(0,num_models):
            Yvalues[:,j] = y_axis[i][j]
            Zvalues[:,j] = z_axis[i][j]
        plt.contourf(Xvalues, Yvalues, Zvalues, cmap = plt.cm.gnuplot, norm=norm, levels = levels)
    
        marker_color = ""
        marker_symbol = ""
        if (min_s1_h1[i] > 0.0001 and min_s2_h1[i] > 0.0001):
            #contact
            marker_color = "k"
            marker_symbol = "s"
        elif (min_s2_h1[i] > min_s1_h1[i]):
            #donor out of MS
            marker_color = "r"
            marker_symbol = "^"
        else:
            #accretor out of MS
            marker_color = "b"
            marker_symbol = "v"
        if y_label == "mass":
            plt.plot(period,min(y_axis[i]),marker_color+marker_symbol,markersize=10)
        elif y_label == "age":
            plt.plot(period,max(y_axis[i]),marker_color+marker_symbol,markersize=10)
    
    #create legend for contact info
    #contArtist = plt.Line2D((0,0),(0,0), color='k', marker='s',linestyle=" ")
    #s1Artist = plt.Line2D((0,0),(0,0), color='r', marker='^',linestyle=" ")
    #s2Artist = plt.Line2D((0,0),(0,0), color='b', marker='v',linestyle=" ")
    #plt.gca().legend([contArtist,s1Artist,s2Artist],['Contact', 'Donor off MS', 'Accretor off MS'],\
    #        numpoints=1,prop={'size':7},loc="lower left",bbox_to_anchor=(-0.20, 1.0),fancybox=True, shadow=True)
    
    plt.clim(min_z, max_z)
    if z_label=="total_mass_lost":
        cb=plt.colorbar(cmap = plt.cm.gnuplot, norm = norm,format='%.2f')
        cb.ax.set_ylabel("$(M_1+M_2)-(M_1+M_2)_i$ [$\\mathrm{M_\\odot}$]",rotation=270)
    elif z_label=="omega_crit":
        cb=plt.colorbar(cmap = plt.cm.gnuplot, norm = norm,format='%.2f')
        cb.ax.set_ylabel("$\\Omega/\\Omega_{crit}\\;accretor$",rotation=270)
    elif z_label=="lg_mstar_dot_1":
        cb=plt.colorbar(cmap = plt.cm.gnuplot, norm = norm,format='%.2f')
        cb.ax.set_ylabel("$\\dot{M}_1$ [$\\mathrm{M_\\odot\;yr^{-1}}$]",rotation=270)
    elif z_label=="lg_mstar_dot_2":
        cb=plt.colorbar(cmap = plt.cm.gnuplot, norm = norm,format='%.2f')
        cb.ax.set_ylabel("$\\dot{M}_2$ [$\\mathrm{M_\\odot\;yr^{-1}}$]",rotation=270)
    elif z_label=="n14":
        cb=plt.colorbar(cmap = plt.cm.gnuplot, norm = norm,format='%.2e')
        cb.ax.set_ylabel("Accretor $^{14}N$ mass fraction",rotation=270)
    
    plt.xlabel("$\mathrm{P_i\;[d]}$", labelpad = -0.01)
    if y_label == "mass":
        plt.ylabel("$\mathrm{M_1/M\odot}$")
    elif y_label == "age":
        plt.ylabel("$\mathrm{age\;[yr]}$")
        
    plt.title("$M_1 = "+mass_donor+"M_\odot$, $M_2 = "+mass_accretor+"M_\odot$")
    plt.savefig(mass_donor+"_"+mass_accretor+"_"+z_label+"_VS_"+y_label+".pdf", dpi=None, facecolor='w', edgecolor='w',
                    orientation='portrait', papertype=None, format=None,
                            transparent=False)
    print "saved output to: "+mass_donor+"_"+mass_accretor+"_"+z_label+"_VS_"+y_label+".pdf"
    
    plt.close()
