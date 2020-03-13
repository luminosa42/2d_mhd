# Python script to read in hydro/mhd solver data and plot contours
# Provided: density, pressure, magnetic field strength
# Yinghe Lu 4/8/16

import os
import numpy as np
import matplotlib.pyplot as plt

for num in range(1,22):

    fname="output_mhd_128_"+str(num).zfill(4)
    print "reading file", num

    xarray = []
    yarray = []
    rhoarray = []
    parray = []
    Barray = []
    x = []
    y = []
    rho = []
    p = []
    B = []

    f = open(fname, 'r')
    lines = f.readlines()

    Nx = int(lines[1].split()[2])
    Ny = int(lines[2].split()[2])
    dx = float(lines[3].split()[2])
    dy = float(lines[4].split()[2])

    for line in lines:
        line = line.split()	
        if len(line) == 0:
            xarray.append(x)
	    yarray.append(y)
	    rhoarray.append(rho)
	    parray.append(p)
	    Barray.append(B)
	    x = []
	    y = []
	    rho = []
	    p = []
	    B = []

        elif len(line) == 9:
	    xtemp = (float(line[0]) - 0.5 * Nx) * dx
            ytemp = (float(line[1]) - 0.5 * Ny) * dy
	    x.append(xtemp)
	    y.append(ytemp)
	    rho.append(float(line[2]))
	    p.append(float(line[3]))
	    B_x = float(line[7])
            B_y = float(line[8])
	    B.append((B_x*B_x+B_y*B_y) ** 0.5)


    f.close()

    #print len(xarray[0]), len(yarray[0]), len(rhoarray[0])
    #X, Y = np.meshgrid(xarray, yarray)

    pres = np.ma.array(parray)
    CS = plt.contourf(xarray,yarray,pres, 15, vmax = 10.,
                  alpha=0.5,cmap=plt.cm.hsv)
    plt.colorbar()
    plt.title('Pressure contour: # '+str(num))
    plt.savefig("contour_pres"+str(num).zfill(2)+".png")
    plt.clf()

    dens = np.ma.array(rhoarray)
    plt.contourf(xarray,yarray,dens, 15, vmin=0.25, vmax=1.65,
                  alpha=0.5,cmap=plt.cm.hsv)
    plt.colorbar()
    plt.title('Density contour: # '+str(num))
    plt.savefig("contour_dens"+str(num).zfill(2)+".png")
    plt.clf()

    mag = np.ma.array(Barray)
    plt.contourf(xarray,yarray,mag, 15, xmin = 7., xmax = 13.,
                  alpha=0.5,cmap=plt.cm.hsv)
    plt.colorbar()
    plt.title('Magnetic strength contour: # '+str(num))
    plt.savefig("contour_mag"+str(num).zfill(2)+".png")
    plt.clf()
