import matplotlib.pyplot as plt
import numpy as np

method = "mhd"
grid = "128"

rho = []
P = []
r = []
vr = []
vt = []  #transverse velocity
B = []

rho_a = []
P_a = []
r_a = []
vr_a = []
vt_a = []

simfile = open("output_"+method+"_"+grid+"_0011")
lines = simfile.readlines()

#time = float(lines[0].split()[3])
Nx = int(lines[1].split()[2])
Ny = int(lines[2].split()[2])
dx = float(lines[3].split()[2])
dy = float(lines[4].split()[2])

for line in lines:
    line = line.split()
    if len(line) == 9:
        x = (float(line[0]) - 0.5 * Nx) * dx
        y = (float(line[1]) - 0.5 * Ny) * dy
        radius = (x*x + y*y) ** 0.5
        r.append(radius)
        rho.append(float(line[2]))        
        P.append(float(line[3]))        
        try:
            v_x = float(line[5])
        except:
            v_x = 0.
        try:
            v_y = float(line[6])
        except:
            v_y = 0.
        if radius == 0.:
            vr.append(0.)
            vt.append(0.)
	else:
            vr.append((v_x * x + v_y * y) / radius)
            vt.append((-v_x * y + v_y * x) / radius)

        B_x = float(line[7])
        B_y = float(line[8])
	B.append((B_x*B_x+B_y*B_y) ** 0.5)
   

simfile.close()

anafile	= open("sedov_analytic_512.out")
lines = anafile.readlines()

Nx = int(lines[1].split()[2])
Ny = int(lines[2].split()[2])
dx = float(lines[3].split()[2])
dy = float(lines[4].split()[2])

for line in lines:
    line = line.split()
    if len(line) == 7:
        x = (float(line[0]) - 0.5 * Nx) * dx
        y = (float(line[1]) - 0.5 * Ny) * dy
       	radius = (x*x + y*y) ** 0.5
       	r_a.append(radius)
        rho_a.append(float(line[2]))
        P_a.append(float(line[3]))
        v_x = float(line[5])
        v_y = float(line[6])
        if radius == 0.:
            vr_a.append(0.)
            vt_a.append(0.)
        else:
            vr_a.append((v_x * x + v_y * y) / radius)
            vt_a.append((-v_x * y + v_y * x) / radius)

anafile.close()

plt.plot(r, rho, "b.", markersize = 1, label = "numerical")
plt.plot(r_a, rho_a, "r.", markersize = 1, label = "analytical")
plt.xlabel("r")
plt.ylabel(r"$\rho$")
plt.legend()
#plt.ylim(0,0.1)
plt.xlim(0,1.0)
#plt.show()
plt.savefig("density_"+method+grid+"old.png")
plt.clf()

plt.plot(r, P, "b.", markersize = 1, label = "numerical")
plt.plot(r_a, P_a, "r.", markersize = 1, label = "analytical")
plt.xlabel("r")
plt.ylabel(r"P")
plt.legend()
#plt.show()
plt.savefig("pressure_"+method+grid+"old.png")
plt.clf()

plt.plot(r, vr, "b.", markersize = 1, label = "numerical")
plt.plot(r_a, vr_a, "r.", markersize = 1, label = "analytical")
plt.xlabel("r")
plt.ylabel(r"v$_{rad}$")
plt.legend()
#plt.show()
plt.savefig("radVel_"+method+grid+"old.png")
plt.clf()

plt.plot(r, vt, "b.", markersize = 1, label = "numerical")
plt.plot(r_a, vt_a, "r.", markersize = 1, label = "analytical")
plt.xlabel("r")
plt.ylabel(r"v$_{trans}$")
plt.legend()
#plt.show()
plt.savefig("transVel_"+method+grid+"old.png")
plt.clf()

plt.plot(r, B, "b.", markersize = 1, label = "numerical")
plt.xlabel("r")
plt.ylabel(r"B")
plt.legend()
#plt.show()
plt.savefig("magStr_"+method+grid+"old.png")
plt.clf()
