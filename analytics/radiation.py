import matplotlib.pyplot as plt
import numpy as np

err = []
nums = np.array([8,7,6,5,4,3])
bins = 10**nums
bin2 = bins**(-0.5)

for num in nums:

    files = 'out1e'+str(num)+'.dat'
    print "reading", files
    data = np.loadtxt(files, unpack = True)
    theta = data[0]
    iten = data[1]

    if(num==8) :
        theta_an = theta
        iten_an = iten
    
#    plt.plot(theta, iten, linewidth=2)
#    plt.xlim(0,6.28)
#    plt.xlabel(r"$\theta$", fontsize = 16)
#    plt.ylabel(r"$\rm I_{\nu}/B$", fontsize = 16)
#    plt.savefig("plot8.png")
#    plt.clf()

    diff = abs(iten-iten_an)
#    print diff
    err.append(np.max(diff))

 
fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(12,5))
ax0, ax1 = axes.flat

ax0.plot(bins,err)
ax0.set_xscale('log')
ax0.set_xlabel('log N')
ax0.set_ylabel(r"$\epsilon$")

ax1.plot(bin2, err)
ax1.set_xlabel(r"$N^{1/2}$")

#plt.show()
plt.tight_layout()
plt.savefig('scale.png')

