import numpy as np
from numpy.linalg import inv
from matplotlib import pyplot as plt

xmin=1.0
xmax=20.0
npoints=12
sigma=0.2
lx=np.zeros(npoints)
ly=np.zeros(npoints)
ley=np.zeros(npoints)
pars=[0.5,1.3,0.5]

from math import log
def f(x,par):
    return par[0]+par[1]*log(x)+par[2]*log(x)*log(x) #a + b*log(x) + c*(log(x))^2

from random import gauss
def getX(x):  # x = array-like
    step=(xmax-xmin)/npoints
    for i in range(npoints):
        x[i]=xmin+i*step #fills in array with x_i = x_0 + i*delta(x). for bins i think?
        
def getY(x,y,ey):  # x,y,ey = array-like
    for i in range(npoints):
        y[i]=f(x[i],pars)+gauss(0,sigma) #f(x) + gaus(0, sigma). y is f with a gaus smearing 
        ey[i]=sigma

# get a random sampling of the (x,y) data points, rerun to generate different data sets for the plot below
'''
getX(lx)
getY(lx,ly,ley)

fig, ax = plt.subplots()
ax.errorbar(lx, ly, yerr=ley)
ax.set_title("Pseudoexperiment")
fig.show() #I don't run interactive sessions so my worst enemy is .show() lol
'''

# *** modify and add your code here ***
nexperiments = 1000

par_a = np.zeros(nexperiments)
par_b = np.zeros(nexperiments)
par_c = np.zeros(nexperiments)
chi2 = np.zeros(nexperiments)
chi2_reduced = np.zeros(nexperiments)

for iexp in range(nexperiments):
    x = np.zeros(npoints)
    y = np.zeros(npoints)
    ey = np.zeros(npoints)
    getX(x)
    getY(x, y, ey)

    A = np.zeros((npoints, 3))
    for i in range(npoints):
        lx_i = log(x[i])
        A[i, 0] = 1.0
        A[i, 1] = lx_i
        A[i, 2] = lx_i**2

    W = np.diag(1.0 / (ey**2))
    cov = inv(A.T @ W @ A)
    best_fit = cov @ (A.T @ W @ y)

    par_a[iexp], par_b[iexp], par_c[iexp] = best_fit

    yfit = A @ best_fit
    chi2_val = np.sum(((y - yfit) / ey)**2)
    chi2[iexp] = chi2_val
    chi2_reduced[iexp] = chi2_val / (npoints - len(best_fit))

#I almost missed that I needed to make 2 plots, so I did reduced chi sq first
fig, axs = plt.subplots(2, 2, figsize=(8, 8))
plt.tight_layout()

axs[0, 0].hist2d(par_a, par_b, bins=60)
axs[0, 0].set_title('Parameter b vs a')
axs[0, 0].set_xlabel('a')
axs[0, 0].set_ylabel('b')

axs[0, 1].hist2d(par_a, par_c, bins=60)
axs[0, 1].set_title('Parameter c vs a')
axs[0, 1].set_xlabel('a')
axs[0, 1].set_ylabel('c')

axs[1, 0].hist2d(par_b, par_c, bins=60)
axs[1, 0].set_title('Parameter c vs b')
axs[1, 0].set_xlabel('b')
axs[1, 0].set_ylabel('c')

axs[1, 1].hist(chi2_reduced, bins=80)
axs[1, 1].set_title('Reduced $\chi^2$ distribution')
axs[1, 1].set_xlabel('$\chi^2_{red}$')
axs[1, 1].set_ylabel('Frequency')

fig.tight_layout()
fig.savefig("./redChi_py.png")

#now the other guys
fig2, axs2 = plt.subplots(2, 2, figsize=(8, 8))
plt.tight_layout()

axs2[0, 0].hist(par_a, bins=60)
axs2[0, 0].set_title('Parameter a distribution')
axs2[0, 0].set_xlabel('a')
axs2[0, 0].set_ylabel('Count')

axs2[0, 1].hist(par_b, bins=60)
axs2[0, 1].set_title('Parameter b distribution')
axs2[0, 1].set_xlabel('b')
axs2[0, 1].set_ylabel('Count')

axs2[1, 0].hist(par_c, bins=60)
axs2[1, 0].set_title('Parameter c distribution')
axs2[1, 0].set_xlabel('c')
axs2[1, 0].set_ylabel('Count')

axs2[1, 1].hist(chi2, bins=80)
axs2[1, 1].set_title('Normal $\chi^2$ distribution')
axs2[1, 1].set_xlabel('$\chi^2$')
axs2[1, 1].set_ylabel('Frequency')

fig2.tight_layout()
fig2.savefig("./chi_py.png")

# -----------------------------
input("hit Enter to exit")
