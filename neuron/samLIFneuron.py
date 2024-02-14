import numpy as np
import matplotlib.pyplot as plt

# FUNCTIONS
#######################################
def getWidth(xpos):
    if xpos > 0 and xpos < L:
        return (xpos*((w2-w1)/L) + w1)
    else:
        print('\nError: domain wall outside bounds\n')
        return
    
def gaussian(x,mu,sig):
    return np.exp(-np.power(x - mu,2.)/(2*np.power(sig,2.)))

def getVJ(J):
    return -J*(g*mu_B*P)/(2*e*Ms*(1+alpha**2))

def plotNotchPos():
    notchline = np.arange(0,L+(1e-9),notch_L)
    plt.hlines(notchline*1e9,0,maxt*1e9,colors='r',linestyles='dotted',linewidth=0.9)
    
def getNotchH(xpos,index):
    notchline = np.arange(0,L+(1e-9),notch_L)
    tempH = 0
    k_notch = np.array([3.5,7,9,11])*1e4
    for i in range(0,len(notchline)):
        if xpos < notchline[i]:
            tempH += k_notch[index]*gaussian(xpos,notchline[i],notch_diam/2)
        else:
            tempH -= k_notch[index]*gaussian(xpos,notchline[i],notch_diam/2)
    return tempH 
#######################################

# geometry
L = 250e-9
w1 = 25e-9
w2 = 100e-9
d = 1.5e-9
wfixed = 10e-9
deltw = (w2-w1)/L

# flags
notch = 1
if w1 == w2:
    sloped = 0
else:
    sloped = 1

# material parameters
Ms = 8e5 # saturation magnetization
g = 2 # Lande factor
e = 1.602e-19 # electron charge
mu_B = 9.274e-24 # Bohr magneton
mu_0 = 4e-7*np.pi # vacuum permeability
GammaE = 1.7595e11 # electron gyromagnetic ratio
GammaLL = mu_0*GammaE # Gilbert gyromagnetic ratio
alpha = 0.05 # damping parameter
P = 0.7 # polarization
Xi = 0.05 
Aex = 1.3e-11 # exchange parameter
Ku = 5e5 # uniaxial anisotropy
Keff = Ku - (0.5*0.9*mu_0*(Ms**2)) # effective anisotropy
Delt = np.sqrt(Aex/Keff) # DW width param (Bloch)
DW_w = np.pi*Delt # actual DW width
if sloped == 1:
    mu_H = GammaLL*Delt*alpha # mobility for Walker breakdown
else:
    mu_H = GammaLL*Delt/alpha # mobility for Neel
xmin = wfixed + DW_w*0.5 # can choose to hard code for better fitting
xmax = L - wfixed - DW_w*0.4

# notch geometry
notch_L = 25e-9
notch_diam = 10e-9

# current input
I = np.array([-15,-25,-35,-45])*1e-6

# velocity biases (aka fitting variables)
Hslope = -(6.000171568956471e4*deltw**2 + 1.084982104400131e+04*deltw + 1.228766907095338e+03)
k_repel = 1000

# run
tsteps = 10000
maxt = 35e-9
t = np.linspace(0,maxt,tsteps+1)
tdelta = maxt/tsteps
x = np.zeros((len(t),len(I)))
x[0,:] = xmin
v = np.zeros_like(x)
for i in range(0,len(I)):
    for j in range(1,tsteps+1):
        notch_pos = np.floor(x[j-1,i]/notch_L)
        w = getWidth(x[j-1,i])
        J = I[i]/(d*w)
        v_J = getVJ(J)
        if notch == 1:            
            Hnotch = getNotchH(x[j-1,i],i)
        else:
            Hnotch = 0
        Heff = Hslope + Hnotch
        v[j,i] = v_J + mu_H*(Heff)
        x[j,i] = x[j-1,i] + v[j,i]*tdelta # next xpos
        if x[j,i] > xmax:
            x[j,i] = xmax
            v[j,i] = 0
    # for j in range(int(tsteps/2),len(t)):
    #     Heff = Hslope
    #     v[j,i] = mu_H*(Heff)
    #     x[j,i] = x[j-1,i] + v[j,i]*tdelta # next xpos
        # if x[j,i] < notch_pos*notch_L:
        #     x[j,i] = notch_pos*notch_L
        
    
# plotting
if notch == 0:
    tmm = np.genfromtxt('t15.csv',delimiter='\n')
    x15 = np.genfromtxt('x15.csv',delimiter='\n')
    x25 = np.genfromtxt('x25.csv',delimiter='\n')
    x35 = np.genfromtxt('x35.csv',delimiter='\n')
    x45 = np.genfromtxt('x45.csv',delimiter='\n')
    xmm = np.vstack((x15,x25,x35,x45)).T
    plt.plot(tmm,xmm)
    plt.title(r'Smooth constant current propagation')
else:
    t15notch = np.genfromtxt('t15notch.csv',delimiter='\n')
    t25notch = np.genfromtxt('t25notch.csv',delimiter='\n')
    t35notch = np.genfromtxt('t35notch.csv',delimiter='\n')
    t45notch = np.genfromtxt('t45notch.csv',delimiter='\n')
    x15notch = np.genfromtxt('x15notch.csv',delimiter='\n')
    x25notch = np.genfromtxt('x25notch.csv',delimiter='\n')
    x35notch = np.genfromtxt('x35notch.csv',delimiter='\n')
    x45notch = np.genfromtxt('x45notch.csv',delimiter='\n')
    plt.plot(t15notch,x15notch)
    plt.plot(t25notch,x25notch)
    plt.plot(t35notch,x35notch)
    plt.plot(t45notch,x45notch)
    plt.title(r'Notched constant current propagation')
label_ana = [r'$15$ $\mu A$',r'$25$ $\mu A$',r'$35$ $\mu A$',r'$45$ $\mu A$']
plt.legend(label_ana)
plt.plot(t*1e9,x*1e9,'k--')
if notch == 1:
    plotNotchPos()
plt.grid()
plt.ylim((0.250))
plt.xlim((0,35))
plt.xlabel(r'Time ($ns$)')
plt.ylabel(r'DW position ($nm$)')
plt.figure()
plt.plot(t[1:]*1e9,v[1:,:])
plt.grid()
plt.xlabel(r'Time ($ns$)')
plt.ylabel(r'DW velocity ($m/s$)')
plt.title(r'DW velocity :: $w_1=25$ $nm$ $w_2=100$ $nm$')
plt.legend(label_ana)