from scipy.special import comb
import numpy as np
import xyz_py as xyz
import math
import pickle
import matplotlib.pyplot as plt
from scipy.constants import e, epsilon_0, k

import matplotlib.font_manager
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from matplotlib import cm

#using data from an xyz file for a layered oxide material, we can find the free energy and chemical potential data 

plt.rcParams.update({'font.sans-serif': 'Arial'})


epsilon = 2.770


def epsilon(c, v, x, y, z):
    alpha = 0.03*c + 0.266*(x-v) + 0.544*y + 0.508*z + 2*2.75*1
    #molecular mass
    MM = 6.941*c + 58.693*(x-v) + 54.938*y + 58.933*z + 1*2*15.999
    #number density
    #density of NMC~ 2e6 per A^3
    prefac = 6.02e23*2.3e6/(1e30)
    numdensity = prefac/MM
    eps = (-3-8*alpha*numdensity*math.pi)/(-3+4*alpha*numdensity*math.pi)
    return eps

def NMC532_Colclasure20(y):
    OCV = lambda x : (5.314735633000300E+00 +
    -3.640117692001490E+03*x**14.0 + 1.317657544484270E+04*x**13.0
    - 1.455742062291360E+04*x**12.0 - 1.571094264365090E+03*x**11.0
    + 1.265630978512400E+04*x**10.0 - 2.057808873526350E+03*x**9.0
    - 1.074374333186190E+04*x**8.0 + 8.698112755348720E+03*x**7.0
    - 8.297904604107030E+02*x**6.0 - 2.073765547574810E+03*x**5.0
    + 1.190223421193310E+03*x**4.0 - 2.724851668445780E+02*x**3.0
    + 2.723409218042130E+01*x**2.0 - 4.158276603609060E+00*x +
    -5.573191762723310E-04*np.exp(6.560240842659690E+00*x**4.148209275061330E+01)
    )
    muR_ref = 0
    muR = (OCV(y)-muR_ref)
    return muR


# load data, set cutoff coordinates

labels, coords = xyz.load_xyz('../../LiNiO2_-10_10.xyz')

#center atom
j = 19590

#set NMC ratio
#x = 0.33333
#y = 0.33333
#z = 0.33333
x = 0.5
y = 0.3
z = 0.2
#x = 0.8
#y = 0.1
#z = 0.1
#x = 0.6
#y = 0.2
#z = 0.2

#define cutoff distance
cutoff_dist = 20

#define discretization of calculations
n_grid = 400

#initialize concentration
c = np.reshape(np.linspace(0,1,n_grid), (1,n_grid))
#drop edges
v = np.reshape(np.linspace(0,x,n_grid), (n_grid,1))
#this is the percentage of the concentration that


E_0 = np.zeros((n_grid,n_grid)) #this is multiple particle interactions
# first axis is vacancy concentration, while second axis is 

#find the neighboring points around the atom number j and j+1 midpoint
ref_atom = np.squeeze((coords[j]+coords[j+1])/2)
dip_direction = (coords[j+1]-coords[j])/np.linalg.norm(coords[j+1]-coords[j])

#Cm = 2.928e-10 #distance between dipoles
Cm = 2.103e-10 #distance between dipoles
Cm2 = 3.602e-10 #distance between dipoles
q = -1*e #TM charge

dip_a = (1.88*(x-2*v)+1.75*y+1.84*z-0.912*c)*e*Cm #dipole strength
#dip_a = (1.88*(x-2*v)+1.75*y+1.84*z-0.912)*e*Cm #dipole strength
#just coulomb charges
#dip_a = (c-2*v)*e*Cm
vec_a = np.multiply(np.expand_dims(dip_a,2),np.expand_dims(np.expand_dims(dip_direction, 0), 0))
#vec_a = dip_a*dip_direction
#process atoms

neighbor_coords_Li = []
neighbor_coords_Ni = []
for i in range(len(labels)):
    if i != j and i != j+1: #not ref atoms
        #split into original Li-Ni bonds or Ni-X bonds
        rij = (coords[i]-ref_atom)*1e-10
        if labels[i] == "Li":
            # "odd" LAYERS
            #drop oxygen atoms
            if np.linalg.norm(rij) <= cutoff_dist*1e-10:
                #add to neighbor list
                neighbor_coords_Li.append(coords[i])
                #original Li-Ni bonds broken and turned into Ni-NMC. also even plane interactions
                E_0 = E_0 - (1/(4*np.pi*epsilon(c,v,x,y,z)*epsilon_0)*np.dot(vec_a,rij)*q*(c-v)/np.linalg.norm(rij)**3/e)
 #               E_0 = E_0 - (1/(4*np.pi*epsilon_0)*np.dot(vec_a,rij)*q/np.linalg.norm(rij)**3/e)
        if labels[i] == "Ni":
            # "even" LAYERS
            #drop oxygen atoms
            if np.linalg.norm(rij) <= cutoff_dist*1e-10:
                #add to neighbor list
                neighbor_coords_Ni.append(coords[i])
                #original NMC bonds broken and turned into Li-Ni. also odd plane interactions
                E_0 = E_0 + (1/(4*np.pi*epsilon(c,v,x,y,z)*epsilon_0)*np.dot(vec_a,rij)*q*(c-v)/np.linalg.norm(rij)**3/e)
 #               E_0 = E_0 + (1/(4*np.pi*epsilon_0)*np.dot(vec_a,rij)*q/np.linalg.norm(rij)**3/e)
        #tack on extra terms

#E_0 = E_0 + (c-v)*(1-2*v+c)*(q**2/(16*np.pi*epsilon(c,v,x,y,z)*epsilon_0*e*Cm))
#E_0 = E_0 + (1-2*v+c)*(q**2/(16*np.pi*epsilon(c,v,x,y,z)*epsilon_0*e*Cm))
E_0 = E_0 + (q**2/(12*np.pi*epsilon(c,0,x,y,z)*epsilon_0*e))*(1/Cm - 1/Cm2)
#total enthalpic interactions in units of eV


# define plotting function
def plot_func(max_value, cmap):
    """Defines a clean plotting plot.
    max_value is the nickel concentration. cmap is the cmap to be used
    Returns the plt"""

    plt.style.use('classic')
    fig, ax = plt.subplots(figsize=(5.5,6))
        
    ax = plt.gca()
    for axis in ['top','bottom','left','right']:
        ax.spines[axis].set_linewidth(2.0)
    ax.tick_params(axis='x', direction='in')
    ax.tick_params(axis='y', direction='in')
    ax.xaxis.set_tick_params(width=1.5)
    ax.yaxis.set_tick_params(width=1.5)
    # Change the size of the fontsize
    axfontsize = 22
    
    ax.xaxis.set_ticks(np.arange(0, 1.1, 0.125))
    #ax.yaxis.set_ticks(np.arange(0, 11, 1.25))
    if max_value == 0.5:
        ax.yaxis.set_ticks(np.arange(0.5, -0.01, -0.125))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    elif max_value == 0.33333:
        ax_ticks = np.arange(0.3, -0.01, -0.075)
        ax_ticks[-1] = 0
        ax.yaxis.set_ticks(ax_ticks)
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    elif max_value == 0:
        ax.yaxis.set_ticks(np.arange(0, 5.25, 0.5))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%d'))
    elif max_value == -1:
        ax.yaxis.set_ticks(np.arange(2, 8.1, 0.5))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%d')) 
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    
    plt.set_cmap(cmap)
    
    for tick in ax.xaxis.get_major_ticks():
          tick.label1.set_fontsize(axfontsize)
    for tick in ax.yaxis.get_major_ticks():
          tick.label1.set_fontsize(axfontsize)
    
    for tick in ax.get_xticklabels():
          tick.set_fontname("Arial")
    for tick in ax.get_yticklabels():
          tick.set_fontname("Arial")
    
    every_nth_y = 2
    for n, label in enumerate(ax.yaxis.get_ticklabels()):
        if (n) % every_nth_y != 0:
            label.set_visible(False)
    
    every_nth_x = 2
    for n, label in enumerate(ax.xaxis.get_ticklabels()):
        if (n) % every_nth_x != 0:
            label.set_visible(False)
    
    plt.tight_layout()

    return plt

def cbar_plot(cbar):

    cbar.ax.set_frame_on(True)

    cbar.outline.set_linewidth(2)

    cbar.ax.tick_params(axis='y', direction='in')
    cbar.ax.yaxis.set_tick_params(width=1.5)
    
    axfontsize = 22
    
    for tick in cbar.ax.yaxis.get_major_ticks():
          tick.label2.set_fontsize(axfontsize)
    
    for tick in cbar.ax.get_yticklabels():
          tick.set_fontname("Arial")
    
    return

plt = plot_func(x, 'viridis')

plt.imshow(E_0,extent = [0,1,x,0])
#plt.ylabel("Disorder % of the System")
#plt.xlabel("Concentration")
#plt.title('Enthalpy (eV)')
plt.colorbar()
#plt.clim([0,4])
plt.savefig('Ni'+str(x)+'Mn'+str(y)+'Co'+str(z)+"_enthalpy.png", dpi = 300, bbox_inches = 'tight')
plt.show()

#entropy in units of kBT
entropy2 = -c*np.log(c)-(1-c-v)*np.log(1-c-v)-2*v*np.log(v)+x*np.log(x)-(x-v)*np.log(x-v)

plt = plot_func(x, 'plasma')
plt.imshow(0.0259*entropy2,extent = [0,1,x,0])
cb = plt.colorbar(ticks = np.arange(-0.001, 0.037, 0.009), fraction=0.024, pad=0.06)
cbar_plot(cb)
#plt.ylabel("Disorder % of the System")
#plt.xlabel("Concentration")
#plt.title('Entropy (eV/K)')
#plt.colorbar()
plt.savefig('Ni'+str(x)+'Mn'+str(y)+'Co'+str(z)+"_entropy.png", dpi = 300, bbox_inches = 'tight')
plt.show()

#total free energy in eV
H = E_0 -0.0257*entropy2
##enthalpic interactions in eV
H1 = E_0
plt = plot_func(x, 'viridis')
plt.imshow(H,extent = [0,1,x,0])
cb = plt.colorbar(ticks = np.arange(0,4.1,1),fraction=0.024, pad=0.06)
cbar_plot(cb)
#plt.ylabel("Disorder % of the System")
#plt.xlabel("Concentration")
#plt.title('Effective Hamiltonian (eV)')
#plt.clim([0,4])
plt.savefig('Ni'+str(x)+'Mn'+str(y)+'Co'+str(z)+"_FE.png", dpi = 300, bbox_inches = 'tight')
plt.show()


mu_v = (H[1:,:]-H[0:-1,:])/(v[1:,:]-v[0:-1,:])
##second term is only from the enthalpic forces
plt = plot_func(x, 'viridis')
mu_v_ent = (H1[1:,:]-H1[0:-1,:])/(v[1:,:]-v[0:-1,:])
plt.imshow(mu_v,extent = [0,1,x,0])
plt.colorbar()
#cb = plt.colorbar(ticks = np.arange(-3.6, -0.9, 0.6), fraction=0.024, pad=0.06)
#cbar_plot(cb)
#plt.clim([-1,-3.9])
#plt.ylabel("Disorder % of the System")
#plt.xlabel("Concentration")
#plt.title('Chemical Potential for Disorder (eV)')
#plt.colorbar()
plt.savefig('Ni'+str(x)+'Mn'+str(y)+'Co'+str(z)+"_muv.png", dpi = 300, bbox_inches = 'tight')
plt.show()


mu_c = (H[:,1:] - H[:,0:-1])/(c[:,1:] - c[:,0:-1])
plt = plot_func(x, 'viridis')
##second term is only from the enthalpic forces
mu_c_ent = (H1[:,1:] - H1[:,0:-1])/(c[:,1:] - c[:,0:-1])
plt.imshow(mu_c,extent = [0,1,x,0])
plt.colorbar()
#cb = plt.colorbar(ticks = np.arange(-0.5, 4.6, 1), fraction=0.024, pad=0.06)
#cbar_plot(cb)
#plt.clim([-0.5,4.5])
#plt.ylabel("Disorder % of the System")
#plt.xlabel("Concentration")
#plt.title('Chemical Potential for Intercalation (eV)')
#plt.colorbar()
plt.savefig('Ni'+str(x)+'Mn'+str(y)+'Co'+str(z)+"_muc.png", dpi = 300, bbox_inches = 'tight')
plt.show()


c_av = np.squeeze((c[:,:-1]+c[:,1:])/2)

#set reference potential
murRef = np.nanmax(mu_c[10,:])

plt = plot_func(0, 'viridis')
plt.plot(c_av, -(mu_c[10,:] - murRef), 'g', linewidth = 4, label = 'v ~ 0')
plt.plot(c_av, -(mu_c_ent[10,:] - murRef), 'g--', linewidth = 4, label = 'Enthalpy Contribution')
plt.plot(c_av, -(mu_c[10,:]-mu_c_ent[10,:]), 'g', linestyle = 'dotted', linewidth = 4, label = 'Entropy Contribution')
plt.plot(np.linspace(0,1,500), (NMC532_Colclasure20(np.linspace(0,1,500))), 'r', linewidth = 4, label = 'Colcalsure 2020 Experimental OCV')
#plt.xlabel('Concentration')
#plt.ylabel('Chemical Potential for Intercalation (eV)')
plt.legend(framealpha = 0, frameon = False, fontsize = 15, loc = 'lower center')
plt.savefig('Ni'+str(x)+'Mn'+str(y)+'Co'+str(z)+"_muc_single.png", dpi = 300, bbox_inches = 'tight')
plt.show()


v_av = (v[:-1]+v[1:])/2

max_ind = np.nanargmax(mu_v[:,-10])

plt.plot(v_av[:max_ind], -(mu_v[:max_ind,-10]-murRef), 'b', label = 'c ~ 1')
plt.plot(v_av[:max_ind], -(mu_v_ent[:max_ind,-10]-murRef), 'b--', label = 'Enthalpy Contribution')
#plt.plot(v_av[:max_ind], mu_v[:max_ind,-10]-mu_v_ent[:max_ind,-10], 'b', linestyle = 'dotted', label = 'Entropy Contribution')
plt.legend()
plt.xlabel('Vacancy Concentration')
plt.ylabel('Chemical Potential for Disorder (eV)')
plt.savefig('Ni'+str(x)+'Mn'+str(y)+'Co'+str(z)+"_muv_single.png", dpi = 300, bbox_inches = 'tight')
plt.show()


plt = plot_func(-1, 'viridis')
plt.plot(np.linspace(0,1,n_grid), -np.squeeze(mu_v[10,:]), 'b', linewidth = 4, label = 'v ~ 0')
plt.plot(np.linspace(0,1,n_grid), -np.squeeze(mu_v_ent[10,:]), 'b--', linewidth = 4, label = 'Enthalpy Contribution')
plt.plot(np.linspace(0,1,n_grid), -np.squeeze(mu_v[10,:]+mu_v_ent[10,:]), 'b', linestyle = 'dotted', linewidth = 4, label = 'Entropy Contribution')
plt.legend()
plt.legend(framealpha = 0, frameon = False, fontsize = 15, loc = 'upper center')
#plt.xlabel('Concentration')
#plt.ylabel('Chemical Potential for Disorder (eV)')
plt.savefig('Ni'+str(x)+'Mn'+str(y)+'Co'+str(z)+"_muv_single2.png", dpi = 300, bbox_inches = 'tight')
plt.show()

np.savetxt('Ni'+str(x)+'Mn'+str(y)+'Co'+str(z)+'_muc.txt', (mu_c[1:,:]+mu_c[:-1,:])/2)
np.savetxt('Ni'+str(x)+'Mn'+str(y)+'Co'+str(z)+'_muv.txt', (mu_v[:,1:]+mu_v[:,:-1])/2)
np.savetxt('Ni'+str(x)+'Mn'+str(y)+'Co'+str(z)+'_muc_ent.txt', (mu_c_ent[1:,:]+mu_c_ent[:-1,:])/2)
np.savetxt('Ni'+str(x)+'Mn'+str(y)+'Co'+str(z)+'_muv_ent.txt', (mu_v_ent[:,1:]+mu_v_ent[:,:-1])/2)
np.savetxt('Ni'+str(x)+'Mn'+str(y)+'Co'+str(z)+'_c_v.txt', np.hstack((np.reshape(c_av, (-1,1)), v_av)))

with open('neighbor_coords_Ni.p', 'wb') as fp:
    pickle.dump(neighbor_coords_Ni, fp)
with open('neighbor_coords_Li.p', 'wb') as fp:
    pickle.dump(neighbor_coords_Li, fp)
