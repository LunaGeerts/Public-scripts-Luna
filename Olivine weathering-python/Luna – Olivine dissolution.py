# -*- coding: utf-8 -*-
"""
Created on Wed Apr 15 12:42:36 2020

@author: Luna Geerts
"""


#%% Importing liberaries

import numpy as np
from numpy import matlib as ml
from scipy import constants
import matplotlib.pyplot as plt
from scipy.integrate import odeint
import pandas as pd


#making a nutrient limited chemostat from the reactran excercises.

#model domain: undefined / no saturation occurs
#time frame: months 


##############################################
# Functions
##############################################

#---------------------------------------------
# Conversion of rate constant depending on temperature (supplm. Montserat)/ahrrenius equation
# Turned out to be unneeded as I used an empirical relationship which includes T, left in as reference
#---------------------------------------------

Ea           = 79.5 *1000    # [ j mol^-1 ]
constants.R                  # [ J K^-1 mol^-1]
Ref_T1       = constants.convert_temperature(25,'C','K' ) #[K]
bulk_dens    = 1750000       # [ gram m^-3]

#%% function that calculates new k values
def new_k (k2, New_T2 ) :
    K_T2 = constants.convert_temperature( New_T2,'C','K' )
    exp1 = Ea / (K_T2 * constants.R)
    exp2 = Ea / (Ref_T1 * constants.R)
    k1   = k2 / (ml.exp(exp1)/ml.exp(exp2) )
    return k1

new_k(1,15) #so the reaction (r) is 0.32 than that of at 25 degrees

T = np.linspace(0,30,20)
rate = new_k(1,T)
new_k(1,10)
plt.grid()
plt.plot(T,rate,'r')
plt.xlabel("Temperature (Celsius)")
plt.ylabel("Rate")
title = "Rate increase with temperature"

plt.title(title)
plt.show()
#%%

#Function that checks if input is numeric

def is_num(f):
    try:
        float(f)
        return True
    except ValueError:
        return False







##############################################
#%% Parameters
##############################################



#---------------------------------------------
# Sea temperatures per month , averages over past 10 yr
#---------------------------------------------

T_hawaii = np.array([25.,25.,25.,25.,25.,26.,26.,26.,27.,27.,26.,25.]) # Kahaluu-Keauhou
T_Ostend = np.array([7.,6.,7.,9.,12.,16.,18.,19.,17.,15.,12.,9.])

months   = np.linspace([0.],[1.],12) #converting months on a year basis


poly_months   = np.linspace([0.],[1.],100) #converting months on a year basis

#making a function/ curve from this data

months_5 = np.linspace([0.],[5.],12*5)
curve_fit_Ostend = np.polyfit(months[:,0],T_Ostend, 5)
curve_fit_Hawaii = np.polyfit(months[:,0],T_hawaii, 6)




def T_funct_Ostend (t):
    if t <= 1:
        return(np.polyval(curve_fit_Ostend,t))
    elif t > 1:
        new_t   = t-np.trunc(t) #I only got values for t between 0 and 1
        return(np.polyval(curve_fit_Ostend,new_t))        

def T_funct_Hawaii (t):
    if t <= 1:
        return(np.polyval(curve_fit_Hawaii,t))
    elif t > 1:
        new_t   = t-np.trunc(t) #I only got values for t between 0 and 1
        return(np.polyval(curve_fit_Hawaii,new_t))        

T_funct_vec_Ostend = np.vectorize(T_funct_Ostend)
T_funct_vec_Hawaii = np.vectorize(T_funct_Hawaii)


plt.plot(months,T_Ostend, '.')
plt.plot(months,T_funct_vec_Ostend(months), 'b', label='Ostend')
plt.plot(months,T_hawaii, 'r.')
plt.plot(months,T_funct_vec_Hawaii(months), 'r', label='Hawaii')


plt.xlabel("Years")
plt.ylabel("Sea-Temperature")
title = "Yearly variation of temperatures"
plt.grid()
plt.title(title)
plt.legend(loc='lower right')
plt.show()

#Bingo! :)

#---------------------------------------------
# Olivine constants
#---------------------------------------------

Mol_Mass = 140.693 # [ gram mol^-1] #From Rimstidt et al.,
Mol_Dens_OLI = 4.365*10**-5 #[ m^3 mol^-1 ] also Rimstidt
dia      = 100 *10**-6     # [ meter]
pH       = np.linspace([0.],[12.])


C_OLI    = 10 / Mol_Mass   #lets say per cubic meter of seafloor we want 10 grams?
gamma_eff= 4 #we assume no secondary reactions occur so the gamma efficiency is 4 aka 4 mols CO2 capture per mol olivine

#---------------------------------------------
#%% Olivine dissolution function in function of pH and temperature
#---------------------------------------------


def r_OLI (pH,T_K): #the pH and Kelvin temperatures are needed
    #for pH ranges <5.6
    if pH <= 5.6 :
        log_r_out = 6.05-0.46*pH-3683.0*(1/T_K)
    #for pH ranges >5.6   
    elif pH > 5.6:
        log_r_out = 4.07-0.256*pH-3465*(1/T_K)
        
    r_out = log_r_out
    return r_out #r in units [ mol m^2 year^-1] 



#---------------------------------------------
# Vectorizing the function
#---------------------------------------------

r_OLI2 = np.vectorize(r_OLI)

#---------------------------------------------
#%% Surface area and volume sphere for shrinking particle
#---------------------------------------------

#Assuming our olivine particles are perfectly round, knowing the radius of each particle
# we can calculate the surface area with A_sphere

def A_sphere(diameter):
    return(4 * ml.pi * (diameter/2)**2)

def V_sphere(diameter):
    return(4/3 * ml.pi * (diameter/2)**3)

#%%
#---------------------------------------------
# Making an identical plot as in the reference paper -> DID NOT USE IN THE END
#---------------------------------------------
pH= np.linspace(3,12,10)

plt.plot(pH,r_OLI2(pH,298.15)) #compared to hangx and spiers paper its way too high...
plt.show()

plt.plot(pH,ml.log(ml.exp(r_OLI2(pH,298.15))/10000))  #making a plot like in Oelkers, I do / 10000 to convert
# m^2 to cm^2

plt.plot(pH,(ml.log(2)/r_OLI2(pH,298.15)) ) #compared to hangx and spiers paper its way too high...


#rendering
plt.xlabel("pH")
plt.ylabel("ln 2 / log r")
title = "Dissolution in function of pH – DID NOT USE"
plt.grid()
plt.title(title)

plt.show()


#---------------------------------------------
# CO2 capture in function of Olivine dissolution
#---------------------------------------------


#MAIN FUNCTION

def decr_OLI(y,t,r,nr_grains,Temp_func,uncertainty,purity):
    
    Vol_grain = y[1]/ nr_grains
    Vol_grain_up = y[3]/ nr_grains
    Vol_grain_down = y[5]/ nr_grains
    
    diameter  =(Vol_grain/((4/3) * ml.pi))**(1/3)*2
    diameter_up=(Vol_grain_up/((4/3) * ml.pi))**(1/3)*2
    diameter_down=(Vol_grain_down/((4/3) * ml.pi))**(1/3)*2
    
    A_grain   = A_sphere(diameter) #representative for one grain but...
    A_grain_up= A_sphere(diameter_up)
    A_grain_down=A_sphere(diameter_down)
    #we need something in units [mu_m^2 g^-1]
    
    A     =  A_grain / Vol_grain / (bulk_dens *10**-18)    # [mu_m^2 g^-1]
    A_up  =  A_grain_up / Vol_grain_up / (bulk_dens *10**-18)    # [mu_m^2 g^-1]
    A_down =  A_grain_down / Vol_grain_down / (bulk_dens *10**-18)    # [mu_m^2 g^-1]

    
    
    #Part of code that finds the correct temperature for the timestep:
        
    if Temp_func == "Ostend":
        Temp = T_funct_vec_Ostend(t)
    if Temp_func == "Hawaii":
        Temp = T_funct_Hawaii(t)
    if is_num(Temp_func):
        Temp = Temp_func

    #since the Rimstidt function did not work as intended we will use values of
    #Hangx and Spiers and convert these values using the Ahrenius equation
    #Detailed as in Hangx and Spiers:
        
    r_adj =        (r*3600*24*365.25)        *new_k(1,Temp) #adjust rate for temperature
    r_adj_up   =   (r + np.abs(uncertainty)) *3600*24*365.25 * new_k(1,Temp) #adjust rate for temperature
    r_adj_down =   (r - np.abs(uncertainty)) *3600*24*365.25 * new_k(1,Temp) #adjust rate for temperature
    
    r_new      =  (r_adj*10**-12)  * Mol_Mass*10**-6 # [ mol mu_m^-2 time^-1 * g mol^-1] = [ g mu_m^-2 year^-1]
    r_new_up   =  (r_adj_up*10**-12)  * Mol_Mass*10**-6
    r_new_down =  (r_adj_down*10**-12)  * Mol_Mass*10**-6
    
    
    dCdt      = (-(A*y[0]*r_new)*10**6 )*purity   # [ A = mu_m^2 g^-1, y[0]= g m^-2(seafloor), r= g mu_m^-2 year^-1]
                                   # [dCdt = gram m^-2 (seafloor) year^-1]
    
    dCdt_Up   = (-(A_up*y[2]*r_new_up)*10**6 )*purity 
    dCdt_Down = (-(A_down*y[4]*r_new_down)*10**6 )*purity 
    
    #converting gram olivine dissolved to volume change
    
    if y[1]  < 1.05*10**11 : #its a messy solution... but a solution nontheless!
         dVoldt = 0  
    else :
        dVoldt = dCdt   / (bulk_dens*10**-18) # [ mu_m^3 m^-2(seafloor) year^-1 ]
    
    if y[3]  < 1.05*10**11 : #its a messy solution... but a solution nontheless!
         dVoldt_Up = 0  
    else :
        dVoldt_Up = dCdt_Up   / (bulk_dens*10**-18) # [ mu_m^3 m^-2(seafloor) year^-1 ]

    if y[5]  < 1.05*10**11 : #its a messy solution... but a solution nontheless!
         dVoldt_Down = 0  
    else :
        dVoldt_Down = dCdt_Down   / (bulk_dens*10**-18) # [ mu_m^3 m^-2(seafloor) year^-1 ]


    
    #Since there is an uncertainty of a magnitude I want to plot these too
    
    #dVoldt  = 0
        
    dydt     = [dCdt, dVoldt, dCdt_Up, dVoldt_Up ,dCdt_Down,dVoldt_Down]
    return dydt


#WRAPPER FUNCTION FOR MAIN FUNCTION, handles the input !

def OLI_model2 (t,r,gram_oli,dia,Temp_func=25,uncertainty=0,purity=1):


        
    muVol     = (gram_oli / bulk_dens) *10**18 # [mu_m^3 m^-2] volume per m^2 sea floor
   
    nr_grains  = muVol / V_sphere(dia) #total volume per m^2 sea floor divided by volume of one grain
                                     #gives nr of grains
    y = [gram_oli , muVol , gram_oli, muVol, gram_oli, muVol ]
    
    sol_tot  = odeint(decr_OLI, y, t, args=(r,nr_grains,Temp_func,uncertainty,purity))

    return sol_tot


#==============================================================================
# Trying to give a ballpark estimate of C_OLI decrease
#==============================================================================

r = ml.exp(r_OLI(pH=8.2,T_K=298.15))# [mol m^-2 s^-1]
#but this rate is too fast I think? since Hangx and spiers (other paper uses other values...)
#but the rate calculated seem to be consistent with the papers tho... -> DID NOT USE

r = 1.58*10**-10 #mol m^-2 s^-1 #but this rate by Hangx and Spiers is used
error_r = 1.40*10**-10

acid_r = 1.80*10**-10 #mol m^-2 s^-1, dissolution at pH 7.8

#The time windows:

t1 = np.linspace(0,1000,1000)
t2 = np.linspace(0,1000,1000)
t3 = np.linspace(0,1000,1000)


#==============================================================================
# Computing different grain sizes...
#==============================================================================

sol1 = OLI_model2(t=t1,r=r,gram_oli=10,dia=10,uncertainty=error_r,purity=1)
sol2 = OLI_model2(t=t2,r=r,gram_oli=10,dia=100,uncertainty=error_r,purity=1)
sol3 = OLI_model2(t=t3,r=r,gram_oli=10,dia=1000,uncertainty=error_r,purity=1)

df_sol1 = pd.DataFrame(sol1[:,(2,0,4)],columns=['lwr','Rt','upr'],index=t1)
df_sol2 = pd.DataFrame(sol2[:,(2,0,4)],columns=['lwr','Rt','upr'],index=t2)
df_sol3 = pd.DataFrame(sol3[:,(2,0,4)],columns=['lwr','Rt','upr'],index=t3)


plt.plot(df_sol1.index, df_sol1.Rt  ,'b',label='10 micrometer diameter')
plt.fill_between(df_sol1.index,df_sol1.lwr,df_sol1.upr,alpha=0.1,facecolor='blue')

plt.plot(df_sol2.index, df_sol2.Rt  ,'g',label='100 micrometer diameter')
plt.fill_between(df_sol2.index,df_sol2.lwr,df_sol2.upr,alpha=0.1,facecolor='green')

plt.plot(df_sol3.index, df_sol3.Rt  ,'r',label='1000 micrometer diameter')
plt.fill_between(df_sol3.index,df_sol3.lwr,df_sol3.upr,alpha=0.1,facecolor='red')


plt.legend(loc='best')

plt.grid()

#rendering
plt.xlabel("Years")
plt.ylabel("Gram Olivine per square meter seafloor")
title = "Dissolution with different starting diameters"
plt.title(title)

plt.show()


#trying to make a graph such as we find in Hangx and spiers:


plt.plot(t1, 100-sol1[:,0]/10*100, 'b', label='10 micrometer diameter')
plt.plot(t2, 100-sol2[:,0]/10*100, 'g', label='100 micrometer diameter')
plt.plot(t3, 100-sol3[:,0]/10*100, 'r', label='1000 micrometer diameter')

plt.legend(loc='best')
plt.grid()
plt.xlabel("Years")
plt.ylabel("% Olivine dissolved")
title = "Dissolution with different starting diameters"
plt.title(title)

plt.show()
    

#==============================================================================
# For different pH scenarios, with different diameters
#==============================================================================

t1_acid = np.linspace(0,100,1000)

sol1_acid = OLI_model2(t=t1_acid,r=acid_r,gram_oli=10,dia=10,uncertainty=error_r,purity=1)
sol2_acid = OLI_model2(t=t2,r=acid_r,gram_oli=10,dia=100,uncertainty=error_r,purity=1)
sol3_acid = OLI_model2(t=t3,r=acid_r,gram_oli=10,dia=1000,uncertainty=error_r,purity=1)

df_sol1_acid = pd.DataFrame(sol1_acid[:,(2,0,4)],columns=['lwr','Rt','upr'],index=t1_acid)
df_sol2_acid = pd.DataFrame(sol2_acid[:,(2,0,4)],columns=['lwr','Rt','upr'],index=t2)
df_sol3_acid = pd.DataFrame(sol3_acid[:,(2,0,4)],columns=['lwr','Rt','upr'],index=t3)

#plotting:






plt.plot(df_sol1.index, df_sol2.Rt  ,'b--')
plt.plot(df_sol2.index, df_sol2.Rt  ,'g--')
plt.plot(df_sol3.index, df_sol3.Rt  ,'r--')
    


plt.plot(df_sol1_acid.index, df_sol2_acid.Rt  ,'b',label='10 micrometer diameter - pH 7.8')
plt.fill_between(df_sol1_acid.index,df_sol2_acid.lwr,df_sol1_acid.upr,alpha=0.1,facecolor='blue')


plt.plot(df_sol2_acid.index, df_sol2_acid.Rt  ,'g',label='100 micrometer diameter - pH 7.8')
plt.fill_between(df_sol2_acid.index,df_sol2_acid.lwr,df_sol2_acid.upr,alpha=0.1,facecolor='green')


plt.plot(df_sol3_acid.index, df_sol3_acid.Rt  ,'r',label='1000 micrometer diameter - pH 7.8')
plt.fill_between(df_sol3_acid.index,df_sol3_acid.lwr,df_sol3_acid.upr,alpha=0.1,facecolor='red')
plt.plot(0,0,'k--',label='pH 8.2 - reference')
plt.legend(loc='best')

plt.grid()

#rendering
plt.xlabel("Years")
plt.ylabel("Gram Olivine per square meter seafloor")
title = "Dissolution in an acidified ocean"
plt.title(title)

plt.show()


#==============================================================================
# Trying out the function for the hawaii and Ostend scenario
#==============================================================================

sol_Ostend = OLI_model2(t=t2,r=r,gram_oli=10,dia=100,Temp_func="Ostend",uncertainty=error_r,purity=1)
sol_Hawaii =  OLI_model2(t=t2,r=r,gram_oli=10,dia=100,Temp_func="Hawaii",uncertainty=error_r,purity=1)

df_solOst = pd.DataFrame(sol_Ostend[:,(2,0,4)],columns=['lwr','Rt','upr'],index=t2)
df_solHai = pd.DataFrame(sol_Hawaii[:,(2,0,4)],columns=['lwr','Rt','upr'],index=t2)



plt.plot(df_solOst.index, df_solOst.Rt  ,'b',label='Ostend')
plt.fill_between(df_solOst.index,df_solOst.lwr,df_solOst.upr,alpha=0.1,facecolor='blue')

plt.plot(df_solHai.index, df_solHai.Rt  ,'r',label='Hawaii')
plt.fill_between(df_solHai.index,df_solHai.lwr,df_solHai.upr,alpha=0.1,facecolor='red')


plt.legend(loc='best')

plt.grid()


#rendering
plt.xlabel("Years")
plt.ylabel("Gram Olivine per square meter seafloor")
title = "Dissolution of olivine (d=100 micrometer) at different locations"
plt.title(title)

plt.show()












