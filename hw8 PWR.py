import numpy as np
import matplotlib.pyplot as plt
from pyfluids import Fluid, FluidsList, Input

###variables###

H = 4 # m
He = 4.3 # m - extrapolated height
D_rod = 0.95 * 0.01 #m diameter
pitch = 1.26 * 0.01 #m distance between fuel pellets
D_fuel = 0.82 * 0.01#m diameter 
gap_thick = 0.006 * 0.01 #m space between fuel and clad
k_gap = 0.25 # W / m degC 
k_fuel = 3.6 # W / m degC
k_clad = 21.5 # W / m degC

r_fuel = D_fuel/2
r_clad_i = r_fuel + gap_thick
r_clad_o = D_rod/2

G = 4000 #kg/ m2 s or density times velocity
q0 = 380*100 # W/m linear initial heat generation
P_low = 15*(10**6) #MPa at z=-H/2
T_f_low = 277 #degC at z=-H/2
g = 9.8 #m/s2 acceleration due to grabity

wet_perim = np.pi * D_rod #m
area = pitch**2 - (np.pi*((D_rod/2)**2))
D_equiv = 4*area / wet_perim

water = Fluid(FluidsList.Water).unspecify_phase()
############################################

def tempQualProp(temp, qual, prop):
    """
    Parameters
    -----------
    temp: int
        Degrees celsius
    qual: int
        Number between 0 and 100 for flow quality
    prop: str
        one of the properties in the following list["vol", "intEnrg", "dynVisc", "enth", "rho"]
    """
    water_state = water.with_state(Input.temperature(temp), Input.quality(qual))
    if prop == "vol":
        return water_state.specific_volume #m3/kg
    elif prop == "intEnrg":
        return water_state.internal_energy #J/kg
    elif prop == "dynVisc":
        return water_state.dynamic_viscosity #Pa*s
    elif prop == "enth":
        return water_state.enthalpy #J/kg
    elif prop == "rho":
        return water_state.density #kg / m3

def presQualProp(pres, qual, prop):
    """
    Parameters
    -----------
    pres: int
        pressure in pascals
    qual: int
        Number between 0 and 100 for flow quality
    prop: str
        one of the properties in the following list["vol", "intEnrg", "dynVisc", "enth", "rho"]
    """
    water_state = water.with_state(Input.pressure(pres), Input.quality(qual))

    if prop == "vol":
        return water_state.specific_volume #m3/kg
    elif prop == "intEnrg":
        return water_state.internal_energy #J/kg
    elif prop == "dynVisc":
        return water_state.dynamic_viscosity #Pa*s
    elif prop == "enth":
        return water_state.enthalpy # J/kg
    elif prop =="rho":
        return water_state.density #kg / m3

def tempPresProp(temp, pres, prop):
    """
    Parameters
    -----------
    pres: int
        pressure in pascals
    temp: int
        Degrees celsius
    prop: str
        one of the properties in the following list["vol", "intEnrg", "dynVisc", "enth", "rho"]
    """
    water_state = water.with_state(Input.temperature(temp), Input.pressure(pres))

    if prop == "vol":
        return water_state.specific_volume #m3/kg
    elif prop == "intEnrg":
        return water_state.internal_energy #J/kg
    elif prop == "dynVisc":
        return water_state.dynamic_viscosity #Pa*s
    elif prop == "enth":
        return water_state.enthalpy # J/kg
    elif prop =="rho":
        return water_state.density #kg / m3

def calchfg(temp=0, pres=0):
    if temp != 0:
        return tempQualProp(temp, 100, "enth") - tempQualProp(temp, 0, "enth") #J/kg
    else:
        return presQualProp(pres, 100, "enth") - presQualProp(pres, 0, "enth") #J/kg

def reynolds(visc):
    return (G*D_equiv/visc) # Unitless

def frictionfactor(re):
    return 0.316 * re**(-0.25)

def deltaP(rho, f, drho):
    var = 0.5*f*(G**2/rho)*(wet_perim/area) + rho*g
    if drho != 0:
        return -1*(var + (G**2)*drho)
    else:
        return -1*var

def calcXe(pres, temp):
    h = tempPresProp(pres=pres, temp=temp, prop="enth")
    hfg = calchfg(pres = pres)
    hf = presQualProp(pres=pres, qual=0, prop="enth")
    #print(h, hfg, hf, type(hf))
    return (h-hf)/hfg

def calcNewTemp(Xe, pressure):
    hf = presQualProp(pres=pressure, qual=0, prop="enth")
    hfg = calchfg(pres=pressure)
    #print(Xe)
    #print(Xe*hfg)
    h = Xe*hfg + hf
    #print(Xe, h, hf, hfg)
    
    corresTemp = water.with_state(Input.enthalpy(h), Input.pressure(pressure)).temperature
    
    return corresTemp, h # degrees Celsius
    

#############################



###############################

option = 1 # 1 : PWR, 2: BWR

fig,ax = plt.subplots()

if option == 1 :

    for numOfPoints in [200]:#[5, 10, 20, 50,200]:
        
        ###########################################################

        Ps = [P_low]
        T_f_s = [T_f_low]
        Xes = [calcXe(P_low, T_f_low)]
        #print(Xes[0], calcNewTemp(Xes[0], Ps[0]))

        #print(f"P_i = {P_low}, T_f_i = {T_f_low}, Xe_i={Xes[0]}, hfg_i = {calchfg(pres=Ps[0])}")

        rhos = [tempPresProp(T_f_low, P_low, "rho")]
        mus = [tempPresProp(T_f_low, P_low, "dynVisc")]
        #print(f"rho_i = {rhos[0]}, mu_i = {mus[0]}")

        Res = [reynolds(mus[0])]
        frics = [frictionfactor(Res[0])]

        hs = [tempPresProp(T_f_low, P_low, "enth")]

        ##########################################################

        zs, deltaz = np.linspace(-H/2, H/2, numOfPoints, retstep=True)

        qlins = q0*np.cos(np.pi*(zs/He)) # W / m
        qdoubles = qlins/(np.pi*D_fuel)
        #print(qlins[:10])

        for i in range(1, len(zs)):
            #print(i)
            
        
            if i != 1 :
                drho = ((1/(rhos[i-1]))-(1/rhos[i-2])) * (1/(deltaz))
            else:
                drho = 0
            
            #print(f"drho = {drho}")
            

            deltaPCurrent = deltaP(rhos[0], frics[0], 0) #delta P / delta z using constant properties
            #print(f"delta P = {deltaPCurrent}")

            newP = deltaPCurrent*deltaz + Ps[i-1]
            #print(f"P = {newP}")
            Ps.append(newP)

            newhfg = calchfg(pres=Ps[i])
            #print(f"newhfg = {newhfg}")
            qdouble = qdoubles[i]

            #print(qdouble, qdouble*wet_perim/G*area)
            newXe = ((qdouble*wet_perim*deltaz)/(area*G*newhfg))+Xes[i-1]

            #print(f"hfg = {newhfg}, Xe = {newXe}")
            Xes.append(newXe)

            
            
            newT_f, newh = calcNewTemp(Xes[i], pressure=Ps[i])
            T_f_s.append(newT_f)
            hs.append(newh)

            rhos.append(tempPresProp(pres=Ps[i], temp=T_f_s[i], prop="rho"))
            #print(rhos[i])
            mus.append(tempPresProp(pres=Ps[i], temp=T_f_s[i], prop="dynVisc"))
            Res.append(reynolds(mus[i]))
            frics.append(frictionfactor(Res[i]))
            
            #print(f"P_{i} = {Ps[i]}, T_f_{i} = {T_f_s[i]}, Xe_{i}={Xes[i]}")
            #print(f"rho_{i} = {rhos[i]}, mu_{i} = {mus[i]}")

        #########################################
        ###Temperature inside the fuel rod#######
        

        T_c_in = []
        T_c_out = []
        T_f_c = []

        for j in range(len(zs)):
            curh = hs[j]
            curT_f = T_f_s[j]
            qlin = qlins[j]

            c_3 = (-qlin)/(2*np.pi*k_gap)
            c_5 = k_gap * c_3 / k_clad
            c_6 = (c_5*((-k_clad/(curh*r_clad_o)) - np.log(r_clad_o))) + curT_f
            c_4 = (c_3 *np.log(r_clad_i)*((k_gap/k_clad)-1))+c_6
            c_2 = (c_3*np.log(r_fuel)) + c_4 + (qlin/(4*np.pi*k_fuel))

            
            T_f_c.append(c_2)
            
            T_c_in.append(c_3*np.log(r_clad_i) + c_4)
            T_c_out.append(c_5*np.log(r_clad_o) + c_6)


        Tsats = []
        for press in Ps:
            Tsats.append(water.dew_point_at_pressure(press).temperature)

        #ax.plot(qlins, zs,  label="q\'")
        #ax.plot(np.array(Ps)/(10**6), zs, label="Fluid Pressure")

        ax.plot(T_c_out, zs, label="Outer Clad Surface Temperature")
        ax.plot(T_f_s, zs, label="Fluid Temperature")
        #ax.plot(T_f_c, zs, label="Centerline Temperature")        
        #ax.plot(Tsats, zs, label="Saturation Temperature")
        #ax.plot(T_c_in, zs, label="Inner Clad Surface Temperature")

        #ax.plot(Xes, zs, label=f"Equilibrium Quality")

        #ax.plot(rhos, zs, label="Densities")
        #ax.plot(mus, zs, label="Dynamic Viscosities")
        #ax.plot(Res, zs, label="reynolds")
        #ax.plot(frics, zs,  label = "frics")



#ax.set_xscale('log')
ax.set_ylabel(f"Vertical height [m]")
ax.set_xlabel(f"Temperature $^\circ$C")
ax.grid()
ax.legend()
plt.tight_layout()

plt.show()


