import numpy as np
import matplotlib.pyplot as plt
from pyfluids import Fluid, FluidsList, Input
import scipy.integrate as spi


MM_water = 18.02 #g/mol

water = Fluid(FluidsList.Water).unspecify_phase()
critP = water.critical_pressure

def tempQualProp(temp, qual, prop):
    """
    Parameters
    -----------
    temp: int
        Degrees celsius
    qual: int
        Number between 0 and 100 for flow quality
    prop: str
        one of the properties in the following list["vol", "intEnrg", "dynVisc", "enth", "rho", "pr", "k"]
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
    elif prop == "pr":
        return water_state.prandtl
    elif prop == "k":
        return water_state.conductivity
    
def presQualProp(pres, qual, prop):
    """
    Parameters
    -----------
    pres: int
        pressure in pascals
    qual: int
        Number between 0 and 100 for flow quality
    prop: str
        one of the properties in the following list["vol", "intEnrg", "dynVisc", "enth", "rho", "pr", "k"]
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
    elif prop == "pr":
        return water_state.prandtl
    elif prop == "k":
        return water_state.conductivity

def tempPresProp(temp, pres, prop):
    """
    Parameters
    -----------
    pres: int
        pressure in pascals
    temp: int
        Degrees celsius
    prop: str
        one of the properties in the following list["vol", "intEnrg", "dynVisc", "enth", "rho", "pr", "k"]
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
    elif prop == "pr":
        return water_state.prandtl
    elif prop == "k":
        return water_state.conductivity

def calchfg(temp=0, pres=0):
    if temp != 0:
        return tempQualProp(temp, 100, "enth") - tempQualProp(temp, 0, "enth") #J/kg
    else:
        return presQualProp(pres, 100, "enth") - presQualProp(pres, 0, "enth") #J/kg

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
    
def calcvfg(temp=0, pres=0):
    if temp != 0:
        return tempQualProp(temp, 100, "vol") - tempQualProp(temp, 0, "vol") #J/kg
    else:
        return presQualProp(pres, 100, "vol") - presQualProp(pres, 0, "vol") #J/kg

def qCritUniform(pressure, xe, G, hf, hin, Dh):
    pressure = pressure/(10**6) #convert to MPa
    hf = hf/(10**3) #convert to kJ/kg
    hin = hin/(10**3) #convert to kJ/kg
    term1 = ((2.022 -0.06238*pressure)+(0.1722-0.01427*pressure)*(np.e**((18.177-0.5987*pressure)*xe)))
    term2 = ((0.1484-1.596*xe + 0.1729*xe*np.abs(xe))*2.326*G + 3271)
    term3 = (0.2664 + 0.8357*(np.e**(-124.1*Dh)))*(0.8258+0.0003414*(hf-hin))
    return term1*term2*term3*1000 #W/m2


option = 2 # 1 : PWR, 2: BWR
plotted = ["Pressure"]

# Options for plotted: 
# "Linear Heat Generation", "Pressure", "Wall Temperature", "Fluid Temperature", 
# "Centerline Temperature", "Saturation Temperature", "Clad Inner Temperature", 
# "Equilibrium Quality", "Density", "Viscosity",  "Reynolds Number", "Friction Factor",
# "Heat Flux"


versions = [300]#[5, 10, 20, 100, 400] #number of points to calculate at over entire height

fig,ax = plt.subplots()


if option == 1 :


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

elif option == 2:
    ###variables###

    H = 4.1 # m
    He = 4.4 # m - extrapolated height
    D_rod = 1.227 * 0.01 #m diameter
    pitch = 1.62 * 0.01 #m distance between fuel pellets
    D_fuel = 1.04 * 0.01#m diameter 
    gap_thick = 0.010 * 0.01 #m space between fuel and clad
    k_gap = 0.25 # W / m degC 
    k_fuel = 3.6 # W / m degC
    k_clad = 21.5 # W / m degC

    r_fuel = D_fuel/2
    r_clad_i = r_fuel + gap_thick
    r_clad_o = D_rod/2

    G = 2350 #kg/ m2 s or density times velocity
    q0 = 605*100 # W/m linear initial heat generation
    P_low = 7.5*(10**6) #MPa at z=-H/2
    T_f_low = 272 #degC at z=-H/2
    g = 9.8 #m/s2 acceleration due to grabity

    wet_perim = np.pi * D_rod #m
    area = pitch**2 - (np.pi*((D_rod/2)**2))
    D_equiv = 4*area / wet_perim



def reynolds(visc):
    return (G*D_equiv/visc) # Unitless

def frictionfactor1phase(re):
    return 0.316 * re**(-0.25)

def frictionfactor2phase(re, pressure, chi):
    mug = presQualProp(pressure, 100, "dynVisc")
    muf = presQualProp(pressure, 0, "dynVisc")
    #mux = ((muf*chi + mug - chi*mug)/(muf*mug))**-1
    mux = chi*(mug-muf) + muf
    #print(re, muf/mux)
    return (0.046 * (re**(-0.2)) * (muf/mux)**(-0.2)), mux

def dittusBoelter(k, re, pr):
    return 0.023 * re**(0.8) * pr**0.4 * (k/D_equiv)

def liuWinterton(S, F, htc, hnb, tf, tsat, tw):
    temporary = (F*htc*(tw-tf))**2 + (S*hnb*(tw-tsat))**2
    return np.sqrt(temporary)
    
Ps = [P_low]
T_f_s = [T_f_low]
Xes = [calcXe(pres=Ps[0], temp=T_f_s[0])]
chis = [0]
mus = [tempPresProp(T_f_low, P_low, 'dynVisc')]
rhos = [tempPresProp(T_f_low, P_low, 'rho')]

Res = [reynolds(mus[0])]
frics = [frictionfactor1phase(Res[0])]

hs = [tempPresProp(T_f_low, P_low, "enth")]
hfs = [presQualProp(P_low, 0, "enth")]
hfgs = [calchfg(pres=P_low)]
Prs = [tempPresProp(T_f_low, P_low, "pr")]
ks = [tempPresProp(T_f_low, P_low, "k")]
htcs = [dittusBoelter(ks[0], Res[0], Prs[0])]

Tsats = [water.dew_point_at_pressure(P_low).temperature]


for numOfPoints in versions:

    ############################################
    
    zs, deltaz = np.linspace(-H/2, H/2, numOfPoints, retstep=True)

    qlins = q0*np.cos(np.pi*(zs/He)) # W / m
    qdoubles = qlins/(np.pi*D_fuel)
    
    ############################################

    c_3 = (-qlins[0])/(2*np.pi*k_gap)
    c_5 = k_gap * c_3 / k_clad
    c_6 = (c_5*((-k_clad/(htcs[0]*r_clad_o)) - np.log(r_clad_o))) + T_f_s[0]
    c_4 = (c_3 *np.log(r_clad_i)*((k_gap/k_clad)-1))+c_6
    c_2 = (c_3*np.log(r_fuel)) + c_4 + (qlins[0]/(4*np.pi*k_fuel))

    T_f_c = [c_2]
    T_c_in = [c_3*np.log(r_clad_i) + c_4]
    T_ws = [c_5*np.log(r_clad_o) + c_6]
    
    ############################################
    jump = 0
    for i in range(1, len(zs)):
        print(f"{i}, {Ps[-1]/(10**6):0.2e}MPa, {Tsats[-1]:0.2f}C, {T_ws[-1]:0.2f}C, {T_f_s[-1]:0.2f}C,  {mus[-1]}, {Xes[-1]:0.2e}, {chis[-1]:0.2e}, {frics[-1]:0.4e}")

    
        if Xes[-1] < 0:
            f1phase = frics[0]
            frics.append(f1phase)
            #rhos.append(tempPresProp(T_f_s[-1], Ps[-1], 'rho'))
            deltaPcurrent = -1*((2*f1phase*(G**2)/(rhos[-1]*D_equiv)) + (rhos[-1]*g))
            
        else:

            if jump==0 :
                jump = i

            Res.append(reynolds(mus[-1]))
            f2phase, mux = frictionfactor2phase(Res[-1], Ps[-1], chis[-1])
            mus.append(mux)
            frics.append(f2phase)
            vfg = calcvfg(pres=Ps[-1])
            vf  = presQualProp(Ps[-1], 0, "vol")
            deltaPcurrent = -1*(((G**2)*vfg*(chis[-1]-chis[-2])/deltaz) + (f2phase*(G**2)*(vfg*chis[-1] + vf)*2/D_equiv) + (g/(vfg*chis[-1]+vf)))
        
        P_new = deltaPcurrent*deltaz + Ps[-1]
        Ps.append(P_new)

        #print(P_new)
        Tsats.append(water.dew_point_at_pressure(Ps[-1]).temperature)
        
        hfgs.append(calchfg(pres=Ps[-1]))
        hfs.append(presQualProp(Ps[-1], 0, "enth"))
        newXe = (((qdoubles[i]*4*deltaz/(D_equiv*G))-(hfs[-1]-hfs[-2])-(Xes[-1]*(hfgs[-2]-hfgs[-1])))/hfgs[-1])+Xes[-1]
        Xes.append(newXe)
        newT_f, newh = calcNewTemp(Xes[-1], pressure=Ps[-1])
        
        if Xes[-1] < 0:
            if newT_f < Tsats[-1]:
                T_f_s.append(newT_f)
            else:
                T_f_s.append(Tsats[-1])
            hs.append(newh)
            htcs.append(dittusBoelter(ks[-1], Res[-1], Prs[-1]))
            c_3 = (-qlins[i])/(2*np.pi*k_gap)
            c_5 = k_gap * c_3 / k_clad
            if T_ws[i-1] <= Tsats[i-1]:
                c_6 = (c_5*((-k_clad/(htcs[i]*r_clad_o)) - np.log(r_clad_o))) + T_f_s[i]
            else:
                curP = Ps[i]
                F = (1+(chis[-1]*Prs[-1]*((presQualProp(curP, 0, "rho")/presQualProp(curP, 100, "rho"))-1)))**0.35
                S = (1+(0.055*(F**0.1)*(Res[0]**0.16)))**(-1)
                hnb = 55*((curP/critP)**0.12)*(qdoubles[i]**(2/3))*((-1*np.log10(curP/critP))**-0.55)*(MM_water**-0.5)
                '''tw = Tsats[i]
                #print(S, F, hnb, tw, T_f_s[i], Tsats[i], liuWinterton(S, F, htcs[i], hnb, T_f_s[i], Tsats[i], tw))
                
                while (qdoubles[i] - liuWinterton(S, F, htcs[i], hnb, T_f_s[i], Tsats[i], tw) >=0.0001):
                    tw = tw*1.001'''

                A1 = (F**2) * (htcs[i]**2)
                A2 = (S**2) * (hnb**2) #If F and S are flipped there is a discontinuity in the graph
                L1 = A1+A2
                L2 = -2*(A1*T_f_s[i] + A2*Tsats[i])
                L3 = A1*(T_f_s[i]**2) + A2*(Tsats[i]**2) - ((k_clad*c_5/r_clad_o)**2)
                #print(f"F = {F}, S = {S:0.2f}, hnb={hnb:0.2f}, A1 = {A1:0.2f}, A2 = {A2:0.2f}, L1 = {L1:0.2f}, L2 = {L2:0.2f}, L3 = {L3:0.2f}")
                #T_w_min = ((-L2-((L2**2 - (4*L1*L3))**0.5))/(2*L1))
                T_w_plus = ((-L2 +((L2**2 - (4*L1*L3))**0.5))/(2*L1))
                #print("Wall temp", T_w_plus)
                c_6 = T_w_plus - (c_5*np.log(r_clad_o))



            c_4 = (c_3 *np.log(r_clad_i)*((k_gap/k_clad)-1))+c_6
            c_2 = (c_3*np.log(r_fuel)) + c_4 + (qlins[i]/(4*np.pi*k_fuel))
        else:
            T_f_s.append(Tsats[-1])
            hs.append(newh)
            htcs.append(dittusBoelter(ks[-1], Res[-1], Prs[-1]))
            c_3 = (-qlins[i])/(2*np.pi*k_gap)
            c_5 = k_gap * c_3 / k_clad
            curP = Ps[i]

            F = (1+(chis[-1]*Prs[-1]*((presQualProp(curP, 0, "rho")/presQualProp(curP, 100, "rho"))-1)))**0.35
            S = (1+(0.055*(F**0.1)*(Res[0]**0.16)))**(-1)
            hnb = 55*((curP/critP)**0.12)*(qdoubles[i]**(2/3))*((-1*np.log10(curP/critP))**-0.55)*(MM_water**-0.5)
            '''tw = Tsats[i]
            #print(S, F, hnb, tw, T_f_s[i], Tsats[i], liuWinterton(S, F, htcs[i], hnb, T_f_s[i], Tsats[i], tw))
            
            while (qdoubles[i] - liuWinterton(S, F, htcs[i], hnb, T_f_s[i], Tsats[i], tw) >=0.0001):
                tw = tw*1.001'''

            A1 = (F**2) * (htcs[i]**2)
            A2 = (S**2) * (hnb**2) #If F and S are flipped there is a discontinuity in the graph
            L1 = A1+A2
            L2 = -2*(A1*T_f_s[i] + A2*Tsats[i])
            L3 = A1*(T_f_s[i]**2) + A2*(Tsats[i]**2) - ((k_clad*c_5/r_clad_o)**2)
            #print(f"F = {F}, S = {S:0.2f}, hnb={hnb:0.2f}, A1 = {A1:0.2f}, A2 = {A2:0.2f}, L1 = {L1:0.2f}, L2 = {L2:0.2f}, L3 = {L3:0.2f}")
            #T_w_min = ((-L2-((L2**2 - (4*L1*L3))**0.5))/(2*L1))
            T_w_plus = ((-L2 +((L2**2 - (4*L1*L3))**0.5))/(2*L1))
            #print("Wall temp", T_w_plus)
            c_6 = T_w_plus - (c_5*np.log(r_clad_o))

            c_4 = (c_3 *np.log(r_clad_i)*((k_gap/k_clad)-1))+c_6
            c_2 = (c_3*np.log(r_fuel)) + c_4 + (qlins[i]/(4*np.pi*k_fuel))
        
        T_f_c.append(c_2)
        T_c_in.append(c_3*np.log(r_clad_i) + c_4)
        T_ws.append(c_5*np.log(r_clad_o) + c_6)

        if Xes[-1] < 0:
            chis.append(0)
        elif Xes[-1] > 1:
            chis.append(1)
        else:
            chis.append(Xes[-1])

    qcrits = []
    dnbrs = []
    intercept = None
    for i in range(len(zs)):
        qCrit = (qCritUniform(Ps[i], Xes[i], G, hfs[i], hs[0], D_equiv))
        F = 1
        qcrits.append(qCrit/F)
        dnbrs.append(qCrit/(qdoubles[i]*F))
    



if "Linear Heat Generation" in plotted:   
    ax.plot(qlins, zs,  label="q\'")
if "Heat Flux" in plotted:   
    ax.plot(qdoubles, zs,  label="q\'\'")
if "Pressure" in plotted:
    ax.plot(np.array(Ps)/(10**6), zs, label="Fluid Pressure")
if "Wall Temperature" in plotted:
    ax.plot(T_ws, zs, label="Outer Clad (Wall) Surface Temperature")
if "Fluid Temperature" in plotted:
    ax.plot(T_f_s, zs, label="Fluid Temperature")
if "Centerline Temperature" in plotted:
    ax.plot(T_f_c, zs, label="Centerline Temperature")   
if "Saturation Temperature" in plotted:
    ax.plot(Tsats, zs, label="Saturation Temperature")
if "Clad Inner Temperature" in plotted:
    ax.plot(T_c_in, zs, label="Inner Clad Surface Temperature")
if "Equilibrium Quality" in plotted:
    ax.plot(Xes, zs, label=f"$\Delta z$ = {deltaz:0.2f} m")

if "Density" in plotted:
    ax.plot(rhos, zs, label="Densities")
if "Viscosity" in plotted:
    ax.plot(mus, zs, label="Dynamic Viscosities")
if "Reynolds Number" in plotted:
    ax.plot(Res, zs, label="Reynolds Number")
if "Friction Factor" in plotted:
    ax.plot(frics, zs,  label = "Friction Factor")
if "Critical Heat Flux" in plotted:
    ax.plot(qcrits, zs,  label = "Critical Heat Flux")
if "DNBR" in plotted:
    ax.plot(dnbrs, zs,  label = "DNBR")


#ax.set_xscale('log')
#ax.axhline(zs[jump], 0, 1, color='k', linestyle = '--', linewidth= 1)
ax.set_ylabel(f"Vertical height [m]")
ax.set_xlabel(f"{plotted[-1]}")
ax.grid()
ax.legend()
plt.tight_layout()

plt.show()
