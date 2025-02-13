import numpy as np

#Altitude conversions
def hPa_to_hft(press):
    return 1453.6645 * (1 - (press/1013.25)**0.190284)

def hPa_to_km(press):
    return hPa_to_hft(press)*0.1/3.28084

def hft_to_hPa(alt):
    return 1013.25*((1-(alt/1453.6645))**(1/0.190284))

def km_to_hPa(alt):
    return hft_to_hPa(alt*3.28084/0.1)

# Saturation vapour pressure
def psat_liq(T):
    return 100*np.exp(-6096.9386/T+16.635794-0.02711198*T+1.673953e-5*T**2+2.433502*np.log(T))

def dpsat_liq(T):
    return psat_liq(T)*(-6096.9386*(-1/T**2)-0.02711198+1.673953e-5*2*T+2.433502/T)

def psat_ice(T):
    return 100*np.exp(-6024.5282/T+24.7219+0.010613868*T-1.3198825e-5*T**2-0.49382577*np.log(T))

def dpsat_ice(T):
    return psat_ice(T)*(-6024.5282*(-1/T**2)+0.010613868-1.3198825e-5*2*T-0.49382577/T)
    
def psat_ice2(T):
    return np.exp(43.494-6545.8/((T-273.15)+278))/((T-273.15)+868)**2

def rhl_to_rhi(rhl, T):
    sh = rhl*psat_liq(T)
    return sh/psat_ice(T)

def rhi_to_rhl(rhi, T):
    sh = rhi*psat_ice(T)
    return sh/psat_liq(T)


# SAC calcualtions
def sac_temp_slow(RH, p, eta=0.35):
    '''RH as fraction over liquid
    p in hPa'''
    cp = 1004
    eps = 0.622
    EI = 1.223 # WV eission index for kerosene
    Q = 43.2e6 # Lower heating value
    G = EI * p*100 * cp / (eps * Q * (1-eta))

    #Tf = -44.16 + 9.43*np.log(G-0.053) + 0.72*(np.log(G-0.053))**2
    Tf = -46.46 + 9.08*np.log(G-0.053) + 0.72*(np.log(G-0.053))**2
    Tf = 273.15+Tf

    Tt = Tf-0.001
    evalRH = (psat_liq(Tf)-G*(Tf-Tt))/psat_liq(Tt)
    while evalRH>RH:
        Tt -= 0.01
        evalRH = (psat_liq(Tf)-G*(Tf-Tt))/psat_liq(Tt)
    return Tt


def sac_temp(RH, p, eta=0.35):
    '''RH as fraction
    p in hPa'''
    cp = 1004
    eps = 0.622
    EI = 1.223 # WV eission index for kerosene
    Q = 43.2e6 # Lower heating value
    G = EI * p*100 * cp / (eps * Q * (1-eta))

    #Tf = -44.16 + 9.43*np.log(G-0.053) + 0.72*(np.log(G-0.053))**2
    Tf = -46.46 + 9.08*np.log(G-0.053) + 0.72*(np.log(G-0.053))**2
    Tf = 273.15+Tf

    Tt = Tf-0.001
    
    for i in range(6):
        Tt_new = Tt - (Tt-Tf + (psat_liq(Tf)-RH*psat_liq(Tt))/G)/(1-(RH/G)*dpsat_liq(Tt))
        Tt = Tt_new
    return Tt
