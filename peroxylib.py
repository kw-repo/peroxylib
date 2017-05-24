"This is the only actual version. Anything other is not the workfile."
"Peroxylib v.0.4.001"
"Rearrange shit, especially vap_comp_rT, which, obviously, sucks fucking ass at beforelast points of tables"

from math import *
__name__ = 'peroxylib'

############## Math functions ###############

def parabola(X_, Y_):
    
    if (X_[0] == X_[1]) | (X_[1] == X_[2]) | (X_[0] == X_[2]):
        print "Error calculating parabolic coefficients: curve is not a function"
        return [0, 0, 0]
    
    X = [0.,0.,0.]
    Y = [0.,0.,0.]    
    
    X[0] = float(X_[0])
    X[1] = float(X_[1])
    X[2] = float(X_[2])
    
    Y[0] = float(Y_[0])
    Y[1] = float(Y_[1])
    Y[2] = float(Y_[2])
            
    
    
    A = Y[2] - Y[1] - (Y[1] - Y[0]) * (X[2] - X[1]) / (X[1] - X[0])
    A = A / ((X[2] - X[1]) * (X[2] - X[0]))

    B = A * (X[0]*X[0] - X[1]*X[1]) + Y[1] - Y[0]
    B = B / (X[1] - X[0])
    
    C = Y[0] - A * X[0]*X[0] - B * X[0]
    
    return (A, B, C)
    
    
def sqf(coeff, x):
    
    return coeff[0]*x*x + coeff[1]*x + coeff[2]


def parabolic(points, x):
    
    X_ = points[0]
    Y_ = points[1]
    X = X_
    Y = Y_    

    l = len(X_)
    N = l - 1
    
    for i in range(0, l):
        X[i] = float(X_[i])
        Y[i] = float(Y_[i])
    
    p = 1
    
    while x > X[p]:
        if p + 1 > N:
            break
        p = p + 1
    p = p - 1
    
    
    B = [0, 0, 0]
    F = [0, 0, 0]
    
    if p < N-2:
        input_X = [X[p], X[p+1], X[p+2]]
        input_Y = [Y[p], Y[p+1], Y[p+2]]
        F = parabola(input_X, input_Y)
        
    if p > 0:
        input_X = [X[p-1], X[p], X[p+1]]
        input_Y = [Y[p-1], Y[p], Y[p+1]]
        B = parabola(input_X, input_Y)
        
    if B == [0, 0, 0]:
        B = F
        
    if F == [0, 0, 0]:
        F = B
    
    valB = sqf(B, x)
    valF = sqf(F, x)
    
    h = X[p+1] - X[p]
    
    wB = (x - X[p]) / h
    wF = (X[p+1] - x) / h
    
    return sqf(B, x) * wB + sqf(F, x) * wF
   

def rk4(x, y, h, f):
        
    k1 = h * f(x, y)
    k2 = h * f(x + 0.5*h, y + 0.5*k1)
    k3 = h * f(x + 0.5*h, y + 0.5*k2)
    k4 = h * f(x + h, y + k3)
    
    delta = (k1 + 2*k2 + 2*k3 + k4)/6
    
    return y + delta


def bisec(a_, b_, delta, f):
    
    a = a_
    b = b_  
        
    while b - a > delta:
            
        mid = (a + b) / 2.
        
        if f(a) * f(mid) > 0:
            a = mid
        else:
            b = mid
                
    return (a+b)/2


############## Basic functions ############## 

   
def bar():
    
    return 101325
    
    
def F_to_C(F):

    return (F-32)/1.8


def F_to_K(F):

    return F_to_C(F) + 273.15

def K_to_F(K):

    return 1.8*(K - 273.15) + 32

def K_to_C(K):
    
    return K - 273.15


def C_to_K(C):
    
    return C + 273.15
    
    
def C_to_F(C):
    
    return 1.8*C + 32


def antoine(A, B, C, T):
    
    power = A - B/(C+T)
    return pow(10, power)
    
def rev_antoine(A, B, C, P):
    
    return B/(A-log(P, 10.)) - C

def molar_to_mass(x):

    return 17. * x / (8. * x + 9)
    
    
def mass_to_molar(m):
    
    return 9. * m / (17 - 8. * m)
    
def molar_mass(x):
    
    return 34. * x + 18 * (1-x)    
    
def molar_mass_w(m):
    
    x = mass_to_molar(m)
    return molar_mass(x)

def density(x, T):

    c0 = [-0.168143, 0.641843, 0.9998]
    c100 = [-0.098476, 0.505076, 0.9584]

    ro0 = sqf(c0, x)
    ro100 = sqf(c100, x)

    T_C = K_to_C(T)

    return 1000 * (ro0 + T_C * (ro100 - ro0) / 100)

def mix(c1, c2, s1, s2):
    
    return (s1*c1 + s2*c2) / (s1+s2)


############## Physical functions ############## 


def w_vap_press(T):
        
        A = 10.19621
        B = 1730.63 
        C = -39.724
        
        return antoine (A, B, C, T)         
        
        
def hp_vap_press(T_):
    
    T = float(T_)
    e = 8.92536 - 2482.6/T - 24675/(T*T)
    
    P = pow(10, e)
    P = P * 133.322
    
    return P
        
        
def w_vap_temp(P):
    
    A = 10.19621
    B = 1730.63 
    C = 233.426 
    
    return C_to_K(rev_antoine(A, B, C, P))
        
        
def hp_vap_temp(P):
    
    P_ = P/133.322
    
    a = -24675.
    b = -2482.6
    c = 8.92536 - log(P_, 10)
    
    delta = b*b - 4*a*c
    sdelta = sqrt(delta)    
    
    x = 2*a / (-b - sdelta)

    return x
    
        
def ac(x, T):     
      
    T_C = K_to_C(T)
    
    X = [0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.]

    Yw90 = [1., 0.990, 0.960, 0.910, 0.842, 0.770, 0.696, 0.628, 0.563, 0.502, 0.455]
    Yw60 = [1., 0.989, 0.953, 0.895, 0.825, 0.746, 0.665, 0.588, 0.518, 0.456, 0.405]
    Yw30 = [1., 0.987, 0.946, 0.880, 0.803, 0.715, 0.628, 0.546, 0.471, 0.410, 0.355]
    Yw00 = [1., 0.985, 0.940, 0.865, 0.773, 0.675, 0.587, 0.497, 0.420, 0.353, 0.302]

    Yw = [Yw00, Yw30, Yw60, Yw90]

    Yhp90 = [0.365, 0.455, 0.552, 0.646, 0.733, 0.820, 0.890, 0.940, 0.975, 0.995, 1.]
    Yhp60 = [0.320, 0.405, 0.502, 0.606, 0.701, 0.797, 0.875, 0.931, 0.971, 0.993, 1.]
    Yhp30 = [0.272, 0.358, 0.458, 0.558, 0.668, 0.765, 0.850, 0.918, 0.966, 0.992, 1.]
    Yhp00 = [0.223, 0.306, 0.407, 0.515, 0.620, 0.731, 0.827, 0.904, 0.960, 0.990, 1.]

    Yhp = [Yhp00, Yhp30, Yhp60, Yhp90]

    T_tab = [0., 30., 60., 90.]        
    
    yw = [0., 0., 0., 0.]
    yhp = [0., 0., 0., 0.]
    
    for i in range(0, 4):
        yw[i] = parabolic([X, Yw[i]], x)
        yhp[i] = parabolic([X, Yhp[i]], x)
        
    points = [T_tab, yw]
    yw_ = parabolic(points, T_C)

    points = [T_tab, yhp]
    yhp_ = parabolic(points, T_C)
    
    return [yhp_, yw_]
        
        
def sol_vap_press(x_, T):
    
    ac_ = ac(x_, T)
    
    return hp_vap_press(T)*x_*ac_[0] + w_vap_press(T)*(1.-x_)*ac_[1] 
         
        
def sol_vap_temp(x_, P):
        
    a = 270.
    b = 430.  
        
    def testfun(T):
            
        return P - sol_vap_press(x_, T)
        
    while b - a > 0.04:
            
        mid = (a + b) / 2.
        
        if testfun(a) * testfun(mid) > 0:
            a = mid
        else:
            b = mid
                
    return (a+b)/2
        
        
def vap_press(y_, T):

    def testfun(x):
        
        return sol_vap_comp_T(x, T) - y_
    
    x = bisec(0., 1., 0.0005, testfun)
    
    return sol_vap_press(x, T)
        

def vap_temp(y_, P):
    
    def testfun(x):
        
        return sol_vap_comp_P(x, P) - y_
        
    x = bisec(0., 1., 0.0005, testfun)
    
    return sol_vap_temp(x, P)

        
def sol_vap_comp_T(x_, T):
        
    ac_ = ac(x_, T)
    
    return hp_vap_press(T)*x_*ac_[0] / sol_vap_press(x_, T)        
        
        
def sol_vap_comp_P(x_, P):
        
    T = sol_vap_temp(x_, P)
        
    return sol_vap_comp_T(x_, T)


def vap_heat_T(x, T):
    
    y = sol_vap_comp_T(x, T)
    M = molar_mass(y)
    
    return  M * 4.185 * (88*x*x - 300*x + 570 + (1.333 - 0.666*x)*(0.0004762*T*T - 0.6433*T + 156.4))


def vap_heat_P(x, P):
    
    T = sol_vap_temp(x, P)

    return vap_heat_T(x, T)
    

def cond_heat_T(y, T):
    
    x = vap_cond_T(y, T)
    M = molar_mass(y)
    
    return  M * 4.185 * (88*x*x - 300*x + 570 + (1.333 - 0.666*x)*(0.0004762*T*T - 0.6433*T + 156.4))
    

def cond_heat_P(y, P):
    
    T = vap_temp(y, P)    
    
    return cond_heat_T(y, T)


############## Numerical functions ############## 


def vap_cond_T(y_, T):
         
    a = 0.
    b = 1.  
        
    def testfun(x):
            
        return y_ - sol_vap_comp_T(x, T)
        
    while b - a > 0.0001:
            
        mid = (a + b) / 2
            
        if testfun(a) * testfun(mid) > 0:
            a = mid
        else:
            b = mid
                
    return (a+b)/2
        
        
def vap_cond_P(y_, P):
        
    a = 0.
    b = 1.  
        
    def testfun(x):
            
        return y_ - sol_vap_comp_P(x, P)
        
    while b - a > 0.0001:
            
        mid = (a + b) / 2.
            
        if testfun(a) * testfun(mid) > 0:
            a = mid
        else:
            b = mid
                
    return (a+b)/2


def condensate_to_T(y_, P, T):

    def der(eff, y):
        x = vap_cond_P(y, P)
        return -((1-y)*x - y*(1-x))

    part = 1.
    y = y_
    Hc = 0.
    yn = 0.

    while (vap_temp(y, P) > T) & (part > 0.001):
        h = 0.05
        yn = rk4(1, y, h, der)
        Hc = Hc + part* h * cond_heat_T(0.5 * (yn + y), T)
        y = yn
        part = part * (1-h)
    return [y, Hc, 1 - part]

    
def condensate(y_, P, eff_): #sumting wong
    
    def der(eff, y):
        
        return (y - vap_cond_P(y, P)) / (1. - eff)
    
    def rk4(eff, y, h):
        
        k1 = h * der(eff, y)
        k2 = h * der(eff + 0.5*h, y + 0.5*k1)
        k3 = h * der(eff + 0.5*h, y + 0.5*k2)
        k4 = h * der(eff + h, y + k3)
    
        delta = (k1 + 2*k2 + 2*k3 + k4)/6
    
        return y + delta
        
    progress = 0.   
    y = y_;    
    
    while progress < eff_:
        
        h = (1.01 - progress)/60.
        y = rk4(progress, y, h)
        progress  = progress + h
        print sol_vap_temp(vap_cond_P(y, P), P)
        if y < 0:
            return 0

    return y        
        
        ###########


def evaporate(liq_comp_mol, P, eff, output):
    
    def der(eff, liq_comp_mass):
        
        der_liq_comp_mol = mass_to_molar(liq_comp_mass)
        
        sol_vap_mol = sol_vap_comp_P(der_liq_comp_mol, P, output)
        sol_vap_mass = molar_to_mass(sol_vap_mol)
        
        return (liq_comp_mass - sol_vap_mass) / (1. - eff)
    
    def rk4(eff, liq_comp_mass, h):
        
        k1 = h * der(eff, liq_comp_mass)
        k2 = h * der(eff + 0.5*h, liq_comp_mass + 0.5*k1)
        k3 = h * der(eff + 0.5*h, liq_comp_mass + 0.5*k2)
        k4 = h * der(eff + h, liq_comp_mass + k3)
    
        delta = (k1 + 2*k2 + 2*k3 + k4)/6
    
        return liq_comp_mass + delta
        
    liq_comp_mass = molar_to_mass(liq_comp_mol)
    progress = 0.   
    
    while progress < eff:
        
        liq_comp_mass = rk4(progress, liq_comp_mass, eff/100)
        progress  = progress + eff/100
    
        
        
    return mass_to_molar(liq_comp_mass)
        
            
## SIMULATION BLOCKS ##

class boiler: ## to be worked on

    def __init__(self, type):

        self.P = 0.
        self.x = 0.

        self.T = 0.
        self.y_out = 0.


        self.V_out = 1

        if type != "s":
            self.power = 100.
            self.V_out = 1.
            self.volume = 100.
            self.L_in = 0
            self.x_in = 0

    def initiate(self):
        self.T = sol_vap_temp(self.x, self.P)
        self.y_out = sol_vap_comp_T(self.x, self.T)

        if type != "s":
            self.V_out = power / vap_heat_T(self.x, self.T)

class column:

    def __init__(self, N_, fi_, P_,):

        self.N = N_
        self.fi = fi_
        self.trays = [0]
        self.P = P_
        self.P_keep = 1.


    class tray:

        def __init__(self, P_, fi_, x_):

            self.P = P_
            self.fi = fi_
            self.x = x_

            self.V_in = 0.
            self.y_in = 0.

            self.L_in = 0.
            self.x_in = 0.

            self.V_lo = 0.
            self.y_lo = 0.

            self.V_em = 0.
            self.y_em = 0.

            self.V_out = 0.
            self.y_out = 0.

            self.L_out = 0.
            self.x_out = 0.


        def calculate(self):


            T_TL = sol_vap_temp(self.x, self.P)
            vap_heat = vap_heat_T(self.x, T_TL)

            V_in_temp = vap_temp(self.y_in, self.P)

            T_TL = V_in_temp - self.fi * (V_in_temp - T_TL)

            LO = condensate_to_T(self.y_in, self.P, T_TL)

            self.V_lo = 1 - LO[2]
            self.y_lo = LO[0]

            heat_rec = self.V_in * LO[1]

            self.V_em = heat_rec / vap_heat
            self.y_em = sol_vap_comp_P(self.x, self.P)

            self.V_out = self.V_em + self.V_lo
            self.y_out = mix(self.y_lo, self.y_em, self.V_lo, self.V_em)

            self.L_out = self.L_in + self.V_in - self.V_out
            self.x = (self.V_in * self.y_in + self.L_in * self.x_in - self.V_out * self.y_out) / self.L_out
            self.x_out = self.x

        def connect(self, a, b):
            a.x_in = self.x_out

    def initiate(self):

        for i in range (1, self.N + 1):
            P_i = self.P - i * self.P * (1 - self.P_keep) / self.N
            self.trays.append(self.tray(P_i, self.fi, 0.))


class reflux:

    def __init__(self, type_, param):

        self.type = type_
        self.R = param

        self.P = 0.

        self.V_in = 0.
        self.y_in = 0.

        self.L_out = 0.
        self.x_out = 0.

    def calculate(self):

        if self.type == "f":
            self.L_out = self.R * V_in
            self.x_out = y_in

        if self.type == "p":
            T_W = vap_temp(0, self.P)
            T_V = vap_temp(self.y_in, self.P)

            T = T_W + self.R * (T_V - T_W)
            LO = condensate_to_T(self.y_in, self.P, T)

            self.L_out = self.V_in * LO[2]
            self.x_out = (self.y_in - (1 - self.L_out) * LO[0]) / self.L_out

        else:
            print "Warning: in peroxylib.reflux(): wrong type of reflux specified. Specify \"f\" for full, or \"p\" for partial"

