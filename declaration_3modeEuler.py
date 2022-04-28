from functions import *
from parametersPlot import *

#-----------------------
# Arrays
#-----------------------
N      = 64
kG     = int(N/3)
L      = 2.0 * np.pi
dx     = L / N
eps    = 1e-2
T      = 1.0
dt     = 0.001
t_N    = int(T/dt) + 1
tim_N    = int(T/dt) + 1
x_ar   = np.linspace(0,N-1,N) * dx
y_ar   = np.linspace(0,N-1,N) * dx
u_ar   = np.array([[0.0]*N]*N)
v_ar   = np.array([[0.0]*N]*N)
w_ar   = np.array([[0.0]*N]*N)
t_ar   = np.linspace(0,T,t_N)
en_t   = np.array([0.0]*t_N)
es_t   = np.array([0.0]*t_N)
en_k_t = np.array([0.0]*t_N)
en_q_t = np.array([0.0]*t_N)
en_p_t = np.array([0.0]*t_N)
es_k_t = np.array([0.0]*t_N)
es_q_t = np.array([0.0]*t_N)
es_p_t = np.array([0.0]*t_N)
sav_dir= 'results/'
#-----------------------
# Scalar functions - Initial condition
#-----------------------
Ak = (0.0 - 1.0j) * eps / kG
Aq = 1.0 + 2.0j
Ap = 0.0 + 0.0j
Bk = (1.0 + 0.0j) * eps
Bq = 1.0 + 3.0j
Bp = 0.0 + 0.0j
y  = np.array([Ak,Aq,Ap,Bk,Bq,Bp])
#-----------------------
# Modes
#-----------------------
qx = 1.0
qy = 0.0
q2 = qx**2.0 + qy**2.0
kx = -1.0
ky = -float(kG)
k2 = kx**2.0 + ky**2.0
px = 0.0
py = float(kG)
p2 = px**2.0 + py**2.0
de = qx * py
#-----------------------
# Differential Equation
#-----------------------
def evolution_eqn(Ak,Aq,Ap,Bk,Bq,Bp):
    Bk_dot = de * np.conjugate( Ap * Bq - Aq * Bp )
    Bq_dot = de * np.conjugate( Ak * Bp - Ap * Bk )
    Bp_dot = de * np.conjugate( Aq * Bk - Ak * Bq )
    Ak_dot = de * ( p2 - q2 ) * np.conjugate( Aq * Ap ) / k2
    Aq_dot = de * ( k2 - p2 ) * np.conjugate( Ap * Ak ) / q2
    Ap_dot = de * ( q2 - k2 ) * np.conjugate( Ak * Aq ) / p2
    ar_dot = np.array([Ak_dot,Aq_dot,Ap_dot,Bk_dot,Bq_dot,Bp_dot])
    return ar_dot

def velocity_modes(Ak,Aq,Ap):
    uk   = +(0.0+1.0j) * ky * Ak
    uq   = +(0.0+1.0j) * qy * Aq
    up   = +(0.0+1.0j) * py * Ap
    vk   = -(0.0+1.0j) * kx * Ak
    vq   = -(0.0+1.0j) * qx * Aq
    vp   = -(0.0+1.0j) * px * Ap
    ar_v = np.array([uk,uq,up,vk,vq,vp])
    return ar_v

def phase_modes(x,y):
    ph_k = np.exp( (0.0+1.0j) * ( kx * x + ky * y ) )
    ph_q = np.exp( (0.0+1.0j) * ( qx * x + qy * y ) )
    ph_p = np.exp( (0.0+1.0j) * ( px * x + py * y ) )
    phase= np.array([ph_k,ph_q,ph_p])
    return phase

def velocity_field(uk,uq,up,vk,vq,vp,Bk,Bq,Bp):
    global u_ar,v_ar,w_ar
    for i_x in range(0,N):
        for i_y in range(0,N):
            ph=phases(x_ar[i_x],y_ar[i_y])
            u_c           = uk * ph[0] + uq * ph[1] + up * ph[2]
            v_c           = vk * ph[0] + vq * ph[1] + vp * ph[2]
            w_c           = Bk * ph[0] + Bq * ph[1] + Bp * ph[2]
            u_ar[i_x,i_y] = 2.0 * u_c.real
            v_ar[i_x,i_y] = 2.0 * v_c.real
            w_ar[i_x,i_y] = 2.0 * w_c.real

def energy_k(A,B):
    return 0.5 * ( k2*(abs(A)**2) + abs(B)**2 )

def energy_q(A,B):
    return 0.5 * ( q2*(abs(A)**2) + abs(B)**2 )

def energy_p(A,B):
    return 0.5 * ( p2*(abs(A)**2) + abs(B)**2 )

def energy_field(Ak,Aq,Ap,Bk,Bq,Bp):
    return energy_k(Ak,Bk) + energy_q(Aq,Bq) + energy_p(Ap,Bp)

def enstrophy_field(Ak,Aq,Ap,Bk,Bq,Bp):
    return k2*energy_k(Ak,Bk) + q2*energy_q(Aq,Bq) + p2*energy_p(Ap,Bp)
