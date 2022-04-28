from declaration_3modeEuler import *
# ----------------------------------------------------
# RUNGA KUTTA 4th ORDER ALGORITHM
# ----------------------------------------------------
def time_evolution(Ak,Aq,Ap,Bk,Bq,Bp):
    y  = np.array([Ak,Aq,Ap,Bk,Bq,Bp])
    y_t= y
    dy_1 = dt * evolution_eqn( y_t[0], y_t[1], y_t[2], y_t[3], y_t[4], y_t[5] )
    y_t  = y + 0.5 * dy_1
    dy_2 = dt * evolution_eqn( y_t[0], y_t[1], y_t[2], y_t[3], y_t[4], y_t[5] )
    y_t  = y + 0.5 * dy_2
    dy_3 = dt * evolution_eqn( y_t[0], y_t[1], y_t[2], y_t[3], y_t[4], y_t[5] )
    y_t  = y + dy_3
    dy_4 = dt * evolution_eqn( y_t[0], y_t[1], y_t[2], y_t[3], y_t[4], y_t[5] )
    y    = y + ( dy_1 + 2.0 * ( dy_2 + dy_3 ) + dy_4 ) / 6.0
    return y

def solve_system():
    global Ak,Aq,Ap,Bk,Bq,Bp
    for t in range(0,t_N):
        en_k_t[t] = energy_k(Ak,Bk)
        en_q_t[t] = energy_q(Aq,Bq)
        en_p_t[t] = energy_p(Ap,Bp)
        es_k_t[t] = k2 * energy_k(Ak,Bk)
        es_q_t[t] = q2 * energy_q(Aq,Bq)
        es_p_t[t] = p2 * energy_p(Ap,Bp)
        en_t[t] = energy_field(Ak,Aq,Ap,Bk,Bq,Bp)
        es_t[t] = enstrophy_field(Ak,Aq,Ap,Bk,Bq,Bp)
        y  = time_evolution(Ak,Aq,Ap,Bk,Bq,Bp)
        Ak = y[0]
        Aq = y[1]
        Ap = y[2]
        Bk = y[3]
        Bq = y[4]
        Bp = y[5]

def write_data():
    DataOut = np.column_stack((t_ar[:],en_k_t[:],en_q_t[:],en_p_t[:]))
    path    = sav_dir+'en_modes_S'+time_stamp()+'.dat'
    np.savetxt(path,DataOut,fmt=('%8.4f','%34.17f','%34.17f','%34.17f') )
    DataOut = np.column_stack((t_ar[:],es_k_t[:],es_q_t[:],es_p_t[:]))
    path    = sav_dir+'es_modes_S'+time_stamp()+'.dat'
    np.savetxt(path,DataOut, fmt=('%8.4f','%34.17f','%34.17f','%34.17f') )
    DataOut = np.column_stack((t_ar[:],en_t[:],es_t[:]))
    path=sav_dir+'evolution_S'+time_stamp()+'.dat'
    np.savetxt(path,DataOut, fmt=('%8.4f','%34.17f','%34.17f') )

solve_system()
write_data()
