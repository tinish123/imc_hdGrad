import numpy as np

######Energy numbers taken from Table S1 & S3######
E_wl = 0.27
E_xbar1 = 0.008
E_xbar0 = 0.00008
E_tia = 0.26
E_sa = 0.0023
E_mux = 0.00023
E_wta = 0.037
E_dff = 0.023
E_hc = 0.25

######Area numbers taken from Table S1 & S3######
A_wl = 2.73
A_mem = 0.0123
A_tia = 10.5135
A_sa = 8.058
A_dff = 1.915
A_mux = 2.084
A_wta = 29.01
A_hc = 3.48

#EPA stands for Energy Power Area

def getEPA_wsat(N,M,t_clk):
    E_xbar = M*3*E_xbar1+M*(2*N-3)*E_xbar0
    E_fp = 2*N*E_wl+E_xbar+M*E_tia+2*M*E_sa
    E_rs = M*E_dff
    E_bmp = E_wl+E_xbar+N*E_mux+N*E_tia+N*E_sa
    E_bbp = M*E_wl+E_xbar+N*E_mux+3*E_tia+3*E_sa
    E_bwp = E_xbar+3*E_tia+N*E_wta
    E_fl = 2*N*E_dff
    E_avg_wsat = (E_fp+E_rs+E_bmp+E_bbp+E_bwp+E_fl)/3
    P_wsat = (E_avg_wsat*1e-12)/t_clk
    A_xbar = 2*N*M*A_mem
    A_wsat = 2*N*A_wl+A_xbar+M*A_tia+2*M*A_sa+M*A_dff+M*A_wl+A_xbar+N*A_mux+N*A_tia+N*A_sa+N*A_wta+2*N*A_dff

    return E_avg_wsat*1e-12, P_wsat, A_wsat
    
def getEPA_hohnn(N,M,t_clk):
    E_xbar = M*3*E_xbar1+M*(2*N-3)*E_xbar0
    E_fp = 2*N*E_dff+2*N*E_wl+E_xbar+M*E_tia+2*M*E_sa
    E_bp = 2*M*E_dff+M*E_wl+E_xbar+2*N*E_tia+N*E_hc
    E_avg_hohnn = (E_fp+E_bp)
    P_hohnn = (E_avg_hohnn*1e-12)/t_clk
    A_xbar = 2*N*M*A_mem
    A_hohnn = 2*N*A_wl+A_xbar+M*A_tia+2*M*A_sa+2*M*A_dff+M*A_wl+A_xbar+2*N*A_tia+N*A_hc+2*N*A_dff

    return E_avg_hohnn*1e-12, P_hohnn, A_hohnn
    
def getEPA_sohnn(N,M,t_clk):
    E_xbar = (N*20+M*3)*E_xbar1+(2*(N+M)*(N+M)-(N*20+M*3))*E_xbar0
    E_avg_sohnn = (N+M)*E_dff+(N+M)*E_wl+E_xbar+2*(N+M)*E_tia+(N+M)*E_hc
    P_sohnn = (E_avg_sohnn*1e-12)/t_clk
    A_xbar = 2*(N+M)*(N+M)*A_mem
    A_sohnn = (N+M)*A_dff+(N+M)*A_wl+A_xbar+2*(N+M)*A_tia+(N+M)*A_hc
    
    return E_avg_sohnn*1e-12, P_sohnn, A_sohnn