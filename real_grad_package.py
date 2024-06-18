import numpy as np
import random
import math

def parse_poly(filename):

    with open(filename) as f:
        lines = f.readlines()

    line_counter = 0
    for line in lines:
        if line_counter==0:
            N = int(line.split()[2])
            M = int(line.split()[3])
            F_mat = np.zeros((N,M))
            B_mat = np.zeros((M,N))
        else:
            line_temp = line.split()
            if line_temp[0]!='0':
                term_size = len(line_temp)-2
                for i in range(term_size):
                    F_mat[int(line_temp[i])-1,line_counter-1] = 1
                    B_mat[line_counter-1,int(line_temp[i])-1] = int(line_temp[-1])
        line_counter = line_counter + 1
        
    return F_mat, B_mat

def get_IdealRealGradient(B_mat,v_arr):
    N = B_mat.shape[1]
    M = B_mat.shape[0]
    idealgrad_arr = np.zeros((1,N))
    for i in range(N):
        for j in range(M):
            if B_mat[j,i]!=0:
                member_indices = np.where(B_mat[j,]!=0)[0]
                if len(member_indices)>1:
                    member_indices = member_indices[np.where(member_indices!=i)]
                    idealgrad_arr[0,i] = idealgrad_arr[0,i] + B_mat[j,i]*np.prod(v_arr[0,member_indices])
                    
    return idealgrad_arr
    
    
class real_grad:
    
    def __init__(self, N, M, F_mat, B_mat, do_xbar_variation, **kwargs):
        
        self.N = N
        self.M = M
        self.max_var_degree = np.max(F_mat)
        self.do_xbar_variation = do_xbar_variation
        self.xmin = kwargs["xmin"] if "xmin" in kwargs.keys() else 0.028
        self.xmax = kwargs["xmax"] if "xmax" in kwargs.keys() else 1.0
        self.vcg_max = kwargs["vcg_max"] if "vcg_max" in kwargs.keys() else 2.0
        self.vt0 = kwargs["vt0"] if "vt0" in kwargs.keys() else 2.5
        self.vth_min = kwargs["vth_min"] if "vth_min" in kwargs.keys() else 2.0
        self.vth_max = kwargs["vth_max"] if "vth_max" in kwargs.keys() else 3.5
        self.x_dr = self.xmax/self.xmin
        self.vb1 = kwargs["vb1"] if "vb1" in kwargs.keys() else 0.1
        self.vb3 = kwargs["vb3"] if "vb3" in kwargs.keys() else self.vb1
        self.R = kwargs["R"] if "R" in kwargs.keys() else 1e4
        self.Is = kwargs["Is"] if "Is" in kwargs.keys() else 0.03e-15
        self.I0 = kwargs["I0"] if "I0" in kwargs.keys() else 100e-9#30e-12
        self.m1 = kwargs["m1"] if "m1" in kwargs.keys() else 0.02788
        self.m2 = kwargs["m2"] if "m2" in kwargs.keys() else 0.13
        self.delta_vx = math.log(self.x_dr)*self.m1
        self.Av1 = 0.2/self.delta_vx
        self.Gon = kwargs["Gon"] if "Gon" in kwargs.keys() else 100e-6
        self.Goff = kwargs["Goff"] if "Goff" in kwargs.keys() else 1e-6
        self.sigma_Gon = kwargs["sigma_Gon"] if "sigma_Gon" in kwargs.keys() else 3e-6
        self.sigma_Goff = kwargs["sigma_Goff"] if "sigma_Goff" in kwargs.keys() else 0.25e-6
        self.Iflash_precision = kwargs["Iflash_precision"] if "Iflash_precision" in kwargs.keys() else 2
        self.vth_bound = kwargs["vth_bound"] if "vth_bound" in kwargs.keys() else (self.Iflash_precision*self.m2)/50
        self.Rf = kwargs["Rf"] if "Rf" in kwargs.keys() else self.m2/(self.m1*self.Av1*((self.Gon-self.Goff)/self.max_var_degree))
        self.Av2 = kwargs["Av2"] if "Av2" in kwargs.keys() else self.Rf*self.Av1*((self.Gon-self.Goff)/self.max_var_degree)
        if self.xmax <=1:
            self.vb2 = (self.m2*math.log(self.xmax)+(self.m2/self.m1)*self.vb1+self.m2*math.log(self.R*self.Is)-self.vcg_max)/(1+(self.m2/self.m1))
        
        #self.vb2 = kwargs["vb2"] if "vb2" in kwargs.keys() else self.m1*math.log(self.R*self.Is)
        self.v_arr = kwargs["v_arr"] if "v_arr" in kwargs.keys() else (self.vb1+self.xmin)*np.ones((self.N,1))
        #self.vbot = self.vb1+self.m1*math.log(self.R*self.Is)
        self.vbot = (1+self.Av1)*self.vb3 - self.Av1*self.vb1 - self.Av1*self.m1*math.log(self.R*self.Is)
        self.gf_mat = self.get_Gf_mat(F_mat)
        self.vth_mat = self.get_vth_mat(B_mat)
        
    def get_Gf_mat(self,F_mat):
        gf_mat = np.zeros((self.N,2*self.M))
        for i in range(F_mat.shape[0]):
            for j in range(F_mat.shape[1]):
                if self.do_xbar_variation==1:
                    new_Goff = np.random.normal(self.Goff,self.sigma_Goff)
                    new_Gon = np.random.normal(self.Gon,self.sigma_Gon)
                    gf_mat[i,2*j] = new_Goff + F_mat[i,j]*((new_Gon-new_Goff)/self.max_var_degree)
                    gf_mat[i,2*j+1] = new_Goff
                else:    
                    gf_mat[i,2*j] = self.Goff + F_mat[i,j]*((Gon-Goff)/self.max_var_degree)
                    gf_mat[i,2*j+1] = self.Goff
        
        return gf_mat
    
    def get_vth_mat(self,B_mat):
        vth_mat = np.zeros((self.M,self.N))
        #vth_const = self.Rf*((self.Gon-self.Goff)/self.max_var_degree)*self.vb1 - self.vb2
        vth_const = self.vt0
        for i in range(B_mat.shape[0]):
            for j in range(B_mat.shape[1]):
                if B_mat[i,j]==0:
                    vth_mat[i,j] = self.vth_max
                else:
                    desired_vth = vth_const - self.m2*math.log(B_mat[i,j])
                    if self.do_xbar_variation==1:
                        vth_mat[i,j] = np.random.uniform(desired_vth-self.vth_bound,desired_vth+self.vth_bound)
                    else:
                        vth_mat[i,j] = desired_vth
                
        return vth_mat

    def get_Vx(self,v_arr):
        return (1+self.Av1)*self.vb3 - self.Av1*self.vb1 + self.Av1*self.m1*np.log((v_arr-self.vb1)/(self.R*self.Is))
        #return self.vb1-self.m1*np.log((v_arr-self.vb1)/(self.R*self.Is))
    
    def get_Vc(self,vx_arr):
        vmem_arr = vx_arr - self.vbot
        fa_currents_diff = np.dot(vmem_arr,self.gf_mat)
        tia_voltages_diff = self.vbot - self.Rf*fa_currents_diff
        net_voltages = tia_voltages_diff[0,1::2]-tia_voltages_diff[0,0::2]
        return net_voltages.reshape(1,vx_arr.shape[1])

    def get_Vs(self,v_arr):
        return (1+self.Av2)*self.vb2 - self.Av2*(self.vb1-self.m1*np.log((v_arr-self.vb1)/(self.R*self.Is)))
    
    def get_Igrad(self,vc_arr,vs_arr):
        Igrad_arr = np.zeros((1,self.N))
        for i in range(self.N):
            Igrad_arr[0,i] = np.sum(self.I0*np.exp((vc_arr-vs_arr[0,i]-self.vth_mat[:,i].T)/self.m2))
        
        normalizing_constant = self.I0*math.e**((-(1+self.Av2)*self.vb2+self.Av2*self.vb1+self.Av2*self.m1*math.log(self.R*self.Is)-self.vt0)/self.m2)
            
        return Igrad_arr/(normalizing_constant)
    
    def get_realGradient(self,v_arr):
        vx_arr = self.get_Vx(v_arr)
        vs_arr = self.get_Vs(v_arr)
        vc_arr = self.get_Vc(vx_arr)
        Igrad_arr = self.get_Igrad(vc_arr,vs_arr)
        return Igrad_arr