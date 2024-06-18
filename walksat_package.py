import numpy as np
import random
import sys
import math
from multiprocessing import Process
import concurrent.futures
import itertools
import pickle

def get_errorMap(load_dir,sim_type,**kwargs):
    
    error_map = np.zeros((17,2))
    
    if sim_type==1:
        if 'N' in kwargs:
            N = kwargs['N']
        else:
            N = 20
        with open(load_dir+'bv_vs_nvar_errormap_N'+str(N)+'.csv','rb') as f:
            lst = pickle.load(f)
        error_map = lst[0]
    elif sim_type==2:
        if 'N' in kwargs:
            N = kwargs['N']
        else:
            N = 20
        if 'Gh_sigma' in kwargs:
            Gh_sigma = kwargs['Gh_sigma']
        else:
            Gh_sigma = 3
        if 'Gl_sigma' in kwargs:
            Gl_sigma = kwargs['Gl_sigma']
        else:
            Gl_sigma = 0.25
        with open(load_dir+'bv_vs_sigma_errormap_N'+str(N)+'_Gh_sigma_'+str(Gh_sigma)+'_Gl_sigma_'+str(Gl_sigma)+'.csv','rb') as f:
            lst = pickle.load(f)
        error_map = lst[0]
    else:
        error_map[:,0] = np.arange(0,17,1)
        error_map[:,1] = 0

    return error_map
    
def sat_eval(var_val,clause_mat):
    
    nvar = len(var_val)
    nclause = clause_mat.shape[0]
    clause_stat_mat = np.zeros((nclause,1),dtype=bool)
    clause_f_mat = np.zeros((nclause,1))
    sat_result = 1
    
    for i in range(nclause):
        
        var_ind = np.where(clause_mat[i,]!=0)[0]
        clause_stat_mat[i], clause_f_mat[i] = clause_eval(var_val[var_ind],clause_mat[i,var_ind])
        sat_result = sat_result and clause_stat_mat[i]
        
    return clause_stat_mat, sat_result, clause_f_mat
        
def clause_eval(var_val,var_sign):
    
    nvar = len(var_val)
    clause_sum = 0
    clause_f_count = 0
    for i in range(nvar):
        if var_sign[i]==1:
            clause_sum = clause_sum or var_val[i]
            clause_f_count = clause_f_count + (1 if var_val[i] else 0)
        elif var_sign[i]==-1:
            clause_sum = clause_sum or ~var_val[i]
            clause_f_count = clause_f_count + (1 if ~var_val[i] else 0)
        else:
            clause_sum = clause_sum
            
    return clause_sum, clause_f_count

def rand_init(nvar):
    
    var_val_mat = np.random.rand(nvar,1)
    var_val_mat = [1 if (i>0.5) else 0 for i in var_val_mat]
    var_val_mat = np.array(var_val_mat,dtype=bool)
    
    return var_val_mat

def number_of_sat_clause(clause_stat_mat):
    
    return len(np.where(clause_stat_mat==1)[0])


def get_stats(walksat_result):
        
    prob_s = np.mean(walksat_result[:,1])
    success_ind = np.where(walksat_result[:,1]==1)[0]
    avg_flips = np.mean(walksat_result[success_ind,0])
    std_flips = np.std(walksat_result[success_ind,0])
    
    max_restarts = np.shape(walksat_result)[0]
    p_targ = 0.99
    sorted_arr = walksat_result[walksat_result[:,0].argsort()]
    if prob_s>=p_targ:
        ind_tts = math.ceil(p_targ*max_restarts)-1
        tts_99 = sorted_arr[ind_tts,0]
    else:
        if prob_s == 0:
            tts_99 = 0
        else:
            tts_99 = sorted_arr[max_restarts-1,0]*(math.log((1-p_targ),10)/math.log((1-prob_s),10))
    
    return avg_flips, prob_s, std_flips, tts_99

    
def SKC_heuristics(p,var_val,clause_mat,cnf_mat,clause_f_mat,clause_ind,clause_stat_mat,err_map):
    
    nvar = len(var_val)
    nclause = clause_mat.shape[0]
    var_list = np.where(clause_mat[clause_ind,]!=0)[0]
    num_sat_clause_init = number_of_sat_clause(clause_stat_mat)
    clause_stat_change_mat = np.zeros((len(var_list),1))
    var_break_value_mat = np.zeros((len(var_list),1))
    var_counter = 0
    xx = random.random()
    
    s_vec = np.zeros((nclause,1))
    z_vec = np.zeros((nclause,1))
    
    s_vec[np.where(clause_f_mat==0)[0],0] = 1
    z_vec[np.where(clause_f_mat==1)[0],0] = 1
     
    a_vec = np.transpose(cnf_mat)@s_vec
    b_vec = np.transpose(cnf_mat)@z_vec
    
    var1_indices = np.where(var_val==1)[0]
    var0_indices = np.where(var_val==0)[0]
    lit1_indices = np.sort(np.concatenate((2*var1_indices,2*var0_indices+1)))
    lit0_indices = np.sort(np.concatenate((2*var0_indices,2*var1_indices+1)))
    sw_bv_vec = b_vec[lit1_indices]
    mv_vec = a_vec[lit0_indices]
    gain_vec = mv_vec - sw_bv_vec
    sw_bv_vec = sw_bv_vec.astype(int)
    winning_bv = np.min(sw_bv_vec[var_list])
    hw_bv_vec = np.random.normal(np.array(list(map(lambda x: err_map[x,0], list(sw_bv_vec))),dtype=float),np.array(list(map(lambda x: err_map[x,1], list(sw_bv_vec)))))
    bv_0_1_thresh = (err_map[0,0]+err_map[1,0])/2
    hw_bv_vec[np.where(hw_bv_vec<=bv_0_1_thresh)] = 0
    
    var_break_value_mat = hw_bv_vec[var_list]
    clause_stat_change_mat = gain_vec[var_list]
    
    break_value_zero_ind = np.where(var_break_value_mat==0)[0]
    
    if len(break_value_zero_ind)!=0:
        chosen_var = var_list[np.random.choice(break_value_zero_ind,size=1)[0]]
    elif xx <= p:
        chosen_var = np.random.choice(var_list[np.where(var_break_value_mat==np.min(var_break_value_mat))[0]],size=1)[0]
    else:
        chosen_var = np.random.choice(var_list,size=1)[0]
        
    
    var_val[chosen_var] = ~var_val[chosen_var]
    clause_stat_mat_new, sat_result_new, clause_f_mat_new = get_clause_stat(chosen_var,clause_ind,clause_mat,clause_stat_mat,clause_f_mat)
    
    return var_val, clause_stat_mat_new, sat_result_new, clause_f_mat_new, sw_bv_vec, hw_bv_vec, winning_bv


def get_clause_stat(i,clause_ind,clause_mat,clause_stat_mat,clause_f_mat):
    var_sign = clause_mat[clause_ind,i]
    A_sign = var_sign
    B_sign = -var_sign
    A_clause_ind = np.where(clause_mat[:,i]==A_sign)[0]
    A_clause_ind_pa = np.intersect1d(A_clause_ind,np.where(clause_stat_mat==0)[0])
    A_clause_ind_pb = np.intersect1d(A_clause_ind,np.where(clause_stat_mat==1)[0])
    B_clause_ind = np.where(clause_mat[:,i]==B_sign)[0]
    B_clause_ind_pa = np.intersect1d(B_clause_ind,np.where(clause_f_mat==1)[0])
    B_clause_ind_pb = np.intersect1d(B_clause_ind,np.where(clause_f_mat!=1)[0])
    clause_stat_mat_new = np.array(clause_stat_mat)
    clause_f_mat_new = np.array(clause_f_mat)
    clause_stat_mat_new[A_clause_ind_pa] = True
    clause_stat_mat_new[B_clause_ind_pa] = False
    clause_f_mat_new[A_clause_ind] = clause_f_mat_new[A_clause_ind] + 1
    clause_f_mat_new[B_clause_ind] = clause_f_mat_new[B_clause_ind] - 1
    sat_result_new = 1 if (np.sum(clause_stat_mat_new)==len(clause_stat_mat_new)) else 0
    
    return clause_stat_mat_new, sat_result_new, clause_f_mat_new

def clausemat2cnfmat(clause_mat):
    
    nvar = clause_mat.shape[1]
    nclause = clause_mat.shape[0]
    nlit = 2*nvar
    cnfmat = np.zeros((nclause,nlit))
    
    for i in range(nvar):
        pos_lit_ind = np.where(clause_mat[:,i]==1)[0]
        if len(pos_lit_ind)!=0:
            cnfmat[pos_lit_ind,2*i] = 1
            
        neg_lit_ind = np.where(clause_mat[:,i]==-1)[0]
        if len(neg_lit_ind)!=0:
            cnfmat[neg_lit_ind,2*i+1] = 1
            
    return cnfmat


def WalkSat_Solver(p,clause_mat,max_iter,max_flips,err_map):
    
    nvar = clause_mat.shape[1]
    nclause = clause_mat.shape[0]
    cnf_mat = clausemat2cnfmat(clause_mat)
    
    walksat_result = np.zeros((max_iter,2))

    for r in range(max_iter):

        var_val_mat = rand_init(nvar)
        clause_stat_mat, sat_result, clause_f_mat = sat_eval(var_val_mat,clause_mat)
        
        for f in range(max_flips):

            if sat_result == 1:
                walksat_result[r,0] = f
                walksat_result[r,1] = 1
                break

            unsat_clause_list = np.where(clause_stat_mat==0)[0]
                
            chosen_clause_ind = np.random.choice(unsat_clause_list,size=1)[0]
                      
            var_val_mat, clause_stat_mat, sat_result, clause_f_mat, sw_bv_vec, hw_bv_vec, winning_bv = SKC_heuristics(p,var_val_mat,clause_mat,cnf_mat,clause_f_mat,chosen_clause_ind,clause_stat_mat,err_map)
            
        if walksat_result[r,1] == 0:
            clause_stat_mat, sat_result, clause_f_mat = sat_eval(var_val_mat,clause_mat)
            walksat_result[r,0] = max_flips
            if sat_result == 1:
                walksat_result[r,1] = 1  

    avg_flips, prob_s, std_flips, tts_99 = get_stats(walksat_result)
    
    return avg_flips, prob_s, std_flips, tts_99, walksat_result


def WalkSat_Solver_full(nvar,i,params):
    
    max_iter = params['max_iter']
    max_flips = params['max_flip']
    p = params['p']
    err_map = params['err_map']
    sat_dir = params['sat_dir']
    
    clause_mat = get_kSATprob(sat_dir,nvar,i)
    
    avg_flips, prob_s, std_flips, tts_99, walksat_result = WalkSat_Solver(p,clause_mat,max_iter,max_flips,err_map)
    
    return avg_flips, prob_s, std_flips, i, tts_99, walksat_result


def get_kSATprob(sat_dir,nvar,instance):
    
    if nvar==20:
        instance_name = 'uf20-0'+str(instance)
        line_start = 8
        k = 3
        nclause = 91
        if instance > 1000:
            print("Instance Number not present in Instance Directory!!")
            sys.exit()
    elif nvar==14:
        instance_name = 'uf14-0'+str(instance)
        line_start = 1
        k = 3
        nclause = 64
        if instance > 1000:
            print("Instance Number not present in Instance Directory!!")
            sys.exit()
    elif nvar==50:
        instance_name = 'uf50-0'+str(instance)
        line_start = 8
        k = 3
        nclause = 218
        if instance > 1000:
            print("Instance Number not present in Instance Directory!!")
            sys.exit()
    elif nvar==75:
        instance_name = 'uf75-0'+str(instance)
        line_start = 8
        k = 3
        nclause = 325
        if instance > 100:
            print("Instance Number not present in Instance Directory!!")
            sys.exit()
    elif nvar==100:
        instance_name = 'uf100-0'+str(instance)
        line_start = 8
        k = 3
        nclause = 430
        if instance > 1000:
            print("Instance Number not present in Instance Directory!!")
            sys.exit()
    elif nvar==150:
        instance_name = 'uf150-0'+str(instance)
        line_start = 8
        k = 3
        nclause = 645
        if instance > 100:
            print("Instance Number not present in Instance Directory!!")
            sys.exit()
    elif nvar==200:
        instance_name = 'uf200-0'+str(instance)
        line_start = 8
        k = 3
        nclause = 860
        if instance > 100:
            print("Instance Number not present in Instance Directory!!")
            sys.exit()
    elif nvar==225:
        instance_name = 'uf225-0'+str(instance)
        line_start = 8
        k = 3
        nclause = 960
        if instance > 100:
            print("Instance Number not present in Instance Directory!!")
            sys.exit()
    elif nvar==250:
        instance_name = 'uf250-0'+str(instance)
        line_start = 8
        k = 3
        nclause = 1065
        if instance > 100:
            print("Instance Number not present in Instance Directory!!")
            sys.exit()
    else:
        instance_name = 'uf20-0'+str(instance)
        line_start = 8
        k = 3
        nclause = 91
        if instance > 1000:
            print("Instance Number not present in Instance Directory!!")
            sys.exit()
    
    
    #file_addr = sat_dir+dir_name+'/'+instance_name+'.cnf'
    file_addr = sat_dir+'/'+instance_name+'.cnf'
    
    with open(file_addr) as f:
        lines = f.readlines()
    f.close()
    
    line_end = line_start + nclause
    clause_mat = np.zeros((nclause,nvar))

    for i in np.arange(line_start,line_end,1):
        cl_temp = np.array([int(kk) for kk in lines[i].split()])
        clause_mat[i-line_start,np.absolute(cl_temp[0:k])-1] = np.sign(cl_temp[0:k])
        
    return clause_mat

def gen_kSATprob(nvar,nclause,k):
    
    rng = np.random.default_rng()
    clause_mat = np.zeros((nclause,nvar))
    literal_arr = np.arange(1,2*nvar+1,1)
    for i in range(nclause):
        pass_flag = 0
        while pass_flag==0:
            
            literal_cand = rng.choice(literal_arr,k,replace=False)
            inv_literal_cand = np.array(literal_cand)
            clause_entry = np.zeros((1,nvar))
            neg_lit_ind = np.where(literal_cand>nvar)[0]
            pos_lit_ind = np.where(literal_cand<nvar+1)[0]
            if len(neg_lit_ind)!=0:
                inv_literal_cand[neg_lit_ind] = inv_literal_cand[neg_lit_ind] - nvar
                clause_entry[0,inv_literal_cand[neg_lit_ind]-1] = -1
            if len(pos_lit_ind)!=0:
                inv_literal_cand[pos_lit_ind] = inv_literal_cand[pos_lit_ind] + nvar
                clause_entry[0,literal_cand[pos_lit_ind]-1] = 1
            
            tautological_arr = np.intersect1d(literal_cand,inv_literal_cand)
            if len(tautological_arr)==0:
                if i!=0:
                    clause_entry_rep = np.repeat(clause_entry,i,axis=0)
                    clause_mat_curr = clause_mat[0:i,:]
                    temp1 = np.absolute(clause_mat_curr - clause_entry_rep)
                    temp2 = np.sum(temp1,axis=1)
                    if len(np.where(temp2==0)[0])==0:
                        pass_flag = 1
                        clause_mat[i,] = clause_entry
                    else:
                        pass_flag = 0
                else:
                    pass_flag = 1
                    clause_mat[i,] = clause_entry
            else:
                pass_flag = 0
                
    return clause_mat
                
def clausemat2CNF(dir_name,clause_mat,instance_id):
    nvar = clause_mat.shape[1]
    nclause = clause_mat.shape[0]
    k = len(np.where(clause_mat[0,]!=0)[0])
    filename = dir_name+'uf'+str(nvar)+'-0'+str(instance_id)+'.cnf'
    file = open(filename, 'w')
    file.write('p cnf '+str(nvar)+' '+str(nclause)+'\n')
    for i in range(nclause):
        var_indices = np.where(clause_mat[i,]!=0)[0]
        line = ''
        for j in range(len(var_indices)):
            if clause_mat[i,var_indices[j]]>0:
                line = line+str(var_indices[j]+1)+' '
            else:
                line = line+str(-var_indices[j]-1)+' '
        
        line = line + '0\n'
        file.write(line)
    
    file.close()