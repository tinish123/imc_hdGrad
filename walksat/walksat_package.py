import numpy as np
import random
import sys
import math
from multiprocessing import Process
import concurrent.futures
import itertools

    
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

def G_heuristics(p,var_val,clause_mat,clause_f_mat,clause_ind,clause_stat_mat):
    
    nvar = len(var_val)
    nclause = clause_mat.shape[0]
    var_list = np.where(clause_mat[clause_ind,]!=0)[0]
    num_sat_clause_init = number_of_sat_clause(clause_stat_mat)
    clause_stat_change_mat = np.zeros((len(var_list),1))
    var_counter = 0
    xx = random.random()
    if xx <= p:
        var_candidate = var_list
    else:
        for i in var_list:

            clause_stat_change_mat[var_counter], _ = delta_E(i,clause_ind,clause_mat,clause_stat_mat,clause_f_mat)
            var_counter = var_counter + 1

        var_candidate = var_list[np.where(clause_stat_change_mat>0)[0]]
        if len(var_candidate)==0:
            var_candidate = var_list
        #else:
         #   var_candidate = var_list[np.where(clause_stat_change_mat==np.max(clause_stat_change_mat))[0]]
    
    chosen_var = np.random.choice(var_candidate,size=1)[0]
    var_val[chosen_var] = ~var_val[chosen_var]
    clause_stat_mat_new, sat_result_new, clause_f_mat_new = get_clause_stat(chosen_var,clause_ind,clause_mat,clause_stat_mat,clause_f_mat)
    
    return var_val, clause_stat_mat_new, sat_result_new, clause_f_mat_new

def Shonning_heuristics(p,var_val,clause_mat,clause_f_mat,clause_ind,clause_stat_mat):
    
    nvar = len(var_val)
    nclause = clause_mat.shape[0]
    var_list = np.where(clause_mat[clause_ind,]!=0)[0]    
    chosen_var = np.random.choice(var_list,size=1)[0]
    var_val[chosen_var] = ~var_val[chosen_var]
    clause_stat_mat_new, sat_result_new, clause_f_mat_new = get_clause_stat(chosen_var,clause_ind,clause_mat,clause_stat_mat,clause_f_mat)
    
    return var_val, clause_stat_mat_new, sat_result_new, clause_f_mat_new

#def get_stats(walksat_result, solvability_stat):
def get_stats(walksat_result):
        
    prob_s = np.mean(walksat_result[:,1])
    success_ind = np.where(walksat_result[:,1]==1)[0]
    avg_flips = np.mean(walksat_result[success_ind,0])
    std_flips = np.std(walksat_result[success_ind,0])
    
    max_restarts = np.shape(walksat_result)[0]
    p_targ = 0.99
    sorted_arr = walksat_result[walksat_result[:,0].argsort()]
    #print(sorted_arr)
    #prob_arr = sorted_arr[:,1]*np.arange(1,max_restarts+1,1)/max_restarts
    #prob_arr[np.where(prob_arr==0)] = np.max(prob_arr)
    #print(prob_arr)
    #tts_arr = sorted_arr[:,0]*(math.log((1-p_targ),10)/np.log10(1-prob_arr))
    #print(tts_arr)
    if prob_s>=p_targ:
        ind_tts = math.ceil(p_targ*max_restarts)-1
        tts_99 = sorted_arr[ind_tts,0]
    else:
        if prob_s == 0:
            tts_99 = 0
        else:
            tts_99 = sorted_arr[max_restarts-1,0]*(math.log((1-p_targ),10)/math.log((1-prob_s),10))
    
    #return avg_flips, prob_s, std_flips, tts_99, np.mean(solvability_stat[:,0])
    return avg_flips, prob_s, std_flips, tts_99

def SKC_heuristics(p,var_val,clause_mat,cnf_mat,clause_f_mat,clause_ind,clause_stat_mat):
    
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
    bv_vec = b_vec[lit1_indices]
   
    for i in var_list:

        clause_stat_change_mat[var_counter], var_break_value_mat[var_counter] = delta_E(i,clause_ind,clause_mat,clause_stat_mat,clause_f_mat)
        var_counter = var_counter + 1

    break_value_zero_ind = np.where(var_break_value_mat==0)[0]
    
    if len(break_value_zero_ind)!=0:
        chosen_var = var_list[np.random.choice(break_value_zero_ind,size=1)[0]]
    elif xx <= p:
        chosen_var = np.random.choice(var_list[np.where(clause_stat_change_mat==np.max(clause_stat_change_mat))[0]],size=1)[0]
    else:
        chosen_var = np.random.choice(var_list,size=1)[0]
    
    var_val[chosen_var] = ~var_val[chosen_var]
    clause_stat_mat_new, sat_result_new, clause_f_mat_new = get_clause_stat(chosen_var,clause_ind,clause_mat,clause_stat_mat,clause_f_mat)
    
    return var_val, clause_stat_mat_new, sat_result_new, clause_f_mat_new, bv_vec
    
def B_heuristics(p,var_val,clause_mat,cnf_mat,clause_f_mat,clause_ind,clause_stat_mat,lin_map,wta_bias_flag,biased_select_flag):
    
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
    #print(bv_vec)
    #hw_bv_vec = np.random.normal(sw_bv_vec, s)
    #print(sw_bv_vec.dtype)
    #print(np.array(list(map(lambda x: lin_map[x,0], list(np.int(sw_bv_vec))))))
    sw_bv_vec = sw_bv_vec.astype(int)
    #print(sw_bv_vec)
    #print(np.where(sw_bv_vec==0)[0])
    winning_bv = np.min(sw_bv_vec[var_list])
    hw_bv_vec = np.random.normal(np.array(list(map(lambda x: lin_map[x,0], list(sw_bv_vec))),dtype=float),np.array(list(map(lambda x: lin_map[x,1], list(sw_bv_vec)))))
    #print(hw_bv_vec)
    #print(hw_bv_vec[np.where(sw_bv_vec==0)[0]])
    #print(bv_vec)
    #print(hw_bv_vec)
    bv_0_1_thresh = (lin_map[0,0]+lin_map[1,0])/2
    hw_bv_vec[np.where(hw_bv_vec<=bv_0_1_thresh)] = 0
    #print(hw_bv_vec)
    
    var_break_value_mat = hw_bv_vec[var_list]
    #print("Potential BVs: ",var_break_value_mat)
    clause_stat_change_mat = gain_vec[var_list]
    
    # for i in var_list:

        # clause_stat_change_mat[var_counter], var_break_value_mat[var_counter] = delta_E(i,clause_ind,clause_mat,clause_stat_mat,clause_f_mat)
        # var_counter = var_counter + 1
 
    #print("Break-Value array from Old method: ",var_break_value_mat)
    #print("Break-Value array from New method: ",bv_vec[var_list])
    
    break_value_zero_ind = np.where(var_break_value_mat==0)[0]
    
    if len(break_value_zero_ind)!=0:
        if biased_select_flag==1:
            chosen_var = randVarSelect(var_list[break_value_zero_ind],nvar)
        else:
            chosen_var = var_list[np.random.choice(break_value_zero_ind,size=1)[0]]
    elif xx <= p:
        l11 = var_list[np.where(var_break_value_mat==np.min(var_break_value_mat))[0]]
        if wta_bias_flag==0:
            chosen_var = np.random.choice(var_list[np.where(var_break_value_mat==np.min(var_break_value_mat))[0]],size=1)[0]
        else:
            chosen_var = l11[0]
        #print("List of Candidate Variables when WTA: ",l11)
        #print("Chosen Variables when WTA: ",chosen_var)
    else:
        if biased_select_flag==1:
            chosen_var = chosen_var = randVarSelect(var_list,nvar)
        else:
            chosen_var = np.random.choice(var_list,size=1)[0]
        
    
    var_val[chosen_var] = ~var_val[chosen_var]
    clause_stat_mat_new, sat_result_new, clause_f_mat_new = get_clause_stat(chosen_var,clause_ind,clause_mat,clause_stat_mat,clause_f_mat)
    
    return var_val, clause_stat_mat_new, sat_result_new, clause_f_mat_new, sw_bv_vec, hw_bv_vec, winning_bv
    
def SKC2_heuristics(p,var_val,clause_mat,clause_f_mat,clause_ind,clause_stat_mat):
    
    nvar = len(var_val)
    nclause = clause_mat.shape[0]
    var_list = np.where(clause_mat[clause_ind,]!=0)[0]
    num_sat_clause_init = number_of_sat_clause(clause_stat_mat)
    clause_stat_change_mat = np.zeros((len(var_list),1))
    var_break_value_mat = np.zeros((len(var_list),1))
    var_counter = 0
    xx = random.random()
   
    for i in var_list:

        clause_stat_change_mat[var_counter], var_break_value_mat[var_counter] = delta_E(i,clause_ind,clause_mat,clause_stat_mat,clause_f_mat)
        var_counter = var_counter + 1

    break_value_zero_ind = np.where(var_break_value_mat==0)[0]
    
    if len(break_value_zero_ind)!=0:
        chosen_var = var_list[np.random.choice(break_value_zero_ind,size=1)[0]]
    elif xx <= p:
        chosen_var = np.random.choice(var_list[np.where(var_break_value_mat==np.min(var_break_value_mat))[0]],size=1)[0]
    else:
        chosen_var = np.random.choice(var_list,size=1)[0]
    
    var_val[chosen_var] = ~var_val[chosen_var]
    clause_stat_mat_new, sat_result_new, clause_f_mat_new = get_clause_stat(chosen_var,clause_ind,clause_mat,clause_stat_mat,clause_f_mat)
    
    return var_val, clause_stat_mat_new, sat_result_new, clause_f_mat_new
    
def gn_heuristics(p,var_val,cnf_mat,clause_mat,clause_f_mat,clause_stat_mat):
    
    nvar = len(var_val)
    nclause = cnf_mat.shape[0]
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
    gain_vec = a_vec[lit0_indices] - b_vec[lit1_indices]
    rand_vec = p*np.random.normal(0,1,(nvar,1))
    #rand_vec = p*np.random.uniform(0,1,(nvar,1))
    nonzero_ind = np.where(a_vec[lit0_indices]!=0)[0]
    #print(a_vec[lit0_indices].T)
    #print(gain_vec.T)
    #print(nonzero_ind)
    noisy_gain_vec = gain_vec + rand_vec
    noisy_gain_vec_v2 = gain_vec[nonzero_ind] + rand_vec[nonzero_ind]
    #print(noisy_gain_vec.T)
    #print(noisy_gain_vec_v2.T)
    #chosen_var = np.argmax(noisy_gain_vec)
    chosen_var = np.argmax(noisy_gain_vec_v2)
    chosen_var = nonzero_ind[chosen_var]
    
    var_init = var_val[chosen_var]
    var_val[chosen_var] = ~var_val[chosen_var]
    clause_stat_mat_new, sat_result_new, clause_f_mat_new = get_clause_stat_v2(chosen_var,var_init,clause_mat,clause_stat_mat,clause_f_mat)
    
    return var_val, clause_stat_mat_new, sat_result_new, clause_f_mat_new
    
def delta_E(i,clause_ind,clause_mat,clause_stat_mat,clause_f_mat):
    
    var_sign = clause_mat[clause_ind,i]
    A_sign = var_sign
    B_sign = -1*var_sign
    A_clause_ind1 = np.where(clause_mat[:,i]==A_sign)[0]
    A = len(np.where(clause_stat_mat[A_clause_ind1]==0)[0])
    B_clause_ind1 = np.where(clause_mat[:,i]==B_sign)[0]
    B = len(np.where(clause_f_mat[B_clause_ind1]==1)[0])
    gain = A-B
    
    return gain, B
    
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

def randVarSelect(var_list,nvar):
    num_candVar = len(var_list)
    filler_array = np.zeros((64-nvar,1))
    candVar_reg = np.zeros((nvar,1))
    candVar_reg[var_list,0] = 1
    candVar_l0 = np.array(np.concatenate((candVar_reg,filler_array),axis=0),dtype=bool)
    prngbits_randVar = np.array(np.random.randint(2, size=6),dtype=bool)
    #prngbits_randVar[5] = False
    
    candVar_l1 = np.zeros((32,1),dtype=bool)
    candVar_l2 = np.zeros((16,1),dtype=bool)
    candVar_l3 = np.zeros((8,1),dtype=bool)
    candVar_l4 = np.zeros((4,1),dtype=bool)
    candVar_l5 = np.zeros((2,1),dtype=bool)
    randVar_sel_l5 = np.zeros((2,1),dtype=bool)
    randVar_sel_l4 = np.zeros((4,1),dtype=bool)
    randVar_sel_l3 = np.zeros((8,1),dtype=bool)
    randVar_sel_l2 = np.zeros((16,1),dtype=bool)
    randVar_sel_l1 = np.zeros((32,1),dtype=bool)
    randVar_sel_l0 = np.zeros((64,1),dtype=bool)
    
    for j in range(32):
        candVar_l1[j,0] = candVar_l0[2*j,0] or candVar_l0[2*j+1,0]
    
    for j in range(16):
        candVar_l2[j,0] = candVar_l1[2*j,0] or candVar_l1[2*j+1,0]
    
    for j in range(8):
        candVar_l3[j,0] = candVar_l2[2*j,0] or candVar_l2[2*j+1,0]
    
    for j in range(4):
        candVar_l4[j,0] = candVar_l3[2*j,0] or candVar_l3[2*j+1,0]
    
    for j in range(2):
        candVar_l5[j,0] = candVar_l4[2*j,0] or candVar_l4[2*j+1,0]
    
    candVar_present = candVar_l5[0,0] or candVar_l5[1,0]
    randVar_sel_l5[0,0] = candVar_present and ((candVar_l5[0,0] and not(candVar_l5[1,0])) or (not(prngbits_randVar[5]) and candVar_l5[0,0] and candVar_l5[1,0]))
    randVar_sel_l5[1,0] = candVar_present and ((not(candVar_l5[0,0]) and candVar_l5[1,0]) or (prngbits_randVar[5] and candVar_l5[0,0] and candVar_l5[1,0]))
    
    for j in range(2):
        randVar_sel_l4[2*j,0] = randVar_sel_l5[j,0] and ((candVar_l4[2*j,0] and not(candVar_l4[2*j+1,0])) or (not(prngbits_randVar[4]) and candVar_l4[2*j,0] and candVar_l4[2*j+1,0]))
        randVar_sel_l4[2*j+1,0] = randVar_sel_l5[j,0] and ((not(candVar_l4[2*j,0]) and candVar_l4[2*j+1,0]) or (prngbits_randVar[4] and candVar_l4[2*j,0] and candVar_l4[2*j+1,0]))
    
    for j in range(4):
        randVar_sel_l3[2*j,0] = randVar_sel_l4[j,0] and ((candVar_l3[2*j,0] and not(candVar_l3[2*j+1,0])) or (not(prngbits_randVar[3]) and candVar_l3[2*j,0] and candVar_l3[2*j+1,0]))
        randVar_sel_l3[2*j+1,0] = randVar_sel_l4[j,0] and ((not(candVar_l3[2*j,0]) and candVar_l3[2*j+1,0]) or (prngbits_randVar[3] and candVar_l3[2*j,0] and candVar_l3[2*j+1,0]))
    
    for j in range(8):
        randVar_sel_l2[2*j,0] = randVar_sel_l3[j,0] and ((candVar_l2[2*j,0] and not(candVar_l2[2*j+1,0])) or (not(prngbits_randVar[2]) and candVar_l2[2*j,0] and candVar_l2[2*j+1,0]))
        randVar_sel_l2[2*j+1,0] = randVar_sel_l3[j,0] and ((not(candVar_l2[2*j,0]) and candVar_l2[2*j+1,0]) or (prngbits_randVar[2] and candVar_l2[2*j,0] and candVar_l2[2*j+1,0]))
    
    for j in range(16):
        randVar_sel_l1[2*j,0] = randVar_sel_l2[j,0] and ((candVar_l1[2*j,0] and not(candVar_l1[2*j+1,0])) or (not(prngbits_randVar[1]) and candVar_l1[2*j,0] and candVar_l1[2*j+1,0]))
        randVar_sel_l1[2*j+1,0] = randVar_sel_l2[j,0] and ((not(candVar_l1[2*j,0]) and candVar_l1[2*j+1,0]) or (prngbits_randVar[1] and candVar_l1[2*j,0] and candVar_l1[2*j+1,0]))
    
    for j in range(32):
        randVar_sel_l0[2*j,0] = randVar_sel_l1[j,0] and ((candVar_l0[2*j,0] and not(candVar_l0[2*j+1,0])) or (not(prngbits_randVar[0]) and candVar_l0[2*j,0] and candVar_l0[2*j+1,0]))
        randVar_sel_l0[2*j+1,0] = randVar_sel_l1[j,0] and ((not(candVar_l0[2*j,0]) and candVar_l0[2*j+1,0]) or (prngbits_randVar[0] and candVar_l0[2*j,0] and candVar_l0[2*j+1,0]))
    
    chosen_var = np.where(randVar_sel_l0==1)[0]
    return chosen_var[0]

def randClauseSelect(clause_stat_mat):
    nclause = clause_stat_mat.shape[0]
    filler_array = np.zeros((256-nclause,1))
    true0_reg = np.zeros((nclause,1))
    true0_reg[np.where(clause_stat_mat==0)[0],0] = 1
    true0_l0 = np.array(np.concatenate((true0_reg,filler_array),axis=0),dtype=bool)
    prngbits_randClause = np.array(np.random.randint(2, size=8),dtype=bool)
    
    true0_l1 = np.zeros((128,1),dtype=bool)
    true0_l2 = np.zeros((64,1),dtype=bool)
    true0_l3 = np.zeros((32,1),dtype=bool)
    true0_l4 = np.zeros((16,1),dtype=bool)
    true0_l5 = np.zeros((8,1),dtype=bool)
    true0_l6 = np.zeros((4,1),dtype=bool)
    true0_l7 = np.zeros((2,1),dtype=bool)
    randClause_sel_l7 = np.zeros((2,1),dtype=bool)
    randClause_sel_l6 = np.zeros((4,1),dtype=bool)
    randClause_sel_l5 = np.zeros((8,1),dtype=bool)
    randClause_sel_l4 = np.zeros((16,1),dtype=bool)
    randClause_sel_l3 = np.zeros((32,1),dtype=bool)
    randClause_sel_l2 = np.zeros((64,1),dtype=bool)
    randClause_sel_l1 = np.zeros((128,1),dtype=bool)
    randClause_sel_l0 = np.zeros((256,1),dtype=bool)
    
    for j in range(128):
        true0_l1[j,0] = true0_l0[2*j,0] or true0_l0[2*j+1,0]
    
    for j in range(64):
        true0_l2[j,0] = true0_l1[2*j,0] or true0_l1[2*j+1,0]
    
    for j in range(32):
        true0_l3[j,0] = true0_l2[2*j,0] or true0_l2[2*j+1,0]
    
    for j in range(16):
        true0_l4[j,0] = true0_l3[2*j,0] or true0_l3[2*j+1,0]
    
    for j in range(8):
        true0_l5[j,0] = true0_l4[2*j,0] or true0_l4[2*j+1,0]
    
    for j in range(4):
        true0_l6[j,0] = true0_l5[2*j,0] or true0_l5[2*j+1,0]
    
    for j in range(2):
        true0_l7[j,0] = true0_l6[2*j,0] or true0_l6[2*j+1,0]
    
    unsat = true0_l7[0,0] or true0_l7[1,0]
    randClause_sel_l7[0,0] = unsat and ((true0_l7[0,0] and not(true0_l7[1,0])) or (not(prngbits_randClause[7]) and true0_l7[0,0] and true0_l7[1,0]))
    randClause_sel_l7[1,0] = unsat and ((not(true0_l7[0,0]) and true0_l7[1,0]) or (prngbits_randClause[7] and true0_l7[0,0] and true0_l7[1,0]))
    
    for j in range(2):
        randClause_sel_l6[2*j,0] = randClause_sel_l7[j,0] and ((true0_l6[2*j,0] and not(true0_l6[2*j+1,0])) or (not(prngbits_randClause[6]) and true0_l6[2*j,0] and true0_l6[2*j+1,0]))
        randClause_sel_l6[2*j+1,0] = randClause_sel_l7[j,0] and ((not(true0_l6[2*j,0]) and true0_l6[2*j+1,0]) or (prngbits_randClause[6] and true0_l6[2*j,0] and true0_l6[2*j+1,0]))
    
    for j in range(4):
        randClause_sel_l5[2*j,0] = randClause_sel_l6[j,0] and ((true0_l5[2*j,0] and not(true0_l5[2*j+1,0])) or (not(prngbits_randClause[5]) and true0_l5[2*j,0] and true0_l5[2*j+1,0]))
        randClause_sel_l5[2*j+1,0] = randClause_sel_l6[j,0] and ((not(true0_l5[2*j,0]) and true0_l5[2*j+1,0]) or (prngbits_randClause[5] and true0_l5[2*j,0] and true0_l5[2*j+1,0]))
    
    for j in range(8):
        randClause_sel_l4[2*j,0] = randClause_sel_l5[j,0] and ((true0_l4[2*j,0] and not(true0_l4[2*j+1,0])) or (not(prngbits_randClause[4]) and true0_l4[2*j,0] and true0_l4[2*j+1,0]))
        randClause_sel_l4[2*j+1,0] = randClause_sel_l5[j,0] and ((not(true0_l4[2*j,0]) and true0_l4[2*j+1,0]) or (prngbits_randClause[4] and true0_l4[2*j,0] and true0_l4[2*j+1,0]))
    
    for j in range(16):
        randClause_sel_l3[2*j,0] = randClause_sel_l4[j,0] and ((true0_l3[2*j,0] and not(true0_l3[2*j+1,0])) or (not(prngbits_randClause[3]) and true0_l3[2*j,0] and true0_l3[2*j+1,0]))
        randClause_sel_l3[2*j+1,0] = randClause_sel_l4[j,0] and ((not(true0_l3[2*j,0]) and true0_l3[2*j+1,0]) or (prngbits_randClause[3] and true0_l3[2*j,0] and true0_l3[2*j+1,0]))
    
    for j in range(32):
        randClause_sel_l2[2*j,0] = randClause_sel_l3[j,0] and ((true0_l2[2*j,0] and not(true0_l2[2*j+1,0])) or (not(prngbits_randClause[2]) and true0_l2[2*j,0] and true0_l2[2*j+1,0]))
        randClause_sel_l2[2*j+1,0] = randClause_sel_l3[j,0] and ((not(true0_l2[2*j,0]) and true0_l2[2*j+1,0]) or (prngbits_randClause[2] and true0_l2[2*j,0] and true0_l2[2*j+1,0]))
    
    for j in range(64):
        randClause_sel_l1[2*j,0] = randClause_sel_l2[j,0] and ((true0_l1[2*j,0] and not(true0_l1[2*j+1,0])) or (not(prngbits_randClause[1]) and true0_l1[2*j,0] and true0_l1[2*j+1,0]))
        randClause_sel_l1[2*j+1,0] = randClause_sel_l2[j,0] and ((not(true0_l1[2*j,0]) and true0_l1[2*j+1,0]) or (prngbits_randClause[1] and true0_l1[2*j,0] and true0_l1[2*j+1,0]))
    
    for j in range(128):
        randClause_sel_l0[2*j,0] = randClause_sel_l1[j,0] and ((true0_l0[2*j,0] and not(true0_l0[2*j+1,0])) or (not(prngbits_randClause[0]) and true0_l0[2*j,0] and true0_l0[2*j+1,0]))
        randClause_sel_l0[2*j+1,0] = randClause_sel_l1[j,0] and ((not(true0_l0[2*j,0]) and true0_l0[2*j+1,0]) or (prngbits_randClause[0] and true0_l0[2*j,0] and true0_l0[2*j+1,0]))

    #print("Clause Status Mat: ",clause_stat_mat)
    
    chosen_clause = np.where(randClause_sel_l0==1)[0]
    #print("Chosen Clause: ",chosen_clause)
    return chosen_clause[0]

def WalkSat_Solver_sameInit(p,clause_mat,max_restarts,max_flips,heuristics,lin_map):

    nvar = clause_mat.shape[1]
    nclause = clause_mat.shape[0]
    cnf_mat = clausemat2cnfmat(clause_mat)
    sw_bv_vec_total = []
    hw_bv_vec_total = []
    winning_bv_total = []
    
    walksat_result1 = np.zeros((max_restarts,2))
    walksat_result2 = np.zeros((max_restarts,2))
    walksat_result3 = np.zeros((max_restarts,2))
    walksat_result4 = np.zeros((max_restarts,2))

    for r in range(max_restarts):

        var_init = rand_init(nvar)
        
        var_val_mat = var_init.copy()
        clause_stat_mat, sat_result, clause_f_mat = sat_eval(var_val_mat,clause_mat)
        for f in range(max_flips):
            if sat_result == 1:
                walksat_result1[r,0] = f
                walksat_result1[r,1] = 1
                break
            unsat_clause_list = np.where(clause_stat_mat==0)[0]
            chosen_clause_ind = np.random.choice(unsat_clause_list,size=1)[0]
            var_val_mat, clause_stat_mat, sat_result, clause_f_mat, _, _, _ = B_heuristics(0.5,var_val_mat,clause_mat,cnf_mat,clause_f_mat,chosen_clause_ind,clause_stat_mat,lin_map,0,0)

        if walksat_result1[r,1] == 0:
            clause_stat_mat, sat_result, clause_f_mat = sat_eval(var_val_mat,clause_mat)
            walksat_result1[r,0] = max_flips
            if sat_result == 1:
                walksat_result1[r,1] = 1  

        var_val_mat = var_init.copy()
        clause_stat_mat, sat_result, clause_f_mat = sat_eval(var_val_mat,clause_mat)
        for f in range(max_flips):
            if sat_result == 1:
                walksat_result2[r,0] = f
                walksat_result2[r,1] = 1
                break
            unsat_clause_list = np.where(clause_stat_mat==0)[0]
            chosen_clause_ind = np.random.choice(unsat_clause_list,size=1)[0]
            var_val_mat, clause_stat_mat, sat_result, clause_f_mat, _, _, _ = B_heuristics(p,var_val_mat,clause_mat,cnf_mat,clause_f_mat,chosen_clause_ind,clause_stat_mat,lin_map,1,0)

        if walksat_result2[r,1] == 0:
            clause_stat_mat, sat_result, clause_f_mat = sat_eval(var_val_mat,clause_mat)
            walksat_result2[r,0] = max_flips
            if sat_result == 1:
                walksat_result2[r,1] = 1  


        var_val_mat = var_init.copy()
        clause_stat_mat, sat_result, clause_f_mat = sat_eval(var_val_mat,clause_mat)
        for f in range(max_flips):
            if sat_result == 1:
                walksat_result3[r,0] = f
                walksat_result3[r,1] = 1
                break
            unsat_clause_list = np.where(clause_stat_mat==0)[0]
            chosen_clause_ind = randClauseSelect(clause_stat_mat)
            var_val_mat, clause_stat_mat, sat_result, clause_f_mat, _, _, _ = B_heuristics(p,var_val_mat,clause_mat,cnf_mat,clause_f_mat,chosen_clause_ind,clause_stat_mat,lin_map,0,1)

        if walksat_result3[r,1] == 0:
            clause_stat_mat, sat_result, clause_f_mat = sat_eval(var_val_mat,clause_mat)
            walksat_result3[r,0] = max_flips
            if sat_result == 1:
                walksat_result3[r,1] = 1


        var_val_mat = var_init.copy()
        clause_stat_mat, sat_result, clause_f_mat = sat_eval(var_val_mat,clause_mat)
        for f in range(max_flips):
            if sat_result == 1:
                walksat_result4[r,0] = f
                walksat_result4[r,1] = 1
                break
            unsat_clause_list = np.where(clause_stat_mat==0)[0]
            chosen_clause_ind = randClauseSelect(clause_stat_mat)
            var_val_mat, clause_stat_mat, sat_result, clause_f_mat, _, _, _ = B_heuristics(p,var_val_mat,clause_mat,cnf_mat,clause_f_mat,chosen_clause_ind,clause_stat_mat,lin_map,1,1)

        if walksat_result4[r,1] == 0:
            clause_stat_mat, sat_result, clause_f_mat = sat_eval(var_val_mat,clause_mat)
            walksat_result4[r,0] = max_flips
            if sat_result == 1:
                walksat_result4[r,1] = 1
    
    return walksat_result1, walksat_result2, walksat_result3, walksat_result4 
    

def WalkSat_Solver(p,clause_mat,max_restarts,max_flips,heuristics,lin_map,wta_bias_flag,biased_select_flag):
    
    nvar = clause_mat.shape[1]
    nclause = clause_mat.shape[0]
    cnf_mat = clausemat2cnfmat(clause_mat)
    sw_bv_vec_total = []
    hw_bv_vec_total = []
    winning_bv_total = []
    
    walksat_result = np.zeros((max_restarts,2))
    clause_evol = np.zeros((max_restarts,max_flips))

    for r in range(max_restarts):

        var_val_mat = rand_init(nvar)
        clause_stat_mat, sat_result, clause_f_mat = sat_eval(var_val_mat,clause_mat)
        
        for f in range(max_flips):

            if sat_result == 1:
                walksat_result[r,0] = f
                walksat_result[r,1] = 1
                break

            unsat_clause_list = np.where(clause_stat_mat==0)[0]
            if biased_select_flag==1:
                chosen_clause_ind = randClauseSelect(clause_stat_mat)
            else:
                chosen_clause_ind = np.random.choice(unsat_clause_list,size=1)[0]
            
            clause_evol[r,f] = len(unsat_clause_list)
            if heuristics == "SKC":
                var_val_mat, clause_stat_mat, sat_result, clause_f_mat, bv_vec = SKC_heuristics(p,var_val_mat,clause_mat,cnf_mat,clause_f_mat,chosen_clause_ind,clause_stat_mat)
            elif heuristics == "B":
                var_val_mat, clause_stat_mat, sat_result, clause_f_mat, sw_bv_vec, hw_bv_vec, winning_bv = B_heuristics(p,var_val_mat,clause_mat,cnf_mat,clause_f_mat,chosen_clause_ind,clause_stat_mat,lin_map,wta_bias_flag,biased_select_flag)
            elif heuristics == "SKC2":
                var_val_mat, clause_stat_mat, sat_result, clause_f_mat = SKC2_heuristics(p,var_val_mat,clause_mat,clause_f_mat,chosen_clause_ind,clause_stat_mat)
            elif heuristics == "G":
                var_val_mat, clause_stat_mat, sat_result, clause_f_mat = G_heuristics(p,var_val_mat,clause_mat,clause_f_mat,chosen_clause_ind,clause_stat_mat)
            else:
                var_val_mat, clause_stat_mat, sat_result, clause_f_mat = Shonning_heuristics(p,var_val_mat,clause_mat,clause_f_mat,chosen_clause_ind,clause_stat_mat)

            #print(bv_vec.shape)
            #print(list(bv_vec[:,0]))
            sw_bv_vec_total.extend(list(sw_bv_vec[:,0]))
            hw_bv_vec_total.extend(list(hw_bv_vec[:,0]))
            winning_bv_total.append(winning_bv.item())
            #print(hw_bv_vec_total)
            
        if walksat_result[r,1] == 0:
            clause_stat_mat, sat_result, clause_f_mat = sat_eval(var_val_mat,clause_mat)
            walksat_result[r,0] = max_flips
            if sat_result == 1:
                walksat_result[r,1] = 1  

    avg_flips, prob_s, std_flips, tts_99 = get_stats(walksat_result)
    
    return avg_flips, prob_s, std_flips, tts_99, walksat_result, sw_bv_vec_total, hw_bv_vec_total, winning_bv_total, clause_evol, 0

def WalkSat_Solver_lowMF(p,clause_mat,max_restarts1,max_restarts2,max_flips,heuristics,lin_map,wta_bias_flag,biased_select_flag):
    
    nvar = clause_mat.shape[1]
    nclause = clause_mat.shape[0]
    cnf_mat = clausemat2cnfmat(clause_mat)
    
    walksat_result1 = np.zeros((max_restarts1,2))
    
    for mm in range(max_restarts1):
        
        walksat_result2 = np.zeros((max_restarts2,2))
        
        for r in range(max_restarts2):

            var_val_mat = rand_init(nvar)
            clause_stat_mat, sat_result, clause_f_mat = sat_eval(var_val_mat,clause_mat)

            for f in range(max_flips):

                if sat_result == 1:
                    walksat_result2[r,0] = f
                    walksat_result2[r,1] = 1
                    break

                unsat_clause_list = np.where(clause_stat_mat==0)[0]
                if biased_select_flag==1:
                    chosen_clause_ind = randClauseSelect(clause_stat_mat)
                else:
                    chosen_clause_ind = np.random.choice(unsat_clause_list,size=1)[0]

                if heuristics == "SKC":
                    var_val_mat, clause_stat_mat, sat_result, clause_f_mat, bv_vec = SKC_heuristics(p,var_val_mat,clause_mat,cnf_mat,clause_f_mat,chosen_clause_ind,clause_stat_mat)
                elif heuristics == "B":
                    var_val_mat, clause_stat_mat, sat_result, clause_f_mat, _, _, _ = B_heuristics(p,var_val_mat,clause_mat,cnf_mat,clause_f_mat,chosen_clause_ind,clause_stat_mat,lin_map,wta_bias_flag,biased_select_flag)
                elif heuristics == "SKC2":
                    var_val_mat, clause_stat_mat, sat_result, clause_f_mat = SKC2_heuristics(p,var_val_mat,clause_mat,clause_f_mat,chosen_clause_ind,clause_stat_mat)
                elif heuristics == "G":
                    var_val_mat, clause_stat_mat, sat_result, clause_f_mat = G_heuristics(p,var_val_mat,clause_mat,clause_f_mat,chosen_clause_ind,clause_stat_mat)
                else:
                    var_val_mat, clause_stat_mat, sat_result, clause_f_mat = Shonning_heuristics(p,var_val_mat,clause_mat,clause_f_mat,chosen_clause_ind,clause_stat_mat)


            if walksat_result2[r,1] == 0:
                clause_stat_mat, sat_result, clause_f_mat = sat_eval(var_val_mat,clause_mat)
                walksat_result2[r,0] = max_flips
                if sat_result == 1:
                    walksat_result2[r,1] = 1
                if r==max_restarts2-1:
                    walksat_result1[mm,0] = np.sum(walksat_result2[:,0])
                    if np.sum(walksat_result2[:,1])>0:
                        walksat_result1[mm,1] = 1
            elif walksat_result2[r,1] == 1:
                walksat_result1[mm,0] = np.sum(walksat_result2[:,0])
                walksat_result1[mm,1] = 1
                break
                
        #print("Inner result array: ",walksat_result2)
        #print("Outer result array: ",walksat_result1)

    avg_flips, prob_s, std_flips, tts_99 = get_stats(walksat_result1)
    
    return avg_flips, prob_s, std_flips, tts_99, walksat_result1

def WalkSat_Solver_multiTrial(p,clause_mat,max_restarts1,max_restarts2,max_flips,heuristics,lin_map,wta_bias_flag,biased_select_flag):
    
    nvar = clause_mat.shape[1]
    nclause = clause_mat.shape[0]
    cnf_mat = clausemat2cnfmat(clause_mat)
    
    walksat_result1 = np.zeros((max_restarts1,2))
    
    for mm in range(max_restarts1):
        
        walksat_result2 = np.zeros((max_restarts2,2))
        
        for r in range(max_restarts2):

            var_val_mat = rand_init(nvar)
            clause_stat_mat, sat_result, clause_f_mat = sat_eval(var_val_mat,clause_mat)

            for f in range(max_flips):

                if sat_result == 1:
                    walksat_result2[r,0] = f
                    walksat_result2[r,1] = 1
                    break

                unsat_clause_list = np.where(clause_stat_mat==0)[0]
                if biased_select_flag==1:
                    chosen_clause_ind = randClauseSelect(clause_stat_mat)
                else:
                    chosen_clause_ind = np.random.choice(unsat_clause_list,size=1)[0]

                if heuristics == "SKC":
                    var_val_mat, clause_stat_mat, sat_result, clause_f_mat, bv_vec = SKC_heuristics(p,var_val_mat,clause_mat,cnf_mat,clause_f_mat,chosen_clause_ind,clause_stat_mat)
                elif heuristics == "B":
                    var_val_mat, clause_stat_mat, sat_result, clause_f_mat, _, _, _ = B_heuristics(p,var_val_mat,clause_mat,cnf_mat,clause_f_mat,chosen_clause_ind,clause_stat_mat,lin_map,wta_bias_flag,biased_select_flag)
                elif heuristics == "SKC2":
                    var_val_mat, clause_stat_mat, sat_result, clause_f_mat = SKC2_heuristics(p,var_val_mat,clause_mat,clause_f_mat,chosen_clause_ind,clause_stat_mat)
                elif heuristics == "G":
                    var_val_mat, clause_stat_mat, sat_result, clause_f_mat = G_heuristics(p,var_val_mat,clause_mat,clause_f_mat,chosen_clause_ind,clause_stat_mat)
                else:
                    var_val_mat, clause_stat_mat, sat_result, clause_f_mat = Shonning_heuristics(p,var_val_mat,clause_mat,clause_f_mat,chosen_clause_ind,clause_stat_mat)


            if walksat_result2[r,1] == 0:
                clause_stat_mat, sat_result, clause_f_mat = sat_eval(var_val_mat,clause_mat)
                walksat_result2[r,0] = max_flips
                if sat_result == 1:
                    walksat_result2[r,1] = 1
                        
        walksat_result1[mm,0] = np.min(walksat_result2[np.where(walksat_result2[:,1]==1)[0],0])
        if np.sum(walksat_result2[:,1])>0:
            walksat_result1[mm,1] = 1

    avg_flips, prob_s, std_flips, tts_99 = get_stats(walksat_result1)
    
    return avg_flips, prob_s, std_flips, tts_99, walksat_result1


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

def WalkSat_Solver_full(nvar,i,params):
    
    max_restarts = params['max_res']
    max_flips = params['max_flip']
    p = params['p']
    var_heur = params['var_heur']
    lin_map = params['lin_map']
    wta_bias_flag = params['wta_bias_flag']
    biased_select_flag = params['biased_select_flag']
    clause_mat = get_kSATprob(nvar,i)
    avg_flips, prob_s, std_flips, tts_99, walksat_result, _ , _, _, _, ttsv_98 = WalkSat_Solver(p,clause_mat,max_restarts,max_flips,var_heur,lin_map,wta_bias_flag,biased_select_flag)
    return avg_flips, prob_s, std_flips, i, tts_99, ttsv_98, walksat_result

def WalkSat_Solver_full_lowMF(nvar,i,params):
    
    max_restarts1 = params['max_res1']
    max_restarts2 = params['max_res2']
    max_flips = params['max_flip']
    p = params['p']
    var_heur = params['var_heur']
    lin_map = params['lin_map']
    wta_bias_flag = params['wta_bias_flag']
    biased_select_flag = params['biased_select_flag']
    clause_mat = get_kSATprob(nvar,i)
    avg_flips, prob_s, std_flips, tts_99, walksat_result = WalkSat_Solver_lowMF(p,clause_mat,max_restarts1,max_restarts2,max_flips,var_heur,lin_map,wta_bias_flag,biased_select_flag)
    return avg_flips, prob_s, std_flips, i, tts_99, walksat_result

def WalkSat_Solver_full_multiTrial(nvar,i,params):
    
    max_restarts1 = params['max_res1']
    max_restarts2 = params['max_res2']
    max_flips = params['max_flip']
    p = params['p']
    var_heur = params['var_heur']
    lin_map = params['lin_map']
    wta_bias_flag = params['wta_bias_flag']
    biased_select_flag = params['biased_select_flag']
    clause_mat = get_kSATprob(nvar,i)
    avg_flips, prob_s, std_flips, tts_99, walksat_result = WalkSat_Solver_multiTrial(p,clause_mat,max_restarts1,max_restarts2,max_flips,var_heur,lin_map,wta_bias_flag,biased_select_flag)
    return avg_flips, prob_s, std_flips, i, tts_99, walksat_result
    
def WalkSat_Solver_full_v3(nvar,clause_mat,params):
    
    max_restarts = params['max_res']
    max_flips = params['max_flip']
    p = params['p']
    var_heur = params['var_heur']
    avg_flips, prob_s, std_flips, tts_99 = WalkSat_Solver(p,clause_mat,max_restarts,max_flips,var_heur)
    return avg_flips, prob_s, std_flips, tts_99  
        
    
# def WalkSat_Solver_full_v2(nvar,i,params):
    
    # max_restarts = params['max_res']
    # max_flips = params['max_flip']
    # p = params['p']
    # var_heur = params['var_heur']
    # clause_mat = get_kSATprob(nvar,i)
    # avg_flips, prob_s, std_flips, tts_99 = WalkSat_Solver_v2(p,clause_mat,max_restarts,max_flips,var_heur)
    # return avg_flips, prob_s, std_flips, i, tts_99
    
# def WalkSat_Solver_full_v2(nvar,i,params):
    
    # max_restarts = params['max_res']
    # max_flips = params['max_flip']
    # p = params['p']
    # var_heur = params['var_heur']
    # clause_mat = get_kSATprob(nvar,i)
    # avg_flips, prob_s, std_flips, tts_99 = WalkSat_Solver_iterated(p,clause_mat,max_restarts,max_flips,var_heur)
    # return avg_flips, prob_s, std_flips, tts_99

def get_kSATprob(nvar,instance):
    
    if nvar==20:
        dir_name = 'uf20-91'
        instance_name = 'uf20-0'+str(instance)
        if instance > 1000:
            print("Instance Number not present in Instance Directory!!")
            sys.exit()
    elif nvar==50:
        dir_name = 'uf50-218'
        instance_name = 'uf50-0'+str(instance)
        if instance > 1000:
            print("Instance Number not present in Instance Directory!!")
            sys.exit()
    elif nvar==75:
        dir_name = 'uf75-325'
        instance_name = 'uf75-0'+str(instance)
        if instance > 100:
            print("Instance Number not present in Instance Directory!!")
            sys.exit()
    elif nvar==100:
        dir_name = 'uf100-430'
        instance_name = 'uf100-0'+str(instance)
        if instance > 1000:
            print("Instance Number not present in Instance Directory!!")
            sys.exit()
    elif nvar==150:
        dir_name = 'uf150-645'
        instance_name = 'uf150-0'+str(instance)
        if instance > 100:
            print("Instance Number not present in Instance Directory!!")
            sys.exit()
    elif nvar==200:
        dir_name = 'uf200-860'
        instance_name = 'uf200-0'+str(instance)
        if instance > 100:
            print("Instance Number not present in Instance Directory!!")
            sys.exit()
    elif nvar==225:
        dir_name = 'uf225-960'
        instance_name = 'uf225-0'+str(instance)
        if instance > 100:
            print("Instance Number not present in Instance Directory!!")
            sys.exit()
    elif nvar==250:
        dir_name = 'uf250-1065'
        instance_name = 'uf250-0'+str(instance)
        if instance > 100:
            print("Instance Number not present in Instance Directory!!")
            sys.exit()
    else:
        dir_name = 'uf20-91'
        instance_name = 'uf20-0'+str(instance)
        if instance > 1000:
            print("Instance Number not present in Instance Directory!!")
            sys.exit()
        
    file_addr = '/home/tinish/SATLIB_Data/'+dir_name+'/'+instance_name+'.cnf'
    
    with open(file_addr) as f:
        lines = f.readlines()
    f.close()
    k = int(lines[5].split(' ')[7])
    nvar = int(lines[7].split(' ')[2])
    nclause = int(lines[7].split(' ')[4])
    line_start = 8
    line_end = line_start + nclause
    clause_mat = np.zeros((nclause,nvar))

    for i in np.arange(line_start,line_end,1):
        cl_temp = np.array([int(kk) for kk in lines[i].split()])
        clause_mat[i-line_start,np.absolute(cl_temp[0:k])-1] = np.sign(cl_temp[0:k])
        
    return clause_mat
    
    # p_targ = 0.99
    # sorted_array = walksat_result[walksat_result.argsort()]
    # if SR >= p_targ:
        # index_tts = math.ceil(p_targ*MAX_TRIES)-1
        # tts_99 = sorted_array[index_tts]
    # else:
        # tts_99 = sorted_array[MAX_TRIES-1]*(math.log((1-p_targ),10)/math.log((1-SR),10))
        
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
            
            tautological_arr =  np.intersect1d(literal_cand,inv_literal_cand)
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
    
def get_satInstance(nvar,nclause,k):
    sat_flag = 0
    num_trials = 0
    while sat_flag==0:
        clause_mat = gen_kSATprob(nvar,nclause,k)
        avg_flips, prob_s, std_flips, tts_99 = WalkSat_Solver(0.5,clause_mat,2,10000,"SKC")
        num_trials = num_trials + 1
        if prob_s!=0:
            sat_flag = 1
            
    return clause_mat
    
def get_someInstances(nvar,nclause,k,num):
    gen_result = np.zeros((num,2))
    for i in range(num):
        clause_mat = gen_kSATprob(nvar,nclause,k)
        avg_flips, prob_s, std_flips, tts_99 = WalkSat_Solver(0.5,clause_mat,1,500000,"SKC")
        if prob_s!=0:
            gen_result[i,0] = 1
            gen_result[i,1] = avg_flips
        else:
            gen_result[i,0] = 0
            gen_result[i,1] = 0
            
    return gen_result