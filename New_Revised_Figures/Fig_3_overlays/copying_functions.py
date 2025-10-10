import numpy as np

factor = 1.2 # in the main text f=1/factor=1/1.2
b = -2 # for anticonformist, b=2 for conformist

def compute_D(greater_vals, b, factor):
    if b>=1:
        y_max = np.max(greater_vals)
        D_max = (1/y_max - 1)*sum(greater_vals)
        return D_max/factor
    elif b<1:
        D_min = -sum(greater_vals)
        return D_min/factor
    

def switching_probabilities(y, b, copying_type=1):
    """ computes switching probabilities based on the type of mate copying: 
    1: Type I
    2: Type II 
    3: Type II (mixed)
    All types are conformist if b>1, anticonformist if b<-1
    Output: array probabilities of switching to a given morph"""

    if copying_type==1:
        probs = (y**b)/(y**b + (1-y)**b)
        
    elif copying_type==2:
        if b==1: probs = y # freq proportional copying
        else: 
            majority = 1/len(y)
            vec = np.array(y)
            
            # partition the frequencies
            greater_idx = np.where(vec > majority)[0]
            equal_idx = np.where(vec == majority)[0]
            less_idx = np.where(vec < majority)[0]

            greater_vals = vec[greater_idx]
            equal_vals = vec[equal_idx]
            less_vals = vec[less_idx]

            # compute value of D
            if len(greater_vals)==0: D = 0
            else: D = compute_D(greater_vals, b, factor)

            # convert to copying probabilities
            # greater than majority
            greater_vals = greater_vals + (greater_vals/sum(greater_vals))*D
            # less than majority
            inv_y = 1/less_vals
            factors = inv_y / np.sum(inv_y)
            less_vals = less_vals - factors*D
            
            # construct vector of copying probabilities
            probs = np.empty_like(vec)
            probs[greater_idx] = greater_vals
            probs[equal_idx] = equal_vals
            probs[less_idx] = less_vals
        #else: 
        #    print("Invalid value of copying parameter b, assuming neutral copying")
        #    probs = y
    elif copying_type==3:
        majority = 1/len(y)
        vec = np.array(y)
    
        # partition the frequencies 
        greater_idx = np.where(vec > majority)[0]
        equal_idx = np.where(vec == majority)[0]
        less_idx = np.where(vec < majority)[0]

        greater_vals = vec[greater_idx]
        equal_vals = vec[equal_idx]
        less_vals = vec[less_idx]

        # compute value of D
        if len(greater_vals)==0: 
            D_conf = 0
            D_anti = 0
        else: 
            D_conf = compute_D(greater_vals, b=2, factor = factor)
            D_anti = compute_D(greater_vals, b=-2, factor = factor)

        if len(y)==2:
            threshold = 0.7
        else: # for more than 2 morphs
            threshold = 0.5

        if len(greater_vals)!=0: 
            #mask_conf = greater_vals >= threshold
            #mask_anti = greater_vals < threshold
            mask_conf = greater_vals < threshold
            mask_anti = greater_vals >= threshold

            # convert to copying probabilities
            # greater than majority (adjust according to threshold)
            greater_vals[mask_conf] += (greater_vals[mask_conf] / np.sum(greater_vals)) * D_conf
            greater_vals[mask_anti] += (greater_vals[mask_anti] / np.sum(greater_vals)) * D_anti
            # less than majority
            inv_y = np.zeros_like(less_vals, dtype=float)
            mask = np.abs(less_vals) > 1e-5
            inv_y[mask] = 1 / less_vals[mask]

            if np.sum(inv_y) == 0:
                factors = np.zeros_like(inv_y)
            else:
                factors = inv_y / np.sum(inv_y)
        
            less_vals = less_vals - factors * D_conf
            
        # construct vector of copying probabilities
        probs = np.empty_like(vec)
        probs[greater_idx] = greater_vals
        probs[equal_idx] = equal_vals
        probs[less_idx] = less_vals

        probs = np.where(probs > 1e-5, probs, 0)
        #print(probs)           
        #print(y_t)

        #normalize
        prob_sum = sum(probs)
        probs = probs/prob_sum
              
    return probs

result = switching_probabilities(np.array([0.6, 0.3, 0.1]), b, copying_type=2)
print(result)