# Linear combination of a series of values
# For values [v_0, v_1,... v_k] returns:
# v_0 + challenge * v_1 + ... + challenge^k  * v_k
def Multiset_lc(values, challenge):
    kth_val = values.elements[-1]
    for val in reversed(values.elements[:-1]):
        for i in range(len(kth_val)):
            kth_val[i] = kth_val[i].mul(challenge)
            kth_val[i] = kth_val[i].add(val[i])
    return kth_val

def lc(values:list, challenge):
    kth_val = values[-1]
    for val in reversed(values[:-1]):
        kth_val = kth_val.mul(challenge)
        kth_val = kth_val.add(val)
    return kth_val

