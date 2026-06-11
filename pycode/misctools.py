def ext_tup(t,x):
    return tuple(list(t)+[x])

def rotate_tuple(v:tuple,k:int):
    k = k % len(v)
    if k == 0:
        return v
    return tuple(list(v[k:])+list(v[:k]))

def sort_tuple(v:tuple[int],l_scores:dict={}):
    if len(v)<2:
        return v
    lv = list(v)
    if len(l_scores)==0:
        lv.sort()
    else:
        lv.sort(key=lambda l:l_scores[l])
    return tuple(lv)

def merge_dicts(diclist:list[dict])->dict:
    if len(diclist)==0:
        return {}
    dic = diclist[0]
    for dic1 in diclist[1:]:
        for x in dic1:
            if x in dic:
                assert dic[x] == dic1[x]
            else:
                dic[x]=dic1[x]
    return dic

def invert_dict(dict):
    return {dict[x]:x for x in dict}

# Given two dictionaries with the same keys, returns
# a dictionary whose keys are keys of dict, and values are keys of second dict
def compdiv_dics(dic1:dict,dic2:dict)->dict:
    if len(dic1)!=len(dic2):
        raise ValueError('Dictionaries must have same size')
    common_elts = {d for d in dic1 if d in dic2}
    if len(common_elts)<len(dic1):
        raise ValueError('Dictionaries should have same keys')
    dic1inv = {dic1[x]:x for x in dic1}
    return {y:dic2[dic1inv[y]] for y in dic1inv}