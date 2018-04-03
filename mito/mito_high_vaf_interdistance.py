
# coding: utf-8

# In[14]:

import matplotlib.pyplot as plt
import math 
get_ipython().magic('matplotlib inline')


# In[21]:

high_vaf_pos = [72,152,204,248,263,302,310,316,515,567,750,1438,2487,3109,3492,6419,7028,10277,10306,12684,12705,12825,13062,13095,13105,13650,15326,16129,16183,16189,16218,16230,16249,16259,16263,16264,16274,16278,16284,16288,16293,16301,16311,16355,16356,16368,16390,16399,16427,16444,16496,16519,16527]
def rainfall(poslist):
    dist_list = [0]
    for i in range(0, len(poslist)):
        if i != len(poslist) -1:
            dist = poslist[i+1] - poslist[i]
            dist_list.append(math.log10(dist))
        
    return dist_list


# In[25]:

plt.scatter(high_vaf_pos, dist_list)

for i in range(0, len(high_vaf_pos)):
    print(f'hsMT\t{high_vaf_pos[i]}\t{high_vaf_pos[i]}\t{dist_list[i]}')


# In[4]:

pos = [int (i) for i in '''302,
 352,
 432,
 493,
 533,
 540,
 556,
 567,
 636,
 801,
 866,
 955,
 1082,
 1222,
 1374,
 1530,
 1836,
 1900,
 2487,
 3565,
 3584,
 3892,
 4136,
 4317,
 4390,
 4435,
 4833,
 5436,
 5600,
 5843,
 5894,
 6419,
 6565,
 6583,
 7231,
 7465,
 7490,
 8271,
 8280,
 8405,
 8931,
 9254,
 9531,
 9570,
 10161,
 10277,
 10935,
 10946,
 10978,
 11511,
 11535,
 11672,
 11866,
 12108,
 12193,
 12431,
 12967,
 13025,
 13057,
 13127,
 13676,
 13691,
 14280,
 14339,
 14388,
 14417,
 14424,
 14587,
 14619,
 14628,
 14750,
 14769,
 14813,
 15366,
 15525,
 15536,
 16183,
 16258,
 16374'''.split(',')]


# In[23]:

distance = rainfall(pos)
plt.scatter(pos, distance)


