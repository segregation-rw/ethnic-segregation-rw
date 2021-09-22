"""
This program is free software. You can use, copy, distribute, modify and/or redistribute
it under the terms of the MIT/Expat License. See the file LICENSE for more details.

(c) 2021 Sandro Sousa

If you use this software please consider citing the original paper:

S. Sousa, V. Nicosia
"Quantifying ethnic segregation in cities through random walks"


-----------------------------------------------------------------------------

Compute the class coverage time CCT as a function of c fractions of classes
Values correspond to the number of unique classes seen
at time t averaged over <num> iterations. The time series Y(t) shows
the time evolution of the classes seen for each node.
<edges> reported as DIRECTED, no self-loops, no multiple edges.
<C> number of classes to build the MxC matrix
<num> number of walk repetions for computing the average for each node
<idx> range of nodes <0-10>, full graph <g> or node id
<mode> select the spatial pattern of the population distribution

> output:
    node
    list: [avg trajectory (n classes seens) at each fraction c]

-----------------------------------------------------------------------------
"""

import numpy as np
import random
from itertools import zip_longest
import sys
# from tqdm import tqdm


if len(sys.argv) < 6:
    print("Usage: %s <edges> <C> <num> <idx> <mode>\n" % sys.argv[0])
    exit(1)


def null_model_grid(M, C):
    # null model w 8 cells per class placed at random
    pool = list(range(M))
    mx = np.zeros((M, C), dtype=int)
    for i in range(C):
        for _ in range(int(M/C)):
            loc = np.random.choice(pool)
            pool.remove(loc)
            mx[loc, i] = i+1
    return mx


def null_model_5_classes(M, C):
    # null model w 4x63 cells per class and 1x4 cells, tail graph
    pool = list(range(M))
    mx = np.zeros((M, C), dtype=int)
    for _ in range(63):
        for i in range(4):
            loc = np.random.choice(pool)
            pool.remove(loc)
            mx[loc, i] = 1
    mx[pool, 4] = 1
    return mx


def null_model_5_classes_16c(M, C):
    # null model w 4x63 cells per class and 1x4 cells, tail graph
    pool = list(range(M))
    mx = np.zeros((M, C), dtype=int)
    # populate the 4 larger colours
    for _ in range(60):
        for i in range(4):
            loc = np.random.choice(pool)
            pool.remove(loc)
            mx[loc, i] = 1
    # populate the 5th colour
    mx[pool, 4] = 1
    return mx


def cluster_8_cells(M, C):
    # 1x8 1 block of 8 cells for each class, 32 classes
    L = int(np.sqrt(M))
    mx = np.zeros((M, C), dtype=int)
    nodes = [i*L+j for j in range(0,L,2) for i in range(0,L,4)]
    pool = [k for k in range(C)]
    for col in pool:
        i = np.random.choice(nodes)
        nodes.remove(i)
        mx[i:i+2, col] = col+1
        mx[i+L:i+L+2, col] = col+1
        mx[i+L*2:i+L*2+2, col] = col+1
        mx[i+L*3:i+L*3+2, col] = col+1
    return mx


def cluster_4_cells(M, C):
    # 2x4 2 blocks of 4 cells each class, 32 classes
    L = int(np.sqrt(M))
    mx = np.zeros((M, C), dtype=int)
    nodes = [i*L+j for j in range(0,L,2) for i in range(0,L,2)]
    pool = [k for k in range(C)] * 2
    for col in pool:
        i = np.random.choice(nodes)
        nodes.remove(i)
        mx[i:i+2, col] = col+1
        mx[i+L:i+L+2, col] = col+1
    return mx


def stripe_4_cells(M, C):
    # 2x4x1 2 stripes of 4 cells each class, 32 classes
    L = int(np.sqrt(M))
    mx = np.zeros((M, C), dtype=int)
    nodes = [i*L+j for j in range(0,L) for i in range(0,L,4)]
    pool = [k for k in range(C)] * 2
    for col in pool:
        i = np.random.choice(nodes)
        nodes.remove(i)
        mx[i, col] = col+1
        mx[i+L, col] = col+1
        mx[i+L*2, col] = col+1
        mx[i+L*3, col] = col+1
    return mx


def stripe_2_cells(M, C):
    # 4x2 4 stripes of 2 cells each class, 32 classes
    L = int(np.sqrt(M))
    mx = np.zeros((M, C), dtype=int)
    nodes = [i*L+j for j in range(0,L) for i in range(0,L,2)]
    pool = [k for k in range(C)] * 4
    for col in pool:
        i = np.random.choice(nodes)
        nodes.remove(i)
        mx[i, col] = col+1
        mx[i+L, col] = col+1
    return mx


def cluster_8_cells_tail(M, C):
    # 1x8 1 block of 8 cells each class, tail graph
    L = 14
    mx = np.zeros((M, C), dtype=int)
    nodes = [i*L+j for j in range(0,L,2) for i in range(0,16,4)] + \
            [224, 226, 228, 230]
    pool = [k for k in range(C)]
    for col in pool:
        i = np.random.choice(nodes)
        nodes.remove(i)
        L = 14
        if i >= 224:
            L = 8
        mx[i:i+2, col] = col+1
        mx[i+L:i+L+2, col] = col+1
        mx[i+L*2:i+L*2+2, col] = col+1
        mx[i+L*3:i+L*3+2, col] = col+1
    return mx


def cluster_4_cells_tail(M, C):
    # 2x34 2 blocks of 4 cells each class, tail graph
    L = 14
    mx = np.zeros((M, C), dtype=int)
    pool = [k for k in range(C)] * 2
    nodes = [i*L+j for j in range(0,L,2) for i in range(0,16,2)] + \
            [224, 226, 228, 230, 240, 242, 244, 246]
    for col in pool:
        i = np.random.choice(nodes)
        nodes.remove(i)
        L = 14
        if i >= 224:
            L = 8
        mx[i:i+2, col] = col+1
        mx[i+L:i+L+2, col] = col+1
    return mx


def stripe_4_cells_tail(M, C):
    # 2x4x1 2 stripes of 4 cells each class, tail graph
    L = 14
    mx = np.zeros((M, C), dtype=int)
    nodes = [i*L+j for j in range(0,L) for i in range(0,16,4)] + \
            list(range(224,232))
    pool = [k for k in range(C)] * 2
    for col in pool:
        i = np.random.choice(nodes)
        nodes.remove(i)
        L = 14
        if i >= 224:
            L = 8
        mx[i, col] = col+1
        mx[i+L, col] = col+1
        mx[i+L*2, col] = col+1
        mx[i+L*3, col] = col+1
    return mx


def stripe_2_cells_tail(M, C):
    # 4x2 4 stripes of 2 cells each class, tail graph
    L = 14
    mx = np.zeros((M, C), dtype=int)
    pool = [k for k in range(C)] * 4
    nodes = [i*L+j for j in range(0,L) for i in range(0,16,2)]  + \
            list(range(224,232)) + list(range(240,248))
    for col in pool:
        i = np.random.choice(nodes)
        nodes.remove(i)
        L = 14
        if i >= 224:
            L = 8
        mx[i, col] = col+1
        mx[i+L, col] = col+1
    return mx


def cluster_63_4_cells(M, C, mask):
    # 4x63 4 blocks 63 cells and a class w 4 cells
    mx = np.zeros((M, C), dtype=int)
    L = int(np.sqrt(M))
    for i in range(L):
        if i*L < 127:
            # fill top left and right squares
            mx[L*i:L*i+8, 0] = 1
            mx[L*i+8:L*i+16, 1] = 2
        else:
            # fill bottom left and right squares
            mx[L*i:L*i+8, 2] = 3
            mx[L*i+8:L*i+16, 3] = 4
    # reposition block of 4
    mx[mask, 4] = 5
    mx[mask, :4] = 0
    if 14 in mask or 42 in mask:
        mx[[119,135,136], :] = 0
        mx[[119,135,136], 1] = 2
    return mx


def cluster_60_16_cells(M, C, mask):
    # 4x60 4 blocks 60 cells and a class w 16 cells
    mx = np.zeros((M, C), dtype=int)
    L = int(np.sqrt(M))
    for i in range(L):
        if i*L < 127:
            # fill top left and right squares
            mx[L*i:L*i+8, 0] = 1
            mx[L*i+8:L*i+16, 1] = 2
        else:
             # fill bottom left and right squares
            mx[L*i:L*i+8, 2] = 3
            mx[L*i+8:L*i+16, 3] = 4
    # reposition block of 16
    temp = []
    for v in mask:
        temp += [v,v+1,v+L,v+L+1]
    mx[temp, :] = 0
    mx[temp, 4] = 5
    # adjust 2nd colour size for uniformity on 1 quadrant cases
    if 12 in mask or 25 in mask:
        recolour = []
        for v in [102,134,136]:
            recolour += [v,v+1,v+L,v+L+1]
        mx[recolour, :] = 0
        mx[recolour, 1] = 2
    return mx


def cluster_2_cells_scaled(M, C):
    # 1x2 1 cluster of 2 cells each class, 32 classes
    L = int(np.sqrt(M))
    mx = np.zeros((M, C), dtype=int)
    nodes = [i*L+j for j in range(0,L) for i in range(0,L,2)]
    pool = [k for k in range(C)]
    for col in pool:
        i = np.random.choice(nodes)
        nodes.remove(i)
        mx[i, col] = col+1
        mx[i+L, col] = col+1
    return mx


def cluster_32_cells_scaled(M, C):
    # 1x32 1 block of 32 cells for each class, 32 classes
    L = int(np.sqrt(M))
    mx = np.zeros((M, C), dtype=int)
    nodes = [i*L+j for j in range(0,L,4) for i in range(0,L,8)]
    pool = [k for k in range(C)]
    for col in pool:
        i = np.random.choice(nodes)
        nodes.remove(i)
        for j in range(i,i+L*8,32):
            mx[j:j+4, col] = col+1
    return mx


def scale_up_matrix(mx):
    L = int(np.sqrt(mx.shape[0])) # side
    # duplicate fist half, mirror to right
    mx_new = []
    for i in range(L):
        tmp = np.concatenate((mx[L*i:(L*i)+L,:], mx[L*i:(L*i)+L,:]), axis=0)
        if len(mx_new) == 0:
            mx_new = tmp
        else:
            mx_new = np.concatenate((mx_new, tmp), axis=0)
    scaled = np.concatenate((mx_new, mx_new), axis=0)
    return scaled


def compute_neighbours(edges):
    # dictionary with node neighbours
    E = dict()
    for i, j in edges:
        if i in E:
            E[i].append(j)
        else:
            E[i] = [j]
    return E


def walk(i, obs):
    nc = []  # number of classes seen
    while 0 in obs:
        # selects one of i neighbours
        j = random.choice(neigh_dict[i])
        obs = obs + classes[j]
        # get current number of classes seen
        non_z = int(np.count_nonzero(obs))
        nc.append(non_z)
        i = j
    return nc


# initialize dict with neighbours
neigh_dict = compute_neighbours(np.loadtxt(sys.argv[1], dtype='int'))
C = int(sys.argv[2])  # number of classes
num = int(sys.argv[3])  # number of iterations
M = max(neigh_dict) + 1  # number of nodes
mode = int(sys.argv[5]) # numeric code for model to run

# process range or single node from param
if "-" in sys.argv[4]:
    start, end = [int(x) for x in sys.argv[4].split('-')]
if sys.argv[4] == 'g':
    start, end = (0, M)
if sys.argv[4].isnumeric():
    start, end = [int(sys.argv[4]), int(sys.argv[4])+1]


## --------------------------------------##
if mode == 100:
    # 100 = null model 8 cells per class, placed randomly
    classes = null_model_grid(M, C)
if mode == 200:
    # 200 = null model w 4x63 cells per class and one w 4 cells
    classes = null_model_5_classes(M, C)
if mode == 300:
    # 300 = null model reading nodes file
    classes = np.loadtxt(sys.argv[6], dtype='int')[:, 1:]
if mode == 400:
    # 400 = null model w 4x60 cells per class and one w 16 cells
    classes = null_model_5_classes_16c(M, C)
## --------------------------------------##
if mode == 118:
    # 118 = 1x8 1 block of 8 cells each class
    classes = cluster_8_cells(M, C)
if mode == 124:
    # 124 = 2x4 2 blocks of 4 cells each class
    classes = cluster_4_cells(M, C)
if mode == 141:
    # 141 = 2x4 2 stripes of 4 cells each class
    classes = stripe_4_cells(M, C)
if mode == 142:
    # 142 = 4x2 4 stripes of 2 cells each class
    classes = stripe_2_cells(M, C)
## --------------------------------------##
if mode == 218:
    # 218 = 1x8 1 block of 8 cells each class tail graph
    classes = cluster_8_cells_tail(M, C)
if mode == 224:
    # 224 = 2x4 2 block of 4 cells each class tail graph
    classes = cluster_4_cells_tail(M, C)
if mode == 241:
    # 241 = 2x4x1 2 stripes of 4 cells each class tail graph
    classes = stripe_4_cells_tail(M, C)
if mode == 242:
    # 242 = 4x2 4 stripes of 2 cells each class tail graph
    classes = stripe_2_cells_tail(M, C)
## --------------------------------------##
if mode == 341:
    # 363 = 4x63 4 blocks 63 cells each, a center class w 4 cells
    classes = cluster_63_4_cells(M,C, [119,120,135,136])
if mode == 342:
    # 363 = 4x63 4 blocks 63 cells each, class w 4 cells each corner
    classes = cluster_63_4_cells(M,C, [0,15,240,255])
if mode == 343:
    # 363 = 4x63 4 blocks 63 cells each, class w 4 cells at corner block
    classes = cluster_63_4_cells(M,C, [14,15,30,31])
if mode == 344:
    # 363 = 4x63 4 blocks 63 cells each, class w 4 cells spread quadrant
    classes = cluster_63_4_cells(M,C, [42,61,90,109])
## --------------------------------------##
if mode == 441:
    ## 1x2 1 cluster of 2 cells each class, 32 classes
    classes = cluster_2_cells_scaled(M, C)
if mode == 442:
    ## 1x32 1 block of 32 cells for each class, 32 classes
    classes = cluster_32_cells_scaled(M, C)
## --------------------------------------##
if mode == 541:
    # 541 = 4x60 4 blocks 60 cells each, a class w 16-cells cluster at center
    classes = cluster_60_16_cells(M,C, [102,104,134,136])
if mode == 542:
    # 542 = 4x60 4 blocks 60 ells each, a class w 16-cells 4 spread all corners
    classes = cluster_60_16_cells(M,C, [0,14,224,238])
if mode == 543:
    # 543 = 4x60 4 blocks 60 cells each, class w 16 cells cluster in one corner
    classes = cluster_60_16_cells(M,C, [12,14,44,46])
if mode == 544:
    # 544 = 4x60 4 blocks 60 cells each, class w 16 cells blocks spread quadrant
    classes = cluster_60_16_cells(M,C, [25,29,89,93])
## --------------------------------------##



# loop over all nodes given on input
for node in range(start, end):
# for node in tqdm(range(start, end)):
    obs = classes[node]
    avg_series = []
    # repeat "num" iterations for each node
    for n in range(num):
        series = walk(node, obs)
        # check if is first trajectory
        if len(avg_series) == 0:
            avg_series = series
            continue
        len_avg = len(avg_series)
        len_sr = len(series)
        # concatenate series with avg and fill missing values
        if len_avg < len_sr:
            summed = [sum(x) for x in zip_longest(avg_series, series, fillvalue=C*n)]
        elif len_avg > len_sr:
            summed = [sum(x) for x in zip_longest(avg_series, series, fillvalue=C)]
        else:
            summed = [sum(x) for x in zip(avg_series, series)]
        # update avg_series with trajectories summed
        avg_series = summed

    # compute avg of classes seen at time t over all trajectories
    avg_series = [nc/num for nc in avg_series]

    # compute the avg time to find each ratio of classes in [0,1]
    th_values = []
    for t in range(1, 101):
        # th = int(np.ceil((t*C)/100))
        th = (t*C) // 100  # threshold value
        time = 0
        # check only series w/ first elem > th
        for s in avg_series:
            if s >= th:
                time += 1
                break
            else:
                time += 1
        # time for each th per node (series)
        th_values.append(time)
    # dump node id and class coverage profile (time for each fraction)
    print(node, *th_values)
