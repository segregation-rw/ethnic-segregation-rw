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
<prop> count of each group per area, at MxC matrix format
<num> number of walk repetions for computing the average for each node
<idx> range of nodes or single id to compute, e.g.: 1, 0-10

> output:
    node
    list: [avg trajectory (n classes seens) at each fraction c]

-----------------------------------------------------------------------------
"""

import numpy as np
import random
from itertools import zip_longest
import sys


if len(sys.argv) < 5:
    print("Usage: %s <edges> <prop> <num> <idx> \n" % sys.argv[0])
    exit(1)


def compute_neighbours(edges):
    # dictionary with node neighbours
    E = dict()
    for i, j in edges:
        if i in E:
            E[i].append(j)
        else:
            E[i] = [j]
    return E


def fill_empty_classes(i, idxs):
    # find groups not represented in city, update w/ 1
    obs = classes[i]
    for n in idxs:
        # set group count to 1
        obs[n] = 1
    return obs


def walk(i, obs):
    nc = []  # number of classes seen
    while 0 in obs:
        # selects one of i neighbours
        j = random.choice(neigh_dict[i])
        obs = obs + classes[j]
        # get current number of classes seen
        non_z = int(np.count_nonzero(obs) - len(idxs))
        nc.append(non_z)
        i = j
    return nc

# %%
# initialize dict with neighbours
neigh_dict = compute_neighbours(np.loadtxt(sys.argv[1], dtype='int'))

# load population data and set variables
prop = np.loadtxt(sys.argv[2], dtype='int')
classes = prop[:, 1:]

# total of class at city level
classes_sum = classes.sum(axis=0)
# total number of classes in the system
C = (classes_sum >= 1).sum()

# index of missing groups at city level
idxs = np.where(classes_sum == 0)[0]

# process range or single node from param
if "-" in sys.argv[4]:
    start, end = [int(x) for x in sys.argv[4].split('-')]
else:
    start, end = [int(sys.argv[4]), int(sys.argv[4])+1]
num = int(sys.argv[3])


# loop over all nodes given on input
for node in range(start, end):
    # fill empty classes
    if len(idxs) > 0:
        obs = fill_empty_classes(node, idxs)
    else:
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

    # compute avg over all trajectories
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
