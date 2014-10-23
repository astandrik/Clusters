__author__ = 'anton'


import numpy as np
import math
import Sorts
import sys
import scipy
import scipy.cluster.hierarchy as hier
import scipy.spatial.distance as dist
import matplotlib.pyplot as plt
from hcluster import pdist, linkage, dendrogram
import getopt



def squared_criterion():
    criterion = 0
    for c in Cluster.clusters:
        for x in c.elements:
            criterion += ((x - c.get_center()) ** 2).sum()
    return criterion

def avg_dist(element, clust):
    distance = 0
    for x in clust.elements:
        distance += euqlid_distance(element, x)
    return (distance / len(clust.elements))


def silhouette_criterion_el(element, in_claster):
    clusters_to_check = Cluster.clusters[:]
    a_i = avg_dist(element, in_claster);
    clusters_to_check.remove(in_claster)
    b_list = []
    for x in clusters_to_check:
        b_list.append(avg_dist(element,x))
    b_i = min(b_list)
    return (b_i - a_i) / (max(a_i, b_i))

def silhouette_criterion():
    all_dots_num = 0
    all_silhouette = 0
    for c in Cluster.clusters:
        for el in c.elements:
            all_dots_num += 1
            all_silhouette += silhouette_criterion_el(el, c)
    return (all_silhouette/all_dots_num)

def euqlid_distance(element1, element2):
    return np.sqrt((((element2 - element1)**2).sum()))

def count_criteriums(element1, element2, J):
    counted = J(element1, element2)
    Cluster.counted_criteriums[(element1, element2)] = counted;


def find_best_merge(J):
    min_criterium = sys.float_info.max
    for i in range(len(Cluster.clusters)):
        for j in range(i + 1, len(Cluster.clusters)):
            x = Cluster.clusters[i]
            y = Cluster.clusters[j]

            if (x,y) in Cluster.counted_criteriums:
                criterium = Cluster.counted_criteriums[(x,y)]
            else:
                count_criteriums(x,y, J)
                criterium = Cluster.counted_criteriums[(x,y)]

            if criterium < min_criterium:
                min_set = (x, y)
                min_criterium = criterium
    else:
        return min_set


class Cluster:
    clusters = []
    counted_criteriums = {}
    squared_criterion_values = []
    silhouette_criterion_values = []
    def __init__(self, element):
        self.elements = np.array([element])
        self.num_elements = 1

    def get_center(self):
        sum = self.elements[0] - self.elements[0]
        for x in self.elements:
            sum += x
        sum /= self.num_elements
        return sum

    def add_element(self, element):
        self.elements = np.array(self.elements.tolist() + element)
        self.num_elements += 1

    def merge(self, c2):
        self.elements = np.array(self.elements.tolist() + c2.elements.tolist())
        self.num_elements += len(c2.elements)
        Cluster.clusters.remove(c2)



def d_max(c1, c2):
    max_dist = 0
    distances = []
    for x in c1.elements:
        for y in c2.elements:
            distances.append(euqlid_distance(x,y))
    return max(distances)


def d_e(c1,c2):
    n1 = len(c1.elements)
    n2 = len(c2.elements)
    return math.sqrt(n1*n1 / (n1 + n2)) * euqlid_distance(c1.get_center(), c2.get_center())

def process_and_merge_clusters():
    for x in Cluster.clusters:
        pass

def swo(K, J, crit_func):
    while len(Cluster.clusters) > K:
        if len(Cluster.clusters) < 30:
            if(crit_func == silhouette_criterion):
                Cluster.silhouette_criterion_values.append(silhouette_criterion())
            else:
                Cluster.squared_criterion_values.append(squared_criterion())
        print "Clusters now: ", len(Cluster.clusters)
        C_a,C_b = find_best_merge(J)
        C_a.merge(C_b)



def main(argv):

    print argv
    if(len(argv) > 0):
        params = argv[::2]
        param_values = argv[1::2]

    crit_func = silhouette_criterion
    merge_func = d_e
    for i in range(0,len(argv),2):
        if params[i] == "--criterium":
            if param_values[i+1] == "silhoette":
                crit_func = silhouette_criterion
            elif param_values[i+1] == "squared":
                crit_func = squared_criterion
            else:
                crit_func = silhouette_criterion
        elif params[i] == "--merge":
            if param_values[i+1] == "de":
                merge_func = d_e
            elif param_values[i+1] == "dmax":
                merge_func = d_max
            else:
                merge_func = d_e

    Cluster.clusters = []
    Cluster.squared_criterion_values = []
    Cluster.silhouette_criterion_values = []
    Cluster.counted_criteriums.clear()

    my_data = np.genfromtxt('/home/anton/workspace/Proj/data.csv', delimiter=',', dtype=float)
    #Make only clusterization params in array
    data_list = my_data[1:].tolist()
    maximum = 0
    for i in range(len(data_list)):
        data_list[i] = data_list[i][2:]
        if max(data_list[i]) > maximum:
            maximum = max(data_list[i])

    #normalize all lists:
    data_list = np.array(data_list)

    for i in range(len(data_list[0])):
        data_list[:,i] /= data_list[:,i].max()

    data_list = data_list[:100]
    #Make each element = 1 cluster


    Y = dist.pdist(data_list)
    Z = linkage(Y)
    print Z
    plt.subplot(121)
    dendrogram(Z)

    for x in data_list:
        Cluster.clusters.append(Cluster(x))

    print(len(Cluster.clusters))
    swo(1, merge_func, crit_func)
    squared_criterion_values = Cluster.squared_criterion_values[::-1]
    silhouette_criterion_values = Cluster.silhouette_criterion_values[::-1]
    plt.subplot(122)
    if(crit_func == silhouette_criterion) :
        plt.plot(range(len(silhouette_criterion_values)), silhouette_criterion_values, 'ro')
        plt.axis([0,10,min(silhouette_criterion_values),max(silhouette_criterion_values)])
    else:
        plt.plot(range(len(squared_criterion_values)), squared_criterion_values, 'ro')
        plt.axis([0,10,min(squared_criterion_values),max(squared_criterion_values)])

    plt.show()



    for x in Cluster.clusters:
        print len(x.elements), '\n'


if __name__ == "__main__":
    if(len(sys.argv) > 1):
        main(sys.argv[1:])
    else:
        main([])