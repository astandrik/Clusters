__author__ = 'anton'


import numpy as np
import math
import sys
import scipy.spatial.distance as dist
import matplotlib.pyplot as plt
from hcluster import pdist, linkage, dendrogram



def squared_criterion():
    criterion = 0
    for c in Cluster.clusters:
        for x in c.elements:
            criterion += ((x - c.get_center()) ** 2).sum()
    return criterion




def avg_dist(element, clust):
    distance = 0
    for x in clust.elements:
        distance += squared_euqlid_distance(element, x)
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

def get_distance(element1, element2):
    return Cluster.counted_distances[(tuple(element1), tuple(element2))]

def euqlid_distance(element1, element2):
    return np.sqrt((((element2 - element1)**2).sum()))


def hexic_euqlid_distance(element1, element2):
    return (((element2 - element1)**2).sum())**15

def squared_euqlid_distance(element1, element2):
    return (((element2 - element1)**2).sum())

def sqrt_distance(element1, element2):
    return np.sqrt(np.sqrt((((element2 - element1)**2).sum())))

def quadric_euqlid_distance(element1, element2):
    return (((element2 - element1)**2).sum())**2

def manhattan_distance(element1, element2):
    return (((element2 - element1)).sum())

def chebishev_distance(element1, element2):
    distances = []
    for i in range(len(element1)):
        for j in range(len(element2)):
            distances.append(abs(element1[i] - element2[j]))
    return max(distances)

def count_criteriums(element1, element2, J):
    counted = J(element1, element2)
    return counted


def find_best_merge(J):
    min_criterium = sys.float_info.max
    for i in range(len(Cluster.clusters)):
        for j in range(i + 1, len(Cluster.clusters)):
            x = Cluster.clusters[i]
            y = Cluster.clusters[j]


            criterium = count_criteriums(x,y, J)

            if criterium < min_criterium:
                min_set = (x, y)
                min_criterium = criterium
    else:
        return min_set


def purity():
    pass



class Cluster:
    clusters = []
    counted_distances = {}
    squared_criterion_values = []
    silhouette_criterion_values = []
    cluster_num = 0
    merge_history = np.array([np.array([0,0,0,0])])
    etalon_clasters = {}

    def get_class(self):
        unique = set(self.etalon_map)
        cl = 0
        for u in unique:
            m = 0
            counter = [x for x in self.etalon_map if x == u]
            if len(counter) > m:
                m = len(counter)
                cl = u
        return cl





    def __init__(self, element):
        self.elements = np.array([element])
        self.clust = Cluster.cluster_num
        Cluster.cluster_num += 1

    def etalon_to_current_mapping(self):
        self.etalon_map = []
        for x in self.elements:
            self.etalon_map.append(Cluster.etalon_clasters[tuple(x[2:])])


    def get_center(self):
        sum = self.elements[0] - self.elements[0]
        for x in self.elements:
            sum += x
        sum /= len(self.elements)
        return sum



    def merge(self, c2, J):
        self.elements = np.array(self.elements.tolist() + c2.elements.tolist())
        hist_el = np.array([np.array([self.clust, c2.clust, count_criteriums(self, c2, J), len(self.elements)])])
        Cluster.merge_history = np.array(Cluster.merge_history.tolist() + hist_el.tolist())
        Cluster.clusters.remove(c2)


def d_max(c1, c2):
    distances = []
    for x in c1.elements:
        for y in c2.elements:
            distances.append(get_distance(x,y))
    return max(distances)

def d_avg_pairs(c1, c2):
    num_pairs = 0
    summ = 0
    for x in c1.elements:
        for y in c2.elements:
            summ += get_distance(x,y)
            num_pairs += 1
    return (summ / float(num_pairs))

def d_min(c1, c2):
    distances = []
    for x in c1.elements:
        for y in c2.elements:
            distances.append(get_distance(x,y))
    return min(distances)

def d_centroid(c1,c2):
    return get_distance(c1.get_center(), c2.get_center())

def d_e(c1,c2):
    n1 = float(len(c1.elements))
    n2 = float(len(c2.elements))
    return math.sqrt(n1*n1 / (n1 + n2)) * squared_euqlid_distance(c1.get_center(), c2.get_center())



def swo(K, J, crit_func):
    while len(Cluster.clusters) > K:
        if len(Cluster.clusters) < 30:
            if(crit_func == silhouette_criterion):
                Cluster.silhouette_criterion_values.append(silhouette_criterion())
            else:
                Cluster.squared_criterion_values.append(squared_criterion())
        print "Clusters now: ", len(Cluster.clusters)
        C_a,C_b  = find_best_merge(J)
        C_a.merge(C_b, J)
        print '!!!!'
        for c in Cluster.clusters:
            if(len(c.elements) > 1):
                print len(c.elements)
        print "________"



def main(argv):

    print argv
    if(len(argv) > 0):
        params = argv[::2]
        param_values = argv[1::2]

    crit_func = squared_criterion
    merge_func = d_min
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

    my_data = np.genfromtxt('/home/anton/workspace/Proj/data.csv', delimiter=',', dtype=float)
    #Make only clusterization params in array
    data_list = my_data[1:].tolist()
    maximum = 0




    data_list = data_list[:]
    etalon = data_list[:]

    for i in range(len(data_list)):
        data_list[i] = data_list[i][2:]

    #normalize all lists:
    data_list = np.array(data_list)


    #count all distances
    print "Precounting distances"
    for i in range(len(data_list)):
        for j in range(len(data_list)):
            print ".",
            Cluster.counted_distances[(tuple(data_list[i]), tuple(data_list[j]))] = hexic_euqlid_distance(data_list[i], data_list[j])

    print "Distances Counted"

    for i in range(len(data_list)):
        Cluster.etalon_clasters[tuple(data_list[i][2:])] = etalon[i][1]

    print Cluster.etalon_clasters.values()
    #Make each element = 1 cluster




    for x in data_list:
        Cluster.clusters.append(Cluster(x))

    print(len(Cluster.clusters))
    K_num = 1
    swo(K_num, merge_func, crit_func)
    Y = Cluster.merge_history[1:]
    Z = linkage(Y)
    plt.subplot(121)
    dendrogram(Z, labels=range(len(data_list)))
    squared_criterion_values = Cluster.squared_criterion_values[::-1]
    silhouette_criterion_values = Cluster.silhouette_criterion_values[::-1]
    plt.subplot(122)
    if(crit_func == silhouette_criterion) :
        plt.plot(range(len(silhouette_criterion_values)), silhouette_criterion_values)
        plt.axis([K_num,30,min(silhouette_criterion_values),max(silhouette_criterion_values)])
    else:
        plt.plot(range(len(squared_criterion_values)), squared_criterion_values)
        plt.axis([K_num,30,min(squared_criterion_values),max(squared_criterion_values)])

    plt.show()

    for x in Cluster.clusters:
        x.etalon_to_current_mapping()
        print x.etalon_map


   # for x in Cluster.clusters:
    #    print len(x.elements), '\n'


if __name__ == "__main__":
    if(len(sys.argv) > 1):
        main(sys.argv[1:])
    else:
        main([])