__author__ = 'Debha'

import numpy as np
import random
from operator import itemgetter
from itertools import combinations
import scipy.spatial.distance as dist

def make_cluster_mat(data, gene_input):
    cluster_mat = np.zeros((len(data), len(gene_input)))
    cluster_mat_nonbinary = np.zeros((len(data), len(gene_input)))
    subject_list = []
    gene_list = []
    gene_length_map = {}
    
    sub_num = 0
    for sub in data:
        subject_list.append(sub)
        for m in data[sub]:
            if m.gene in gene_input:
                cluster_mat[sub_num, gene_input.index(m.gene)] = 1
                cluster_mat_nonbinary[sub_num, gene_input.index(m.gene)] += 1

                # map of gene lengths
                gene_length_map[m.gene] = m.length

                if m.gene not in gene_list: gene_list.append(m.gene)
        sub_num += 1

    f = open("cluster_mat.csv", mode="w")
    g = open("cluster_mat_nonbinary.csv", mode="w")
    for row in range(np.shape(cluster_mat)[0]):
        for col in range(np.shape(cluster_mat)[1]):
            f.write(str(cluster_mat[row, col]) + ", ")
            g.write(str(cluster_mat_nonbinary[row, col]) + ", ")
        f.write("\n")
        g.write("\n")
    f.close()
    g.close()

    np.save("cluster_mat", cluster_mat)
    np.save("cluster_mat_nonbinary", cluster_mat_nonbinary)
    np.save("subjects", subject_list); np.save("genes", gene_list)

    return [cluster_mat, cluster_mat_nonbinary, subject_list, gene_list, gene_length_map]

def make_sub_mut_map(cluster_mat, subjects, genes):
    sub_mut_map = {}
    for row in range(np.shape(cluster_mat)[0]):
        for col in range(np.shape(cluster_mat)[1]):
            if cluster_mat[row, col] == 1:
                try:
                    sub_mut_map[subjects[row]].append(genes[col])
                except KeyError:
                    sub_mut_map[subjects[row]] = [genes[col]]
    return sub_mut_map

def make_GO_cluster_mat(sub_mut_map, GO_map, subjects):
    GO_terms = GO_map.keys()
    np.save("GO_terms", list(GO_terms))
    GO_cluster_mat = np.zeros([len(subjects), len(GO_terms)])

    f = open("GO_cluster_mat.csv", mode="w")
    row = 0
    for s in subjects:
        sub_muts = sub_mut_map[s]
        col = 0
        for g in GO_terms:
            GO_term_muts = GO_map[g]
            overlap = list(set(sub_muts) & set(GO_term_muts))
            frac_overlap = len(overlap)/len(GO_term_muts)

            GO_cluster_mat[row, col] = frac_overlap
            f.write(str(frac_overlap) + ", ")
            col += 1
        f.write("\n")
        row += 1

    f.close()
    np.save("GO_cluster_mat", GO_cluster_mat)
    return GO_cluster_mat

def clean(GO_cluster_mat, GO_terms):
    # Remove singleton ontologies
    removed_GO = []
    indices = []

    for i in range(np.shape(GO_cluster_mat)[1]):
        if np.count_nonzero(GO_cluster_mat[:, i]) <= 1:
            removed_GO.append(GO_terms[i])
            indices.append(i)
    c_GO_cluster_mat = np.delete(GO_cluster_mat, indices, axis=1)
    c_GO_terms = np.delete(GO_terms, indices)

    # scale each subject to max ontology value
    for i in range(np.shape(c_GO_cluster_mat)[0]):
        c_GO_cluster_mat[i, :] = c_GO_cluster_mat[i, :] / float(max(c_GO_cluster_mat[i, :]))



    np.save("../intermediates/c_GO_terms", c_GO_terms)
    np.save("../intermediates/c_GO_cluster_mat", c_GO_cluster_mat)
    print(len(c_GO_terms), " GO terms used")
    return [c_GO_cluster_mat, c_GO_terms]

def cluster_points(X,mu,subjects):
    clusters = {}
    subject_clusters = {}
   #for x in X: #Assign each point to the centroid it is closest to
    for i in range(len(X)):
        norms = []
        x = np.array(X[i])
        for m in enumerate(mu):
            cent = np.array(m[1])
           #norms.append(np.linalg.norm(x-cent))
            norms.append(dist.cosine(x, cent))
        bestmukey = min(enumerate(norms), key=itemgetter(1))[0]

        try:#Store vectors
            clusters[bestmukey].append(x)
        except KeyError:
            clusters[bestmukey] = [x]

        try:#Store subject identifiers. AS WRITTEN, ONLY STORES
            #INDEX, USE subjects[i] instead of i to store real names
           subject_clusters[bestmukey].append(i)
        except KeyError:
           subject_clusters[bestmukey] = [i]

    return [clusters, subject_clusters]

def reevaluate_centers(mu,clusters):
    newmu = []
    keys = sorted(clusters.keys())
    for k in keys:
        newmu.append(np.mean(clusters[k], axis=0))
    return newmu

def has_converged(mu,oldmu):
    return set([tuple(a) for a in mu]) == set([tuple(a) for a in oldmu])

def find_distortion(mu, clusters):
    distortion = 0
    for c in clusters:
        cent = mu[c]
        for x in clusters[c]:
            distortion += np.linalg.norm(x-cent)
    return distortion

def find_centers(X,K,subjects):
    oldmu = random.sample(X,K)
    mu = random.sample(X,K)
    while not has_converged(mu,oldmu):
        oldmu = mu
        [clusters, subject_clusters] = cluster_points(X,mu,subjects)
        mu = reevaluate_centers(oldmu,clusters)
    distortion = find_distortion(mu,clusters)
    return [mu,clusters,distortion,subject_clusters]

def cluster_main(X,K,n,subjects):
    d = 0
    best_mu = 0
    best_clusters = {}
    best_subject_clusters = {}
    for i in range(n):
        [mu,clusters,distortion,subject_clusters] = find_centers(X,K,subjects)
        if d == 0: d = distortion
        elif distortion < d:
            d = distortion
            best_mu = mu
            best_clusters = clusters
            best_subject_clusters = subject_clusters
    return [best_mu,best_clusters,distortion,best_subject_clusters]

def pairwiseDist(cluster_mat):
    # Create list of pairwise distances
    pairD = []
    for c in combinations(cluster_mat,2):
        # d = np.linalg.norm(np.array(c[0])-np.array(c[1]))
        d = dist.cosine(np.array(c[0]), np.array(c[1]))
        pairD.append(d)

    f = open("../outputs/pairwise_distance.txt", mode="w")
    for i in range(len(pairD)):
        f.write(str(pairD[i]) + "\n")
    f.close()

def sim_matrix(cluster_mat):
    # Create similarity matrix
    n = len(cluster_mat)
    f = open("../outputs/similarity_matrix.csv", mode="w")
    for i in range(n):
        for j in range(n):
            if i == j: d = 0
            else: d = dist.cosine(np.array(cluster_mat[i,:]), np.array(cluster_mat[j,:]))
            f.write(str(d))
            if j+1 != n: f.write(",")
        f.write("\n")
    f.close()   

def clustered_GOs(clusters, c_GO_terms, mu):
    f = open("../outputs/centroids.txt", mode="w")

    for c in range(len(mu)):
        top_GO = mu[c].argsort()[::-1] #indices of top ontologies in descending order
        for i in range(10):
            print(c, " ", c_GO_terms[top_GO[i]])
        print()
        #write centroid
        for i in range(len(mu[c])):
            f.write(str(mu[c][i]) + "\t")
        f.write("\n")
