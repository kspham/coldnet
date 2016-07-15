__author__ = 'Debha'

from read_data import read_data, read_map, gene_muts, make_GO_map
from cluster import make_cluster_mat, make_sub_mut_map, make_GO_cluster_mat, cluster_main, clean, pairwiseDist, clustered_GOs
import numpy as np

data_file = "../inputs/somatic.csv"
map_file = "../inputs/mart_export.csv"
GO_file = "../inputs/c5.all.v5.0.symbols.gmt"

# data = read_data(data_file)
# gene_map = read_map(map_file)
# gene_mut = gene_muts(data, gene_map)
# GO_map = make_GO_map(GO_file)

# [cluster_mat, cluster_mat_nonbinary, subjects, genes, gene_length_map] = make_cluster_mat(data, list(gene_mut.keys()))

cluster_mat = np.load("../intermediates/cluster_mat.npy")
subjects = np.load("../intermediates/subjects.npy"); genes = np.load("../intermediates/genes.npy")
GO_terms = np.load("../intermediates/GO_terms.npy")

# sub_mut_map = make_sub_mut_map(cluster_mat, subjects, genes)
# GO_cluster_mat = make_GO_cluster_mat(sub_mut_map, GO_map, subjects)

GO_cluster_mat = np.load("../intermediates/GO_cluster_mat.npy")
[c_GO_cluster_mat, c_GO_terms] = clean(GO_cluster_mat, GO_terms)

k = 3
iterations = 20
[mu, clusters, distortion] = cluster_main(c_GO_cluster_mat.tolist(), k, iterations)

for key in clusters:
    print(key, len(clusters[key]))

clustered_GOs(clusters, c_GO_terms, mu)
pairwiseDist(c_GO_cluster_mat)


