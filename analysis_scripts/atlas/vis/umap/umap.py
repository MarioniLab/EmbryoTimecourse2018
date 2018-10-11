import numpy as np
import scanpy.api as sc

#All cells
pca = np.loadtxt("/nfs/research1/marioni/jonny/embryos/scripts/vis/umap/scanpy_input.tab")
obj = sc.AnnData(np.random.random((pca.shape[0], 100)))
obj.obsm["X_pca"] = pca

sc.pp.neighbors(obj, n_neighbors = 20, use_rep = "X_pca", random_state = 42)
sc.tl.umap(obj, min_dist = 0.7, random_state = 42)
out = obj.obsm["X_umap"]
np.savetxt("/nfs/research1/marioni/jonny/embryos/scripts/vis/umap/umap.tab", out, delimiter="\t")

# Tests
# for i in np.arange(0.1, 0.8, 0.1):
# 	sc.tl.umap(obj, min_dist = i, random_state = 42)
# 	out = obj.obsm["X_umap"]
# 	np.savetxt("/nfs/research1/marioni/jonny/embryos/scripts/vis/umap/mindist_" + str(i) + ".tab", out, delimiter="\t")

# Embryonic cells
pca = np.loadtxt("/nfs/research1/marioni/jonny/embryos/scripts/vis/umap/scanpy_input_embryonic.tab")
obj = sc.AnnData(np.random.random((pca.shape[0], 100)))
obj.obsm["X_pca"] = pca

sc.pp.neighbors(obj, n_neighbors = 20, use_rep = "X_pca", random_state = 42)
sc.tl.umap(obj, min_dist = 0.7, random_state = 42)
out = obj.obsm["X_umap"]
np.savetxt("/nfs/research1/marioni/jonny/embryos/scripts/vis/umap/umap_embryonic.tab", out, delimiter="\t")

