import numpy as np
import scanpy.api as sc

pca = np.loadtxt("/nfs/research1/marioni/jonny/embryos/scripts/endoderm/corrected.tab")

obj = sc.AnnData(np.random.random((pca.shape[0], 100)))

obj.obsm["X_pca"] = pca



sc.pp.neighbors(obj, n_neighbors = 20, use_rep = "X_pca")




sc.tl.umap(obj, min_dist = 0.5)
out = obj.obsm["X_umap"]
np.savetxt("/nfs/research1/marioni/jonny/embryos/scripts/endoderm/umap.tab", out, delimiter="\t")

