import numpy as np
import pandas as pd
import scipy
import scipy.io
from scipy.sparse import csr_matrix, csc_matrix
import anndata
import scanpy
import scanpy.api as sc

fileprefix = 'export_'

adata = sc.read_mtx(fileprefix + 'data.txt')
adata.var_names = np.genfromtxt(fileprefix + 'rownames.txt', dtype=str)
metadata = pd.read_table(fileprefix + 'metadata.txt', index_col=0, header=0)
adata.obs = metadata

adata.write('adata_' + fileprefix + '.h5ad')
