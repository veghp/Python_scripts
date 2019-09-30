from scipy import sparse, io
import pandas as pd
import numpy as np
# See https://github.com/veghp/R_scripts/blob/master/export_cellphonedb.R
# counts.txt : writeMM(seuratobject@raw.data[, cells], file = "counts.txt")
# colnames.txt : write(colnames(counts), file = "colnames.txt")
# rownames.txt : write(ens.rownames, file = "rownames.txt")
counts = io.mmread('counts.txt')
c_dense = counts.toarray()

# Count normalization
c_dense = c_dense / c_dense.sum(axis=0, keepdims=True)
c_dense = c_dense * 10000

var_names = np.genfromtxt('rownames.txt', dtype=str)
col_names = np.genfromtxt('colnames.txt', dtype=str)
df = pd.DataFrame(c_dense, columns=col_names, index=var_names)

df.to_csv('export_counts.txt', sep='\t', header=True, index=True, index_label='Gene')
