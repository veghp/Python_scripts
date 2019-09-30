adatafile = '.h5ad'
n_cell = 5000 # how many cells to show in 3D plot
annotation_col = ''
###############################################################################
# 3D UMAP with scanpy
import numpy as np
import pandas as pd
import scanpy as sc # v1.4.3
sc.settings.set_figure_params(dpi=150)
adata = sc.read_h5ad(adatafile)
adata3d = sc.tl.umap(adata, n_components=3, copy=True)
adata3d.obsm['X_umap']
# tsne can't be made 3d in scanpy
data3d = adata3d.obsm['X_umap']
np.save('data3d.npy', data3d)
adata.obsm['X_umap3d'] = data3d
del adata3d
#############
# 2D plotting
sc.pl.umap(adata)
# built-in plotters return error:
#~ sc.pl.umap(adata, projection='3d')
#~ sc.pl.scatter(adata, color='tissue', basis='umap3d', projection='3d')
#############
# 3D plotting
coord3d = np.split(data3d, indices_or_sections=3, axis=1)
X, Y, Z = coord3d

from mpl_toolkits.mplot3d import axes3d
import matplotlib.pyplot as plt
%matplotlib qt

div = round(X.shape[0] / n_cell)
A = X[::div] # for a fast plot
B = Y[::div]
C = Z[::div]

# Colours
cmap = plt.get_cmap('Dark2') # tab10
colour_names = adata.obs[annotation_col].tolist()
np.linspace(0, 1, len(colour_names))
colours = cmap(np.linspace(0, 1, len(colour_names)))
colours = colours[::div]

# Plot
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(A, B, C, c=colours)

# Rotate the axes and update
for angle in range(0, 360):
    ax.view_init(30, angle)
    plt.draw()
    plt.pause(.001)

plt.close()
