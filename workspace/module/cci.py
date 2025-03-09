import argparse
import json
from pathlib import Path

import commot as ct
import matplotlib.pyplot as plt
import networkx as nx
import numpy as np
import scanpy as sc
from anndata import AnnData
from networkx.readwrite.gexf import GEXFWriter
from scipy.spatial import distance_matrix
from scipy.stats import norm
from sklearn.neighbors import NearestNeighbors


def pre_process(adata: AnnData):
    adata.var_names_make_unique()
    adata.raw = adata
    sc.pp.normalize_total(adata, inplace=True)
    sc.pp.log1p(adata)


def extract_cell_comm_grid_data(adata: AnnData, database_name: str, summary: str) -> dict:
    vf_name = database_name + '-total-total'
    pos_idx: np.ndarray = np.array([0, 1], int)
    normalize_v_quantile: float = 0.95
    grid_density = 1
    grid_thresh = 1.0
    grid_scale = 1.0
    X = adata.obsm["spatial"][:, pos_idx]
    V = adata.obsm[f'commot_{summary}_vf-{vf_name}'][:, pos_idx]
    V = V / np.quantile(np.linalg.norm(V, axis=1), normalize_v_quantile)
    # Get a rectangular grid
    xl, xr = np.min(X[:, 0]), np.max(X[:, 0])
    epsilon = 0.02 * (xr - xl)
    xl -= epsilon
    xr += epsilon
    yl, yr = np.min(X[:, 1]), np.max(X[:, 1])
    epsilon = 0.02 * (yr - yl)
    yl -= epsilon
    yr += epsilon
    ngrid_x = int(50 * grid_density)
    gridsize = (xr - xl) / float(ngrid_x)
    ngrid_y = int((yr - yl) / gridsize)
    meshgrid = np.meshgrid(np.linspace(xl, xr, ngrid_x), np.linspace(yl, yr, ngrid_y))
    grid_pts = np.concatenate((meshgrid[0].reshape(-1, 1), meshgrid[1].reshape(-1, 1)), axis=1)
    grid_knn = int(X.shape[0] / 50)
    nn_mdl = NearestNeighbors()
    nn_mdl.fit(X)
    dis, nbs = nn_mdl.kneighbors(grid_pts, n_neighbors=grid_knn)
    w = norm.pdf(x=dis, scale=gridsize * grid_scale)
    w_sum = w.sum(axis=1)
    V_grid = (V[nbs] * w[:, :, None]).sum(axis=1)
    V_grid /= np.maximum(1, w_sum)[:, None]
    grid_thresh *= np.percentile(w_sum, 99) / 100
    grid_pts, V_grid = grid_pts[w_sum > grid_thresh], V_grid[w_sum > grid_thresh]
    spatial_mapping = adata.uns.get("spatial", {})
    library_id = list(spatial_mapping.keys())[0]
    spatial_data = spatial_mapping[library_id]
    sf = spatial_data['scalefactors']['tissue_hires_scalef']
    # ax.quiver(grid_pts[:,0]*sf, grid_pts[:,1]*sf, V_grid[:,0]*sf, V_grid[:,1]*sf, scale=scale, scale_units='x',
    # width=grid_width, color=arrow_color)
    return {
        'X': list(grid_pts[:, 0] * sf),
        'Y': list(grid_pts[:, 1] * sf),
        'U': list(V_grid[:, 0] * sf),
        'V': list(V_grid[:, 1] * sf),
    }


def cluster_center(adata, clustering, method="representative_point") -> None:
    X = adata.obsm['spatial']
    cluster_pos = {}
    clusternames = list(adata.obs[clustering].unique())
    clusternames.sort()
    clusterid = np.array(adata.obs[clustering], str)
    for name in clusternames:
        tmp_idx = np.where(clusterid == str(name))[0]
        tmp_X = X[tmp_idx]
        if method == "geometric_mean":
            X_mean = np.mean(tmp_X, axis=0)
        else:  # representative_point
            tmp_D = distance_matrix(tmp_X, tmp_X)
            tmp_D = tmp_D ** 2
            X_mean = tmp_X[np.argmin(tmp_D.sum(axis=1)), :]
        cluster_pos[str(name)] = X_mean
    adata.uns['cluster_pos-' + clustering] = cluster_pos


def linear_clamp_value(x, lower_bound, upper_bound, out_min, out_max):
    if x <= lower_bound:
        y = out_min
    elif x >= upper_bound:
        y = out_max
    else:
        y = out_min + (x - lower_bound) / (upper_bound - lower_bound) * (out_max - out_min)
    return y


def plot_cluster_communication_network(adata):
    uns_names = ['commot_cluster-cluster-cellchat-total-total']
    quantile_cutoff: float = 0.99
    p_value_cutoff: float = 0.05
    self_communication_off: bool = False
    clustering: str = 'cluster'
    nx_pos_idx: np.ndarray = np.array([0, 1], int)
    edge_width_lb_quantile = 0.05
    edge_width_ub_quantile = 0.95
    node_size = 0.2
    edge_width_min = 1
    edge_width_max = 4
    X_tmp = adata.uns[uns_names[0]]['communication_matrix'].copy()
    labels = list(X_tmp.columns.values)
    X = np.zeros_like(X_tmp.values, float)
    for i in range(len(uns_names)):
        X_tmp = adata.uns[uns_names[i]]['communication_matrix'].values.copy()
        p_values_tmp = adata.uns[uns_names[i]]['communication_pvalue'].values.copy()
        if not quantile_cutoff is None:
            cutoff = np.quantile(X_tmp.reshape(-1), quantile_cutoff)
        else:
            cutoff = np.inf
        tmp_mask = (X_tmp < cutoff) * (p_values_tmp > p_value_cutoff)
        X_tmp[tmp_mask] = 0
        X = X + X_tmp
    X = X / len(uns_names)
    if self_communication_off:
        for i in range(X.shape[0]):
            X[i, i] = 0
    node_pos = [adata.uns["cluster_pos-" + clustering][labels[i]] for i in range(len(labels))]
    node_pos = np.array(node_pos)
    node_pos = node_pos[:, nx_pos_idx]
    lx = np.max(node_pos[:, 0]) - np.min(node_pos[:, 0])
    ly = np.max(node_pos[:, 1]) - np.min(node_pos[:, 1])
    pos_scale = max(lx, ly)
    node_pos = node_pos / pos_scale * 8.0
    background_pos = adata.obsm["spatial"][:, nx_pos_idx]
    background_pos = background_pos / pos_scale * 8.0
    G = nx.MultiDiGraph()
    edge_width_lb = np.quantile(X.reshape(-1), edge_width_lb_quantile)
    edge_width_ub = np.quantile(X.reshape(-1), edge_width_ub_quantile)
    cluster_colors = adata.obs[['cluster', 'color']].drop_duplicates().set_index('cluster')['color'].to_dict()
    node_cmap = [cluster_colors[int(label)] for label in labels]
    for i in range(len(labels)):
        G.add_node(labels[i], shape="point", fillcolor=node_cmap[i], color=node_cmap[i])
        G.nodes[labels[i]]["pos"] = "%f,%f!" % (node_pos[i, 0], node_pos[i, 1])
        G.nodes[labels[i]]["width"] = str(node_size)
    # Draw the edges
    for i in range(X.shape[0]):
        for j in range(X.shape[1]):
            if X[i, j] > 0:
                G.add_edge(labels[i], labels[j], splines="curved")
                G[labels[i]][labels[j]][0]["penwidth"] = str(
                    linear_clamp_value(X[i, j], edge_width_lb, edge_width_ub, edge_width_min, edge_width_max))
    # Draw the network
    # graphviz needs to be installed, I don't want to
    # A = to_agraph(G)
    # A.layout()
    # return A
    return G


def main(adata_path, output_path, species):
    adata_spatial_comm_path = output_path / "adata_spatial_comm.h5ad"
    if adata_spatial_comm_path.exists():
        adata: AnnData = sc.read(adata_spatial_comm_path)
    else:
        adata: AnnData = sc.read_visium(adata_path)
        pre_process(adata)

        # species = 'mouse'
        df_cellchat = ct.pp.ligand_receptor_database(species=species, signaling_type='Secreted Signaling',
                                                     database='CellChat')
        # Filter the LR pairs to keep only the pairs with both ligand and receptor expressed in at least 5% of the
        # spots.
        df_cellchat_filtered = ct.pp.filter_lr_database(df_cellchat, adata, min_cell_pct=0.05)
        print(df_cellchat_filtered.head())

        ct.tl.spatial_communication(adata, database_name='cellchat', df_ligrec=df_cellchat_filtered,
                                    dis_thr=500, heteromeric=True, pathway_sum=True)

        # Determine the spatial direction of a signaling pathway, for example, the PSAP pathway.
        # The interpolated signaling directions for where the signals are sent by the spots
        # and where the signals received by the spots are from are stored in
        # adata.obsm['commot_sender_vf-cellchat-PSAP'] and adata.obsm['commot_receiver_vf-cellchat-PSAP']
        ct.tl.communication_direction(adata, database_name='cellchat', k=5)
        # Summarize the signaling to the cluster. The results are stored in adata_dis500.uns[
        # 'commot_cluster-leiden-cellchat-PSAP']
        ct.tl.cluster_communication(adata, database_name='cellchat', clustering='cluster', n_permutations=100)
        # Calculate the cluster center positions. The results are stored in adata.uns['cluster_pos-cluster']
        cluster_center(adata, clustering='cluster')
        adata.write(output_path / "adata_spatial_comm.h5ad")

    # Sender
    ct.pl.plot_cell_communication(
        adata, database_name='cellchat', plot_method='cell', background_legend=True,
        scale=0.003, ndsize=8, grid_density=0.4, summary='sender', background='image', clustering='cluster',
        cmap='Alphabet', normalize_v=True, normalize_v_quantile=0.995, arrow_color=adata.obs['color'])
    plt.gca().invert_yaxis()
    plt.savefig(output_path / 'cell_communication_sender_cell_plot.svg', format='svg')
    plt.close()
    ct.pl.plot_cell_communication(
        adata, database_name='cellchat', plot_method='stream', background_legend=True,
        scale=0.003, ndsize=8, grid_density=0.4, summary='sender', background='image', clustering='cluster',
        cmap='Alphabet', normalize_v=True, normalize_v_quantile=0.995)
    plt.gca().invert_yaxis()
    plt.savefig(output_path / 'cell_communication_sender_stream_plot.svg', format='svg')
    plt.close()

    # Receiver
    ct.pl.plot_cell_communication(
        adata, database_name='cellchat', plot_method='cell', background_legend=True,
        scale=0.003, ndsize=8, grid_density=0.4, summary='receiver', background='image', clustering='cluster',
        cmap='Alphabet', normalize_v=True, normalize_v_quantile=0.995, arrow_color=adata.obs['color'])
    plt.gca().invert_yaxis()
    plt.savefig(output_path / 'cell_communication_receiver_cell_plot.svg', format='svg')
    plt.close()
    ct.pl.plot_cell_communication(
        adata, database_name='cellchat', plot_method='stream', background_legend=True,
        scale=0.003, ndsize=8, grid_density=0.4, summary='receiver', background='image', clustering='cluster',
        cmap='Alphabet', normalize_v=True, normalize_v_quantile=0.995)
    plt.gca().invert_yaxis()
    plt.savefig(output_path / 'cell_communication_receiver_stream_plot.svg', format='svg')
    plt.close()

    cell_communication_plot_sender_json = extract_cell_comm_grid_data(adata, 'cellchat', 'sender')
    with open(output_path / 'cell_communication_plot_sender.json', 'w') as f:
        json.dump(cell_communication_plot_sender_json, f)
    cell_communication_plot_receiver_json = extract_cell_comm_grid_data(adata, 'cellchat', 'receiver')
    with open(output_path / 'cell_communication_plot_receiver.json', 'w') as f:
        json.dump(cell_communication_plot_receiver_json, f)

    # Plot cluster communication on the spatial grid.
    networkx_obj = plot_cluster_communication_network(adata)
    nx.draw(networkx_obj, node_size=1)
    plt.savefig(output_path / 'cci_cluster_spatial.png')
    plt.close()
    writer = GEXFWriter(encoding='utf-8', prettyprint=True, version='1.2draft')
    writer.add_graph(networkx_obj)
    writer.write(output_path / 'cci_cluster_spatial.gexf')


if __name__ == '__main__':
    # read args here
    parser = argparse.ArgumentParser()
    parser.add_argument('--adata_path', help='Visium path', required=True)
    parser.add_argument('--output_path', help='Output path', required=True)
    parser.add_argument('--species', help='Output path', required=True)
    args = parser.parse_args()

    adata_path = Path(args.adata_path)
    assert adata_path.exists()
    output_path = Path(args.output_path)
    assert output_path.exists()
    species = args.species
    main(adata_path, output_path, species)
