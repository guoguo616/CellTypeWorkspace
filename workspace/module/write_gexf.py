from pathlib import Path
from xml.etree.ElementTree import Element, SubElement

import networkx as nx
import pandas as pd
import scanpy as sc
import squidpy as sq
import scipy as sp
from anndata import AnnData
from matplotlib import pyplot as plt
from networkx.readwrite.gexf import GEXFWriter
import argparse


def read_visium(visium_path: Path) -> AnnData:
    adata: AnnData = sc.read_visium(visium_path)
    return adata


def read_aknno_cluster_result(adata: AnnData, cluster_file_path: Path, is_overridden: bool) -> None:
    file_name = f"aKNNO_cluster_result{'_overridden' if is_overridden else ''}.csv"
    cluster_result = pd.read_csv(cluster_file_path / file_name, header=0, index_col=0)
    adata.obs['cluster'] = cluster_result['Cluster']
    adata.obs['color'] = cluster_result['Color']


def calc_connectivity(adata: AnnData, filter_by_corr: bool) -> None:
    sq.gr.spatial_neighbors(adata, coord_type='grid', spatial_key='spatial', n_neighs=6, n_rings=1)
    if filter_by_corr:
        sc.pp.neighbors(adata, metric='correlation', key_added='corr')
        adata.obsp['spatial_corr_connectivities'] = adata.obsp['corr_connectivities'].multiply(
            adata.obsp['spatial_connectivities'])

        mask = sp.csr_matrix(adata.obsp['spatial_corr_connectivities'] > 0, dtype=bool)
        adata.obsp['spatial_corr_distances'] = adata.obsp['corr_distances'].multiply(mask)


def create_network(adata: AnnData, cluster_obs_name: str, connectivities_key: str, output_dir: Path, is_overridden:bool) -> nx.Graph:
    nx_graph: nx.Graph = nx.Graph(adata.obsp[connectivities_key])
    cell_coord_dict = {idx: (int(row[0]), int(-row[1])) for idx, row in enumerate(adata.obsm['spatial'])}

    name_dict = {idx: cell_name for idx, cell_name in enumerate(adata.obs.index)}
    nx.set_node_attributes(nx_graph, name_dict, 'label')

    x_dict = {idx: row[0] for idx, row in cell_coord_dict.items()}
    nx.set_node_attributes(nx_graph, x_dict, 'x')

    y_dict = {idx: -row[1] for idx, row in cell_coord_dict.items()}
    nx.set_node_attributes(nx_graph, y_dict, 'y')

    cluster_dict = {idx: cluster_name for idx, cluster_name in enumerate(adata.obs[cluster_obs_name].array)}
    nx.set_node_attributes(nx_graph, cluster_dict, 'cluster')

    color_dict = {idx: color for idx, color in enumerate(adata.obs['color'].array)}
    nx.set_node_attributes(nx_graph, color_dict, 'color')

    nx.draw(nx_graph, node_size=1, pos=cell_coord_dict)
    file_name = f"network_spatial{'_overridden' if is_overridden else ''}.png"
    plt.savefig(output_dir / file_name)
    plt.close()

    return nx_graph


def write_gexf(graph_obj: nx.Graph, adata: AnnData, cluster_obs_name: str, file_path: Path) -> None:
    cluster_count_dict = adata.obs[cluster_obs_name].value_counts().to_dict()
    max_cluster_count = max(cluster_count_dict.values())
    cluster_name_list = list(cluster_count_dict.keys())

    writer = GEXFWriter(encoding='utf-8', prettyprint=True, version='1.2draft')
    writer.add_graph(graph_obj)

    info_element: Element = Element('info')
    num_cluster: int = len(cluster_name_list)
    info_element.set('num_cluster', str(num_cluster))
    for cluster_name in cluster_name_list:
        cluster_element = SubElement(info_element, 'cluster')
        cluster_element.set('name', str(cluster_name))
        cluster_element.set('freq',
                            '%.6f' % (cluster_count_dict.get(cluster_name) / max_cluster_count * 100))
        cluster_element.set('color',
                            str(adata.obs[adata.obs[cluster_obs_name] == cluster_name]['color'].array[0]))
    writer.xml.append(info_element)

    writer.write(file_path)


def main(visium_path, output_path, is_overridden):
    print('Reading visium')
    adata = read_visium(visium_path)
    print('Reading cluster result')
    read_aknno_cluster_result(adata, output_path, is_overridden)
    print('Calculating connectivity')
    calc_connectivity(adata, False)
    print('Saving adata')
    adata.write(output_path / 'adata.h5ad')
    print('Creating network')
    graph_obj = create_network(adata, 'cluster', 'spatial_connectivities', output_path, is_overridden)
    print('Writing gexf')
    file_name = f"network{'_overridden' if is_overridden else ''}.gexf"
    write_gexf(graph_obj, adata, 'cluster', output_path / file_name)
    print('Done')


if __name__ == '__main__':
    # read args here
    parser = argparse.ArgumentParser()
    parser.add_argument('--visium_path', help='Visium path', required=True)
    parser.add_argument('--output_path', help='Output path', required=True)
    parser.add_argument('--is_overridden', help='Overridden', required=False)
    args = parser.parse_args()

    visium_path = Path(args.visium_path)
    assert visium_path.exists()
    output_path = Path(args.output_path)
    assert output_path.exists()
    is_overridden = args.is_overridden and args.is_overridden.lower() == 'true'
    main(visium_path, output_path, is_overridden)
