import networkx as nx

def create_gene_network(edges, disease_genes=None):
    """Create and analyze gene interaction network"""
    G = nx.Graph()
    
    # Add weighted edges
    for gene1, gene2, weight in edges:
        G.add_edge(gene1, gene2, weight=weight)
    
    # Calculate network metrics
    info = {
        'num_nodes': G.number_of_nodes(),
        'num_edges': G.number_of_edges(),
        'avg_degree': sum(dict(G.degree()).values()) / G.number_of_nodes(),
        'density': nx.density(G)
    }
    
    # Add disease gene info if provided
    if disease_genes:
        disease_connections = []
        for gene in disease_genes:
            if gene in G:
                neighbors = list(G.neighbors(gene))
                disease_connections.append({
                    'gene': gene,
                    'degree': len(neighbors),
                    'neighbors': neighbors
                })
        info['disease_genes'] = disease_connections
    
    return G, info

def find_shortest_path(G, start_gene, end_gene):
    """Find shortest path between two genes"""
    if start_gene not in G or end_gene not in G:
        return {'error': 'One or both genes not found in network'}
    
    try:
        path = nx.shortest_path(G, start_gene, end_gene, weight='weight')
        path_edges = list(zip(path[:-1], path[1:]))
        edge_weights = [G[u][v]['weight'] for u, v in path_edges]
        
        return {
            'path': path,
            'length': len(path) - 1,
            'edge_weights': edge_weights
        }
    except nx.NetworkXNoPath:
        return {'error': 'No path found between genes'}