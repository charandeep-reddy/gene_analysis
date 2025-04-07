import networkx as nx

# Disease gene database
DISEASE_DATABASE = {
    "BRCA1": "Breast & Ovarian Cancer",
    "TP53": "Multiple Cancers (including Breast & Ovarian)",
    "EGFR": "Lung Cancer",
    "MDM2": "Sarcoma, Breast Cancer",
    "CDK2": "Melanoma, Leukemia",
    "MTHFR": "Cardiovascular disease, Neural tube defects",
    # Add more as needed
}

def create_gene_network(edges, disease_genes=None):
    """
    Create a gene interaction network and analyze its properties.
    
    Args:
        edges: List of tuples (gene1, gene2, weight)
        disease_genes: List of disease-related genes
        
    Returns:
        NetworkX graph object
    """
    G = nx.Graph()
    for edge in edges:
        G.add_edge(edge[0], edge[1], weight=edge[2])
    
    # Calculate network metrics
    degrees = dict(G.degree())
    avg_degree = sum(degrees.values()) / len(degrees) if degrees else 0
    
    network_info = {
        "nodes": list(G.nodes),
        "edges": [(u, v, d['weight']) for u, v, d in G.edges(data=True)],
        "avg_degree": round(avg_degree, 2),
    }
    
    # Identify disease genes
    if disease_genes:
        detected_disease_genes = set(G.nodes) & set(disease_genes)
        network_info["detected_disease_genes"] = list(detected_disease_genes)
        
        # Identify associated diseases
        diseases = []
        for gene in detected_disease_genes:
            if gene in DISEASE_DATABASE:
                diseases.append({"gene": gene, "disease": DISEASE_DATABASE[gene]})
        network_info["diseases"] = diseases
    
    return G, network_info

def find_shortest_path(G, start_gene, end_gene):
    """
    Find the shortest path between two genes in the network.
    
    Args:
        G: NetworkX graph
        start_gene: Starting gene
        end_gene: Target gene
        
    Returns:
        Dictionary with path information
    """
    if start_gene not in G.nodes or end_gene not in G.nodes:
        return {"error": f"One or both genes are not in the network"}
    
    try:
        path = nx.shortest_path(G, source=start_gene, target=end_gene, weight='weight')
        path_length = nx.shortest_path_length(G, source=start_gene, target=end_gene, weight='weight')
        
        return {
            "path": path,
            "path_length": path_length
        }
    except nx.NetworkXNoPath:
        return {"error": f"No path found between {start_gene} and {end_gene}"}