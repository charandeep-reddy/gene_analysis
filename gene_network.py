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
    # Make sure the DISEASE_DATABASE variable is defined
    # If it's not defined in the file, add it:
    if 'DISEASE_DATABASE' not in globals():
        DISEASE_DATABASE = {
            "BRCA1": "Breast cancer",
            "TP53": "Li-Fraumeni syndrome",
            "EGFR": "Lung cancer",
            # Add more as needed
        }

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
    # Check if both genes exist in the network
    if start_gene not in G.nodes:
        return {"error": f"Start gene '{start_gene}' not found in the network."}
    if end_gene not in G.nodes:
        return {"error": f"End gene '{end_gene}' not found in the network."}
    
    try:
        # Find shortest path using NetworkX
        path = nx.shortest_path(G, source=start_gene, target=end_gene)
        path_length = nx.shortest_path_length(G, source=start_gene, target=end_gene)
        return {"path": path, "path_length": path_length}
    except nx.NetworkXNoPath:
        return {"error": f"No path exists between '{start_gene}' and '{end_gene}'."}
    except Exception as e:
        return {"error": f"An error occurred: {str(e)}"}