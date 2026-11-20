
import networkx as nx


def resnik_similarity(term1, term2, hpo_graph, information_content):
    """
    Calculate Resnik semantic similarity between two HPO terms.

    Parameters
    ----------
    term1 : str
        First HPO term identifier.
    term2 : str
        Second HPO term identifier.
    hpo_graph : networkx.DiGraph
        HPO ontology graph.
    information_content : dict
        Information Content value for each HPO term.

    Returns
    -------
    float
        Resnik semantic similarity score.
    """

    if term1 not in hpo_graph or term2 not in hpo_graph:
        raise ValueError("One or both HPO terms are not in the ontology.")

    # Include each term itself and all its ancestors.
    ancestors1 = nx.descendants(hpo_graph, term1) | {term1}
    ancestors2 = nx.descendants(hpo_graph, term2) | {term2}

    # Find common ancestors.
    common_ancestors = ancestors1.intersection(ancestors2)

    if not common_ancestors:
        return 0.0

    # Find the most informative common ancestor (MICA).
    mica = max(
        common_ancestors,
        key=lambda term: information_content.get(term, 0.0)
    )

    return information_content.get(mica, 0.0)



def bma_similarity(profile1, profile2, hpo_graph, information_content):
    """
    Calculate Best Match Average (BMA) similarity between
    two phenotype profiles using Resnik semantic similarity.

    Parameters
    ----------
    profile1 : set
        HPO terms describing the first phenotype profile.
    profile2 : set
        HPO terms describing the second phenotype profile.
    hpo_graph : networkx.DiGraph
        HPO ontology graph.
    information_content : dict
        Information Content value for each HPO term.

    Returns
    -------
    float
        Best Match Average similarity score.
    """

    if not profile1 or not profile2:
        raise ValueError("Both phenotype profiles must be non-empty.")

    # For each term in profile1, find its best match in profile2.
    best_matches_1 = [
        max(
            resnik_similarity(term1, term2, hpo_graph, information_content)
            for term2 in profile2
        )
        for term1 in profile1
    ]

    # For each term in profile2, find its best match in profile1.
    best_matches_2 = [
        max(
            resnik_similarity(term2, term1, hpo_graph, information_content)
            for term1 in profile1
        )
        for term2 in profile2
    ]

    # Average all best-match scores from both directions.
    return (
        sum(best_matches_1) + sum(best_matches_2)
    ) / (len(profile1) + len(profile2))
