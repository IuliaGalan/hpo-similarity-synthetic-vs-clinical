
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



def bma_similarity(profile1, profile2, hpo_graph, information_content, similarity_measure=resnik_similarity):
    """
    Calculate Best Match Average (BMA) similarity between
    two phenotype profiles using a semantic similarity measure.
    """

    if not profile1 or not profile2:
        raise ValueError("Both phenotype profiles must be non-empty.")

    best_matches_1 = [
        max(
            similarity_measure(term1, term2, hpo_graph, information_content)
            for term2 in profile2
        )
        for term1 in profile1
    ]

    best_matches_2 = [
        max(
            similarity_measure(term2, term1, hpo_graph, information_content)
            for term1 in profile1
        )
        for term2 in profile2
    ]

    return (
        sum(best_matches_1) + sum(best_matches_2)
    ) / (len(profile1) + len(profile2))



def asymmetric_similarity(
    patient_profile,
    disease_profile,
    hpo_graph,
    information_content,
    similarity_measure=resnik_similarity
):
    """
    Calculate asymmetric similarity between a patient and a disease
    using a semantic similarity measure.
    """

    if not patient_profile or not disease_profile:
        raise ValueError("Both patient and disease profiles must be non-empty.")

    best_matches = [
        max(
            similarity_measure(
                patient_term,
                disease_term,
                hpo_graph,
                information_content
            )
            for disease_term in disease_profile
        )
        for patient_term in patient_profile
    ]

    return sum(best_matches) / len(best_matches)



def lin_similarity(term1, term2, hpo_graph, information_content):
    """
    Calculate Lin semantic similarity between two HPO terms.
    """
    if term1 not in hpo_graph or term2 not in hpo_graph:
        raise ValueError("One or both HPO terms are not in the ontology.")

    ancestors1 = nx.descendants(hpo_graph, term1) | {term1}
    ancestors2 = nx.descendants(hpo_graph, term2) | {term2}

    common_ancestors = ancestors1.intersection(ancestors2)

    if not common_ancestors:
        return 0.0

    mica = max(
        common_ancestors,
        key=lambda term: information_content.get(term, 0.0)
    )

    ic_mica = information_content.get(mica, 0.0)
    ic_term1 = information_content.get(term1, 0.0)
    ic_term2 = information_content.get(term2, 0.0)

    denominator = ic_term1 + ic_term2

    if denominator == 0:
        return 0.0

    return (2 * ic_mica) / denominator



