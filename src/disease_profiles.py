
def build_disease_profiles(phenotype_annotations, hpo_graph):
    """
    Build phenotype profiles for each disease.

    Each disease is represented by the set of its
    positively annotated HPO terms.

    Parameters
    ----------
    phenotype_annotations : pandas.DataFrame
        Positive disease-phenotype annotations.
    hpo_graph : networkx.DiGraph
        HPO ontology graph.

    Returns
    -------
    dict
        Dictionary mapping disease identifiers to sets
        of HPO term identifiers.
    """

    required_columns = {"database_id", "hpo_id"}

    if not required_columns.issubset(phenotype_annotations.columns):
        raise ValueError(
            "The annotation DataFrame must contain "
            "'database_id' and 'hpo_id' columns."
        )

    # Keep only the columns needed for disease profiles.
    annotations = phenotype_annotations[
        ["database_id", "hpo_id"]
    ].dropna().drop_duplicates()

    # Keep only HPO terms present in the ontology.
    annotations = annotations[
        annotations["hpo_id"].isin(hpo_graph.nodes)
    ]

    # Group HPO terms by disease.
    disease_profiles = (
        annotations
        .groupby("database_id")["hpo_id"]
        .apply(set)
        .to_dict()
    )

    return disease_profiles
