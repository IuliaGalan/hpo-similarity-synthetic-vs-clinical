

def frequency_to_probability(frequency):
    """
    Convert an HPO frequency annotation into a numeric value.

    Returns
    -------
    float or None
        Estimated frequency between 0 and 1.
        Returns None when the frequency is unknown.
    """

    if frequency is None:
        return None

    frequency = str(frequency).strip()

    if frequency.lower() in {"nan", "none", ""}:
        return None

    # HPO qualitative frequency categories.
    hpo_frequencies = {
        "HP:0040280": 1.00,             # Obligate
        "HP:0040281": (0.80 + 0.99) / 2,  # Very frequent
        "HP:0040282": (0.30 + 0.79) / 2,  # Frequent
        "HP:0040283": (0.05 + 0.29) / 2,  # Occasional
        "HP:0040284": (0.01 + 0.04) / 2   # Very rare
    }

    if frequency in hpo_frequencies:
        return hpo_frequencies[frequency]

    # Handle percentages, e.g. "25%".
    if frequency.endswith("%"):
        value = float(frequency[:-1]) / 100

        if not 0 <= value <= 1:
            raise ValueError(f"Invalid frequency: {frequency}")

        return value

    # Handle fractional frequencies, e.g. "1/2".
    if "/" in frequency:
        numerator, denominator = frequency.split("/")

        numerator = int(numerator)
        denominator = int(denominator)

        if denominator <= 0 or numerator < 0 or numerator > denominator:
            raise ValueError(f"Invalid frequency: {frequency}")

        return numerator / denominator

    raise ValueError(f"Unsupported frequency format: {frequency}")



import pandas as pd


def aggregate_phenotype_frequencies(annotations):
    """
    Aggregate frequency annotations for each disease-phenotype pair.

    Unknown frequencies are excluded from the mean but counted separately.
    The original annotations DataFrame is not modified.
    """

    data = annotations[
        ["database_id", "hpo_id", "frequency"]
    ].copy()

    data["numeric_frequency"] = data["frequency"].apply(
        frequency_to_probability
    )

    aggregated = (
        data.groupby(["database_id", "hpo_id"])
        .agg(
            frequency=("numeric_frequency", "mean"),
            annotation_count=("frequency", "size"),
            unknown_frequency_count=("numeric_frequency",
                                     lambda values: values.isna().sum())
        )
        .reset_index()
    )

    return aggregated


def build_disease_frequency_profiles(aggregated_frequencies):
    """
    Build a mapping between disease IDs and phenotype frequencies.

    Each disease maps HPO terms to their aggregated frequency.
    Unknown frequencies are represented by None.
    """

    disease_frequency_profiles = {}

    for _, row in aggregated_frequencies.iterrows():
        disease_id = row["database_id"]
        hpo_id = row["hpo_id"]
        frequency = row["frequency"]

        if disease_id not in disease_frequency_profiles:
            disease_frequency_profiles[disease_id] = {}

        disease_frequency_profiles[disease_id][hpo_id] = (
            None if pd.isna(frequency) else float(frequency)
        )

    return disease_frequency_profiles
