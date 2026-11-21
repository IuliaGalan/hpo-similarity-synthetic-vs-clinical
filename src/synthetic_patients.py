
import random


def generate_synthetic_patient(
    disease_profiles,
    hpo_graph,
    retention_rate=0.6,
    noise_count=2,
    random_seed=None
):
    """
    Generate one synthetic patient from a randomly selected disease.

    Parameters
    ----------
    disease_profiles : dict
        Mapping between disease IDs and their HPO phenotype sets.
    hpo_graph : networkx.DiGraph
        HPO ontology graph.
    retention_rate : float
        Proportion of disease phenotypes to retain.
    noise_count : int
        Number of unrelated phenotype terms to add.
    random_seed : int or None
        Seed for reproducibility.

    Returns
    -------
    dict
        Synthetic patient information, including the true disease ID,
        retained phenotype terms and added noise terms.
    """

    if not 0 < retention_rate <= 1:
        raise ValueError("retention_rate must be in (0, 1].")

    if noise_count < 0:
        raise ValueError("noise_count must be non-negative.")

    rng = random.Random(random_seed)

    # Keep diseases with at least one valid HPO phenotype.
    valid_profiles = {
        disease_id: set(terms) & set(hpo_graph.nodes)
        for disease_id, terms in disease_profiles.items()
        if set(terms) & set(hpo_graph.nodes)
    }

    if not valid_profiles:
        raise ValueError("No diseases with valid HPO phenotype terms.")

    # Randomly select a disease.
    disease_id = rng.choice(list(valid_profiles.keys()))
    disease_terms = valid_profiles[disease_id]

    # Select a subset of the disease phenotypes.
    retained_count = max(
        1,
        round(len(disease_terms) * retention_rate)
    )

    original_terms = set(
        rng.sample(sorted(disease_terms), retained_count)
    )

    # Build a pool of HPO terms not belonging to the selected disease.
    noise_pool = set(hpo_graph.nodes) - disease_terms

    if noise_count > len(noise_pool):
        raise ValueError("Not enough HPO terms available for noise.")

    noise_terms = set(
        rng.sample(sorted(noise_pool), noise_count)
    )

    # Combine original and noise phenotypes.
    patient_terms = original_terms | noise_terms

    return {
        "true_disease_id": disease_id,
        "phenotype_terms": patient_terms,
        "original_terms": original_terms,
        "noise_terms": noise_terms
    }


def sample_phenotypes_from_frequency_profile(
    frequency_profile,
    unknown_frequency_sample_size=11,
    random_seed=None
):
    """
    Sample phenotype terms from a disease frequency profile.

    Parameters
    ----------
    frequency_profile : dict
        Mapping between HPO terms and their associated frequencies.
        Frequencies can be floats or None.

    unknown_frequency_sample_size : int
        Maximum number of phenotypes to sample when all frequencies
        are unknown.

    random_seed : int or None
        Seed for reproducibility.

    Returns
    -------
    dict
        Selected phenotype terms and the sampling strategy used.
    """

    if unknown_frequency_sample_size < 1:
        raise ValueError(
            "unknown_frequency_sample_size must be at least 1."
        )

    rng = random.Random(random_seed)

    # Validate known frequencies.
    for term, frequency in frequency_profile.items():
        if frequency is not None and not 0 <= frequency <= 1:
            raise ValueError(
                f"Invalid frequency for {term}: {frequency}"
            )

    # Identify phenotypes with known positive frequencies.
    positive_frequency_terms = {
        term: frequency
        for term, frequency in frequency_profile.items()
        if frequency is not None and frequency > 0
    }

    # Case 1: At least one known positive frequency.
    if positive_frequency_terms:

        selected_terms = {
            term
            for term, probability in positive_frequency_terms.items()
            if rng.random() < probability
        }

        return {
            "phenotype_terms": selected_terms,
            "sampling_strategy": "frequency_weighted"
        }

    # Case 2: All frequencies are unknown.
    if frequency_profile and all(
        frequency is None
        for frequency in frequency_profile.values()
    ):

        available_terms = sorted(frequency_profile.keys())

        sample_size = min(
            unknown_frequency_sample_size,
            len(available_terms)
        )

        selected_terms = set(
            rng.sample(available_terms, sample_size)
        )

        return {
            "phenotype_terms": selected_terms,
            "sampling_strategy": "uniform"
        }

    # Case 3: No positive known frequencies.
    return {
        "phenotype_terms": set(),
        "sampling_strategy": "ineligible"
    }




def generate_synthetic_patient_for_disease(
    disease_id,
    disease_frequency_profiles,
    unknown_frequency_sample_size=11,
    random_seed=None,
    max_attempts=100
):
    """
    Generate one synthetic patient for a specific disease.

    Parameters
    ----------
    disease_id : str
        Identifier of the target disease.
    disease_frequency_profiles : dict
        Mapping between disease IDs and their phenotype frequencies.
    unknown_frequency_sample_size : int
        Maximum number of phenotypes sampled when all frequencies
        for the disease are unknown.
    random_seed : int or None
        Seed for reproducibility.
    max_attempts : int
        Maximum number of attempts to generate a valid patient.

    Returns
    -------
    dict
        Synthetic patient information, including the true disease ID,
        selected phenotype terms and sampling strategy.
    """

    if disease_id not in disease_frequency_profiles:
        raise KeyError(
            f"Disease ID not found: {disease_id}"
        )

    if max_attempts < 1:
        raise ValueError("max_attempts must be at least 1")

    frequency_profile = disease_frequency_profiles[disease_id]

    rng = random.Random(random_seed)

    for attempt in range(max_attempts):

        sampling_result = sample_phenotypes_from_frequency_profile(
            frequency_profile=frequency_profile,
            unknown_frequency_sample_size=unknown_frequency_sample_size,
            random_seed=rng.randint(0, 2**32 - 1)
        )

        phenotype_terms = sampling_result["phenotype_terms"]
        sampling_strategy = sampling_result["sampling_strategy"]

        if phenotype_terms:
            return {
                "true_disease_id": disease_id,
                "phenotype_terms": phenotype_terms,
                "sampling_strategy": sampling_strategy,
                "is_valid": True
            }

    return {
        "true_disease_id": disease_id,
        "phenotype_terms": set(),
        "sampling_strategy": sampling_strategy,
        "is_valid": False
    }



def generate_synthetic_patients(
    disease_frequency_profiles,
    n_per_disease=100,
    unknown_frequency_sample_size=11,
    random_seed=None
):
    """
    Generate multiple synthetic patients for every disease.

    Parameters
    ----------
    disease_frequency_profiles : dict
        Mapping between disease IDs and phenotype frequency profiles.

    n_per_disease : int
        Number of patients to generate for each disease.

    unknown_frequency_sample_size : int
        Maximum number of phenotypes sampled when all frequencies
        are unknown.

    random_seed : int or None
        Seed for reproducibility.

    Returns
    -------
    list
        List of generated synthetic patients.
    """

    if n_per_disease < 1:
        raise ValueError("n_per_disease must be at least 1.")

    rng = random.Random(random_seed)

    patients = []

    for disease_id in sorted(disease_frequency_profiles):

        for patient_number in range(n_per_disease):

            patient = generate_synthetic_patient_for_disease(
                disease_id=disease_id,
                disease_frequency_profiles=disease_frequency_profiles,
                unknown_frequency_sample_size=unknown_frequency_sample_size,
                random_seed=rng.randint(0, 2**32 - 1)
            )

            patient["patient_id"] = (
                f"{disease_id}_patient_{patient_number + 1}"
            )

            patients.append(patient)

    return patients

