# HPO Similarity: Synthetic vs. Clinical Patients

## 1. Research Objective

This project investigates how the performance of Human Phenotype Ontology (HPO) semantic similarity methods differs between synthetic patients and real clinical cases reported in the literature.

The main objective is to determine whether performance evaluations based on synthetic patients accurately reflect the performance of phenotype-based disease prioritization methods on clinical data.

## 2. Research Questions

- **RQ1:** How does the performance of the same HPO similarity methods differ between synthetic patients and real clinical cases?
- **RQ2:** To what extent can these differences be explained by measurable factors such as the number and specificity of patient phenotypes and disease rarity?

## 3. Data Sources

The project uses the following data sources:

- **HPO Ontology (`hp.obo`):** Provides the ontology structure, represented as a Directed Acyclic Graph (DAG).
- **HPO Disease Annotations (`phenotype.hpoa`):** Provides phenotype annotations associated with diseases.
- **Phenopacket Store:** Planned source of standardized clinical cases containing HPO terms and confirmed diagnoses.

## 4. Current Implementation

The initial implementation focuses on preparing and exploring the HPO data.

The following components have been implemented:

- Loading and parsing the HPO ontology.
- Exploring the ontology graph structure.
- Loading and filtering disease-phenotype annotations.
- Validating annotated HPO terms against the ontology.
- Calculating Information Content (IC) values.
- Implementing Resnik semantic similarity between individual HPO terms.

Initial data exploration identified:
- 19,904 ontology nodes.
- 24,453 ontology edges.
- 286,651 initial annotations.
- 268,179 positive annotations after filtering.
- 12,866 diseases.
- 11,560 unique annotated HPO terms, all present in the ontology.

The current implementation is available in:

`notebooks/01_hpo_data_exploration.ipynb`

Disease-level ranking and experimental evaluation are not yet implemented.

## 5. Planned Experimental Methodology

The planned evaluation compares nine method variants:

1. Resnik + Best Match Average (BMA).
2. Resnik + asymmetric aggregation.
3. Lin + BMA.
4. Lin + asymmetric aggregation.
5. Jiang-Conrath + BMA.
6. Jiang-Conrath + asymmetric aggregation.
7. Phenomizer.
8. Jaccard similarity.
9. LIRICAL (phenotype-only mode).

The evaluation will use:

- Recall@1, Recall@3, Recall@10 and Recall@50.
- Mean Reciprocal Rank (MRR).

Additional experiments will investigate the effects of phenotype count, information content, ontology depth and disease rarity.

## 6. Project Structure

```text
hpo-similarity-synthetic-vs-clinical/
│
├── data/
│   └── hpo/
│
├── notebooks/
│   └── 01_hpo_data_exploration.ipynb
│
├── src/
│
├── results/
│
├── .gitignore
└── README.md
```

## 7. Planned Work

1. Implement disease-level phenotype similarity and ranking.
2. Generate synthetic patient profiles from HPO disease annotations.
3. Implement evaluation metrics (Recall@k and MRR).
4. Run preliminary experiments on synthetic patients.
5. Introduce noise, phenotype imprecision and phenotype dropout.
6. Retrieve and preprocess clinical cases from Phenopacket Store.
7. Compare synthetic and clinical evaluation results.
8. Analyze factors explaining performance differences.

## 8. Research Reproducibility

The project is developed in Python. Source code, Jupyter notebooks and experimental results are maintained in this repository to support reproducibility and facilitate the evaluation of the proposed methodology.
