INPUTS = [
    {
        "source": "SLE_v1",
        "path": "pipeline/3.1/data/raw/GSE65391_series_matrix.txt.gz",
        "target": "SLE_v1/pipeline/3.1/data/raw/GSE65391_series_matrix.txt.gz",
        "stages": ["bulk_prepare"],
    },
    {
        "source": "SLE_v1",
        "path": "pipeline/3.1/data/raw/GPL10558.annot.gz",
        "target": "SLE_v1/pipeline/3.1/data/raw/GPL10558.annot.gz",
        "stages": ["bulk_prepare"],
    },
    *[
        {
            "source": "SLE",
            "path": f"metadata/source/geo_matrix_headers/{dataset}_series_matrix.txt.header.txt",
            "target": f"SLE/metadata/source/geo_matrix_headers/{dataset}_series_matrix.txt.header.txt",
            "stages": ["bulk_metadata"],
        }
        for dataset in ["GSE65391", "GSE88884", "GSE88885", "GSE88886", "GSE135779", "GSE121239"]
    ],
    *[
        {
            "source": "SLE",
            "path": f"metadata/source/supplementary/{name}",
            "target": f"SLE/metadata/source/supplementary/{name}",
            "stages": ["bulk_metadata"],
        }
        for name in [
            "GSE135779_Nature_Supplementary_Tables.xlsx",
            "GSE135779_Meta_caSLE_processed_08092021_small.csv",
        ]
    ],
    *[
        {
            "source": "SLE",
            "path": f"GSE{number}_pipeline/data/raw/downloads/{name}",
            "target": f"SLE/GSE{number}_pipeline/data/raw/downloads/{name}",
            "stages": [f"bulk_hta_{number}"],
        }
        for number, name in [
            ("88884", "GSE88884_ILLUMINATE1and2_SLEbaselineVsHealthy_preprocessed.txt.gz"),
            ("88885", "GSE88885_ILLUMINATE1_pd_preprocessed.txt.gz"),
            ("88886", "GSE88886_ILLUMINATE2_pd_preprocessed.txt.gz"),
            ("88884", "GSE88884_series_matrix.txt.gz"),
            ("88885", "GSE88885_series_matrix.txt.gz"),
            ("88886", "GSE88886_series_matrix.txt.gz"),
        ]
    ],
    *[
        {
            "source": "SLE_v1",
            "path": f"pipeline/3.2/config/gene_sets/{name}",
            "target": f"SLE_v1/pipeline/3.2/config/gene_sets/{name}",
            "stages": ["bulk_adult", "figure_module_annotation"],
        }
        for name in ["c5.go.symbols.gmt", "c2.cp.kegg.symbols.gmt"]
    ],
]
STAGES = [
    {
        "id": "bulk_prepare",
        "script": "steps/bulk/prepare.py",
        "engine": "python",
        "depends": [],
        "outputs": [
            "SLE_v1/pipeline/3.1/config/programs.tsv",
            "SLE_v1/pipeline/3.2/config/analysis_contract.tsv",
        ],
    },
    {
        "id": "bulk_metadata",
        "script": "steps/bulk/metadata.py",
        "engine": "python",
        "depends": ["bulk_prepare"],
        "outputs": [
            "SLE/metadata/frozen/2026-07-25_v1.0/by_dataset/GSE65391_meta.tsv",
            "SLE/metadata/frozen/2026-07-25_v1.0/by_dataset/GSE135779_meta.tsv",
            "SLE/metadata/frozen/2026-07-25_v1.0/by_dataset/GSE121239_meta.tsv",
            "SLE_v1/pipeline/3.3/data/metadata/GSE88884_meta.tsv",
        ],
    },
    {
        "id": "bulk_child",
        "script": "SLE_v1/pipeline/3.1/run_3_1.R",
        "engine": "r",
        "depends": ["bulk_metadata"],
        "outputs": [
            "SLE_v1/pipeline/3.1/data/processed/01_gene_expression.rds",
            "SLE_v1/pipeline/3.1/data/processed/04_wgcna.rds",
            "SLE_v1/results/3.1/06_wgcna/module_membership.tsv",
        ],
    },
    {
        "id": "bulk_child_programs",
        "script": "steps/bulk/child_programs.R",
        "engine": "r",
        "depends": ["bulk_child"],
        "outputs": [
            "SLE_v1/pipeline/3.2/data/reference/expanded_program_membership.tsv",
            "SLE/GSE200306_pipeline/results/07_cross_dataset_bridge/expanded_program_membership.tsv",
        ],
    },
    *[
        {
            "id": f"bulk_hta_{number}",
            "script": f"SLE/GSE{number}_pipeline/run_pipeline.R",
            "engine": "r",
            "depends": ["bulk_child_programs"],
            "outputs": [
                f"SLE/GSE{number}_pipeline/data/processed/analysis_results.rds",
                f"SLE/GSE{number}_pipeline/results/02_program_scores/sample_program_scores.tsv",
            ],
        }
        for number in ["88884", "88885", "88886"]
    ],
    {
        "id": "bulk_hta_inputs",
        "script": "steps/bulk/bridge.py",
        "engine": "python",
        "args": ["hta"],
        "depends": ["bulk_hta_88884", "bulk_hta_88885", "bulk_hta_88886"],
        "outputs": [
            "SLE_v1/pipeline/3.3/data/processed/GSE88884_analysis_results.rds",
            "SLE_v1/pipeline/3.4/data/processed/ILLUMINATE1_sample_program_scores.tsv",
            "SLE_v1/pipeline/3.2/data/metadata/GSE88884_sample_audit.tsv",
            "SLE_R1/checks/r1_revision_20261116/supplement/frozen_bulk_program_coverage.tsv",
        ],
    },
    {
        "id": "bulk_adult",
        "script": "SLE_v1/pipeline/3.2/run_3_2.R",
        "engine": "r",
        "depends": ["bulk_hta_inputs"],
        "outputs": [
            "SLE_v1/results/3.2/01_expression/full_gene_expression.rds",
            "SLE_v1/results/3.2/03_preservation/cross_cohort_module_membership.tsv",
            "SLE_v1/results/3.2/02_deg_gsea/full_transcriptome_limma.tsv",
        ],
    },
    {
        "id": "bulk_cross_cohort",
        "script": "steps/bulk/bridge.py",
        "engine": "python",
        "args": ["cross_cohort"],
        "depends": ["bulk_adult"],
        "outputs": [
            "SLE_v1/pipeline/3.5/data/processed/cross_cohort_module_membership.tsv",
            "SLE/SLE_mainline_research/02_cross_cohort_module_preservation/results/01_preservation/cross_cohort_module_membership.tsv",
        ],
    },
    {
        "id": "bulk_grouping",
        "script": "SLE_v1/pipeline/3.3/run_3_3.R",
        "engine": "r",
        "depends": ["bulk_hta_inputs"],
        "outputs": ["SLE_v1/results/3.3/01_grouping/frozen_grouping_model.rds"],
    },
    {
        "id": "bulk_trial_inputs",
        "script": "steps/bulk/bridge.py",
        "engine": "python",
        "args": ["grouping"],
        "depends": ["bulk_grouping"],
        "outputs": ["SLE_v1/pipeline/3.4/data/processed/frozen_grouping_model.rds"],
    },
    {
        "id": "bulk_trials",
        "script": "SLE_v1/pipeline/3.4/run_3_4.R",
        "engine": "r",
        "depends": ["bulk_trial_inputs"],
        "outputs": ["SLE_v1/results/3.4/_SUCCESS.rds"],
    },
    {
        "id": "bulk_landscape",
        "script": "SLE_v1/analysis/exploratory_program_landscape.R",
        "engine": "r",
        "depends": ["bulk_trial_inputs"],
        "outputs": [
            "SLE_v1/results/exploratory_program_landscape/01_pairwise_clinical/program_high_vs_immune_low_clinical_contrasts.tsv",
            "SLE_v1/results/exploratory_program_landscape/02_continuous_landscape/clinical_model_comparison.tsv",
        ],
    },
    {
        "id": "bulk_clinical_validation",
        "script": "SLE_v1/analysis/reviewer_clinical_landscape_validation.R",
        "engine": "r",
        "depends": ["bulk_trial_inputs"],
        "outputs": [
            "SLE_v1/results/reviewer_response/02_clinical_landscape_validation/formal_cv_continuous_vs_four_group.tsv",
            "SLE_v1/results/reviewer_response/02_clinical_landscape_validation/heldout_continuous_vs_four_group.tsv",
        ],
    },
    {
        "id": "bulk_trial_ancova",
        "script": "SLE_v1/analysis/reviewer_trial_reporting_ancova_permutation.R",
        "engine": "r",
        "depends": ["bulk_trial_inputs"],
        "outputs": [
            "SLE_v1/results/reviewer_response/08_trial_reporting_ancova_permutation/ancova_coordinate_effects_trial_specific.tsv",
            "SLE_v1/results/reviewer_response/08_trial_reporting_ancova_permutation/ancova_vector_bootstrap_draws.rds",
        ],
    },
    {
        "id": "bulk_abundance",
        "script": "SLE_v1/analysis/reviewer_bcell_abundance_analysis.R",
        "engine": "r",
        "depends": ["bulk_trial_inputs", "sc_legacy_scores"],
        "outputs": [
            "SLE_v1/results/reviewer_response/07_bcell_abundance/abundance_adjusted_treatment_models.tsv",
            "SLE_v1/results/reviewer_response/07_bcell_abundance/b_program_abundance_correlations.tsv",
        ],
    },
    {
        "id": "bulk_module_sensitivity",
        "script": "SLE_R1/Try/R2_02_03_10_12_modules/run_modules.R",
        "engine": "r",
        "depends": ["bulk_child"],
        "outputs": ["SLE_R1/Try/R2_02_03_10_12_modules/participant_split_summary.tsv"],
    },
    {
        "id": "bulk_plasmablast",
        "script": "SLE_R1/Try/R2_02_03_10_12_modules/plasmablast_sensitivity.R",
        "engine": "r",
        "depends": ["bulk_module_sensitivity"],
        "outputs": ["SLE_R1/Try/R2_02_03_10_12_modules/plasmablast_proxy_clinical_sensitivity.tsv"],
    },
    {
        "id": "bulk_adult_network",
        "script": "SLE_R1/Try/R2_05_08_adult/run_adult.R",
        "engine": "r",
        "depends": ["bulk_adult", "bulk_trial_inputs"],
        "outputs": ["SLE_R1/Try/R2_05_08_adult/adult_network.rds"],
    },
]
R_PACKAGES = [
    "data.table",
    "digest",
    "matrixStats",
    "limma",
    "WGCNA",
    "dynamicTreeCut",
    "fastcluster",
    "fgsea",
    "AnnotationDbi",
    "hta20transcriptcluster.db",
    "org.Hs.eg.db",
    "GO.db",
    "ggplot2",
    "ggrepel",
    "patchwork",
    "scales",
    "cluster",
    "ComplexHeatmap",
    "circlize",
    "ggalluvial",
    "lme4",
    "emmeans",
    "sandwich",
    "lmtest",
    "SeuratObject",
    "Matrix",
    "quantiseqr",
    "doParallel",
    "foreach",
    "readr",
]
PY_PACKAGES = ["openpyxl"]
