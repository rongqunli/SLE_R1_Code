from pathlib import Path
import csv
import os

root = Path(os.environ["SLE_RUN_ROOT"])


def write_rows(relative, fields, rows):
    target = root / relative
    target.parent.mkdir(parents=True, exist_ok=True)
    with target.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(fields)
        writer.writerows(rows)


write_rows(
    "SLE/GSE135779_pipeline/config/pbmc_markers.tsv",
    ["cell_type", "marker_genes"],
    [
        ["CD4_T_Naive", "CCR7;IL7R;LTB;MAL;TCF7;LEF1;NOSIP"],
        ["CD4_T_Memory", "IL7R;LTB;MAL;AQP3;IL32;CXCR4;GPR183"],
        ["Treg", "FOXP3;IL2RA;CTLA4;TIGIT;IKZF2;RTKN2"],
        ["CD8_T_Naive", "CCR7;LTB;MAL;TCF7;LEF1;CD8A;CD8B"],
        ["CD8_T_GZMK", "CD8A;CD8B;GZMK;CCL5;TRAC;CD3D;CD3E"],
        ["CD8_T_Cytotoxic", "CD8A;CD8B;NKG7;CCL5;GZMH;GZMB;PRF1;GNLY;CX3CR1"],
        ["NK", "NKG7;GNLY;KLRD1;KLRF1;PRF1;GZMB;TYROBP;FCGR3A"],
        ["B_Naive", "MS4A1;CD79A;CD79B;CD37;CD74;TCL1A;IGHM;IGHD"],
        ["B_Memory", "MS4A1;CD79A;CD37;CD74;CD27;GPR183;AIM2"],
        ["B_ABC", "ITGAX;TBX21;FCRL5;ZEB2;CD86;FCGR2B;MS4A1"],
        ["Plasmablast", "JCHAIN;MZB1;XBP1;DERL3;PRDM1;TNFRSF17;SDC1;IGKC"],
        ["CD14_Monocyte", "LYZ;LST1;S100A8;S100A9;S100A10;CTSD;CTSS;FCN1;VCAN"],
        ["CD16_Monocyte", "LYZ;LST1;FCGR3A;MS4A7;IFITM3;LILRB1;SERPINA1;SAT1"],
        ["cDC", "FCER1A;CD1C;CLEC10A;HLA-DRA;HLA-DPA1;CST3"],
        ["pDC", "GZMB;TCF4;CLEC4C;IL3RA;IRF7;SERPINF1"],
        ["Platelet", "PPBP;PF4;NRGN;GP9;ITGA2B;TUBB1"],
        ["Erythroid", "HBB;HBA1;HBA2;ALAS2;AHSP;GYPA"],
        ["Cycling", "MKI67;TOP2A;STMN1;TYMS;UBE2C;NUSAP1"],
    ],
)
write_rows(
    "SLE/GSE135779_pipeline/config/manual_cluster_annotations.tsv",
    ["cluster", "cell_type_major", "cell_type_analysis"],
    [
        ["0", "CD4_T", "CD4_T_Naive"],
        ["1", "CD4_T", "CD4_T"],
        ["2", "CD8_T", "CD8_T_Naive"],
        ["3", "Monocyte", "CD14_Monocyte"],
        ["4", "B", "B_Naive"],
        ["5", "CD8_T", "CD8_T_Cytotoxic"],
        ["6", "NK", "NK"],
        ["7", "Monocyte", "CD14_Monocyte"],
        ["8", "CD8_T", "CD8_T_GZMK"],
        ["9", "CD4_T", "CD4_T_Naive"],
        ["10", "Monocyte", "CD16_Monocyte"],
        ["11", "CD8_T", "CD8_T_Naive"],
        ["12", "Erythroid", "Erythroid"],
        ["13", "Monocyte", "Monocyte"],
        ["14", "Platelet", "Platelet"],
        ["15", "pDC", "pDC"],
        ["16", "B", "B"],
    ],
)
write_rows(
    "SLE/GSE135779_pipeline/config/sle_programs.tsv",
    ["program", "genes"],
    [
        [
            "IFN_I",
            "IFI6;IFI27;IFI35;IFI44;IFI44L;IFIT1;IFIT2;IFIT3;IFITM1;IFITM3;ISG15;MX1;MX2;OAS1;OAS2;OAS3;OASL;RSAD2;USP18;XAF1",
        ],
        ["IFN_II", "GBP1;GBP2;GBP4;GBP5;IRF1;STAT1;CXCL9;CXCL10;HLA-DRA;HLA-DPA1;HLA-DPB1;TAP1"],
        ["TNF_NFKB", "NFKBIA;NFKBIZ;TNFAIP3;REL;RELA;ICAM1;IL1B;CXCL2;CCL3;CCL4;JUN;FOS;IER3"],
        ["TLR_MYELOID", "TLR7;TLR8;MYD88;IRAK4;IRF7;NFKB1;NFKB2;S100A8;S100A9;FCN1;VCAN;CTSS"],
        [
            "BAFF_BCELL",
            "TNFRSF13B;TNFRSF13C;TNFRSF17;CD79A;CD79B;MS4A1;CD74;BANK1;BLK;CD37;HLA-DRA",
        ],
        ["ABC_BCELL", "ITGAX;TBX21;FCRL5;ZEB2;CD86;FCGR2B;AIM2;TLR7;LGALS3"],
        ["PLASMABLAST", "JCHAIN;MZB1;XBP1;DERL3;PRDM1;TNFRSF17;SDC1;SSR4;FKBP11;IGHG1;IGKC"],
        ["CYTOTOXICITY", "NKG7;GNLY;PRF1;GZMB;GZMH;CCL5;CTSW;KLRD1;KLRF1;FGFBP2"],
        ["T_EXHAUSTION", "PDCD1;LAG3;TIGIT;HAVCR2;CTLA4;TOX;ENTPD1;CXCL13"],
        [
            "ANTIGEN_PRESENTATION",
            "HLA-DRA;HLA-DRB1;HLA-DPA1;HLA-DPB1;HLA-DQA1;HLA-DQB1;CD74;CIITA;TAP1;TAP2",
        ],
        ["COMPLEMENT_MYELOID", "C1QA;C1QB;C1QC;C3;C5AR1;CFD;CFH;SERPING1;CTSC"],
        ["GLYCOLYSIS", "HK1;HK2;PFKP;PFKM;ALDOA;GAPDH;PGK1;ENO1;PKM;LDHA;SLC2A1"],
        [],
    ],
)
write_rows(
    "SLE/SLE_mainline_research/05_single_cell_module_projection/config/analysis_contract.tsv",
    ["key", "value"],
    [
        ["study_id", "05"],
        ["input_object", "../../../GSE135779_pipeline/data/processed/02_annotated_seurat.rds"],
        [
            "module_membership",
            "../../02_cross_cohort_module_preservation/results/01_preservation/cross_cohort_module_membership.tsv",
        ],
        [
            "pseudobulk_results",
            "../../../GSE135779_pipeline/results/05_pseudobulk/pseudobulk_all_genes.rds",
        ],
        ["primary_cohort", "Child"],
        ["supporting_cohort", "Adult"],
        ["reference_group", "HD"],
        ["test_group", "SLE"],
        ["primary_scoring_method", "UCell"],
        ["sensitivity_scoring_method", "log_normalized_mean"],
        [
            "module_map",
            "magenta=MAGENTA_IFN_ANTIVIRAL;cyan=CYAN_NEUTROPHIL_IMMATURE;greenyellow=GREENYELLOW_B_CELL_COMPLEMENT",
        ],
        ["focus_cell_types", "B_Naive;CD8_T_Cytotoxic;CD8_T_GZMK;NK;CD14_Monocyte"],
        ["excluded_clusters_primary", "13"],
        ["minimum_matched_genes", "10"],
        ["minimum_match_fraction", "0.20"],
        ["minimum_cells_per_donor_celltype", "20"],
        ["minimum_donors_per_group", "5"],
        ["minimum_pseudobulk_detected_genes", "5"],
        ["covariates_child", "batch;age;sex"],
        ["covariates_adult", "batch;age"],
        ["fdr_method", "BH"],
        ["fdr_threshold", "0.05"],
        ["confidence_level", "0.95"],
        ["random_seed", "20260730"],
    ],
)
write_rows(
    "SLE_v1/pipeline/3.5/config/analysis_contract.tsv",
    ["key", "value"],
    [
        ["study_id", "3.5_single_cell_module_projection"],
        ["donor_scores", "data/processed/donor_celltype_program_scores.tsv"],
        ["donor_proportions", "data/processed/donor_celltype_proportions.tsv"],
        ["program_audit", "data/processed/frozen_program_gene_matching.tsv"],
        ["module_membership", "data/processed/cross_cohort_module_membership.tsv"],
        ["pseudobulk_results", "data/processed/pseudobulk_all_genes.rds"],
        ["primary_cohort", "Child"],
        ["supporting_cohort", "Adult"],
        ["reference_group", "HD"],
        ["test_group", "SLE"],
        ["focus_cell_types", "B_Naive;CD8_T_Cytotoxic;CD8_T_GZMK;NK;CD14_Monocyte"],
        ["covariates_child", "batch;age;sex"],
        ["covariates_adult", "batch;age"],
        ["minimum_donors_per_group", "5"],
        ["confidence_level", "0.95"],
        ["fdr_method", "BH"],
        ["fdr_threshold", "0.05"],
        ["minimum_pseudobulk_detected_genes", "5"],
        ["random_seed", "20260730"],
    ],
)
