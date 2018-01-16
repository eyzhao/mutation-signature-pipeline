def get_tcga_snv_targets(wildcards):
    df = glob_to_df(
        'data/tcga/{method}/cohorts/{cohort}/snv/*.maf'.format(
            method = wildcards.method,
            cohort = wildcards.cohort
        ),
        r'data\/tcga\/' + wildcards.method + '\/cohorts\/' + wildcards.cohort + '\/snv\/(.*?).maf',
        ['samples']
    )
    return df

target_functions['tcga'] = {
    'snv': get_tcga_snv_targets
}

project_cohorts['tcga'] = [
    'BLCA',
    'BRCA',
    'CESC',
    'COAD',
    'LUAD',
    'LUSC',
    'SKCM',
    'STAD',
    'UCEC'
]

project_genome['tcga'] = 'BSgenome.Hsapiens.UCSC.hg38'
maf_colnames['tcga'] = 'Chromosome,Start_Position,Reference_Allele,Allele'

rule get_tcga_snv_counts:
    input:
        'data/tcga/exome/snv_maf/TCGA.{cohort}.somatic.merged.maf'
    output:
        'data/tcga/exome/snv_maf_summary/TCGA.{cohort}.tsv'
    shell:
        replace_newlines('''
            Rscript -e "
                library(tidyverse);
                read_tsv('{input}') %>%
                    mutate(include = var_caller_count > 1) %>%
                    group_by(include) %>%
                    summarise(snv_count = n()) %>%
                    mutate(cohort = '{wildcards.cohort}') %>%
                    write_tsv('{output}');
            "
        ''')

rule aggregate_tcga_snv_counts:
    input:
        expand('data/tcga/exome/snv_maf_summary/TCGA.{cohort}.tsv', cohort = project_cohorts['tcga'])
    params:
        comma_separated = lambda wildcards, input: ','.join(input)
    output:
        'data/tcga/exome/snv_count_summary.tsv'
    shell:
        'Rscript scripts/pipelines/basic-functions/row_bind_tables.R -p {params.comma_separated} -o {output} --trim-paths'
