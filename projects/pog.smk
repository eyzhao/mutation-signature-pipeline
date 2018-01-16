class CohortError(Exception):
    pass

pog_analysis_groups = pd.read_csv('manual/POG500_analysis_groups.txt', sep='\t')
pog_analysis_groups = pog_analysis_groups[
    (~ pd.isnull(pog_analysis_groups['ID'])) &
    (pog_analysis_groups['Biopsy #'] == 1.0)
]

pog_custom_cohorts = pd.read_csv('manual/pog_custom_cohorts.tsv', sep='\t')

def get_pogs_in_cohort(cohort):
    if cohort in pog_analysis_groups['ANALYSIS COHORT'].tolist():
        pogs = pog_analysis_groups[
            (pog_analysis_groups['ANALYSIS COHORT'] == cohort) &
            (pog_analysis_groups['Biopsy #'] == 1.0)
        ]['ID']
        return pd.DataFrame({'samples': pogs})
    elif cohort in pog_custom_cohorts['cohort'].tolist():
        pogs = pog_custom_cohorts[
            pog_custom_cohorts['cohort'] == cohort
        ]['id']
        return pd.DataFrame({'samples': pogs})
    else:
        raise(CohortError('No POG cohort found: {}'.format(cohort)))

def get_pog_snv_targets(wildcards):
    return get_pogs_in_cohort(wildcards.cohort)

target_functions['pog'] = {
    'snv': get_pog_snv_targets
}

project_cohorts['pog'] = [
    'HNSC'
]

project_genome['pog'] = 'BSgenome.Hsapiens.UCSC.hg19'
maf_colnames['pog'] = 'chr,pos,ref,alt'

pog_project_root = '/projects/POG/POG_data'

def get_pog_snv_files(wildcards, sample_label = 'sample'):
    print(wildcards)
    return(collapse_samples(wildcards.sample))

def is_numeric(x):
    try:
        int(x)
        return True
    except:
        return False

def collapse_samples(patient):
    ''' Given a POG sample with multiple samples, this function automatically selects the preferred
        sample and its highest numbered strelka run.
    '''
    available_samples = os.listdir('/projects/POG/POG_data/{}/wgs'.format(patient))
    sample = [s for s in available_samples if len(s.split('_')) == 6 and s.startswith('biop1')][0]

    (t_prefix, t, t_lib, n_prefix, n, n_lib) = sample.split('_')
    path_template = pog_project_root + '/{patient}/wgs/{sample}/{t_lib}_{n_lib}/strelka/{{strelka_run_id}}/bwa/results/passed.somatic.snvs.vcf'
    path = path_template.format(
        patient = patient,
        sample = sample,
        t_lib = t_lib,
        n_lib = n_lib
    )

    (all_strelka_runs,) = glob_wildcards(path)
    all_strelka_runs = [i for i in all_strelka_runs if is_numeric(i)]
    numeric_strelka_runs = [int(i) for i in all_strelka_runs]
    chosen_run = [str(max(numeric_strelka_runs))]
    return expand(path, strelka_run_id = chosen_run)

rule pog_vcf_to_tsv:
    input:
         get_pog_snv_files
    output:
        'data/pog/genome/cohorts/{cohort}/snv/{sample}.maf'
    shell:
        'Rscript scripts/pipelines/somatic-characterization/vcf_to_tsv.R -v {input} -o {output}'
