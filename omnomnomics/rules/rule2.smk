

## Omnomnomics Snake Rule  ##
import os

def update_semaphore_file(experiment_dir, rule_name):
    semaphore_file = os.path.join(experiment_dir, "omnomnomics.semaphore")
    with open(semaphore_file, "a") as f:
        f.write(f"{rule_name}\n")

rule rule2:
    input:
        "trimmed_FASTQ/rule1_output.{sample}.fastq.gz"
    output:
        "fastqc_reports/rule2_output.{sample}.fastq.gz"
    params:
        experiment_dir=config['EXPERIMENT_DIR']
    threads:
        6
    resources:
        mem_mb = 4000
    shell:
        """
        echo "Rule 2" >> {input}
        cp {input} {output}
        echo "2" >> omnomnomics.semaphore
        """ 
    # run:
    #     update_semaphore_file(params.experiment_dir, 2)

# onsuccess:
#     update_semaphore_file(snakemake.params.experiment_dir,2)