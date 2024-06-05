

## Omnomnomics Snake Rule  ##
rule rule4:
    input:
        "trimmed_FASTQ/rule1_output.{sample}.fastq.gz"
    output:
        "BAM/rule4_output.{sample}.bam"
    threads:
        6
    resources:
        mem_mb = 4000
    shell:
        """
        echo "Rule 4" >> {input}
        cp {input} {output}
        """
 