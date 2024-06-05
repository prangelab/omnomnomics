

## Omnomnomics Snake Rule  ##
rule rule3:
    input:
        "trimmed_FASTQ/rule1_output.{sample}.fastq.gz"
    output:
        "BAM/rule3_output.{sample}.bam"
    threads:
        6
    resources:
        mem_mb = 4000
    shell:
        """
        echo "Rule 3" >> {input}
        cp {input} {output}
        """