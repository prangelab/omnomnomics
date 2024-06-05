

## Omnomnomics Snake Rule  ##
import os
import glob

def update_semaphore_file(experiment_dir, rule_name):
    semaphore_file = os.path.join(experiment_dir, "omnomnomics.semaphore")
    with open(semaphore_file, "a") as f:
        f.write(f"{rule_name}\n")


rule rule1:
    input:
        fastq1 = "FASTQ/{sample}.fastq.gz"
        #glob.glob("FASTQ/*.fastq.gz")
    output:
       # expand("trimmed_FASTQ/rule1_output.{sample}.fastq.gz", sample = samples)
       "trimmed_FASTQ/rule1_output.{sample}.fastq.gz"
    params:
        #experiment_dir=config['EXPERIMENT_DIR']
        con = True
    threads:
        max((master_config['mincores_single_sample_step1_9'][master_config['trim_rule_num'] - 1]
         if 'mincores_single_sample_step1_9' in master_config
         and isinstance(master_config['mincores_single_sample_step1_9'], list)
         and len(master_config['mincores_single_sample_step1_9']) > (master_config['trim_rule_num'] - 1)
         and master_config['mincores_single_sample_step1_9'][master_config['trim_rule_num'] - 1] is not None
         and isinstance(master_config['mincores_single_sample_step1_9'][master_config['trim_rule_num'] - 1], int)
         and master_config['mincores_single_sample_step1_9'][master_config['trim_rule_num'] - 1] > master_config['min_slice_cores']
         else master_config['min_slice_cores']),
         min((master_config['maxcores_single_sample_step1_9'][master_config['trim_rule_num'] - 1]
         if 'maxcores_single_sample_step1_9' in master_config
         and isinstance(master_config['maxcores_single_sample_step1_9'], list)
         and len(master_config['maxcores_single_sample_step1_9']) > (master_config['trim_rule_num'] - 1)
         and master_config['maxcores_single_sample_step1_9'][master_config['trim_rule_num'] - 1] is not None
         and isinstance(master_config['maxcores_single_sample_step1_9'][master_config['trim_rule_num'] - 1], int)
         and master_config['maxcores_single_sample_step1_9'][master_config['trim_rule_num'] - 1] < master_config['cores_per_node']
         else master_config['cores_per_node']), ((master_config['nodes_in_partition']*master_config['cores_per_node'])/num_samples)) )
    resources:
        mem_mb = (master_config['min_mem_mb'][master_config['trim_rule_num'] - 1]
            if 'min_mem_mb' in master_config
            and isinstance(master_config['min_mem_mb'], list)
            and len(master_config['min_mem_mb']) > (master_config['trim_rule_num'] - 1)
            and master_config['min_mem_mb'][master_config['trim_rule_num'] - 1] is not None
            and isinstance(master_config['min_mem_mb'][master_config['trim_rule_num'] - 1], int)
            and master_config['min_mem_mb'][master_config['trim_rule_num'] - 1] > master_config['min_slice_mem']
            else (master_config['min_slice_mem'] 
                if('min_mem_mb' in master_config
                    and isinstance(master_config['min_mem_mb'], list)
                    and len(master_config['min_mem_mb']) > (master_config['trim_rule_num'] - 1)
                    and master_config['min_mem_mb'][master_config['trim_rule_num'] - 1] is not None
                    and isinstance(master_config['min_mem_mb'][master_config['trim_rule_num'] - 1], int)
                    and master_config['min_mem_mb'][master_config['trim_rule_num'] - 1] <= master_config['min_slice_mem'])
                else (master_config['mincores_single_sample_step1_9'][master_config['trim_rule_num'] - 1]
                if 'mincores_single_sample_step1_9' in master_config
                and isinstance(master_config['mincores_single_sample_step1_9'], list)
                and len(master_config['mincores_single_sample_step1_9']) > (master_config['trim_rule_num'] - 1)
                and master_config['mincores_single_sample_step1_9'][master_config['trim_rule_num'] - 1] is not None
                and isinstance(master_config['mincores_single_sample_step1_9'][master_config['trim_rule_num'] - 1], int)
                and master_config['mincores_single_sample_step1_9'][master_config['trim_rule_num'] - 1] > master_config['min_slice_cores']
                else master_config['min_slice_cores'])* master_config['max_mem_per_core_mb']))


    # onsuccess:s
    #     update_semaphore_file(snakemake.params.experiment_dir,1)
    run:
        # counter = 0
        # for file, file2 in zip(input,output):
        #     print(file)
        #     print(file2)
        #     def choice(con, file):
        #         if con: 
        #             print(f"{config['THEGENOME']}")
        #             shell(
        #                 f"""
        #                 echo "Rule 1" >> {file}
        #                 """)
        #             shell(
        #                 f"""cp {file} {file2} """) 
        #         else:
        #             print("HERE")
        #     choice(params.con, file)
        # log_it(logfile, "BLABLABLA")


        def choice(con, file):
            semaphore_path = os.path.join(experiment_dir, "omnomnomics.semaphore")
            if con: 
                print(f"{config['THEGENOME']}")
                shell(
                    f"""
                    echo "Rule 1" >> {file}
                    """)
                shell(
                    f"""cp {file} {output} """) 
                #shell(f"""echo "1" > {omnomnomics.semaphore}""")
                #print(max(master_config['maxcores_single_sample_step1_9'][master_config['trim_rule_num']-1], min(master_config['maxcores_single_sample_step1_9'][master_config['trim_rule_num']-1],((master_config['max_nodes']*master_config['max_cores_per_node'])/1))))
                print(threads)
                print(resources.mem_mb)
            else:
                print("HERE")
        choice(params.con, input.fastq1)

