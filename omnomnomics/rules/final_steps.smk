# Final Rule

## Omnomnomics Snake Rule  ##
rule final_steps:
   output:
       "pipeline_completed.txt"
   threads:
       127
   resources:
       mem_mb = 4000*127
   run:
       import os
       import glob
       #---------------------------------------------------------------------------------------------------------------
       # Run multiqc to gather all stats
       #---------------------------------------------------------------------------------------------------------------
       if config['NO_MULTIQC'] == "0":
           log_it(logfile, "Running multiQC...", "STATS")
           subprocess.run(
           f"micromamba activate multiqc && "
           f"multiqc --filename omnomnomics.run.{run_date}.multiqc_report.html --dirs --export . && "
           f"micromamba deactivate",
           shell=True,
           check=True
           )
      
       #---------------------------------------------------------------------------------------------------------------
       # Clean up: Compress and package the tag directories if needed ####################################can delete this I think
       #---------------------------------------------------------------------------------------------------------------
       if os.path.isdir("HOMER_tagDirs"):
           log_it(logfile, "Compressing tagDirs...", "CLEANUP")
           tag_dirs = glob.glob("HOMER_tagDirs/*tagDir")
          
           for tag_dir in tag_dirs:
               log_it(logfile, f"Making tar ball of {tag_dir}...")
               tar_output = f"HOMER_tagDirs/{os.path.basename(tag_dir)}.tar.gz"
               shell(f"tar czf {tar_output} {tag_dir}")
          
           log_it(logfile, "Waiting for tar to complete...")
          
           for tag_dir in tag_dirs:
               log_it(logfile, f"Deleting uncompressed tagDir {tag_dir}...")
               shell(f"rm -r {tag_dir}")


           log_it(logfile, "Cleanup complete")

        #delete the step 11 random final from rule_all if snakemake allows that. else just leave it
      
       shell("echo 'Pipeline Execution Complete!' > pipeline_completed.txt")
       log_it(logfile, "All done!" "FINAL REMARKS")
       log_it(logfile, "Good luck with your downstream analyses!")