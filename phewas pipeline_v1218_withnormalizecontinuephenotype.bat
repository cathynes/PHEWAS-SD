###new phewas SD with normalize continue phenotypes 

cd $OAK/catherine/phewas_SD_1218/settings

# STEP0: Check that the setup seems OK.
sbatch -p normal,owners,jpriest --qos=normal -N 1 -n 1 -t 00:40:00 -J step0.phewas -o step0.log -e step0.log --open-mode=append /oak/stanford/groups/jpriest/catherine/phwas.script/phewas_step0_check_setup.sh

# ... wait for 20-30 minutes

squeue -u tcheandj

# STEP1: Create the input files with the phenotypes used in the analysis.
sbatch -p normal,owners,jpriest --qos=normal -N 1 -n 1 -t 00:40:00 -J step1.phewas -o step1.log -e step1.log --open-mode=append /oak/stanford/groups/jpriest/catherine/phwas.script/phewas_step1_prepare_phenotype_files.sh

# ... wait for 20-30 minutes



# STEP2: Extract the genotype data for the selected loci.
sbatch -p normal,owners,jpriest --qos=normal -N 1 -n 1 -t 01:14:00 --mem=8000 -J step2.phewas -o step2.log -e step2.log --open-mode=append /oak/stanford/groups/jpriest/catherine/phwas.script/phewas_step2_extract_variants.sh

# ... wait for 1-2 hours
# You will have to re-run this command once the queue is empty to pick up failed jobs.



# STEP3: Run the SKAT tests.
sbatch -p normal,owners,jpriest --qos=normal -N 1 -n 1 -t 01:00:00  --mem=30000 -J step3.phewas -o step3.log -e step3.log --open-mode=append /oak/stanford/groups/jpriest/catherine/phwas.script/phewas_step3_run_multi_variant_test.sh

# ... wait for 1-2 days
# You will have to re-run this command once the queue is empty to pick up failed jobs.


# STEP4: Process the SKAT results.
sbatch -p normal,owners,jpriest --qos=normal -N 1 -n 1 -t 01:00:00 --mem=30000 -J step4.phewas -o step4.log -e step4.log --open-mode=append /oak/stanford/groups/jpriest/catherine/phwas.script/phewas_step4_process_multi_variant_results.sh

# ... wait for a few minutes



# STEP5: Run the single variant tests.
sbatch -p normal,owners,jpriest --qos=normal -N 1 -n 1 -t 01:00:00 -J step5.phewas -o step5.log -e step5.log --open-mode=append /oak/stanford/groups/jpriest/catherine/phwas.script/phewas_step5_run_single_variant_test.sh

# ... wait for a few hours
# You will have to re-run this command once the queue is empty to pick up failed jobs.



# STEP6: Process the single variant results.
sbatch -p normal,owners,jpriest --qos=normal -N 1 -n 1 -t 01:00:00 -J step6.phewas -o step6.log -e step6.log --open-mode=append /oak/stanford/groups/jpriest/catherine/phwas.script/phewas_step6_process_single_variant_results.sh

# ... wait for a few minutes



# STEP7: Create HTML reports with a summary of the analyses.
sbatch -p normal,owners,jpriest --qos=normal -N 1 -n 1 -t 01:00:00 -J step7.phewas -o step7.log -e step7.log --open-mode=append /oak/stanford/groups/jpriest/catherine/phwas.script/phewas_step7_build_report.sh

# ... wait for a few minutes



# STEP8: Finalize and archive the results and logs.
sbatch -p normal,owners,jpriest --qos=normal -N 1 -n 1 -t 01:00:00 -J step8.phewas -o step8.log -e step8.log --open-mode=append /oak/stanford/groups/jpriest/catherine/phwas.script/phewas_step8_finalize.sh

# ... wait for a few minutes



# >>> IMPORTANT: If the final result file is stored in your SCRATCH, please note
# that this is a temporary storage without backup. SCRATCH can at any time be
# cleared of old data. Do not use it for long-time storage. Use OAK instead.

# RS-numbers?
# In case you want to add on rs# to the files, map files are available in the
# folder with genotype data: genotypes/pgen_app1372_hrc/*.incl-snps.txt

### cancel all jobs except one 21657615
scancel `squeue -u tcheandj| grep -v 21657615|tail -n+2 |awk '{print $1}'`

## snp id location ####
/oak/stanford/projects/ukbb/genotypes/pgen_app1372_hrc/ukb_imp_chr10_v2.hrc1.1_sites.incl-snps.txt



