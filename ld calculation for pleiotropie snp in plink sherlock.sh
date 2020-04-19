## lD calculatioin in plink

### FBN1
/home/users/tcheandj/phewas/plink2 \
--pfile /oak/stanford/projects/ukbb/genotypes/pgen_app13721_v3/ukb_imp_chr15_v3.mac1 \
--chr 15 \
--from-bp 48696728 \
--to-bp 48941370 \
--force-intersect \
--extract /oak/stanford/groups/jpriest/catherine/data/pleiotropie.snp.phewas.txt \
--keep /scratch/users/tcheandj/subset.ukbb.for.clumping/subset.ind.UKBB.txt \
--make-bed \
--remove-nosex \
--out /oak/stanford/groups/jpriest/catherine/data/plink.file.subset.ukbb/FBN1.pleiotropie.snp

### MAP2K2
/home/users/tcheandj/phewas/plink2 \
--pfile /oak/stanford/projects/ukbb/genotypes/pgen_app13721_v3/ukb_imp_chr15_v3.mac1 \
--chr 15 \
--from-bp 66678170 \
--to-bp 66783265 \
--extract /oak/stanford/groups/jpriest/catherine/data/pleiotropie.snp.phewas.txt \
--keep /scratch/users/tcheandj/subset.ukbb.for.clumping/subset.ind.UKBB.txt \
--force-intersect \
--make-bed \
--remove-nosex \
--out /oak/stanford/groups/jpriest/catherine/data/plink.file.subset.ukbb/MAP2K2.pleiotropie.snp

### RASA2
/home/users/tcheandj/phewas/plink2 \
--pfile /oak/stanford/projects/ukbb/genotypes/pgen_app13721_v3/ukb_imp_chr3_v3.mac1 \
--extract /oak/stanford/groups/jpriest/catherine/data/pleiotropie.snp.phewas.txt  \
--keep /scratch/users/tcheandj/subset.ukbb.for.clumping/subset.ind.UKBB.txt \
--make-bed \
--remove-nosex \
--out /oak/stanford/groups/jpriest/catherine/data/plink.file.subset.ukbb/RASA2.pleiotropie.snp

### PTPN11
/home/users/tcheandj/phewas/plink2 \
--pfile /oak/stanford/projects/ukbb/genotypes/pgen_app13721_v3/ukb_imp_chr12_v3.mac1 \
--extract /oak/stanford/groups/jpriest/catherine/data/pleiotropie.snp.phewas.txt  \
--keep /scratch/users/tcheandj/subset.ukbb.for.clumping/subset.ind.UKBB.txt \
--make-bed \
--remove-nosex \
--out /oak/stanford/groups/jpriest/catherine/data/plink.file.subset.ukbb/PTPN11.pleiotropie.snp


