#!/bin/sh
# properties = {"type": "single", "rule": "annotate", "local": false, "input": ["/work/FAC/FBM/DMF/smitri/evomicrocomm/genome_size/data/ancestors/oa/assembly.contigs.polypolish.fasta"], "output": ["/work/FAC/FBM/DMF/smitri/evomicrocomm/genome_size/data/ancestors/oa/bakta/assembly.contigs.polypolish.gbff"], "wildcards": {"strain": "oa"}, "params": {"out_dir": "/work/FAC/FBM/DMF/smitri/evomicrocomm/genome_size/data/ancestors/oa/bakta", "genus": "Ochrobactrum", "species": "anthropi"}, "log": [], "threads": 1, "resources": {"tmpdir": "/tmp/24489995"}, "jobid": 0, "cluster": {"mem": "64G", "time": "02:0:0"}}
 cd /users/eulrich/black_queen_hypothesis/scripts/workflows/ancestral && \
/work/FAC/FBM/DMF/smitri/evomicrocomm/eric/miniconda3/envs/bqh/bin/python3.7 \
-m snakemake /work/FAC/FBM/DMF/smitri/evomicrocomm/genome_size/data/ancestors/oa/bakta/assembly.contigs.polypolish.gbff --snakefile /users/eulrich/black_queen_hypothesis/scripts/workflows/ancestral/Snakefile \
--force --cores all --keep-target-files --keep-remote --max-inventory-time 0 \
--wait-for-files /users/eulrich/black_queen_hypothesis/scripts/workflows/ancestral/.snakemake/tmp.bzbe18by /work/FAC/FBM/DMF/smitri/evomicrocomm/genome_size/data/ancestors/oa/assembly.contigs.polypolish.fasta --latency-wait 5 \
 --attempt 1 --force-use-threads --scheduler ilp \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules annotate --nocolor --notemp --no-hooks --nolock --scheduler-solver-path /work/FAC/FBM/DMF/smitri/evomicrocomm/eric/miniconda3/envs/bqh/bin \
--mode 2  --default-resources "tmpdir=system_tmpdir"  && touch /users/eulrich/black_queen_hypothesis/scripts/workflows/ancestral/.snakemake/tmp.bzbe18by/0.jobfinished || (touch /users/eulrich/black_queen_hypothesis/scripts/workflows/ancestral/.snakemake/tmp.bzbe18by/0.jobfailed; exit 1)

