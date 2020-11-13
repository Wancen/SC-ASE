#!/bin/sh
# properties = {"type": "single", "rule": "smurf", "local": false, "input": [], "output": ["csv/n10_cnt100.csv"], "wildcards": {"n": "10", "cnt": "100"}, "params": {}, "log": [], "threads": 1, "resources": {}, "jobid": 2, "cluster": {}}
 cd /proj/milovelab/mu/SC-ASE/simulation && \
/nas/longleaf/apps/python/3.6.6/bin/python3.6 \
-m snakemake csv/n10_cnt100.csv --snakefile /proj/milovelab/mu/SC-ASE/simulation/Snakefile \
--force -j --keep-target-files --keep-remote \
--wait-for-files /proj/milovelab/mu/SC-ASE/simulation/.snakemake/tmp.sdgvc6gt --latency-wait 180 \
 --attempt 1 --force-use-threads \
--wrapper-prefix https://github.com/snakemake/snakemake-wrappers/raw/ \
   --allowed-rules smurf --nocolor --notemp --no-hooks --nolock \
--mode 2  && touch "/proj/milovelab/mu/SC-ASE/simulation/.snakemake/tmp.sdgvc6gt/2.jobfinished" || (touch "/proj/milovelab/mu/SC-ASE/simulation/.snakemake/tmp.sdgvc6gt/2.jobfailed"; exit 1)

