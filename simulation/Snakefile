shell.executable("/bin/bash")
shell.prefix("source ~/.bashrc;")

configfile: "config.json"

rule all:
    input:
        expand("csv/n{n}_cnt{cnt}.csv", n=config["n"], cnt=config["cnt"])

rule smurf:
    output: "csv/n{n}_cnt{cnt}.csv"
    benchmark: "bench/n{n}_cnt{cnt}.bench"
    shell:
        "R CMD BATCH --no-save --no-restore '--args {wildcards.n} {wildcards.cnt} {output}' " 
        " smurf.R log/n{wildcards.n}_cnt{wildcards.cnt}.Rout"
