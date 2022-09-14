configfile: "../config/config.yaml"


include: "dagfile.smk"

rule all:
    input:
        "../results/dag.png"
