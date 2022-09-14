rule dag_file:
    output:
        '../results/dag.png'
    message:
        "Creates a dag file showing what happened."
    shell:
         "snakemake --snakefile main.smk --forceall --dag | dot -Tpng > ../results/dag.png"
