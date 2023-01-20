# Includes the scripts it will need to run
include: 'total.smk'
configfile: "Pipeline/config/config.yaml"

# Main rule to call upon the entire script.
rule execute_pipeline:
    input:
        # The output of any rule can be placed here.
        expand("Postprocess/finished/{sample}_log.file", sample=config['samples'])

# Creates the dag file ONLY if the process is completed successfully.
onsuccess:
    with open("total_dag.txt","w") as f:
        f.writelines(str(workflow.persistence.dag))
    shell("cat total_dag.txt | dot -Tpng > total_dag.png")