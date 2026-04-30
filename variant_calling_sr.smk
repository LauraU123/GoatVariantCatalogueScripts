
SCRATCH_DIR = "$TMPDIR"

rule all:
    input:
        "genotyped_shortread_goat_raw.vcf.gz"


INTERVALS = sorted(glob_wildcards("intervals/{chunk}.interval_list").chunk)
#-M {output.metrics}
rule call_haplo:
    input:
        bam = inputdir + "{name}.md.bam",
        reference = "../../../refs/GCA_040806595.1_T2T-goat1.0_genomic_nice.fa"
    output:
        vcf = "gatkoutput/{name}.vcf.gz"
    resources:
        mem="800G",
        time="65:00:00",
        cpus=30
    threads: 30
    shell:
        """
        module load GATK SAMtools 
        gatk --java-options "-Xmx16G" HaplotypeCaller -R {input.reference}  -I {input.bam} -O {output.vcf} -ERC GVCF --native-pair-hmm-threads {threads}
        """


rule genomicsdb:
    input:
        ref = rules.call_haplo.input.reference,
        vcf = expand(rules.call_haplo.output, name=NAMES),
        interval = "intervals/{chunk}.interval_list"
    output:
        directory("genomicsdb/{chunk}")
        
    params:
        samples =lambda wc: " ".join(f"-V {f}" for f in expand("gatkoutput/{name}.vcf.gz",name=NAMES))
 
    resources:
        mem="250G",
        time="35:00:00",
        cpus=30
    shell:
        """
        #export TMPDIR ="$TMPDIR"
        module load GATK
        gatk --java-options "-Xmx4g -Xms4g" GenomicsDBImport {params.samples} --genomicsdb-workspace-path {output} -L {input.interval}
        """

rule genotype:
    input:
        ref = rules.call_haplo.input.reference,
        db = "genomicsdb/{chunk}",
        interval = "intervals/{chunk}.interval_list"
    output:
        "genotyped/{chunk}.vcf.gz"
    resources:
        mem="250G",
        time="35:00:00",
        cpus=30
    shell:
        """
        module load GATK
        gatk GenotypeGVCFs -R {input.ref} -V gendb://{input.db} -L {input.interval} -O {output}
        """

rule gather:
    input:
        expand("genotyped/{chunk}.vcf.gz", chunk=INTERVALS)
    output:
        "genotyped_shortread_goat_raw.vcf.gz"
    params:
        samples = lambda wc: " ".join(f"-I {f}" for f in expand("genotyped/{chunk}.vcf.gz",chunk=INTERVALS))
    resources:
        mem="250G",
        time="35:00:00",
        cpus=30
    shell:
        """
        module load GATK
        gatk GatherVcfs {params.samples} -O {output}
        """
