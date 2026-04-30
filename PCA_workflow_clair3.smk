
rule all:
    input:
        "clair/pca_clair.eigenval",

rule to_plink:
    input:
        "../../aprilVCFs/output/SNV/filteredDPGQ_no_SVclair.g.vcf.gz"
    output:
        "clair/plinkformat_clair.bed"
    params:
        "clair/plinkformat_clair"
    resources:
        mem="80G",
        time="02:00:00",
        cpus=10
    threads: 30
    shell:
        """
        module load PLINK
        plink --vcf {input} --vcf-half-call missing --chr-set 29  --make-bed --geno 0.1 --maf 0.05 -out {params} --threads {threads}
        """


rule make_PCA_swiss:
    input:
        rules.to_plink.output
    output:
        name = "clair/pca_clair.eigenval"
    resources:
        mem="60G",
        time="02:00:00",
        cpus=10
    threads: 30
    params:
        name = "clair/pca_clair",
        in_ = rules.to_plink.params
    shell:
        """
        module load PLINK
        plink -bfile {params.in_}  --chr-set 29 -pca --make-bed -out {params.name}
        """
