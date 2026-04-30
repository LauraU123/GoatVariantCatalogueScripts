CALLERS = [ "sawfish2", "sniffles2"]

rule all:
    input:
        expand("{caller}/pca_plink_{caller}.eigenval",caller = CALLERS)

rule to_plink:
    input:
        "../../aprilVCFs/output/SV/SV_{caller}_PASS.vcf.gz "
    output:
        "{caller}/plinkformat_{caller}.bed"
    params:
        "{caller}/plinkformat_{caller}"
    resources:
        mem="80G",
        time="02:00:00",
        cpus=10
    threads: 30
    shell:
        """
        module load PLINK
        plink --vcf {input} --vcf-half-call missing --chr-set 29  --make-bed --out {params}
        """


rule make_PCA:
    input:
        rules.to_plink.output
    output:
        name = "{caller}/pca_plink_{caller}.eigenval"
    resources:
        mem="80G",
        time="02:00:00",
        cpus=10
    threads: 30
    params:
        name = "{caller}/pca_plink_{caller}",
        in_ = "{caller}/plinkformat_{caller}"
    shell:
        """
        module load PLINK
        plink -bfile {params.in_} --geno 0.1  --maf 0.05 --chr-set 29 -pca --make-bed --out {params.name}
        """
