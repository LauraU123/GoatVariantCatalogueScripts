
rule all:
    input:
        "deepvariant/pca_dv.eigenval",
        #"SWISS/pca_dv_swiss.eigenval"


rule process:
    input:
        "../../aprilVCFs/output/SNV/filteredDPGQ_no_SVdv.g.vcf.gz"
    output:
        "filtereddv.vcf.gz"
    resources:
        mem="80G",
        time="02:00:00",
        cpus=10
    threads: 30
    shell:
        """
        module load BCFtools
        zgrep '^#' {input} > header.vcf
        #plink has a limit of 16000 length of variant name, thus filtering manually here.
        zgrep -v '^#' {input} | awk 'length($3) <= 16000' > filtered_variants.vcf
        cat header.vcf filtered_variants.vcf > filtereddv.vcf
        bgzip filtereddv.vcf
        tabix -p vcf filtereddv.vcf.gz
       """


rule to_plink:
    input:
        "filtereddv.vcf.gz"
    output:
        "deepvariant/plinkformat_dv.bed"
    params:
        "deepvariant/plinkformat_dv"
    resources:
        mem="150G",
        time="02:00:00",
        cpus=10
    threads: 30
    shell:
        """
        module load PLINK
        apptainer exec plink2.sif plink2 --vcf {input} --vcf-half-call missing --chr-set 29 --make-bed -out {params}  --lax-chrx-import --set-missing-var-ids @_#\$r\$a --threads {threads} --geno 0.1 --maf 0.05
        """


rule make_PCA:
    input:
        "deepvariant/plinkformat_dv.bed"
    output:
        #name = "SWISS/pca_dv.eigenval",
        name = "deepvariant/pca_dv.eigenval"
    resources:
        mem="680G",
        time="02:00:00",
        cpus=30
    threads: 30
    params:
        name = "deepvariant/pca_dv",
        #name = "SWISS/pca_dv_swiss",
        in_ = "deepvariant/plinkformat_dv"
    shell:
        """
        module load PLINK
        plink -bfile {params.in_}  --chr-set 29 -pca --make-bed -out {params.name} --threads 30 #--keep swiss.txt
        """
