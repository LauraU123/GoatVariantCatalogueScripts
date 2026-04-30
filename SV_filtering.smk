CALLERS = ["sniffles2", "sawfish2"]

rule all:
    input:
        expand("output/SV/SV_{caller}_PASS.vcf.gz", caller=CALLERS)


rule remove_SNVs:
    message:
        """Since SVs are defined as longer than 50 bp, remove all that are shorter than 50bp"""
    input:
        "raw_{caller}_annotated_2apz.vcf.gz"
    output:
        "output/SV/SV_{caller}_no_SNVS.vcf.gz"
    resources:
        mem="40G",
        time="10:00:00",
        cpus=30
    shell:
        """
        module load BCFtools
        bcftools view -i '(ABS(INFO/SVLEN) >=50) || (INFO/SVTYPE="BND")' {input} -Oz -o {output}
        """


rule filter_quality_20:
    message:
        """To remove outliers with very low score"""
    input:
        rules.remove_SNVs.output
    output:
        "output/SV/SV_{caller}_no_SNVS_qual20plus.vcf.gz"
    resources:
        mem="40G",
        time="10:00:00",
        cpus=30
    shell:
        """
        module load BCFtools
        bcftools view -i 'QUAL>=20' {input} -Oz -o {output}
        bcftools index {output}
        """


rule PASS_only:
    message:
        """Keeping only variants that PASS"""
    input:
        rules.filter_quality_20.output
    output:
        "output/SV/SV_{caller}_PASS.vcf.gz"
    resources:
        mem="40G",
        time="10:00:00",
        cpus=30
    shell:
        """
        module load BCFtools
        bcftools view -f PASS {input} -Oz -o {output}
        """
