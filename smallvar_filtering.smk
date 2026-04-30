CALLERS = ["dv", "clair"]


rule all:
    input:
        expand("output/SNV/filtered_no_SV_{caller}_qual20plus.vcf.gz", caller=CALLERS),

rule normalize:
    message:
        """Left aligning the variant calls, and splitting multiallelic lines to make them more comparable!"""
    input:
        vcf = "merged_output_annotated_allgoat_{caller}_2apz.g.vcf.gz",
        fasta = "../refs/GCA_040806595.1_T2T-goat1.0_genomic_nice.fa"
    output:
        "output/SNV/{caller}_normalized.g.vcf.gz"
    resources:
        mem="30G",
        time="03:00:00",
        cpus=30
    shell:
        """
        # normalizing the vcf files by left aligning them (ref is given for this purpose), and -m to split multiallelics into biallelics
        module load BCFtools
        bcftools norm -f {input.fasta} {input.vcf} -m -any -Oz -o {output}
        """


        
#instead of both to not remove the SNVs here. 

rule remove_SV:
    input:
        rules.normalize.output
    output:
        "output/SNV/no_SV{caller}.g.vcf.gz"
    resources:
        mem="30G",
        time="03:00:00",
        cpus=30
    threads: 30
    shell:
        """
        module load BCFtools
        bcftools view -e 'strlen(REF) >= 50 || (strlen(ALT[0]) >= 50)' {input} -o {output} -Oz 
        bcftools view -i 'strlen(REF) >= 50 || (strlen(ALT[0]) >= 50)' {input} -o removed_{wildcards.caller}.vcf.gz -Oz
        """


rule filtering:
    input:
        rules.remove_SV.output
    output:
        "output/SNV/filteredDPGQ_no_SV{caller}.g.vcf.gz"
    resources:
        mem="40G",
        time="03:00:00",
        cpus=30
    shell:
        """
        module load BCFtools
        bcftools filter -i 'DP>=8 && GQ>=20' {input} -Oz -o {output}
        """


rule filter_quality_20:
    message:
        """To remove outliers with very low score"""
    input:
        rules.filtering.output
    output:
        "output/SNV/filtered_no_SV_{caller}_qual20plus.vcf.gz"
    resources:
        mem="40G",
        time="02:00:00",
        cpus=30
    shell:
        """
        module load BCFtools
        bcftools view -i 'QUAL>=20' {input} -Oz -o {output}
        bcftools index {output}
        """
