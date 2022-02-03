rule all:
    input:
        "SRR2584857_1.ecoli-rel606.sens.vcf", # sensitive, no filtering
        "SRR2584857_1.ecoli-rel606.spec.vcf", # specific, high-q filtering
        "SRR25844O3_1.ecoli-rel606.sens.vcf",
	"SRR25844O3_1.ecoli-rel606.spec.vcf",
	"SRR25844O4_1.ecoli-rel606.sens.vcf",
	"SRR25844O4_1.ecoli-rel606.spec.vcf",
	"SRR25844O5_1.ecoli-rel606.sens.vcf",
	"SRR25844O5_1.ecoli-rel606.spec.vcf",

rule download_data:
    conda: "env-wget.yml"
    output: "SRR2584857_1.fastq.gz"
    shell: """
        wget https://osf.io/4rdza/download -O {output}
    """

rule download_data_2:
    conda: "env-wget.yml"
    output: "SRR25844O3_1.fastq.gz"
    shell: """
	wget https://osf.io/6jx7z/download -O {output}
    """

rule download_data_3:
    conda: "env-wget.yml"
    output: "SRR25844O4_1.fastq.gz"
    shell: """
	wget https://osf.io/s24ky/download -O {output}
    """

rule download_data_4:
    conda: "env-wget.yml"
    output: "SRR25844O5_1.fastq.gz"
    shell: """
	wget https://osf.io/7qek6/download -O {output}
    """


rule download_genome:
    conda: "env-wget.yml"
    output: "ecoli-rel606.fa.gz"
    shell:
        "wget https://osf.io/8sm92/download -O {output}"

rule map_reads:
    conda: "env-minimap.yml"
    input: ref="ecoli-rel606.fa.gz", reads="{sample}.fastq.gz"
    output: "{sample}.ecoli-rel606.sam"
    shell: """
        minimap2 -ax sr {input.ref} {input.reads} > {output}
    """

rule sam_to_bam:
    conda: "env-minimap.yml"
    input: "{sample}.ecoli-rel606.sam"
    output: "{sample}.ecoli-rel606.bam"
    shell: """
        samtools view -b -F 4 {input} > {output}
     """

rule sort_bam:
    conda: "env-minimap.yml"
    input: "{sample}.ecoli-rel606.bam"
    output: "{sample}.ecoli-rel606.bam.sorted"
    shell: """
        samtools sort {input} > {output}
    """

rule gunzip_fa:
    input:
        "ecoli-rel606.fa.gz"
    output:
        "ecoli-rel606.fa"
    shell: """  
        gunzip -c {input} > {output}
    """

rule call_variants:
    conda: "env-bcftools.yml"
    input:
        ref="ecoli-rel606.fa",
        bamsort="{sample}.ecoli-rel606.bam.sorted"
    output:
        pileup="{sample}.ecoli-rel606.pileup",
        bcf="{sample}.ecoli-rel606.bcf",
        vcf="{sample}.ecoli-rel606.sens.vcf", # sensitive VCF
    shell: """
        bcftools mpileup -Ou -f {input.ref} {input.bamsort} > {output.pileup}
        bcftools call -mv -Ob {output.pileup} -o {output.bcf}
        bcftools view {output.bcf} > {output.vcf}
    """

rule call_variants_spec:
    conda: "env-bcftools.yml"
    input:
        "{sample}.ecoli-rel606.sens.vcf"
    output:
        "{sample}.ecoli-rel606.spec.vcf"
    shell: """
        bcftools filter -Ov {input} > {output}
    """
