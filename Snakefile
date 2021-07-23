import os
from glob import glob
configfile: "config.yml"

rule all:
	input:
		expand("results/{sample}/plots/bam/", sample=config['SAMPLES']),
		expand("results/{sample}/plots/vcf_sv/", sample=config['SAMPLES']),
		expand("results/{sample}/plots/vcf_snp/", sample=config['SAMPLES'])

rule index_bam:
	input: 
		lambda wildcards: glob("%s/%s/*.bam"%(config['SAMPLE_DIR'], wildcards.sample))
	output: 
		dynamic("data/samples/{sample}/{FILE}.xml")
	shell:
		"""
		IFS=' '
		read -ra bamfiles <<< "{input}"
		for f in "${{bamfiles[@]}}"; do 
			name=$(basename $f .bam)
			ln -s $f data/samples/{wildcards.sample}/
			dataset create --generateIndices $PWD/data/samples/{wildcards.sample}/$name.xml $f
		done
		"""		

rule map_read:
	input:
		dynamic("data/samples/{sample}/{FILE}.xml")
	output:
		"results/{sample}/mapped/{sample}.bam"
	params:
		config['REFERENCE']
	log:
		"logs/{sample}/mapping.log"
	shell:
		"""
		echo "Started by " $USER " at" > {log}
		echo "smrtlink pbmm2 version " $(pbmm2 --version) >> {log}
		ls $PWD/{input}/../*.subreads.bam > tmp/bamfiles.fofn
		pbmm2 align {params} tmp/bamfiles.fofn results/{wildcards.sample}/mapped/{wildcards.sample}.bam --sort
		echo "pbmm2 align" {params} "tmp/bamfiles.fofn results/"{wildcards.sample}"/mapped/"{wildcards.sample}".bam --sort" >> {log}
		echo "Ended at "$(date) >> {log}
		"""

rule call_snp:
	input:
		genome=config['REFERENCE'],
		bam="results/{sample}/mapped/{sample}.bam"
	output:
		"results/{sample}/called/{sample}_snp.vcf"
	log:
		"logs/{sample}/snp_calling.log"
	shell:
		"""
		echo "Started by " $USER " at" > {log}
		echo "smrtlink gcpp version" $(gcpp --version) >> {log}
		gcpp -j$(nproc) --algorithm=arrow -r {input.genome} -o {output} {input.bam}
		echo "gcpp -j "$(nproc)"--algorithm=arrow -r "{input.genome}" -o "{output} {input.bam}" >> {log}
		echo "Ended at "$(date) >> {log}
		"""
		
rule call_sv:
	input:
		genome=config['REFERENCE'],
		bam="results/{sample}/mapped/{sample}.bam"
	output:
		svsig="tmp/{sample}.svsig.gz",
		vcf="results/{sample}/called/{sample}_sv.vcf"
	log:
		"logs/{sample}/sv_calling.log"
	shell:
		"""
		echo "Started by " $USER " at" > {log}
		echo "smrtlink pbsv version" $(pbsv --version) >> {log}
		pbsv discover {input.bam} {output.svsig}
		pbsv call {input.genome} {output.svsig} {output.vcf}
		echo "pbsv discover" {input.bam} {output.svsig} >> {log}
		echo "pbsv call" {input.genome} {output.svsig} {output.vcf} >> {log}
		echo "Ended at "$(date) >> {log}
		"""

rule bam_stats:
	input:
		"results/{sample}/mapped/{sample}.bam"
	output:
		"results/{sample}/plots/{sample}.bamstats"
	conda: 
		"envs/stats.yml"
	shell:
		"""
		samtools stats {input} > {output}
		"""
rule plot_bam_stats:
	input:
		"results/{sample}/plots/{sample}.bamstats"
	output:
		directory("results/{sample}/plots/bam/")
	conda: 
		"envs/stats.yml"
	shell:
		"""
		plot-bamstats -p {output} {input}
		"""

rule snp_stats:
	input:
		"results/{sample}/called/{sample}_snp.vcf"
	output:
		"results/{sample}/plots/{sample}.snpstats"
	conda: 
		"envs/stats.yml"
	shell:
		"""
		bcftools stats {input} > {output}
		"""

rule plot_snp_stats:
	input:
		"results/{sample}/plots/{sample}.snpstats"
	output:
		directory("results/{sample}/plots/vcf_snp/")
	conda: 
		"envs/stats.yml"
	shell:
		"""
		plot-vcfstats -P -p {output} {input}
		"""
		
rule sv_stats:
	input:
		"results/{sample}/called/{sample}_sv.vcf"
	output:
		"results/{sample}/plots/{sample}.svstats"
	conda: 
		"envs/stats.yml"
	shell:
		"""
		bcftools stats {input} > {output}
		"""

rule plot_sv_stats:
	input:
		"results/{sample}/plots/{sample}.svstats"
	output:
		directory("results/{sample}/plots/vcf_sv/")
	conda: 
		"envs/stats.yml"
	shell:
		"""
		plot-vcfstats -P -p {output} {input}
		"""