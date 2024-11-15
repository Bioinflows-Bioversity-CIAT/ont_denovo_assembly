
rule make_blast_db:
	input:
		curated_draft = f"{base_dir}/purge/purge_haplotigs/{{sample}}/{{sample}}.fasta"
	output:
		db = multiext( f"{base_dir}/scaffolding/{{sample}}/blastdb/{{sample}}",
            ".ndb",
            ".nhr",
            ".nin",
            ".not",
            ".nsq",
            ".ntf",
            ".nto"
        )
	conda:
		"tools2"
	shell:
		"""
		makeblastdb -in {input.curated_draft} \
		-dbtype nucl \
		-input_type fasta -blastdb_version 5 -parse_seqids \
		-out {base_dir}/scaffolding/{wildcards.sample}/blastdb/{wildcards.sample}
		"""

rule map_markers_to_assembly:
	input:
		markers_sequence = config['linkage_map']['markers'],
		db = multiext( f"{base_dir}/scaffolding/{{sample}}/blastdb/{{sample}}",
            ".ndb",
            ".nhr",
            ".nin",
            ".not",
            ".nsq",
            ".ntf",
            ".nto"
        ),
	output:
		f"{base_dir}/scaffolding/{{sample}}/linkage_map/{{sample}}_markers.blast"
	conda:
		"tools2"
	threads:
		10
	log:
		'logs/{sample}/scaffolding/{sample}_mapmarkers.log'
	shell:
		"""
		blastn -db {base_dir}/scaffolding/{wildcards.sample}/blastdb/{wildcards.sample} \
		-query {input.markers_sequence} \
		-evalue 1e-50 \
		-max_target_seqs 1\
		-outfmt 6 -out {output} 2> {log}
		"""

rule get_genetic_map_table:
	input:
		genetic_map = config['linkage_map']['map'],
		blast = f"{base_dir}/scaffolding/{{sample}}/linkage_map/{{sample}}_markers.blast"
	output:
		f"{base_dir}/scaffolding/{{sample}}/linkage_map/{{sample}}_map.bed"
	params:
		script = "/home/scruz/workflows/denovo_assembly/workflow/scripts/get_genetic_map.R"
	conda:
		"R_base"
	shell:
		"""
		Rscript {params.script} -map {input.genetic_map} \
		-blast {input.blast} -out {output}
		"""

rule scaffolding_job_allmaps:
	input:
		g_map = f"{base_dir}/scaffolding/{{sample}}/linkage_map/{{sample}}_map.bed",
		curated_draft = f"{base_dir}/purge/purge_haplotigs/{{sample}}/{{sample}}.fasta"
	output:
		multiext(f"{base_dir}/scaffolding/{{sample}}/linkage_map/{{sample}}",".fasta", ".agp", ".chain")
	conda:
		"../envs/scaffolding.yaml"
	log:
		'logs/{sample}/scaffolding/{sample}_allmaps.log'
	shell:
		"""
		python -m jcvi.assembly.allmaps path {input.g_map} {input.curated_draft}
		"""