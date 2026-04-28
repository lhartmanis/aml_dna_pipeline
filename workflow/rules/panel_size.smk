rule panel_size:
    input:
        bed=config["reference"]["interval_bed"]
    output:
        "results/panel_metrics/panel_size.txt"
    log:
        "logs/snakemake/panel_size.log"
    shell:
        r"""
        mkdir -p results/panel_metrics logs/snakemake
        awk 'BEGIN{{bp=0}} !/^#/ && NF>=3 {{bp += ($3 - $2)}} END {{printf "panel_bp\t%d\npanel_mb\t%.6f\n", bp, bp/1000000}}' {input.bed} > {output} 2> {log}
        """
