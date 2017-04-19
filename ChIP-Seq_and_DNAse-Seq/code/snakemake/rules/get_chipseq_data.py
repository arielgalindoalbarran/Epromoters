rule get_chipseq_data:
    """
    Get ChIP-seq data for better control selection
    """
    output: h3k27ac_k562='DATA/CHIP-SEQ/BW/ChIP-seq_H3K27ac_K562.bigWig', \
            h3k4me3_k562='DATA/CHIP-SEQ/BW/ChIP-seq_H3K4me3_K562.bigWig', \
            rad21_k562='DATA/CHIP-SEQ/BW/ChIP-seq_RAD21_K562.bigWig', \
            polr2a_k562='DATA/CHIP-SEQ/BW/ChIP-seq_POLR2A_K562.bigWig'
    params: ppn="nodes=1:ppn=1"
    shell:"""
    wget https://www.encodeproject.org/files/ENCFF000BWY/@@download/ENCFF000BWY.bigWig
    mv ENCFF000BWY.bigWig {output.h3k27ac_k562}
    
    wget https://www.encodeproject.org/files/ENCFF000VDQ/@@download/ENCFF000VDQ.bigWig
    mv ENCFF000VDQ.bigWig {output.h3k4me3_k562}
    
    wget https://www.encodeproject.org/files/ENCFF000YXZ/@@download/ENCFF000YXZ.bigWig
    mv ENCFF000YXZ.bigWig {output.rad21_k562}

    wget https://www.encodeproject.org/files/ENCFF000YWS/@@download/ENCFF000YWS.bigWig
    mv ENCFF000YWS.bigWig {output.polr2a_k562}
    """
