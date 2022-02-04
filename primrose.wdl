version 1.0


workflow PrimroseAlign {
    input {
        Array[File] bamFiles
        File refFasta
        String sample
        String refName
        String dockerImage
    }

    scatter(bam in bamFiles) {
        call Extracthifi {
            input:
                inputBam = bam,
                dockerImage = dockerImage
        }

        call Primrose {
            input:
                inputBam = Extracthifi.outputBam,
                dockerImage = dockerImage
        }
    }

    call Pbmm2 {
        input:
            inputBams = Primrose.outputBam,
            fasta = refFasta,
            refName = refName,
            sample = sample,
            dockerImage = dockerImage
    }

    call BamSort {
        input:
            inputBam = Pbmm2.outputBam,
            dockerImage= dockerImage
    }

    call BamIndex {
        input:
            inputBam = BamSort.outputBam,
            dockerImage= dockerImage
    }

    output {
        Array[File] hifi = Extracthifi.outputBam
        Array[File] cpg = Primrose.outputBam
        # File map = Pbmm2.outputBam
        File sort = BamSort.outputBam
        File bai = BamIndex.outputBai
    }
}


task Extracthifi {
    input {
        File inputBam
        String dockerImage
    }

    String hifiBam = basename(inputBam, ".reads.bam") + ".hifi.bam"

    command {
        set -e
        extracthifi ${inputBam} ${hifiBam}
    }

    output {
        File outputBam = "${hifiBam}"
    }

    runtime {
        docker: dockerImage
    }
}

task Primrose {
    input {
        File inputBam
        String dockerImage
    }

    String cpgBam = basename(inputBam, ".bam") + ".cpg.bam"

    command {
        set -e
        primrose \
        --keep-kinetics \
        --store-debug-tags \
        --log-file HG002.hifi.primrose.log \
        ${inputBam} \
        ${cpgBam}
    }

    output {
        File outputBam = "${cpgBam}"
    }

    runtime {
        docker: dockerImage
    }
}

task Pbmm2 {
    input {
        String sample
        String refName
        Array[File] inputBams
        File fasta
        String dockerImage
    }

    String outputBam = "${sample}" + ".hifi.cpg." +  "${refName}" + ".bam"
    command {
        set -e
        printf '%s\n' ${sep=' ' inputBams} > bam.fofn
        pbmm2 align --preset CCS \
        ${fasta} \
        bam.fofn \
        ${outputBam}
    }
    output {
        File outputBam = "${outputBam}"
    }
    runtime {
        docker: dockerImage
    }
}

task BamSort {
    input {
        File inputBam
        String dockerImage
    }
    String sortBam = basename(inputBam, ".bam") + ".sort.bam"
    command {
        set -e
        samtools sort \
        -o ${sortBam} \
        ${inputBam}
    }
    output {
        File outputBam = "${sortBam}"
    }
    runtime {
        docker: dockerImage
    }
}

task BamIndex {
    input {
        File inputBam
        String dockerImage
    }

    String bai = basename(inputBam) + ".bai"

    command {
        set -e
        samtools index \
        ${inputBam} \
        ${bai}
    }
    output {
        File outputBai = "${bai}"
    }
    runtime {
        docker: dockerImage
    }
}