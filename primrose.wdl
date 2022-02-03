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
}


task Extracthifi {
    input {
        File inputBam
        String dockerImage
    }

    String hifiBam = basename(inputBam, ".reads.bam") + ".hifi.bam"

    command {
        set -e
        extracthifi ${inputBam} ${outputBam}
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
        ${outputBam}
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
    command {
        set -e
        echo ${sep='\n' inputBams} > fofn
        pbmm2 align --preset CCS \
        ${fasta} \
        fofn \
        ${outputBam}
    }
    output {
        File outputBam = "${sample}" + ".hifi.cpg." +  "${refName}" + ".bam"
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
        -o ${outputBam} \
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

    command {
        set -e
        samtools index \
        ${inputBam}
    }

    runtime {
        docker: dockerImage
    }
}