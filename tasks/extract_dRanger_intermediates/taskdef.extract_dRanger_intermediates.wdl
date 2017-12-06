task extract_dRanger_intermediates {

    #Inputs and constants defined here
    File dRanger_intermediates
    String id 
    String preemptible_limit="4"
    String output_disk_gb = "100"
    String boot_disk_gb = "50"
    String ram_gb = "3"
    String cpu_cores = "1"
    command {
python_cmd="
import subprocess,os
def run(cmd):
    subprocess.check_call(cmd,shell=True)

run('ln -sT `pwd` /opt/execution')
run('ln -sT `pwd`/../inputs /opt/inputs')
run('/opt/src/algutil/monitor_start.py')

# start task-specific calls
##########################

run('mkdir tar')
os.chdir('tar')
run('pwd')
run('cp  ${dRanger_intermediates} .')
DRANGER_INTERMEDIATES_TAR = os.path.basename('${dRanger_intermediates}')    
run('tar Pxvf '+ DRANGER_INTERMEDIATES_TAR )
run('ls -latrh ')
os.chdir('..')
run('pwd')
run('ls -latrh ')

run('gunzip pipette_jobs/links_for_broad/dRangerPreProcess_Normal_sg_gather/sample.all.isz.gz  ')
run('cat  pipette_jobs/links_for_broad/dRangerPreProcess_Normal_sg_gather/sample.all.isz > ${id}.normal.isz  ')
run('gunzip pipette_jobs/links_for_broad/dRangerPreProcess_Tumor_sg_gather/sample.all.isz.gz  ')
run('cat pipette_jobs/links_for_broad/dRangerPreProcess_Tumor_sg_gather/sample.all.isz > ${id}.tumor.isz  ')
run('gunzip pipette_jobs/links_for_broad/dRangerRun/sample.dRanger_results.mat.gz ')
run('cp pipette_jobs/links_for_broad/dRangerRun/sample.dRanger_results.mat  ${id}.dRanger_results.mat  ')
run('gunzip pipette_jobs/links_for_broad/dRangerRun/sample.dRanger_results.forBP.txt.gz ')
run('cp pipette_jobs/links_for_broad/dRangerRun/sample.dRanger_results.forBP.txt  ${id}.dRanger_results.forBP.txt  ')

run('ls -latrh ')


#########################
# end task-specific calls
run('/opt/src/algutil/monitor_stop.py')
"
        echo "$python_cmd"
        python -c "$python_cmd"

    }

    output {
        File dRanger_tumor_isz="${id}.tumor.isz"
        File dRanger_normal_isz="${id}.normal.isz"
        File dRanger_forBP_mat="${id}.dRanger_results.mat"
        File dRanger_forBP_txt="${id}.dRanger_results.forBP.txt"
    }

    runtime {
        docker : "docker.io/chipstewart/extract_dranger_intermediates:1"
        memory: "${ram_gb}GB"
        cpu: "${cpu_cores}"
        disks: "local-disk ${output_disk_gb} HDD"
        bootDiskSizeGb: "${boot_disk_gb}"
        preemptible: "${preemptible_limit}"
    }


    meta {
        author : "Chip Stewart"
        email : "stewart@broadinstitute.org"
    }

}

workflow extract_dRanger_intermediates_workflow {
    call extract_dRanger_intermediates
}