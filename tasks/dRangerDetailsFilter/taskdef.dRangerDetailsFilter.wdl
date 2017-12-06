task dRangerDetailsFilter {

    #Inputs and constants defined here

    String pair_id
    File dRanger_input_file
    String min_tumor_count_threshold
    String max_normal_count_threshold
    String outputdir
    
    String output_disk_gb = "100"
    String boot_disk_gb = "10"
    String ram_gb = "2"
    String cpu_cores = "1"
    String preemption = "4"


    command {
cat <<EOF > pyscript.py

import subprocess,os,glob
def run(cmd):
    print('about to run')
    print(cmd)
    print('')
    subprocess.check_call(cmd,shell=True)

run('ln -sT `pwd` /opt/execution')
run('ln -sT `pwd`/../inputs /opt/inputs')
run('/opt/src/algutil/monitor_start.py')

# start task-specific calls
##########################



cwd = os.getcwd()

pair_id = '${pair_id}'
dRanger_input_file0  = '${dRanger_input_file}'

run('cp ' + dRanger_input_file0 + '  .')
dpath,dRanger_input_file  = os.path.split(dRanger_input_file0)
dRanger_input_filename, dRanger_input_extension = os.path.splitext(dRanger_input_file)
if dRanger_input_extension == '.gz':
    run('gunzip ' + dRanger_input_file )
    dRanger_input_file=os.path.join(cwd,dRanger_input_filename)


min_tumor_count_threshold='${min_tumor_count_threshold}'
max_normal_count_threshold='${max_normal_count_threshold}'
outputdir='${outputdir}'

run('ls -latr ')
run('ls -latr /opt/src/')
run('ls -latr ' + dRanger_input_file)

cmd = 'python /opt/src/dRangerDetailsFilter.py --pair_id %s  --dRanger_input_file %s \
--min_tumor_count_threshold %s  --max_normal_count_threshold %s \
--output %s '% (pair_id, dRanger_input_file,min_tumor_count_threshold,max_normal_count_threshold,outputdir)


run(cmd)

run('ls -latr ')

import time
#time.sleep(999999999)


#########################
# end task-specific calls
run('/opt/src/algutil/monitor_stop.py')
EOF

        cat pyscript.py 
        python pyscript.py

    }

    output {
        File dRangerDetailsFilter_all="${pair_id}.somatic.sv.detail.all.txt"
        File dRangerDetailsFilter_pass="${pair_id}.somatic.sv.detail.pass.txt"
    }

    runtime {
        docker : "docker.io/chipstewart/drangerdetailsfilter:1"
        memory: "${ram_gb}GB"
        cpu: "${cpu_cores}"
        disks: "local-disk ${output_disk_gb} HDD"
        bootDiskSizeGb: "${boot_disk_gb}"
        preemptible: "${preemption}"
    }


    meta {
        author : "Chip Stewart"
        email : "stewart@broadinstitute.org"
    }

}

workflow dRangerDetailsFilter_workflow {
    call dRangerDetailsFilter
}

