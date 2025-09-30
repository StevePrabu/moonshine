#!/bin/bash
usage()
{
echo "obs_image.sh [-o obsnum] [-s scratchID]
    -o  obsnum      : the observation id
    -s  scratchID   : scratchID
    -c  cal model   : the calibrator model" 1>&2;
exit 1;
}

obsnum=
scratchID=
cal=

while getopts 'o:s:c:' OPTION
do
    case "$OPTION" in
        o)
            obsnum=${OPTARG}
            ;;
        s)
            scratchID=${OPTARG}
            ;;
        c)
            cal=${OPTARG}
            ;;
        ? | : | h)
            usage
            ;;
    esac
done


# if obsid is empty then just pring help
if [[ -z ${obsnum} ]]
then
    echo 'obsnum is missing'
    usage
fi

if [[ -z ${scratchID} ]]
then
    echo 'scratchID is missing'
    usage
fi

if [[ -z ${cal} ]]
then
    echo 'calibrator model is missing'
    usage
fi

## load configurations
source config.txt

## run template script
script="${MYBASE}/queue/auto_${obsnum}.sh"
cat ${base}bin/auto.sh | sed -e "s:OBSNUM:${obsnum}:g" \
                                -e "s:SCRATCHID:${scratchID}:g" \
                                -e "s:CAL:${cal}:g" \
                                -e "s:BASE:${MYBASE}:g"> ${script}

output="${base}queue/logs/auto_${obsnum}.o%A"
error="${base}queue/logs/auto_${obsnum}.e%A"
sub="sbatch --begin=now+15 --output=${output} --error=${error} ${depend} -J auto_${obsnum} -M ${MYCLUSTER} ${script}"
jobid=($(${sub}))
jobid=${jobid[3]}

# rename the err/output files as we now know the jobid
error=`echo ${error} | sed "s/%A/${jobid}/"`
output=`echo ${output} | sed "s/%A/${jobid}/"`

echo "Submitted auto job as ${jobid}"



