#! /bin/bash

##Help message
if [ $# -ne 1 ]
then
	echo " "
	echo "usage pipernas <param_file>"
	echo " "
	echo "param.file: file where  the parameters ares specified"
	exit
fi

##Loading the parameters file
PFILE=$1
echo "Loading the parameters file"

WD=$(grep "working_directory:" $PFILE | awk '{ print $2 }')
FD=$(grep "folder_name:" $PFILE | awk '{ print $2 }')
GD=$(grep "genome:" $PFILE | awk '{ print $2 }')
AN=$(grep "annotation:" $PFILE | awk '{ print $2 }')
NUMSAM=$(grep "number_samples:" $PFILE | awk '{ print $2 }')
FS=$(grep "scripts:" $PFILE | awk '{ print $2 }')
CHIPS=()
CONTROLS=()
i=0
while [ $i -lt $NUMSAM ]
do
	j=$(($i + 1))
	CHIPS[$i]=$( grep "chip_${j}:" $PFILE | awk '{ print $2 }')
 	CONTROLS[$i]=$( grep "control_${j}:" $PFILE | awk '{ print $2 }')
	((i++))
done

echo ""
echo "Parameters has been loaded"
echo ""
echo "Working directory: $WD"
echo "Folder name: $FD"
echo "Genome is at: $GD"
echo "Annotation is at: $AN"
echo ""
echo "Number of samples is: $NUMSAM"
echo "Chips: ${CHIPS[@]}"
echo ""
echo "Controls: ${CONTROLS[@]}"
echo ""

##Creating working directory
cd $WD
mkdir $FD
cd $FD
mkdir genome annotation samples results

##Adding genome.fa and annotation.gtf
cd genome
cp $GD genome.fa.gz
gunzip genome.fa.gz

cd ../annotation
cp $AN annotation.gtf.gz
gunzip annotation.gtf.gz
cd ..

##Making genome index
cd genome
bowtie2-build genome.fa index
cd ..

##Vamos por aqui
##Making samples folder
cd samples

i=0
while [ $i -lt $NUMSAM ]
do
	j=$(($i + 1))
	mkdir sample_$j
	cd sample_$j
	cp ${SAMPLES[$i]} sample_${j}.fq.gz
	cd ..
	((i++))
done
cd ..
##Sample processing
mkdir logs
cd logs
i=0
while [ $i -lt $NUMSAM ]
do
	j=$(($i + 1))
	sbatch $PIPER/sample_processing.sh $WD/$FD/samples/sample_${j} $j $NUMSAM $PIPER
	((i++))
  echo "."
done
echo "All sample processing scripts have been submitted"
