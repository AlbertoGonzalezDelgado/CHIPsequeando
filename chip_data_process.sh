#! /bin/bash

#It is an executable script for ChIP-seq data analysis
#Authors: Maria Ada Comparin, Ricardo Daniel de Arellano Casquete de Prado & Alberto Gonzalez Delgado.
#Date: December 2021

##Help message
if [ $# -ne 1 ]
then
	echo ""
	echo "usage chip_data_process.sh <param_file>"
	echo ""
	echo "param.file: file where the parameters are specified"
	echo ""
	echo ""
	echo ""
	echo "For more details, please check README.md and the examples given in param_input folder."
	echo ""
	echo ""
	echo ""
	echo ""
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
NUMREP=$(grep "number_copy:" $PFILE | awk '{ print $2 }')
FS=$(grep "scripts:" $PFILE | awk '{ print $2 }')
PAIREND=$(grep "number_chain_end:" $PFILE | awk '{ print $2 }')
UPSTREAM=$(grep "upstream_bp:" $PFILE | awk '{ print $2 }')
DOWNSTREAM=$(grep "dowstream_bp:" $PFILE | awk '{ print $2 }')
LENGTH=$(grep "word_length:" $PFILE | awk '{ print $2 }')
SIZE=$(grep "region_size:" $PFILE | awk '{ print $2 }')
PVALUE=$(grep "p_value:" $PFILE | awk '{ print $2 }')
CHIPS=()
CONTROLS=()
i=0
while [ $i -lt $NUMREP ]
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
echo "Number of samples is: $NUMSAM"
echo ""
echo "Chip(s): ${CHIPS[@]}"
echo "Control(s): ${CONTROLS[@]}"
echo ""
echo "    _______________                        |*\_/*|________ "
echo "   |  ___________  |     .-.     .-.      ||_/-\_|______  | "
echo "   | |           | |    .****. .****.     | |           | | "
echo "   | |   0   0   | |    .*****.*****.     | |   0   0   | | "
echo "   | |     -     | |     .*********.      | |     -     | | "
echo "   | |   \___/   | |      .*******.       | |   \___/   | | "
echo "   | |___     ___| |       .*****.        | |___________| | "
echo "   |_____|\_/|_____|        .***.         |_______________| "
echo "     _|__|/ \|_|_.............*.............._|________|_ "
echo "    / ********** \                          / ********** \ "
echo "  /  ************  \                      /  ************  \ "
echo " --------------------                    -------------------- "

##Creating working directory
cd $WD
mkdir $FD
cd $FD
mkdir genome annotation samples results
echo ""
echo ""
echo "Starting to create working directory"
echo ""
echo ""

##Adding genome.fa and annotation.gtf into their folders
cd genome
cp $GD genome.fa.gz
gunzip genome.fa.gz

cd ../annotation
cp $AN annotation.gtf.gz
gunzip annotation.gtf.gz
cd ..

echo ""
echo ""
echo "Now genome.fa and annotation.gtf are in their respective folders."
echo ""
echo ""

##Making genome index
cd genome
bowtie2-build genome.fa index
cd ..

echo ""
echo ""
echo "Genome index built"
echo ""
echo " -. .-.   .-. .-.   .-. .-.   .  "
echo " ||\|||\ /|||\|||\ /|||\|||\ /| "
echo " |/ \|||\|||/ \|||\|||/ \|||\|| "
echo " ~   \__/\_/   \__/\_/   \__/\_ "
echo ""
echo ""
##Making samples folder
cd samples
SE=1
PE=2
if [ $PAIREND -eq $SE ]
then
	echo "The samples are single end"
	echo ""
	i=0
	while [ $i -lt $NUMREP ]
	do
		j=$(($i + 1))
		mkdir control_$j chip_$j
		cd control_$j
		cp ${CONTROLS[$i]} control_${j}.fq.gz
		cd ../chip_$j
		cp ${CHIPS[$i]} chip_${j}.fq.gz
		((i++))
	done
elif [ $PAIREND -eq $PE ]
then
	echo "The samples are paired end"
	i=0
	while [ $i -lt $NUMREP ]
	do
		j=$(($i + 1))
		mkdir control_$j chip_$j
		cd control_$j
		cp ${CONTROLS[$i]} control_${j}_1.fq.gz
		cp ${CONTROLS[$j]} control_${j}_2.fq.gz
		cd ../chip_$j
		cp ${CHIPS[$i]} chip_${j}_1.fq.gz
		cp ${CHIPS[$j]} chip_${j}_2.fq.gz
		((i++))
		((i++))
  	done
else
	echo ""
	echo ""
	echo ""
	echo "Fatal error:"
	echo "If the samples are paired end or not is not specified"
	echo ""
	echo ""
	exit
fi
cd ../..


##Sample processing
mkdir logs
cd logs
i=0
while [ $i -lt $NUMREP ]
do
	j=$(($i + 1))
	sbatch $FS/sample_processing.sh $WD/$FD/samples/control_${j} $j $NUMREP $FS 1 $PAIREND $FD $NUMSAM $UPSTREAM $DOWNSTREAM $PVALUE $LENGTH $SIZE
	sbatch $FS/sample_processing.sh $WD/$FD/samples/chip_${j} $j $NUMREP $FS 2 $PAIREND $FD $NUMSAM $UPSTREAM $DOWNSTREAM $PVALUE $LENGTH $SIZE
	((i++))
done
echo ""
echo "All sample processing scripts have been submitted. Now you can have a cup of coffee or tea."
echo ""
echo "      (  )  (   )  ) "
echo "       ) (   )  (  ( "
echo "       ( )  (    ) ) "
echo "       _____________ "
echo "      <_____________> ___ "
echo "      |             |/ _ \ "
echo "      |               | | | "
echo "      |               |_| | "
echo "   ___|             |\___/ "
echo "  /    \___________/    \ "
echo "  \_____________________/ "
echo " "
echo "Check your queue system to follow the process. All messages will be written in files into log folder."
