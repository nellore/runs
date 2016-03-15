cd /home/student/jfortin/rail/sra
i=1 
for i in {1..25}
do
	qsub -cwd -V -l mem_free=10G,h_vmem=12G ata_ann.sh $i;
	echo $i
	sleep 30
done


cd /home/student/jfortin/rail/sra
i=1 
for i in {1..6}
do
	qsub -cwd -V -l mem_free=10G,h_vmem=12G ata_una.sh $i;
	echo $i
	sleep 30
done













