# For annotated junctions:
cd /users/jfortin/runs/sra/pca_analysis
i=1 
for i in {1..25}
do
	qsub -cwd -V -l mem_free=10G,h_vmem=12G createAta_ann.sh $i;
	echo $i
	sleep 30
done

# For unannotated junctions:
cd /users/jfortin/runs/sra/pca_analysis
i=1 
for i in {1..6}
do
	qsub -cwd -V -l mem_free=10G,h_vmem=12G createAta_una.sh $i;
	echo $i
	sleep 30
done