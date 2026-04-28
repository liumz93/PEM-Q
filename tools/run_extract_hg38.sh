for i in `cat ${1}_posfile_bed.txt`
do
	samtools faidx /Volumes/dataraid/users/mengzl/genome/hg38/hg38.fa $i 
	# samtools faidx /home/ubuntu/genomes/mm9/mm9.fa $i | tail -n 1
done