for sample in $1
do
	echo now processing $1 ...
	awk -F "\t" '!($3 ~ "chrUn" || $3 ~ "random")' ${sample}.tlx > ${sample}_tlxfile
	echo begin calling peaks ...
	tlx2BED-MACS.pl ${sample}_tlxfile ${sample}.bed
	awk '{if ($2 >0) {print $1"\t"int($2)"\t"int($3)"\t"$4"\t"$5"\t"$6} }' ${sample}.bed > ${sample}_1.bed
	bedtools sort -i ${sample}_1.bed > ${sample}_2.bed
	macs2 pileup -i ${sample}_2.bed -o ${sample}_3.bed
	bedtools sort -i ${sample}_3.bed > ${sample}_final.bed
	# macs2 bdgpeakcall -l 20 -g 10 -i ${sample}_final.bed -o ${sample}_peaks
	macs2 callpeak -t ${sample}_final.bed -f BED -g hs --keep-dup all -n ${sample} --nomodel --extsize 50 -q 0.05 --llocal 10000000
	echo extract sequence from bed file ...
	# awk -F "\t" '!($1 ~ "track")' ${sample}_peaks > peaks
	awk -F "\t" '!($1 ~ "track")' ${sample}_peaks.narrowPeak > ${sample}_peaks
	awk '{if($2<30){print $1"\t"$2"\t"$3+30}}' ${sample}_peaks > ${sample}_posfile1.txt
	awk '{if($2>30){print $1"\t"$2-30"\t"$3+30}}' ${sample}_peaks > ${sample}_posfile2.txt
	cat ${sample}_posfile1.txt ${sample}_posfile2.txt > ${sample}_posfile.txt
	cat ${sample}_posfile.txt | awk '{print $1":"$2"-"$3}' > ${sample}_posfile_bed.txt
	run_extract_hg38.sh ${sample} > ${sample}_in.fa
	awk 'BEGIN{FS=" "}{if(!/>/){print toupper($0)}else{print $1}}' ${sample}_in.fa > ${sample}_out.fa
	rm ${sample}_1.bed ${sample}_2.bed ${sample}_3.bed ${sample}_peaks ${sample}_posfile_bed.txt ${sample}_summits.bed
	rm ${sample}_peaks.narrowPeak ${sample}_peaks.xls ${sample}_posfile.txt ${sample}_in.fa ${sample}_posfile1.txt ${sample}_posfile2.txt
	echo detecting off target sites ...
	detectOffTarget.R ${sample}_out.fa $2 $3 $sample
	awk -F : '{print $1}' ${sample}_offtarget-fwd.gff > ${sample}_11.bed
	awk -F : '{print $2}' ${sample}_offtarget-fwd.gff | awk -F - '{print $1}' > ${sample}_22.bed
	awk '{print $4"\t"$5}' ${sample}_offtarget-fwd.gff > ${sample}_33.bed
	paste ${sample}_11.bed ${sample}_22.bed ${sample}_33.bed > ${sample}_44.bed
	awk '{print $1"\t"$2+$3-2"\t"$2+$4-1"\t"".""\t"".""\t""+"}' ${sample}_44.bed > ${sample}_offtarget-fwd.bed
	sed -i '1,2d' ${sample}_offtarget-fwd.bed
	awk -F : '{print $1}' ${sample}_offtarget-rev.gff > ${sample}_111.bed
	awk -F : '{print $2}' ${sample}_offtarget-rev.gff | awk -F - '{print $1}' > ${sample}_222.bed
	awk '{print $4"\t"$5}' ${sample}_offtarget-rev.gff > ${sample}_333.bed
	paste ${sample}_111.bed ${sample}_222.bed ${sample}_333.bed > ${sample}_444.bed
	awk '{print $1"\t"$2+$3-2"\t"$2+$4-1"\t"".""\t"".""\t""-"}' ${sample}_444.bed > ${sample}_offtarget-rev.bed
	sed -i '1,2d' ${sample}_offtarget-rev.bed
	echo Detecting off targets ...
	echo fwd offtarget candidates
	wc -l ${sample}_offtarget-fwd.bed
	echo rev offtarget candidates
	wc -l ${sample}_offtarget-rev.bed
	grabPrey.R ${sample}_offtarget-fwd.bed ${sample}.tlx  ${sample}_ttt | awk '{print $2}' > ${sample}_reads.txt
	wc -l ${sample}_reads.txt
	paste ${sample}_offtarget-fwd.bed ${sample}_reads.txt > ${sample}_OT-fwd.bed
	grabPrey.R ${sample}_offtarget-rev.bed ${sample}.tlx  ${sample}_ttt | awk '{print $2}' > ${sample}_reads.txt
	wc -l ${sample}_reads.txt
	paste ${sample}_offtarget-rev.bed ${sample}_reads.txt > ${sample}_OT-rev.bed
	cat ${sample}_OT-fwd.bed ${sample}_OT-rev.bed > ${sample}_OT.bed
	sort -nr -k7 ${sample}_OT.bed > ${sample}_OT.sort.bed
	awk '{if ($7 >2) {print}}' ${sample}_OT.sort.bed > ${sample}_OffTargets_1.bed
	rm ${sample}_11.bed ${sample}_22.bed ${sample}_33.bed ${sample}_44.bed ${sample}_out.fa ${sample}_tlxfile ${sample}_offtarget-fwd.gff ${sample}_offtarget-rev.gff
	rm ${sample}_111.bed ${sample}_222.bed ${sample}_333.bed ${sample}_444.bed ${sample}_OT-fwd.bed ${sample}_OT-rev.bed ${sample}_reads.txt ${sample}_offtarget-rev.bed ${sample}_offtarget-fwd.bed
	rm ${sample}_OT.sort.bed ${sample}_OT.bed
	rm *ttt* *final.bed ${sample}.bed
	#extract off target sequence
	echo extract off target sequence...
	bedtools getfasta -fi /home/ubuntu/genomes/hg38/hg38.fa -bed ${sample}_OffTargets_1.bed -fo ${sample}_OffTargets -s
	awk '{print toupper($1)}' ${sample}_OffTargets | sed -n 'n;p' > seq.txt
	paste ${sample}_OffTargets_1.bed seq.txt > ${sample}_OffTargets_n.bed
	grep -v N ${sample}_OffTargets_n.bed > ${sample}_OffTargets.bed
	grep -v N ${sample}_OffTargets_n.bed > ${sample}_OffTargets_ne.bed
	grep -v -f ${sample}_OffTargets.bed ${sample}_OffTargets_ne.bed > ${sample}_OffTargets_new.bed
	rm ${sample}_OffTargets_1.bed seq.txt ${sample}_OffTargets
	echo Thanks for/ using program written by Mengzhu, bye~
done