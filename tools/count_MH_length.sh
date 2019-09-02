awk -F "\t" '{if ($10 == "chr19" && length($15) == 0) {print $15}}' YJ051a_transloc_all.tab | wc -l
awk -F "\t" '{if ($10 == "chr19" && length($15) == 1) {print $15}}' YJ051a_transloc_all.tab | wc -l
awk -F "\t" '{if ($10 == "chr19" && length($15) == 2) {print $15}}' YJ051a_transloc_all.tab | wc -l
awk -F "\t" '{if ($10 == "chr19" && length($15) == 3) {print $15}}' YJ051a_transloc_all.tab | wc -l
awk -F "\t" '{if ($10 == "chr19" && length($15) >= 4 && length($15) <= 5) {print $15}}' YJ051a_transloc_all.tab | wc -l
awk -F "\t" '{if ($10 == "chr19" && length($15) >=6 && length($15) <= 10) {print $15}}' YJ051a_transloc_all.tab | wc -l
awk -F "\t" '{if ($10 == "chr19" && length($15) >=11 && length($15) <= 25) {print $15}}' YJ051a_transloc_all.tab | wc -l
awk -F "\t" '{if ($10 == "chr19" && length($15) >25) {print $15}}' YJ051a_transloc_all.tab | wc -l
