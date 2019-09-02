awk -F "\t" '{if ($14 == "A") {print $14}}' YJ058c_transloc_all.tab | wc -l
awk -F "\t" '{if ($14 == "T") {print $14}}' YJ058c_transloc_all.tab | wc -l
awk -F "\t" '{if ($14 == "C") {print $14}}' YJ058c_transloc_all.tab | wc -l
awk -F "\t" '{if ($14 == "G") {print $14}}' YJ058c_transloc_all.tab | wc -l