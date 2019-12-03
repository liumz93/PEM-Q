# echo #########
# cat YJ058c_stats.txt | awk -F ":" '{print $2}'
# echo #########
# awk -F "\t" '{if ($13 == "Confident") {print}}' YJ058c_all_vector.tab | wc -l
# awk -F "\t" '{if ($13 == "Half") {print}}' YJ058c_all_vector.tab | wc -l
# awk -F "\t" '{if ($13 == "Suspected") {print}}' YJ058c_all_vector.tab | wc -l

echo ######insertion type
awk -F "\t" '{if ($13 == "Confident" && $8 == "") {print}}' YJ058c_all_vector.tab | wc -l
awk -F "\t" '{if ($13 == "Confident" && $8 == "chr8") {print}}' YJ058c_all_vector.tab | wc -l
awk -F "\t" '{if ($13 == "Confident" && $8 != "chr8") {print}}' YJ058c_all_vector.tab | wc -l