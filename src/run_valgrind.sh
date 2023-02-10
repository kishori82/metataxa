#valgrind --leak-check=full \
#         --show-leak-kinds=all \
#         --track-origins=yes \
#         --verbose \
#         --log-file=valgrind-out.txt \
#         ./taxdistinct --tax-file ../data/TemR1170802.functional_and_taxonomic_table_pwy0.txt --pwy-file ../data/pwy0_orfs.txt 

./taxdistinct --tax-file ../data/TemR1170802.functional_and_taxonomic_table.txt --pwy-file ../data/pwys_orfs.txt 
