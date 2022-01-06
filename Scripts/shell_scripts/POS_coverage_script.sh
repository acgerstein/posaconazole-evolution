#merge all coverage files
paste G2_Anc.coverage.txt G4_Anc.coverage.txt G10_Ev.coverage.txt G10_Q.coverage.txt G11_Ev.coverage.txt G11_Q.coverage.txt G12_Q.coverage.txt G2_Q.coverage.txt G3_Ev.coverage.txt  G4_Ev.coverage.txt G5_Ev.coverage.txt G6_Ev.coverage.txt G7_Ev.coverage.txt G8_Ev.coverage.txt G9_Ev.coverage.txt | cut -f 1,2,3,6,9,12,15,18,21,24,27,30,33,36,39,42,45 > coverage_sub.txt

# exclude mitochondria DNA

grep -v Ca19-mtDNA coverage_sub.txt > no_mtd_coverage_sub.txt
