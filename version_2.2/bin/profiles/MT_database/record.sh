Please not:

Animal_CDS_protein.fa is a file cat together by following files:

Annelida-segmented-worms_CDS_protein.fa
Arthropoda_CDS_protein.fa
Bryozoa_CDS_protein.fa
Chaetognatha-arrow-worms_CDS_protein.fa
Chordata_CDS_protein.fa
Cnidaria_CDS_protein.fa
Echinodermata_CDS_protein.fa
Mollusca_CDS_protein.fa
Nematoda_CDS_protein.fa
Nemertea-ribbon-worms_CDS_protein.fa
Platyhelminthes-flatworms_CDS_protein.fa
Porifera-sponges_CDS_protein.fa

# change ID format
perl -ne 'chomp; if(/^>/){s/\.1\_/\_/;} print "$_\n";' Animal_CDS_protein.fa > Animal_CDS_protein.fa1

ls *fa | grep -v Anima | while read f; do  perl -ne 'chomp; if(/^>/){s/\.1\_/\_/;} print "$_\n";' $f > $f.1 ;done

ls *.fa.1 | while read f; do sed -i 's#)##g' $f; done

ls *.fa.1 | while read f; do sed -i 's#(##g' $f; done
ls *.fa.1 | perl -ne 'chomp; my $new=$_; $new=~s/\.1$//; `mv $_ $new`;'


# 20180517
# replace following in file 'Animal_CDS_protein.fa' and 'Chordata_CDS_protein.fa'
>gi_NC_KF937873_COX1_Mus_musculus_227_aa
MAYPFQLGLQDATSPIMEELMNFHDHTLMIVFLISSLVLYIISLMLTTKLTHTSTMDAQEVETIWTILPAVILIMIALPSLRILYMMDEINNPVLTVKTMGHQWYWSYEYTDYEDLCFDSYMIPTNDLKPGELRLLEVDNRVVLPMELPIRMLISSEDVLHSWAVPSLGLKTDAIPGRLNQATVTSNRPGLFYGQCSEICGSNHSFMPIVLEMVPLKYFENWSASMI
# with 
>gi_NC_KF937873_COX1_Mus_musculus_514_aa
MFINRWLFSTNHKDIGTLYLLFGAWAGMVGTALSILIRAELGQPGALLGDDQIYNVIVTAHAFVMIFFMVMPMMIGGFGNWLVPLMIGAPDMAFPRMNNMSFWLLPPSFLLLLASSMVEAGAGTGWTVYPPLAGNLAHAGASVDLTIFSLHLAGVSSILGAINFITTIINMKPPAMTQYQTPLFVWSVLITAVLLLLSLPVLAAGITMLLTDRNLNTTFFDPAGGGDPILYQHLFWFFGHPEVYILILPGFGIISHVVTYYSGKKEPFGYMGMVWAMMSIGFLGFIVWAHHMFTVGLDVDTRAYFTSATMIIAIPTGVKVFSWLATLHGGNIKWSPAMLWALGFIFLFTVGGLTGIVLSNSSLDIVLHDTYYVVAHFHYVLSMGAVFAIMAGFVHWFPLFSGFTLDDTWAKAHFAIMFVGVNMTFFPQHFLGLSGMPRRYSDYPDAYTTWNTVSSMGSFISLTAVLIMIFMIWEAFASKREVMSVSYASTNLEWLHGCPPPYHTFEEPTYVKVK
