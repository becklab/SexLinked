#!/bin/csh -f
#Script designed to execute permutation analysis to validate LincsCloud sex-specific drug predictions
#
#Author: Jonathan Ma
#Created: 12/24/2014
#
#Must be run from the /home/user081/ESP_permutation_analysis directory
#CAUTION: run this script one cancer at a time by modifying allcancers.txt to only contain the
#cancers desired. extraneousdeletion.csh must also be running at the same time in order to ensure
#that the server does not run out of space.
#
#A general note on CSH scripts: before running a new script file, you must type the following command:
#chmod +x myscriptname.csh. This changes the permissions of the script, allowing it to be run.
foreach cc (`cat allcancers.txt`)
 #loop through all 1000 signatures
 @ j = 1
 while( $j <= 1000 )

  #set filename through concatenation
  #example: Permuted_Signatures/HNSC/1_HNSC_male_permuted_downInTumors_affy.gmt
  set parentdir="Permuted_Signatures_SAM/${cc}/"

 set maleupfile="$parentdir${j}_${cc}_male_permuted_upInTumors.gmt"
  set femaleupfile="$parentdir${j}_${cc}_female_permuted_upInTumors.gmt"
  set maledownfile="$parentdir${j}_${cc}_male_permuted_downInTumors.gmt"
  set femaledownfile="$parentdir${j}_${cc}_female_permuted_downInTumors.gmt"

  #make directory to put male analysis in
  set maledir="Permuted_Results_SAM/${cc}/male${j}"
  mkdir $maledir

  #make directory to put female analysis in
  set femaledir="Permuted_Results_SAM/${cc}/female${j}"
  mkdir $femaledir

  q sig_quest_tool --uptag $maledownfile --dntag $maleupfile --metric wtcs --out $maledir --mkdir 0
  q sig_quest_tool --uptag $femaledownfile --dntag $femaleupfile --metric wtcs --out $femaledir --mkdir 0

  @ j += 1
 end

end
