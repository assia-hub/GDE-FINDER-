#!/bin/bash
echo -e "                                             ------------------                                               "
echo -e "                                             ------------------                                               "
echo -e "°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°//\\GDE-FINDER//\\°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°"
echo -e "                                             ------------------                                               "
echo -e "                                             ------------------                                               "
echo -e "°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°° LAMME & LBBE °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°"
echo -e "°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°° Carène Rizzon & Emmanuelle Lerat °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°"
echo    "°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°° ||IDENTIFICATION OF GENE DUPLICATION EVENTS|| °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°"
echo -e "°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°° M1 GENIOMHE °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°"
echo    "°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°° DEVELOPER : BENMEHDIA ASSIA °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°"
echo    "°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°° CONTACT : benmehdia.assia@gmail.com °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°"
echo    "°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°° VERSION : 1.0 °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°"
echo -e  "°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°° 2021-2022 °°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°\n"
echo -e "°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°°\n"

echo -e "PLEASE CHOOSE YOUR FILE TYPE\n "

echo -e "1 : FASTA FILE CONTAINING THE PROTEOME "
echo -e "2 : FASTA FILE CONTAINING THE CODING SEQUENCES (CDS)"
echo -e "3 : ALIGNMENT FILE IN PHYLIP FORMAT (NUCLEOTIDE) (FILE NAMES MUST END WITH '.phylip')"
echo -e "4 : ALIGNMENT FILE IN PHYLIP FORMAT (PROTEIN) (FILE NAMES MUST END WITH '.phylip')"
echo -e "5 : PHYML FILES IN NEWICK FORMAT (FILE NAMES MUST END WITH '_phyml_tree.txt' )"
echo -e "6 : QUIT THE PROGRAM\n"

read -p 'ENTER YOUR CHOICE : ' TYPE_FILE



################################################################ CHOICE NUMBER 1 ####################################################################

########## THIS CHOICE ALLOWS US TO DO SOME file processing (PROTEIN FASTA FILE)
########## SET UP SEQUENCE ALIGNMENTS
########## TO CONSTRUCT PHYLOGENETIC TREES
########## REROOT THE TREES
########## IDENTIFICATION OF DUPLICATIONS EVENTS

#####################################################################################################################################################

if (( $TYPE_FILE == 1))
then


    read -p 'PLEASE ENTER THE NAME OF THE PROTEOME FILE : ' FILE_prot
    read -p 'PLEASE ENTER THE NAME OF THE FAMILY GENE FILE : ' FILE_family
    read -p 'PLEASE ENTER THE NAME OF THE DATA BASE FILE : ' FILE_data_base

    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "THIS STEP GENERATE INTERMEDIATE FILES, YOU HAVE THE POSSIBILITY TO SAVE THEM FOR CHECKS, OR ERASE THEM!!! "
    echo -e "1 : SAVE INTERMEDIATE FILES "
    echo -e "2 : DELATE INTERMEDIATE FILES "
    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    read -p 'ENTER YOUR CHOICE : '  INTERMEDIATE_FILES_choice
    echo -e "-------------------------------------------------------------------------------------------------------------\n"

    echo -e "IF THE PROGRAM SUSPECTS UNRESOLVED PHYLOGENETIC TREES IN A GENE FAMILY. YOU HAVE TWO CHOICES :\n"
    echo -e "1 : KEEP ALL DUPLICATION EVENTS FOR THIS FAMILY "
    echo -e "2 : RANDOMLY DRAW ONLY TWO GENES \n"
    read -p 'ENTER YOUR CHOICE : ' USER_CHOICE
    export USER_CHOICE

    echo -e "-------------------------------------------------------------------------------------------------------------"
    echo -e "DESCRIPTION OF INPUT FILES : \n "
    echo "-PROTEOME FILE :" $FILE_prot
    echo "-FAMILY GENE FILE :" $FILE_family
    echo "-DATA BASE FILE :" $FILE_data_base
    echo "-INTERMEDIATE FILE CHOICE :" $INTERMEDIATE_FILES_choice
    echo "-UNRESOLVED PHYLOGENETIC TREE CHOICE : " $USER_CHOICE
    echo -e "-------------------------------------------------------------------------------------------------------------"

    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "CLEAN DIRECTORY "
    echo -e "-------------------------------------------------------------------------------------------------------------\n"

    # remove files

    rm final_output_file
    rm *.phylip_phyml_tree.txt
    rm *.phylip_phyml_stats.txt
    rm *.clade.txt
    rm *.reroot_tree.txt
    rm warning_file.txt
    rm *.phylip
    rm id_*
    rm *.phr
    rm *.pin
    rm *.pog
    rm *.psd
    rm *.psi
    rm *.psq
    rm *.nhr
    rm *.nin
    rm *.nog
    rm *.nsd
    rm *.nsi
    rm *.nsq
    # if the directory exist, delate
    rm -r reroot_tree_intr_file/
    rm -r events_couples_family/
    rm -r phyml_results/
    rm -r preprocessing_file/
    rm -r fasta_files/
    rm -r alignement_file/
    rm -r id_prot_family/




    # pattern of prots_ids
    Patern_Id_Prot=">[[:alnum:]]*"

    # pattern of genes_ids
    Patern_Id_gene="gene:[[:alnum:]]*"

    # Extract prots_ids to file
    grep  -o $Patern_Id_Prot   $FILE_prot | sed  's/>//g' > id_prot_file

    # Extract genes_ids to file
    grep  -o $Patern_Id_gene   $FILE_prot | sed  's/gene://g' > id_gene_file


    # Extract length from file
    awk '/^>/{if (l!="") print l; l=0; next}{l+=length($0)}END{print l}' $FILE_prot > sequence_prot_length

    # Merge the 3 files to have one file with ids.prots length prot and ids.genes
    paste id_prot_file sequence_prot_length id_gene_file  > merge_file

    # extract length of prot in file cont id gene and length of longest proteine
    awk -F'\t' ' a[$3]<$2{a[$3]=$2; b[$3]=$1} END {for (i in a){ print i FS a[i] FS b[i]}}' merge_file | sort -d -k1 > id_gene_id_prot_Max_length


    ####################################################################################################################################################
                                            # At this step we have file with id_gene, length_prot, id_prot(longest)#
    ####################################################################################################################################################



    #### next step is to join FILE_family (id_gene + family name ) with id_gene_id_prot_Max_length

    # sorted_gene_family is file is a sort FILE_family, to joined with id_gene_id_prot_Max_length file
    awk ' $1!= "geneName" && $2!="family" {print$1"\t"$2}' $FILE_family | sort -d -k1   > sorted_gene_family

    # # join files to have general file with id_gene, lenprot, id_prot, Name_family.
    join -11  -21  id_gene_id_prot_Max_length sorted_gene_family | sort -n -k4  > gene_lenp_prot_Nfamily

    # Remove data base file


    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "--------------------------------------------FILES PROCESSING DONE--------------------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"

    # create database (i should activate this step in the end ) ===> takes to much time to be executed

    makeblastdb -in $FILE_data_base -parse_seqids -dbtype prot

    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "--------------------------------------------DATA BASE CREATION DONE--------------------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"



    # loop to file conteini
    while read line
    do
    set $line # cut to colomns
    echo $3 >> id_prot_family_name_$4

    done < gene_lenp_prot_Nfamily

    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "----------------------------------GENE FAMILY DISTRIBUTION FILES DONE----------------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"


    # remove all files with number of lines under 3
    for f in id_prot_family_name_*
    do
    #echo $i
    if [ $(wc -l < $f) -lt 3 ]
        then
        rm $f
    fi
    done




    # fonction allows us to get prot sequences fasta format for all id prot of one family
    get_seq(){
    FICHIER=$1
    DATA_BASE=$2
    while read line
    do
    blastdbcmd -db $2  -entry $line

    done < $FICHIER
}


    # Loop to generate all prot sequences of each family name(number)

    for i in id_prot_family_name_*
    do
    #echo $i
    get_seq $i $FILE_prot >  ${i%}.fa

    done


    # Loup to modifie the id proteine

    for i in id_prot*.fa
    do
        sed  's/>[[:alnum:]]\{4\}/>pp/g' $i >${i%.fa}.fasta
    done


    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "---------------------------------------------FASTA FILES DONE------------------------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"



    ####################################################################################################################################################
                                            # At this step we have fasta files of each family name, we can use them for alignemnt with MAFFT
    ####################################################################################################################################################



    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "----------------------------------------------ALIGNMENT STEP START-------------------------------------------\n"
    echo -e "###  ALIGNMENT PARAMETERS :\n"
    echo -e "- OUTPUT FORMAT ==> Phylip format  / Sorted"
    echo -e "- STRATEGY      ==> --auto\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"


    # Loop to have all multiple alignement of each family name
    for i in id_prot*.fasta
    do
    echo "WORKING ON : " $i
    "/usr/bin/mafft"  --auto --phylipout --reorder $i > ${i%.fasta}.phylip

    done


    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "----------------------------------------CONSTRUCTION OF PHYLOGENETIC TREES-----------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"


    # Loop to have phylogenitics tree for evry family gene
    for i in id_prot*.phylip
    do
    echo "WORKING ON : " $i
    phyml -i $i -d aa

    done

    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "-------------------------------------CONSTRUCTION OF PHYLOGENETIC TREES DONE---------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"

    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "------------------------------------------REROOT PHYLOGENETIC TREES------------------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"


    ######## Reroot trees given by Phyml
    # Make the script executable
    chmod +x reroot_tree.py
    # Run the script
    python3 reroot_tree.py

    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "----------------------------------------REROOT PHYLOGENETIC TREES DONE---------------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"




    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "--------------------------------------IDENTIFICATION OF DUPLICATION EVENTS-----------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"

    ######## Identify the duplication  events
    # Make the script executable
    chmod +x newick.py
    # Run the script
    python3 newick.py

    # This file contains all the pairs of genes/proteins representing the duplication events for each family.
    # The first column: id_1
    # The second column: id_2
    # The third column: family name
    # The file is sorted on the third colomn (name of family / number )

    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "----------------------------------IDENTIFICATION OF DUPLICATION EVENTS DONE ---------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"


    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "---------------------------------------------IDENTIFICATION OF RAKE------------------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"

    chmod +x warning.py
    python3 warning.py


    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "--------------------------------------------IDENTIFICATION OF RAKE DONE--------------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"


    cat id_*.clade.txt | sort -n -k3 | sed s/pp/FBpp/g > final_output_file




    if (( $INTERMEDIATE_FILES_choice == 1 ))
    then


        # Create tmp file of intermediate file of rerooted tree (newick)
        mkdir reroot_tree_intr_file
        mv *.reroot_tree.txt reroot_tree_intr_file/


        # Create tmp file of all files of couples of duplicate genes
        mkdir events_couples_family
        mv *.clade.txt events_couples_family/

        # Create tmp file of all results of PHYML
        mkdir phyml_results
        mv *.phylip_phyml_stats.txt phyml_results/
        mv *.phylip_phyml_tree.txt  phyml_results/

        # Create tmp file of preprocessing steps
        mkdir preprocessing_file
        mv sorted_gene_family preprocessing_file/
        mv merge_file preprocessing_file/
        mv id_gene_file preprocessing_file/
        mv id_prot_file preprocessing_file/
        mv sequence_prot_length preprocessing_file/
        mv id_gene_id_prot_Max_length preprocessing_file/
        mv gene_lenp_prot_Nfamily preprocessing_file/

        # Create tmp file of fasta file
        mkdir fasta_files
        mv id_*.fa  fasta_files/
        mv id_*.fasta fasta_files/

        # Create tmp file of id_CDS_families
        mkdir id_prot_family
        mv id*[[:digit:]] id_prot_family/

        mkdir alignement_file
        mv  *.phylip alignement_file/


        rm *.phr
        rm *.pin
        rm *.pog
        rm *.psd
        rm *.psi
        rm *.psq


        echo -e "-------------------------------------------------------------------------------------------------------------\n"
        echo -e "INTERMEDIATE FILES ARE SAVED INTO FILES ACCORDING TO EACH STEP \n"
        echo -e "RESULTS ARE SAVED TO FILE <<final_output_file>> \n"
        echo -e "-------------------------------------------------------------------------------------------------------------\n"



    elif (( $INTERMEDIATE_FILES_choice == 2 ))
    then
        rm *.reroot_tree.txt
        rm *.clade.txt
        rm *.phylip_phyml_stats.txt
        rm *.phylip_phyml_tree.txt
        rm id_prot_file
        rm sorted_gene_family
        rm merge_file
        rm gene_lenp_prot_Nfamily
        rm CDS_length
        rm gene_lenp_prot_Nfamily
        rm id_gene_file
        rm id_*.fa
        rm id_*.fasta
        rm  id*[[:digit:]]
        rm *.phylip
        rm *.phr
        rm *.pin
        rm *.pog
        rm *.psd
        rm *.psi
        rm *.psq
        echo -e "-------------------------------------------------------------------------------------------------------------\n"
        echo -e "RESULTS ARE SAVED TO FILE <<final_output_file>> \n"
        echo -e "-------------------------------------------------------------------------------------------------------------\n"

    fi
################################################################ CHOICE NUMBER 2 ####################################################################

########## THIS CHOICE ALLOWS US TO DO SOME file processing (NUCLEOTIDE FASTA FILE)
########## SET UP SEQUENCE ALIGNMENTS
########## TO CONSTRUCT PHYLOGENETIC TREES
########## REROOT THE TREES
########## IDENTIFICATION OF DUPLICATIONS EVENTS

#####################################################################################################################################################





elif (( $TYPE_FILE == 2 ))
then

    read -p 'PLEASE ENTER THE NAME OF THE CDS FILE : ' FILE_CDS
    read -p 'PLEASE ENTER THE NAME OF THE FAMILY GENE FILE : ' FILE_family
    read -p 'PLEASE ENTER THE NAME OF THE DATA BASE FILE : ' FILE_data_base

    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "THIS STEP GENERATE INTERMEDIATE FILES, YOU HAVE THE POSSIBILITY TO SAVE THEM FOR CHECKS, OR ERASE THEM!!! "
    echo -e "1 : SAVE INTERMEDIATE FILES "
    echo -e "2 : DELATE INTERMEDIATE FILES "
    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    read -p 'ENTER YOUR CHOICE : '  INTERMEDIATE_FILES_choice
    echo -e "-------------------------------------------------------------------------------------------------------------\n"

    echo -e "IF THE PROGRAM SUSPECTS UNRESOLVED PHYLOGENETIC TREES IN A GENE FAMILY. YOU HAVE TWO CHOICES :\n"
    echo -e "1 : KEEP ALL DUPLICATION EVENTS FOR THIS FAMILY "
    echo -e "2 : RANDOMLY DRAW ONLY TWO GENES \n"
    read -p 'ENTER YOUR CHOICE : ' USER_CHOICE
    export USER_CHOICE


    echo -e "-------------------------------------------------------------------------------------------------------------"
    echo -e "DESCRIPTION OF INPUT FILES : \n "
    echo "-CDS FILE :" $FILE_CDS
    echo "-FAMILY GENE FILE :" $FILE_family
    echo "-DATA BASE FILE :" $FILE_data_base
    echo "-INTERMEDIATE FILE CHOICE : " $INTERMEDIATE_FILES_choice
    echo "-UNRESOLVED PHYLOGENETIC TREE CHOICE : " $USER_CHOICE
    echo -e "-------------------------------------------------------------------------------------------------------------"




    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "CLEAN DIRECTORY "
    echo -e "-------------------------------------------------------------------------------------------------------------\n"

    # remove files
    rm final_output_file
    rm *.phylip_phyml_tree.txt
    rm *.phylip_phyml_stats.txt
    rm *.clade.txt
    rm *.reroot_tree.txt
    rm warning_file.txt
    rm *.phylip
    rm id_*
    rm *.phr
    rm *.pin
    rm *.pog
    rm *.psd
    rm *.psi
    rm *.psq
    rm *.nhr
    rm *.nin
    rm *.nog
    rm *.nsd
    rm *.nsi
    rm *.nsq


    # if the directory exist, delate
    rm -r reroot_tree_intr_file/
    rm -r events_couples_family/
    rm -r phyml_results/
    rm -r preprocessing_file/
    rm -r fasta_files/
    rm -r alignement_file/
    rm -r id_prot_family/

    # pattern of CDS_ids
    Patern_Id_CDS=">[[:alnum:]]*"

    # pattern of proteines_ids
    Patern_Id_prot="FlyBase:FBpp[[:digit:]]*"

    # Pattern of genes_ids

    Patern_Id_gene="=FBgn[[:digit:]]*"

    # Pattern of length CDS

    Patern_length_CDS="length=[[:digit:]]*"



    # Extract CDS_ids to file
    grep  -o $Patern_Id_CDS   $FILE_CDS  | sed "s/>//g" > id_CDS

    grep -o $Patern_Id_prot $FILE_CDS  | sed "s/FlyBase://g" > id_prot

    grep -o $Patern_Id_gene $FILE_CDS  | sed "s/=//g"  > id_gene

    grep -o $Patern_length_CDS $FILE_CDS  | sed "s/length=//g" > CDS_length


    paste -d "\t" id_CDS CDS_length id_prot id_gene  > merge_file


    awk -F'\t' ' a[$4]<$2{a[$4]=$2; b[$4]=$3; c[$4]=$1} END {for (i in a){ print i FS a[i] FS b[i] FS c[i]}}' merge_file | sort -d -k1 > id_gene_Max_length_id_prot_id_CDS


    ##########################################################################################################################################################
                                            # At this step we have file with id_gene, length_CDS, id_CDS(longest)#
    ##########################################################################################################################################################



    awk ' $1!= "geneName" && $2!="family" {print$1"\t"$2}' $FILE_family | sort -d -k1   > sorted_gene_family


    join -11  -21  id_gene_Max_length_id_prot_id_CDS  sorted_gene_family | sort -n -k5  > gene_lenp_prot_CDS_Nfamily

    # Remove data base file


    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "--------------------------------------------FILES PROCESSING DONE--------------------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"

    # Create CDS data base
    makeblastdb -in $FILE_data_base  -parse_seqids -dbtype nucl

    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "--------------------------------------------DATA BASE CREATION DONE--------------------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"



    # loop to file conteini
    while read line
    do
    set $line # cut to colomns
    echo $4 >> id_CDS_family_name_$5

    done < gene_lenp_prot_CDS_Nfamily

    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "----------------------------------GENE FAMILY DISTRIBUTION FILES DONE----------------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"




    # remove all files with number of lines under 3
    for f in id_CDS_family_name_*
    do
    #echo $i
    if [ $(wc -l < $f) -lt 3 ]
        then
        rm $f
    fi
    done


    # # fonction allows us to get prot sequences fasta format for all id prot of one family
    get_seq(){
        FICHIER=$1
        DATA_BASE=$2
        while read line
        do
        blastdbcmd -db $2  -entry $line

        done < $FICHIER
    }

    # # Loop to generate all prot sequences of each family name(number)

    for i in id_CDS_family_name_*
    do
    echo "WORKING ON : " $i
        get_seq $i $FILE_data_base  >  ${i%}.fa

    done

    # Loup to modifie the id proteine

    for i in id_CDS*.fa
    do
        sed  's/>[[:alnum:]]\{5\}/>CDS/g' $i >${i%.fa}.fasta
    done


    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "---------------------------------------------FASTA FILES DONE------------------------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"


    ##########################################################################################################################################################
                                # At this step we have fasta files of each family name, we can use them for alignemnt with MAFFT
    ##########################################################################################################################################################


    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "----------------------------------------------ALIGNMENT STEP START-------------------------------------------\n"
    echo -e "###  ALIGNMENT PARAMETERS :\n"
    echo -e "- OUTPUT FORMAT ==> Phylip format  / Sorted"
    echo -e "- STRATEGY      ==> --auto\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"


    # Loop to have all multiple alignement of each family name
    for i in id_CDS*.fasta
    do
    echo "WORKING ON : " $i
    "/usr/bin/mafft"  --auto --phylipout --reorder $i > ${i%.fasta}.phylip

    done


    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "----------------------------------------CONSTRUCTION OF PHYLOGENETIC TREES-----------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"

    # Loop to have phylogenitics tree for evry family gene
    for i in id_CDS*.phylip
    do
    echo "WORKING ON : " $i
    phyml -i $i -d nt

    done

    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "------------------------------------CONSTRUCTION OF PHYLOGENETIC TREES DONE----------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"


    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "------------------------------------------REROOT PHYLOGENETIC TREES------------------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"


    ######## Reroot trees given by Phyml
    # Make the script executable
    chmod +x reroot_tree.py
    # Run the script
    python3 reroot_tree.py

    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "---------------------------------------REROOT PHYLOGENETIC TREES DONE----------------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"



    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "--------------------------------------IDENTIFICATION OF DUPLICATION EVENTS-----------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"

    ######## Identify the duplication  events
    # Make the script executable
    chmod +x newick.py
    # Run the script
    python3 newick.py

    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "----------------------------------IDENTIFICATION OF DUPLICATION EVENTS DONE----------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"

    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "---------------------------------------------IDENTIFICATION OF RAKE------------------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"

    chmod +x warning.py
    python3 warning.py


    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "--------------------------------------------IDENTIFICATION OF RAKE DONE--------------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"


    cat *.clade.txt |  sort -n -k3  > final_output_file

    if (( $INTERMEDIATE_FILES_choice == 1 ))
    then

        # Create tmp file of intermediate file of rerooted tree (newick)
        mkdir reroot_tree_intr_file
        mv *.reroot_tree.txt reroot_tree_intr_file/


        # Create tmp file of all files of couples of duplicate genes
        mkdir events_couples_family
        mv *.clade.txt events_couples_family/

        # Create tmp file of all results of PHYML
        mkdir phyml_results
        mv *.phylip_phyml_stats.txt phyml_results/
        mv *.phylip_phyml_tree.txt phyml_results/

        # Create tmp file of preprocessing steps
        mkdir preprocessing_file
        mv id_CDS preprocessing_file/
        mv sorted_gene_family preprocessing_file/
        mv merge_file preprocessing_file/
        mv gene_lenp_prot_CDS_Nfamily preprocessing_file/
        mv CDS_length preprocessing_file/
        mv id_gene preprocessing_file/
        mv id_prot preprocessing_file/
        mv id_gene_Max_length_id_prot_id_CDS preprocessing_file/


        # Create tmp file of fasta file
        mkdir fasta_files
        mv id_*.fa fasta_files/
        mv id_*.fasta fasta_files/

        # Create tmp file of id_CDS_families
        mkdir id_CDS_family
        mv id*[[:digit:]] id_CDS_family/

        mkdir alignement_file
        mv  *.phylip alignement_file/

        rm *.nhr
        rm *.nin
        rm *.nog
        rm *.nsd
        rm *.nsi
        rm *.nsq

        echo -e "-------------------------------------------------------------------------------------------------------------\n"
        echo -e "INTERMEDIATE FILES ARE SAVED INTO FILES ACCORDING TO EACH STEP \n"
        echo -e "RESULTS ARE SAVED TO FILE <<final_output_file>> \n"
        echo -e "-------------------------------------------------------------------------------------------------------------\n"



    elif (( $INTERMEDIATE_FILES_choice == 2 ))
    then
        rm *.reroot_tree.txt
        rm *.clade.txt
        rm *.phylip_phyml_stats.txt
        rm *.phylip_phyml_tree.txt

        rm id_CDS
        rm id_prot
        rm sorted_gene_family
        rm merge_file
        rm gene_lenp_prot_CDS_Nfamily
        rm CDS_length
        rm id_gene_file
        rm id_gene
        rm id_*.fa
        rm id_*.fasta
        rm  id*[[:digit:]]
        rm *.phylip
        rm *.nhr
        rm *.nin
        rm *.nog
        rm *.nsd
        rm *.nsi
        rm *.nsq

        echo -e "-------------------------------------------------------------------------------------------------------------\n"
        echo -e "RESULTS ARE SAVED TO FILE <<final_output_file>> \n"
        echo -e "-------------------------------------------------------------------------------------------------------------\n"

    fi

################################################################ CHOICE NUMBER 3 ####################################################################

########## THIS CHOICE ALLOWS US TO CONSTRUCT PHYLOGENETIC TREES FROM NUCLIOTIDE FILE PHYLIP
########## REROOT THE TREES
########## IDENTIFICATION OF DUPLICATIONS EVENTS

#####################################################################################################################################################




elif (( $TYPE_FILE == 3 ))
then

    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "THIS STEP GENERATE INTERMEDIATE FILES, YOU HAVE THE POSSIBILITY TO SAVE THEM FOR CHECKS, OR ERASE THEM!!! "
    echo -e "1 : SAVE INTERMEDIATE FILES "
    echo -e "2 : DELATE INTERMEDIATE FILES "
    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    read -p 'ENTER YOUR CHOICE : '  INTERMEDIATE_FILES_choice
    echo -e "-------------------------------------------------------------------------------------------------------------\n"

    echo -e "IF THE PROGRAM SUSPECTS UNRESOLVED PHYLOGENETIC TREES IN A GENE FAMILY. YOU HAVE TWO CHOICES :\n"
    echo -e "1 : KEEP ALL DUPLICATION EVENTS FOR THIS FAMILY "
    echo -e "2 : RANDOMLY DRAW ONLY TWO GENES \n"
    read -p 'ENTER YOUR CHOICE : ' USER_CHOICE
    export USER_CHOICE


    echo -e "-------------------------------------------------------------------------------------------------------------"
    echo -e "DESCRIPTION OF INPUT FILES : \n "
    echo "-FILE TYPE : NUCLIOTIDE FILE PHYLIP"
    echo "-INTERMEDIATE FILE CHOICE : " $INTERMEDIATE_FILES_choice
    echo "-UNRESOLVED PHYLOGENETIC TREE CHOICE : " $USER_CHOICE
    echo -e "-------------------------------------------------------------------------------------------------------------"


    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "CLEAN DIRECTORY "
    echo -e "-------------------------------------------------------------------------------------------------------------\n"

    rm final_output_file
    rm *.phylip_phyml_tree.txt
    rm *.phylip_phyml_stats.txt
    rm warning_file.txt
    rm *.clade.txt
    rm *.reroot_tree.txt

    # if the directory exist, delate
    rm -r reroot_tree_intr_file/
    rm -r events_couples_family/
    rm -r phyml_results/
    rm -r preprocessing_file/
    rm -r fasta_files/
    rm -r alignement_file/
    rm -r id_prot_family/

    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "------------------------------------CONSTRUCTION OF PHYLOGENETIC TREES---------------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"


    # Loop to have phylogenitics tree for evry family gene
    for i in *.phylip
    do
    echo "WORKING ON THE NUCLEOTIDE FILE : " $i
    phyml -i $i -d nt

    done

    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "------------------------------------CONSTRUCTION OF PHYLOGENETIC TREES DONE ---------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"



    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "------------------------------------------REROOT PHYLOGENETIC TREES------------------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"


    ######## Reroot trees given by Phyml
    # Make the script executable
    chmod +x reroot_tree.py
    # Run the script
    python3 reroot_tree.py

    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "------------------------------------------REROOT PHYLOGENETIC TREES DONE-------------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"



    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "--------------------------------------IDENTIFICATION OF DUPLICATION EVENTS-----------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"

    ######## Identify the duplication  events
    # Make the script executable
    chmod +x newick.py
    # Run the script
    python3 newick.py


    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "------------------------------ ------IDENTIFICATION OF DUPLICATION EVENTS DONE ------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"


    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "---------------------------------------------IDENTIFICATION OF RAKE------------------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"

    chmod +x warning.py
    python3 warning.py


    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "--------------------------------------------IDENTIFICATION OF RAKE DONE--------------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"



    cat *.clade.txt | sort -n -k3 > final_output_file


    if (( $INTERMEDIATE_FILES_choice == 1 ))
    then


        # Create tmp file of intermediate file of rerooted tree (newick)
        mkdir reroot_tree_intr_file
        mv *.reroot_tree.txt reroot_tree_intr_file/


        # Create tmp file of all files of couples of duplicate genes
        mkdir events_couples_family
        mv *.clade.txt events_couples_family/

        # Create tmp file of all results of PHYML
        mkdir phyml_results
        mv *.phylip_phyml_stats.txt phyml_results/
        mv *.phylip_phyml_tree.txt phyml_results/

        echo -e "-------------------------------------------------------------------------------------------------------------\n"
        echo -e "INTERMEDIATE FILES ARE SAVED INTO FILES ACCORDING TO EACH STEP"
        echo -e "-------------------------------------------------------------------------------------------------------------\n"


    elif (( $INTERMEDIATE_FILES_choice == 2 ))
    then
        rm *.reroot_tree.txt
        rm *.clade.txt
        rm *.phylip_phyml_stats.txt
        rm *.phylip_phyml_tree.txt


    fi

################################################################ CHOICE NUMBER 4 ####################################################################

########## THIS CHOICE ALLOWS US TO CONSTRUCT PHYLOGENETIC TREES FROM PROTEIN FILE PHYLIP
########## REROOT THE TREES
########## IDENTIFICATION OF DUPLICATIONS EVENTS

#####################################################################################################################################################


elif (( $TYPE_FILE == 4 ))
then

    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "THIS STEP GENERATE INTERMEDIATE FILES, YOU HAVE THE POSSIBILITY TO SAVE THEM FOR CHECKS, OR ERASE THEM!!! "
    echo -e "1 : SAVE INTERMEDIATE FILES "
    echo -e "2 : DELATE INTERMEDIATE FILES "
    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    read -p 'ENTER YOUR CHOICE : '  INTERMEDIATE_FILES_choice
    echo -e "-------------------------------------------------------------------------------------------------------------\n"

    echo -e "IF THE PROGRAM SUSPECTS UNRESOLVED PHYLOGENETIC TREES IN A GENE FAMILY. YOU HAVE TWO CHOICES :\n"
    echo -e "1 : KEEP ALL DUPLICATION EVENTS FOR THIS FAMILY "
    echo -e "2 : RANDOMLY DRAW ONLY TWO GENES \n"
    read -p 'ENTER YOUR CHOICE : ' USER_CHOICE
    export USER_CHOICE


    echo -e "-------------------------------------------------------------------------------------------------------------"
    echo -e "DESCRIPTION OF INPUT FILES : \n "
    echo "-FILE TYPE : PROTEIN FILE PHYLIP"
    echo "-Intermediate file choice :" $INTERMEDIATE_FILES_choice
    echo "-UNRESOLVED PHYLOGENETIC TREE CHOICE : " $USER_CHOICE
    echo -e "-------------------------------------------------------------------------------------------------------------"

    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "CLEAN DIRECTORY "
    echo -e "-------------------------------------------------------------------------------------------------------------\n"

    rm final_output_file
    rm *.phylip_phyml_tree.txt
    rm *.phylip_phyml_stats.txt
    rm warning_file.txt
    rm *.clade.txt
    rm *.reroot_tree.txt

    # if the directory exist, delate
    rm -r reroot_tree_intr_file/
    rm -r events_couples_family/
    rm -r phyml_results/
    rm -r preprocessing_file/
    rm -r fasta_files/
    rm -r alignement_file/
    rm -r id_prot_family/

    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "------------------------------------------CONSTRUCTION OF PHYLOGENETIC TREES---------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"

    ## Loop to have phylogenitics tree for evry family gene
    for i in *.phylip
    do
    echo "WORKING ON THE PROTEIN FILE : " $i
    phyml -i $i -d aa

    done

    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "------------------------------------CONSTRUCTION OF PHYLOGENETIC TREES DONE ---------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"



    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "------------------------------------------REROOT PHYLOGENETIC TREES------------------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"


    ######## Reroot trees given by Phyml
    # Make the script executable
    chmod +x reroot_tree.py
    # Run the script
    python3 reroot_tree.py

    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "------------------------------------------REROOT PHYLOGENETIC TREES DONE--------------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"




    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "--------------------------------------IDENTIFICATION OF DUPLICATION EVENTS-----------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"

    ######## Identify the duplication  events
    # Make the script executable
    chmod +x newick.py
    # Run the script
    python3 newick.py


    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "------------------------------ ------IDENTIFICATION OF DUPLICATION EVENTS DONE ------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"


    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "---------------------------------------------IDENTIFICATION OF RAKE------------------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"

    chmod +x warning.py
    python3 warning.py


    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "--------------------------------------------IDENTIFICATION OF RAKE DONE--------------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"

    cat *.clade.txt |  sort -n -k3  | sed "s/pp/FBpp/g" > final_output_file

    if (( $INTERMEDIATE_FILES_choice == 1 ))
    then


        # Create tmp file of intermediate file of rerooted tree (newick)
        mkdir reroot_tree_intr_file
        mv *.reroot_tree.txt reroot_tree_intr_file/


        # Create tmp file of all files of couples of duplicate genes
        mkdir events_couples_family
        mv *.clade.txt events_couples_family/

        # Create tmp file of all results of PHYML
        mkdir phyml_results
        mv *.phylip_phyml_stats.txt phyml_results/
        mv *.phylip_phyml_tree.txt phyml_results/

        echo -e "-------------------------------------------------------------------------------------------------------------\n"
        echo -e "INTERMEDIATE FILES ARE SAVED INTO FILES ACCORDING TO EACH STEP \n"
        echo -e "RESULTS ARE SAVED TO FILE <<final_output_file>> \n"
        echo -e "-------------------------------------------------------------------------------------------------------------\n"

    elif (( $INTERMEDIATE_FILES_choice == 2 ))
    then
        rm *.reroot_tree.txt
        rm *.clade.txt
        rm *.phylip_phyml_stats.txt
        rm *.phylip_phyml_tree.txt

        echo -e "-------------------------------------------------------------------------------------------------------------\n"
        echo -e "RESULTS ARE SAVED TO FILE <<final_output_file>> \n"
        echo -e "-------------------------------------------------------------------------------------------------------------\n"
    fi

################################################################ CHOICE NUMBER 5 ####################################################################

########## THIS CHOICE ALLOWS US TO REROOT TREES
########## IDENTIFICATION OF DUPLICATIONS EVENTS

#####################################################################################################################################################

elif (( $TYPE_FILE == 5 ))
then


    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "THIS STEP GENERATE INTERMEDIATE FILES, YOU HAVE THE POSSIBILITY TO SAVE THEM FOR CHECKS, OR ERASE THEM!!! "
    echo -e "1 : SAVE INTERMEDIATE FILES "
    echo -e "2 : DELATE INTERMEDIATE FILES "
    echo -e "-------------------------------------------------------------------------------------------------------------"
    read -p 'ENTER YOUR CHOICE : '  INTERMEDIATE_FILES_choice
    echo -e "-------------------------------------------------------------------------------------------------------------\n"

    echo -e "IF THE PROGRAM SUSPECTS UNRESOLVED PHYLOGENETIC TREES IN A GENE FAMILY. YOU HAVE TWO CHOICES :\n"
    echo -e "1 : KEEP ALL DUPLICATION EVENTS FOR THIS FAMILY "
    echo -e "2 : RANDOMLY DRAW ONLY TWO GENES \n"
    read -p 'ENTER YOUR CHOICE : ' USER_CHOICE
    export USER_CHOICE


    echo -e "-------------------------------------------------------------------------------------------------------------"
    echo -e "DESCRIPTION OF INPUT FILES : \n "
    echo "-FILE TYPE : PHYML FILES IN NEWICK FORMAT"
    echo "-Intermediate file choice :" $INTERMEDIATE_FILES_choice
    echo "-UNRESOLVED PHYLOGENETIC TREE CHOICE : " $USER_CHOICE
    echo -e "-------------------------------------------------------------------------------------------------------------"


    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "CLEAN DIRECTORY "
    echo -e "-------------------------------------------------------------------------------------------------------------\n"

    # remove files to clean the directory
    rm warning_file.txt
    rm *.clade.txt
    rm *.reroot_tree.txt
    rm final_output_file

    # if the directory exist, delate
    rm -r reroot_tree_intr_file/
    rm -r events_couples_family/
    rm -r phyml_results/
    rm -r preprocessing_file/
    rm -r fasta_files/
    rm -r alignement_file/
    rm -r id_prot_family/


    ######## Reroot trees given by Phyml
    # Make the script executable
    chmod +x reroot_tree.py
    # Run the script
    python3 reroot_tree.py

    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "------------------------------------------REROOT PHYLOGENETIC TREES DONE--------------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"




    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "--------------------------------------IDENTIFICATION OF DUPLICATION EVENTS-----------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"

    ######## Identify the duplication  events
    # Make the script executable
    chmod +x newick.py
    # Run the script
    python3 newick.py


    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "------------------------------ ------IDENTIFICATION OF DUPLICATION EVENTS DONE ------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"


    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "---------------------------------------------IDENTIFICATION OF RAKE------------------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"

    chmod +x warning.py
    python3 warning.py


    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "--------------------------------------------IDENTIFICATION OF RAKE DONE--------------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"


    # Create final file of all possible couples
    cat *.clade.txt |  sort -n -k3 | sed 's/pp/FBpp/g'  > final_output_file


    if (( $INTERMEDIATE_FILES_choice == 1 ))
    then


        # Create tmp file of intermediate file of rerooted tree (newick)
        mkdir reroot_tree_intr_file
        mv *.reroot_tree.txt reroot_tree_intr_file/


        # Create tmp file of all files of couples of duplicate genes
        mkdir events_couples_family
        mv *.clade.txt events_couples_family/

        echo -e "-------------------------------------------------------------------------------------------------------------\n"
        echo -e "INTERMEDIATE FILES ARE SAVED INTO FILES ACCORDING TO EACH STEP \n"
        echo -e "RESULTS ARE SAVED TO FILE <<final_output_file>> \n"
        echo -e "-------------------------------------------------------------------------------------------------------------\n"


    elif (( $INTERMEDIATE_FILES_choice == 2 ))
    then
        rm *.reroot_tree.txt
        rm *.clade.txt

        echo -e "-------------------------------------------------------------------------------------------------------------\n"
        echo -e "RESULTS ARE SAVED TO FILE <<final_output_file>> \n"
        echo -e "-------------------------------------------------------------------------------------------------------------\n"
    fi



################################################################ CHOICE NUMBER 6 ####################################################################

########## EXIT THE PROGRM

#####################################################################################################################################################

elif (( $TYPE_FILE == 6 ))
then

    echo -e "******************************************* THANK YOU FOR CHOOSING GDE-FINDER APPLICATION !! *******************************************\n"
    exit 6




################################################################ CHOICE NUMBER 6 ####################################################################

########## IF THE USER DOES NOT CHOOSE ANY OF THE PROPOSED CHOICES, THE PROGRAM STOPS

#####################################################################################################################################################


else
    echo "WRONG CHOICE, PLEASE CHOOSE (1) OR (2) OR (3) OR (4) OR (5) OR (6) ! "
fi
