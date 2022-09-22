# GDE-FINDER-
Gene Duplication Event Finder 



# GDE-FINDER.sh  

GDE-FINDER is a pipeline that allows the identification of all duplication events in a set of gene families. From a family of size 𝑁, it is programmed in such a way to recover (𝑁−1) pairs of genes. The development was done in a machine under the Linux operating system (Ubuntu).

The programming languages used in this pipelines:
- Shell (On Ubuntu)
- Python 3

The program consists of six steps. The first three steps are coded with the Shell language, and the last three are coded with Python.

1. This pipeline offers the user five options.
2. FASTA FILE CONTAINING THE PROTEOME. 
3. FASTA FILE CONTAINING THE CODING SEQUENCES (CDS).
4. ALIGNMENT FILE IN PHYLIP FORMAT (PROTEIN) (FILE NAMES MUST END WITH '.phylip').
5. PHYML FILES IN NEWICK FORMAT (FILE NAMES MUST END WITH '_phyml_tree.txt' ).


# Builder 

- Files:  
dmel-all-CDS-r6.32_id_modi.fasta  
Drosophila_melanogaster.BDGP6.32.pep.all.fa  
Galaxy18-[MCS_on_data_17].tabular  

- Program:  
GDE-FINDER.sh  
reroot_tree.py  
warning.py  
newick.py



Before run the GDE Finder we should make the file executable with this command 

```
chmod +x GDE-FINDER.sh
```
We should also make sure that we have successfully install the following soft wares :

- PHYML 
- Blast +
- mafft 


To run the python program we should install these packages :

- glob
```
pip install glob
```
- ete3

```
pip install ete3
```

- re (regular expressions)
```
pip install re
```
- random
```
pip install random
```
- os 
```
pip install os
```


For options 1 and 2 we should have 2 files :

- proteome file or CDS file.
- genes families names.

For option 3 and 4: 

- files ending with .phylip format.

For option 5 :

- files ending with _phyml_tree.txt .


IF THE PROGRAM SUSPECTS UNRESOLVED PHYLOGENETIC TREES IN A GENE FAMILY. YOU HAVE TWO CHOICES:

- 1 : KEEP ALL DUPLICATION EVENTS FOR THIS FAMILY (even if there are unresolved tree)
- 2 : RANDOMLY DRAW ONLY TWO GENES (only a couple of genes will represent the duplicated gene family)


# results 

- At the and of the process we will have a final file (final_output_file) having all detected duplicates genes. The first and second colomns are genes id proteins (couple of genes duplicated), and the third colomn is gene family name.   
- If we choose the option "SAVE INTERMEDIATE FILES":
we will have the followings intermediate files :
 -- alignement_file (alignment results)  
 -- events_couples_family (detected duplication events, for each family there are all genes couples)  
 -- fasta_files (of each protein)   
 -- id_prot_family  
 -- phyml_results (newick  tree format )  
 -- preprocessing_file   
 -- reroot_tree_intr_file (rerooted tree in newick format)  
 -- warning_file.txt (name of gene families which suppose have a rake or unresolved tree)  
 
 
 
- If we chhose to delate intermediate files:   
we will have just:   
-- final_output_file  
-- warning_file.txt  






##  Shell program File

### Program header 
```
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
# Menu option
echo -e "PLEASE CHOOSE YOUR FILE TYPE\n "
echo -e "1 : FASTA FILE CONTAINING THE PROTEOME "
echo -e "2 : FASTA FILE CONTAINING THE CODING SEQUENCES (CDS)"
echo -e "3 : ALIGNMENT FILE IN PHYLIP FORMAT (NUCLEOTIDE) (FILE NAMES MUST END WITH '.phylip')"
echo -e "4 : ALIGNMENT FILE IN PHYLIP FORMAT (PROTEIN) (FILE NAMES MUST END WITH '.phylip')"
echo -e "5 : PHYML FILES IN NEWICK FORMAT (FILE NAMES MUST END WITH '_phyml_tree.txt' )"
echo -e "6 : QUIT THE PROGRAM\n"
read -p 'ENTER YOUR CHOICE : ' TYPE_FILE
```

### Exemple with option "1"

```
# option 1 description
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
    
    
    # remove files  (To avoid errors in the program)
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
    
    
```

```
    # This step is to detect necessery information
    
    
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
  
```


```
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
    echo -e "--------------------------------------------DATA BASE CREATION DONE------------------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"
```

```
 # loop to file conteini
    while read line
    do
    set $line # cut to colomns
    echo $3 >> id_prot_family_name_$4
    done < gene_lenp_prot_Nfamily
    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "----------------------------------GENE FAMILY DISTRIBUTION FILES DONE----------------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    # remove all files with number of lines under 3 to keep only families size > 2 genes 
    
    for f in id_prot_family_name_*
    do
    #echo $i
    if [ $(wc -l < $f) -lt 3 ]
        then
        rm $f
    fi
    done
```
```
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
```


```
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
```


### Python program (ete package)

    
```
# fonction to reroot trees with dendropy package:
def reroot_tree(path_file, output_file):
    t = ete3.Tree(path_file)
    root_point = t.get_midpoint_outgroup()
    t.set_outgroup(root_point)
    # t.write(outfile=output_file, format=9)
    t.write(outfile=output_file, format=5)
for path_file in glob.glob("*.phylip_phyml_tree.txt"):
    print("working on " ,path_file)
    # replace the old name file by a new one 
    
    output_file = path_file.replace(".phylip_phyml_tree.txt", ".reroot_tree.txt")
    
    # reroot the phylogenitics trees 
    
    reroot_tree(path_file=path_file, output_file=output_file)
    
```

```
    
    
    
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
    
```

## This code allows to get randomly all possible paths for any given Newick file

### Step 0.1: Import Libraries & Packages #####


### STEP 1.0: Functions #####

    # Import Newick format from a file and return simplified label-based format
    def get_newick(file_path):
     """Return a dictionary with the Newick source format, simplified format, label format and a dictionary of leafs vs. labels
      If source file does not exist or empty, it return same dictionary, with None value for each key
       """

    # Check if the source file exists
    if os.path.exists(path_file):
        with open(file_path, "r") as newick_file:
            newick = newick_file.readline()

        # Check if the source file is empty
        if os.stat(file_path).st_size == 0:
            print("Warning: File (" + file_path + ") is empty")

            # Prepare the return dictionary
            return_dict = dict()
            return_dict["newick"] = None
            return_dict["new_format"] = None
            return_dict["label_format"] = None
            return_dict["label_dict"] = None

        else:
            newick = newick.replace(";", "").replace("\n", "").replace("-", "")

            # New format without distances and bootstrap details
            new_format = re.sub("\d+\.\d+", "", newick)
            # new_format = new_format.replace(":", "")
            new_format = re.sub("\:+\w*", "", new_format)
            new_format = new_format.replace("'" , "")

            # Label format by replacing each leaf with Tx where x is an auto-generated index
            terminals = new_format.replace("(", "").replace(")", "").split(",")
            label_format = new_format
            label_dict = dict()
            
            for idx, terminal in enumerate(terminals):
                label_format = label_format.replace(terminal, "T" + str(idx))
                label_dict["T" + str(idx)] = terminal

            # Prepare the return dictionary
            return_dict = dict()
            return_dict["newick"] = newick
            return_dict["new_format"] = new_format
            return_dict["label_format"] = label_format
            return_dict["label_dict"] = label_dict

    else:
        print("Warning: File (" + file_path + ") does not exist")

        # Prepare the return dictionary
        return_dict = dict()
        return_dict["newick"] = None
        return_dict["new_format"] = None
        return_dict["label_format"] = None
        return_dict["label_dict"] = None

    return return_dict

    # Get possible paths for any simplified Newick format
     def get_paths(label_format):
     """Return a dictionary including all possible paths"""

    # Check if the source is empty
    if label_format == None or label_format == "":
        print("Warning: No data to work on!")
        is_empty = True

    else:
        is_empty = False
    
    # Get the neighbors
    scan_idx = 0
    no_neighbors = False
    neighbors_dict = dict()
    all_neighbors = list()

    while not is_empty:
        print("SCAN_No", scan_idx)
        print("WORKING_ON --> ", label_format)

        # Check if direct neighbors exist
        neighbors = re.findall("\(T\d+,T\d+\)|\(T\d+,N\d+\)|\(N\d+,T\d+\)|\(N\d+,N\d+\)", label_format)
        
        if neighbors == []:
            no_neighbors = True
        else:
            print("NEIGHBORS_FOUND =", neighbors)
            all_neighbors.extend(neighbors)

        # Replace existing neighbors
        for idx, neighbor in enumerate(neighbors):
            label_format = label_format.replace(neighbor, "N" + str(scan_idx) + str(idx))

            neighbors_dict["N" + str(scan_idx) + str(idx)] = neighbor.replace("(", "").replace(")", "")

        # If no neighbors check replaced ones - priority for recent events
        if no_neighbors:
            # neighbors = re.findall("T\d+,T\d+", label_format) # DOES NOT CONSIDER LEFT AND RIGHT
            neighbors = re.findall("[A-Z]\d+,[A-Z]\d+", label_format)

            if neighbors == []:
                no_neighbors = True
            else:
                print("NEIGHBORS_FOUND =", neighbors)
                all_neighbors.extend(neighbors)

            # Replace existing neighbors
            for idx, neighbor in enumerate(neighbors):
                label_format = label_format.replace(neighbor, "N" + str(scan_idx) + str(idx))

                neighbors_dict["N" + str(scan_idx) + str(idx)] = neighbor.replace("(", "").replace(")", "")

        # Check if all leafs were checked
        end_test = re.findall(",", label_format)
        if end_test == []:
            print("\nNO_NEIGHBORS!")
            is_empty = True
        
        scan_idx += 1
        print()

    # print("ALL_NEIGHBORS =", all_neighbors)
    return neighbors_dict

    # Translate the paths dictionary
    def translate_paths(neighbors_dict, label_dict):
      """Return all possible paths for each node using real names"""
      all_paths_dict = dict()
       p = 0
        
    for key, value in neighbors_dict.items():
        print("PATH ", key, "----", value)

        paths_list = list()
        
        first_term = re.findall("T\d+,|N\d+,", value)[0].replace(",", "")
        second_term = re.findall(",T\d+|,N\d+", value)[0].replace(",", "")

        first_list = list()
        second_list = list()

        if first_term[0] == "T" and second_term[0] == "T":
            first_list.append(first_term)           
            second_list.append(second_term)

        elif first_term[0] == "T" and second_term[0] == "N":
            first_list.append(first_term)

            tmp_s = neighbors_dict[second_term].replace("(", "").replace(")", "").split(",")
            
            while "N" in " ".join(tmp_s):
                for t in tmp_s:
                    if t[0] == "N":
                        tmp_s.remove(t)
                        tmp_s.extend(neighbors_dict[t].replace("(", "").replace(")", "").split(","))
                
            second_list.extend(tmp_s)

        elif first_term[0] == "N" and second_term[0] == "T":
            tmp_f = neighbors_dict[first_term].replace("(", "").replace(")", "").split(",")
            
            while "N" in " ".join(tmp_f):
                for t in tmp_f:
                    if t[0] == "N":
                        tmp_f.remove(t)
                        tmp_f.extend(neighbors_dict[t].replace("(", "").replace(")", "").split(","))
                
            first_list.extend(tmp_f)

            second_list.append(second_term)

        elif first_term[0] == "N" and second_term[0] == "N":
            tmp_f = neighbors_dict[first_term].replace("(", "").replace(")", "").split(",")
            
            while "N" in " ".join(tmp_f):
                for t in tmp_f:
                    if t[0] == "N":
                        tmp_f.remove(t)
                        tmp_f.extend(neighbors_dict[t].replace("(", "").replace(")", "").split(","))
                
            first_list.extend(tmp_f)

            tmp_s = neighbors_dict[second_term].replace("(", "").replace(")", "").split(",")
            
            while "N" in " ".join(tmp_s):
                for t in tmp_s:
                    if t[0] == "N":
                        tmp_s.remove(t)
                        tmp_s.extend(neighbors_dict[t].replace("(", "").replace(")", "").split(","))
                
            second_list.extend(tmp_s)

        print("F", first_list)
        print("S", second_list) 

        for f in first_list:
            for s in second_list:
                print(label_dict[f] + " <-> " + label_dict[s])

                paths_list.append(label_dict[f] + " <-> " + label_dict[s])

            
        all_paths_dict["P" + str(p)] = paths_list
        p += 1

        print()

    return all_paths_dict

    # Select randomly one path for a given dictionary of paths
    def select_path(all_paths_dict, path_file):
      """Return dictionary contains one random path for each branch and create a file with the results"""
      family_name = re.findall("[0-9]+" ,path_file)
      family_name = family_name[0]
    
    random_paths_dict = dict()
    path_file = path_file.replace(".reroot_tree.txt", ".clade.txt")

    with open(path_file, "w") as path_file:
        print("RANDOM_PATHS")
        for key, value in all_paths_dict.items():
            random_paths_dict[key] = random.choice(value)
            print(key, random_paths_dict[key])

            # Save results to txt file
            path_file.write(random_paths_dict[key].split(" <-> ")[0] + "\t" + random_paths_dict[key].split(" <-> ")[1] + "\t" + family_name + "\n")

    return random_paths_dict

    # Human-friendly display function
    def display_newick(path_file, newick_dict):
       """Better display of the Newick dictionary created using the get_newick function"""

    print("\n########################################")
    print(path_file)
    print("########################################\n")

    print("NEWICK_INITIAL_FORMAT")
    print(newick_dict["newick"], "\n")
    
    print("NEWICK_NEW_FORMAT")
    print(newick_dict["new_format"], "\n")
    
    print("NEWICK_LABEL_FORMAT")
    print(newick_dict["label_format"], "\n")
    
    print("TERMINALS_DESCRIPTION")
    for key, value in newick_dict["label_dict"].items():
        print(key, "\t", value)

    print("\n########################################\n")

    ##### STEP 2.0: Main code #####

    if __name__ == "__main__":
        # Get the newick format
        # path_file = "./data/id_prot_family_name_103.phylip_phyml_tree.txt"

    
    for path_file in glob.glob("*.reroot_tree.txt"):
        
        my_newick = get_newick(path_file)

        display_newick(path_file, my_newick)

        # Get the possible paths
        my_label_format = my_newick["label_format"]
        my_neighbors_dict = get_paths(my_label_format)

        # print(my_neighbors_dict)

        # Translate paths from labels
        my_label_dict = my_newick["label_dict"]
        my_all_paths_dict = translate_paths(my_neighbors_dict, my_label_dict)

        # for key, value in my_all_paths_dict.items():
        #     print(key, value)

        # Select randomly one path for each branch
        my_selected_paths = select_path(my_all_paths_dict, path_file)


```
    
    
    
    
    
    
    # This file contains all the pairs of genes/proteins representing the duplication events for each family.
    # The first column: id_1
    # The second column: id_2
    # The third column: family name
    # The file is sorted on the third colomn (name of family / number )
    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "----------------------------------IDENTIFICATION OF DUPLICATION EVENTS DONE ---------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"
```
```
    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    echo -e "---------------------------------------------IDENTIFICATION OF RAKE------------------------------------------\n"
    echo -e "-------------------------------------------------------------------------------------------------------------\n"
    chmod +x warning.py
    python3 warning.py
    
    def random_clade(file_name):
    """
    This function is used to select randomaly two leafs if rake is detected.
    It replaces also the standard clade.txt file if selected by user.
    """
    
    # Identify the clade file
    file_name = file_name.replace(".reroot_tree.txt", ".clade.txt")
    
    print("WARNING! THE FILE -> " + file_name + " WILL BE REPLACED") 
    # Select randomaly two terminals
    with open(file_name, "r") as init_file:
        lines = init_file.readlines()
    
    random_line = random.choice(lines)
    print("RANDOMALY SELECTED LINE --->", random_line)
    # Replace std clade file with new one
    with open(file_name, "w") as init_file:
        init_file.write(random_line)
  user_choice = os.environ['USER_CHOICE']
  print(user_choice)
  # Open the newick file (rerooted tree)
  for file in glob.glob("id_*.reroot_tree.txt"):
    with open(file, "r") as newick_file:
        # Print the file name 
        print("working on ==>", file)
        # Read the file content 
        newick = newick_file.readline()
        # Pattern to find distances 
        newick_test = re.findall(":\d*\D+\d*| :\d+" , newick)
        #print(newick_test)
        # Save the distances in the liste 
        len_liste = []
        for b_length in newick_test :
            # Replace \: 
            b_length = b_length.replace(":", "")
            print(b_length)
            # b_length = re.sub("\(+[A-Za-z]+", "", b_length)
            len_liste.append(b_length)
        
        # Save the distance modified in data frame 
        df = pd.DataFrame(len_liste)
        #print(df)
        # Data frame with grouped distances 
        df2 = df.groupby(df.columns[0])[df.columns[-1]].count()
        # Loop to parse the distances groupes  
        for i in df2:
            #print(i)
            name_file = "warning_file.txt"
            # If the number of duplicates distances is > 3 , the name file will be save in warning file 
            # And the family tree  may have rake 
            if i >= 3:
                with open(name_file, "a") as file_out:
                    file_out.write(file + "\n")
                
                if user_choice == "1":
                    pass
                elif user_choice == "2":
                    random_clade(file) 
    
    
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
```





Footer
© 2022 GitHub, Inc.
Footer navigation
Terms
Privacy
Security
Status
Docs
Contact GitHub
Pricing
API
Training
Blog
About

