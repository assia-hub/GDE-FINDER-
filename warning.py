
import re 
import pandas as pd 
import glob
import random
import os



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