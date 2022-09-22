"""
This code allows to get randomly all possible paths for any given Newick file
"""

##### Step 0.1: Import Libraries & Packages #####

import os
import glob

import re
import random

##### STEP 1.0: Functions #####

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