import glob
import ete3


# fonction to reroot trees with dendropy package:
def reroot_tree(path_file, output_file):
    t = ete3.Tree(path_file)
    root_point = t.get_midpoint_outgroup()
    t.set_outgroup(root_point)
    # t.write(outfile=output_file, format=9)
    t.write(outfile=output_file, format=5)


for path_file in glob.glob("*.phylip_phyml_tree.txt"):
    print("working on " ,path_file)
    output_file = path_file.replace(".phylip_phyml_tree.txt", ".reroot_tree.txt")
    reroot_tree(path_file=path_file, output_file=output_file)
    
