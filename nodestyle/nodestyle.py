import argparse
import pandas as pd
import os
from ete3 import Tree, TreeStyle, TextFace, NodeStyle

def parse_arguments():
    parser = argparse.ArgumentParser(description='python3 your_script.py /path/to/csv_file.csv /path/to/node_treefile.treefile /path/to/treefile.treefile.')
    parser.add_argument('csv_file', type=str, help='Path to the CSV file')
    parser.add_argument('node_treefile', type=str, help='Path to the node treefile')
    parser.add_argument('treefile', type=str, help='Path to the treefile')
    parser.add_argument('--output_dir', '-o', type=str, default='.', help='Output directory')
    return parser.parse_args()

args = parse_arguments()
#python3 a_read2.py /home/clemens/Downloads/M3560/M3560_2/PF3560.fasta_0.05.satute.csv /home/clemens/Downloads/M3560/M3560_2/PF3560.fasta_0.05.satute.node.treefile /home/clemens/Downloads/M3560/M3560_2/PF3560.fasta.treefile
# Read the CSV file into a pandas DataFrame
df2 = pd.read_csv(args.csv_file)

# Extract the values between the open bracket and the comma and create new columns 'Node_A' and 'Node_B'
df2['Node_A'] = df2['edge'].str.extract(r'\(([^,]+),')
df2['Node_B'] = df2['edge'].str.extract(r',([^,)]+)')

# Read the node treefile
with open(args.node_treefile, 'r') as f:
    s = f.read()

# Read the treefile
with open(args.treefile, 'r') as f:
    b = f.read()

s = s.replace('*', '')

# Print the result
pd.set_option('display.max_rows', None)
pd.set_option('display.max_columns', None)

def parentheses_function(s):
    stack = []
    parentheses = []
    for i, char in enumerate(s):
        if char == '(':
            stack.append(i)
        elif char == ')':
            if stack:
                open_idx = stack.pop()
                parentheses.append((open_idx, i))
    return parentheses


parentheses = parentheses_function(s)

data = []
for i, (open_idx, close_idx) in enumerate(parentheses):
    all_parentheses = s[open_idx:close_idx + 1]
    first_pair_part = all_parentheses

    node = s[close_idx + 1:close_idx + 5]
    last_idx = close_idx + 4
    while last_idx + 1 < len(s) and s[last_idx + 1].isdigit():
        node += s[last_idx + 1]
        last_idx += 1

    # print(first_pair_part, last_idx, s[last_idx+1])
    num_letters = ""
    j = last_idx
    while j < len(s) and s[j].isdigit():  # or s[j] == ".":
        num_letters += s[j]
        j += 1
    s2 = s[j:]

    if s2.find(")") < s2.find("(") or s2.find("(") == -1:
        end_node = s2[s2.find(")") + 1:s2.find(")") + 5]
        for char in s2[s2.find(")") + 5:]:
            if char.isdigit():
                end_node += char
            else:
                break
        data.append([first_pair_part, node, end_node])

    elif s2.find("(") < s2.find(")") and s2.find("(") >= 0:
        stack = []
        parentheses = []

        for i, char in enumerate(s2):
            if char == '(':
                stack.append(i)
            elif char == ')':
                if stack:
                    open_idx = stack.pop()
                    parentheses.append((open_idx, i))

        # Find the smallest open_idx
        smallest_open_idx = min([p[0] for p in parentheses])

        # Print the smallest pair
        for open_idx, close_idx in parentheses:
            if open_idx == smallest_open_idx:
                closed_parenthesis_pair = s2[open_idx:close_idx + 1]

                # Print remaining chars in s
                remaining_chars = s2[close_idx + 1:]
                # print(remaining_chars)
                idx = remaining_chars.find(')') + 1
                output = remaining_chars[idx:idx + 4]
                for char in remaining_chars[idx + 4:]:
                    if char.isdigit():
                        output += char
                    else:
                        break
                end_node = output
                data.append([first_pair_part, node, end_node])
                # print(end_node)
pd.set_option('display.max_colwidth', None)
df1 = pd.DataFrame(data, columns=['leafs_one_side', 'NodeA', 'NodeB'])[:-1]
df2['Node_A-Node_B'] = df2.apply(lambda row: (row['Node_A'] + '-' + row['Node_B']).strip(), axis=1)
df2['Node_A-Node_B'] = df2['Node_A-Node_B'].str.replace('*', '', regex=False)
df2['Node_A-Node_B'] = df2['Node_A-Node_B'].str.replace(r'\s+', '', regex=True)
#print(df2)
df1['NodeA-NodeB'] = df1['NodeA'] + '-' + df1['NodeB']

#print(df2['Node_A-Node_B'])
df2['Node_A'] = df2['Node_A'].str.replace(r'\*', '', regex=True)
df2['Node_B'] = df2['Node_B'].str.replace(r'\*', '', regex=True)




boot_data = []
parentheses = parentheses_function(b)
#print(parentheses)
for i, (open_idx, close_idx) in enumerate(parentheses):
    parentheses_pair = b[open_idx:close_idx + 1]
    last_idx = close_idx
    boot = ""
    while last_idx + 1 < len(b) and b[last_idx + 1].isdigit():
        boot += b[last_idx + 1]
        last_idx += 1
    boot_data.append([boot])
#print(boot_data)
df_boot = pd.DataFrame(boot_data, columns=['boots_value'])[:-1]
#print(df_boot)

df1 = df1.merge(df_boot, left_index=True, right_index=True)

#print(df1)

df1_indices = {row['NodeA-NodeB']: i for i, row in df1.iterrows()}
#
#print(df1_indices)
df2['df1_index'] = df2['Node_A-Node_B'].map(df1_indices)##############################
df2_sorted = df2.sort_values(by='df1_index')#####

# Merge df2_sorted and df1 using the 'df1_index' column
merged_df = df2_sorted.merge(df1, left_on='df1_index', right_index=True, how='left')
#print(merged_df)
# Drop the 'df1_index' column if you no longer need it
#merged_df.drop(columns='df1_index', inplace=True)

# Print the merged DataFrame with columns from both df2_sorted and df1
#print(merged_df[['delta', 'c_s', 'p_value', 'result_test', 'boots_value', 'Node_A', 'Node_B']])

merged_df = merged_df[['delta', 'c_s', 'p_value', 'result_test', 'boots_value', 'Node_A', 'Node_B']]
merged_df.fillna("", inplace=True)
#merged_df = merged_df.dropna(subset=['boots_value'])
#merged_df['boots_value'] = merged_df['boots_value'].str.replace(',', '').str.split('.').str[0].astype(int)
# Reset the index after removing rows
#merged_df.reset_index(drop=True, inplace=True)
#print(merged_df)
newick = s

t = Tree(newick, format=1)


# Create a TreeStyle
ts = TreeStyle()
ts.show_leaf_name = True
ts.show_scale = False  # Disable scale bar (if enabled, it might overlap with labels)

# Render the tree
canvas_width = 800
canvas_height = 600
ts.min_leaf_separation = 5
ts.scale = canvas_width / t.get_distance(t.get_leaves()[0], t.get_leaves()[-1])

t.render("tree_shorter_edges.png", w=canvas_width, h=canvas_height, tree_style=ts)
#ts.mode = "c"




labeled_branches = set()

def add_delta_label(node, delta, fgcolor):
    node.add_face(TextFace(f"{delta:.4f}", fsize=14, fgcolor=fgcolor), column=0, position="branch-top")


def add_boot_label(node, boot, fgcolor):
    node.add_face(TextFace(f"{boot:}", fsize=14, fgcolor=fgcolor), column=0, position="branch-top")

def add_p_label(node, p_value, fgcolor):
    node.add_face(TextFace(f"{p_value:.2e}", fsize=14, fgcolor=fgcolor), column=0, position="branch-top")

for _, row in merged_df.iterrows():
    node_b = row['Node_A']
    #print(node_b)
    delta = float(row['delta'])


    boots_value = row['boots_value']
    boot = '' if boots_value == '' else int(boots_value)  # Replace empty strings with 0

    p_value = float(row['p_value'])

    #print(delta)
    node_b_obj = t.search_nodes(name=node_b)

    if node_b_obj:
        node_b_obj = node_b_obj[0]

        if node_b_obj.name not in labeled_branches:
            add_boot_label(node_b_obj, boot, "black")
            add_delta_label(node_b_obj, delta, "green")
            add_p_label(node_b_obj, p_value, fgcolor="blue")

            labeled_branches.add(node_b_obj.name)




for leaf in t.iter_leaves():
    delta_label = merged_df[merged_df['Node_B'] == leaf.name]['delta'].values




for node in t.traverse():
    node_style = NodeStyle()
    node_style["size"] = 1  # Set node size
    node_style["fgcolor"] = "black"  # Set label color
    node_style = NodeStyle()
    node_style["fgcolor"] = "black"
    node_style["hz_line_width"] = 2  # Set horizontal line width for non-saturated nodes
    node_style["vt_line_width"] = 2  # Set vertical line width for non-saturated nodes
    node.set_style(node_style)

    if node.is_leaf():
        # Enlarge the leaf label
        node.add_face(TextFace(node.name, fsize=20), column=0, position="branch-right")
        # Remove the old leaf label
        node.name = ""

    if not node.is_leaf():
        # Add label to the inner node
        node.add_face(TextFace(node.name, fsize=18), column=0, position="branch-right")

    # Add branch length label
    branch_length = node.dist
    node.add_face(TextFace(f" {branch_length:.4f}", fsize=14, fgcolor="red"), column=1, position="branch-bottom")
t.show(tree_style=ts)


output_directory = args.output_dir

# Prepare output file paths
file_name = os.path.basename(args.csv_file)
output_file_name = os.path.splitext(file_name)[0]
pdf_output_file = os.path.join(output_directory, f"{output_file_name}.pdf")
#svg_output_file = os.path.join(output_directory, f"{output_file_name}.svg")

# Save rendered tree as PDF
t.render(pdf_output_file, tree_style=ts)
# Save rendered tree as SVG (uncomment if needed)
#t.render(svg_output_file, tree_style=ts)





