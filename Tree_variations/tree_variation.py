
from ete3 import Tree
import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import regex as re
from Bio import Phylo


def main(csv_file, print_p_value_output=False):
    directory = os.path.dirname(csv_file)
    b = None
    s = None
    tree = None



    # Traverse through all subdirectories and process tree files
    for root, dirs, files in os.walk(directory):
        for file_name in files:
            if file_name.endswith(".treefile"):
                treefile_path = os.path.join(root, file_name)
                with open(treefile_path, 'r') as f:
                    tree_content = f.read()
                    if "Node" in tree_content:
                        s = tree_content
                        try:
                            tree = Tree(treefile_path, format=1)
                        except Exception as e:
                            print(f"An error occurred while parsing the tree from {treefile_path}: {e}")
                    else:
                        b = tree_content

    try:
        df2 = pd.read_csv(csv_file)
    except FileNotFoundError:
        print(f"No {csv_file} found.")
        return
    except Exception as e:
        print(f"An error occurred while reading the CSV file: {e}")
        return


    nodes1 = set(['t001000', 't001001', 't001010', 't001011', 't001100', 't001101', 't001110', 't001111',
                  't000100', 't000101', 't000110', 't000111', 't000010', 't000011', 't000000', 't000001'])
    nodes2 = set(['t010000', 't010001', 't010010', 't010011', 't010100', 't010101', 't010110', 't010111',
                  't011000', 't011001', 't011010', 't011011', 't011100', 't011101', 't011110', 't011111'])
    nodes3 = set(['t100000', 't100001', 't100010', 't100011', 't100100', 't100101', 't100110', 't100111',
                  't101000', 't101001', 't101010', 't101011', 't101100', 't101101', 't101110', 't101111'])
    nodes4 = set(['t110000', 't110001', 't110010', 't110011', 't110100', 't110101', 't110110', 't110111',
                  't111000', 't111001', 't111010', 't111011', 't111100', 't111101', 't111110', 't111111'])

    all_taxa = nodes1 | nodes2 | nodes3 | nodes4

    # Function to find the center branch and sister clades
    def find_center_and_sister_clades(tree, all_taxa):
        for node in tree.traverse("preorder"):
            if node.is_leaf():
                continue

            left_leaves = set(leaf.name for leaf in node.children[0].iter_leaves()) if len(node.children) > 0 else set()
            right_leaves = set(leaf.name for leaf in node.children[1].iter_leaves()) if len(
                node.children) > 1 else set()

            if len(left_leaves) == 32 and len(right_leaves) == 32:
                left_sister = all_taxa - left_leaves
                right_sister = all_taxa - right_leaves
                return node, node.children[0], node.children[1], left_sister, right_sister

            combined_leaves = left_leaves | right_leaves
            if len(combined_leaves) == 32:
                remaining_leaves = all_taxa - combined_leaves
                if len(remaining_leaves) == 32:
                    left_sister = remaining_leaves & (nodes1 | nodes2)
                    right_sister = remaining_leaves & (nodes3 | nodes4)
                    return node, node.children[0], node.children[1], left_sister, right_sister

        return None, None, None, set(), set()

    parent_map = {}
    for edge in df2['edge']:
        # Remove parentheses and split by comma
        child, parent = edge.strip("()").split(", ")
        parent_map[child] = parent

        df2['Node_A'] = df2['edge'].str.extract(r'\(([^,]+),')
        df2['Node_B'] = df2['edge'].str.extract(r',([^,)]+)')

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

        pd.set_option('display.max_colwidth', None)
        df1 = pd.DataFrame(data, columns=['leafs_one_side', 'NodeA', 'NodeB'])[:-1]
        df1['NodeB'] = df1['NodeB'].replace(';', 'Node1')

        df2 = df2[
            (df2['Node_B'].str.strip().str.startswith("Node")) &
            (df2['Node_A'].str.strip().str.startswith("Node"))
            ]

        df2['Node_A-Node_B'] = df2.apply(lambda row: (row['Node_A'] + '-' + row['Node_B']).strip(), axis=1)
        df2['Node_A-Node_B'] = df2['Node_A-Node_B'].str.replace('*', '', regex=False)
        df2['Node_A-Node_B'] = df2['Node_A-Node_B'].str.replace(r'\s+', '', regex=True)

        df1['NodeA-NodeB'] = df1['NodeA'] + '-' + df1['NodeB']

        df1_indices = {row['NodeA-NodeB']: i for i, row in df1.iterrows()}

        df2['df1_index'] = df2['Node_A-Node_B'].map(df1_indices)  ##############################
        df2_sorted = df2.sort_values(by='df1_index')  #####

        df2_sorted = df2_sorted.drop(columns='df1_index')
        df2_sorted = df2_sorted.reset_index(drop=True)

        merged_df = df1.merge(df2_sorted[['delta', 'c_s', 'branch_length', 'p_value', 'result_test']], left_index=True,
                              right_index=True)

        boot_data = []
        parentheses = parentheses_function(b)

        for i, (open_idx, close_idx) in enumerate(parentheses):
            parentheses_pair = b[open_idx:close_idx + 1]
            last_idx = close_idx
            boot = ""
            while last_idx + 1 < len(b) and b[last_idx + 1].isdigit():
                boot += b[last_idx + 1]
                last_idx += 1
            boot_data.append([boot])

        df_boot = pd.DataFrame(boot_data, columns=['boots_value'])[:-1]

        def convert_to_float(b):
            return int(b) / 100

        df_boot['boots_float'] = df_boot['boots_value'].apply(convert_to_float)
        df_b = merged_df.merge(df_boot[['boots_float']], left_index=True, right_index=True)

        df_delta_node_boot = df_b[
            ['c_s', 'delta', 'p_value', 'branch_length', 'result_test', 'NodeA-NodeB', 'boots_float']]
        pd.set_option('display.max_rows', None)
        pd.set_option('display.max_columns', None)

        df_for_graph = df_delta_node_boot.copy()
        df_for_graph['A-B'] = df_delta_node_boot['NodeA-NodeB'].apply(
            lambda x: f"{int(x.split('-')[0][4:])}-{int(x.split('-')[1][4:])}")

        #print(df_for_graph)

        pd.set_option('display.max_colwidth', None)


    center_node_name, left_subtree, right_subtree, left_sister, right_sister = find_center_and_sister_clades(tree, all_taxa)
    print(f"Center node: {center_node_name.name}")


    center_node = center_node_name.name
    if center_node:
        left_leaves = set(leaf.name for leaf in left_subtree.iter_leaves())
        right_leaves = set(leaf.name for leaf in right_subtree.iter_leaves())
        combined_leaves = left_leaves | right_leaves

        nodes_in_combined = {
            taxon: (taxon in nodes1, taxon in nodes2, taxon in nodes3, taxon in nodes4)
            for taxon in combined_leaves
        }


        new_nodes1 = set()
        new_nodes2 = set()
        new_nodes3 = set()
        new_nodes4 = set()

        #  storing positions of True and False
        true_positions = set()
        false_positions = set(range(4))  # initialize as false

        for memberships in nodes_in_combined.values():
            for i, is_member in enumerate(memberships):
                if is_member:
                    true_positions.add(i)
                    false_positions.discard(i)

        if true_positions:
            true_positions_sorted = sorted(true_positions)
            if len(true_positions_sorted) > 1:
                nodes1_index = true_positions_sorted[0]
                nodes2_index = true_positions_sorted[1]

                if nodes1_index == 0:
                    new_nodes1 = nodes1
                elif nodes1_index == 1:
                    new_nodes1 = nodes2
                elif nodes1_index == 2:
                    new_nodes1 = nodes3
                elif nodes1_index == 3:
                    new_nodes1 = nodes4

                if nodes2_index == 0:
                    new_nodes2 = nodes1
                elif nodes2_index == 1:
                    new_nodes2 = nodes2
                elif nodes2_index == 2:
                    new_nodes2 = nodes3
                elif nodes2_index == 3:
                    new_nodes2 = nodes4

        if false_positions:
            for i, pos in enumerate(false_positions):
                if i == 0:
                    if pos == 0:
                        new_nodes3 |= nodes1
                    elif pos == 1:
                        new_nodes3 |= nodes2
                    elif pos == 2:
                        new_nodes3 |= nodes3
                    elif pos == 3:
                        new_nodes3 |= nodes4
                elif i == 1:
                    if pos == 0:
                        new_nodes4 |= nodes1
                    elif pos == 1:
                        new_nodes4 |= nodes2
                    elif pos == 2:
                        new_nodes4 |= nodes3
                    elif pos == 3:
                        new_nodes4 |= nodes4


        nodes1 = list(new_nodes3)
        nodes2 = list(new_nodes4)
        nodes3 = list(new_nodes1)
        nodes4 = list(new_nodes2)


        def calculate_output2(nodes1, nodes2):
            if nodes1[0][1] == nodes2[0][1]:
                return "yes"
            else:
                return "no"
        output2 = calculate_output2(nodes1, nodes2)
        print(f"Sisters: {output2}")

    def find_associated_nodes(node, parent_map):
        associated_nodes = []
        for child, parent in parent_map.items():
            if parent == node:
                associated_nodes.append(child)
        if node in parent_map:
            associated_nodes.append(parent_map[node])
        return associated_nodes

    def find_all_paths(start_node, parent_map):
        paths = []

        def recursive_find_paths(current_node, current_path, visited):
            if current_node.startswith('t'):
                paths.append(current_path)
                return
            # mark the node as visited
            visited.add(current_node)
            # finding nodes
            associated_nodes = find_associated_nodes(current_node, parent_map)
            for node in associated_nodes:
                if node not in visited:  # avoid loops
                    new_path = current_path + [node]
                    recursive_find_paths(node, new_path, visited.copy())

        recursive_find_paths(start_node, [start_node], set())

        return paths

    initial_node = center_node
    all_paths = find_all_paths(initial_node, parent_map)

    paths_group_1 = [path for path in all_paths if any(node in path for node in nodes1)]
    paths_group_2 = [path for path in all_paths if any(node in path for node in nodes2)]
    paths_group_3 = [path for path in all_paths if any(node in path for node in nodes3)]
    paths_group_4 = [path for path in all_paths if any(node in path for node in nodes4)]





    group_1_output = []
    for idx, group in enumerate([paths_group_1], start=1):
        for path in group:
            formatted_path = ", ".join([f"{path[i]}-{path[i + 1]}" for i in range(len(path) - 1)])
            order = formatted_path.count(",") + 1
            group_1_output.append((formatted_path, order))


    group_2_output = []

    for idx, group in enumerate([paths_group_2], start=1):

        for path in group:
            formatted_path = ", ".join([f"{path[i]}-{path[i + 1]}" for i in range(len(path) - 1)])
            order = formatted_path.count(",") + 1
            group_2_output.append((formatted_path, order))


    group_3_output = []

    for idx, group in enumerate([paths_group_3], start=1):

        for path in group:
            formatted_path = ", ".join([f"{path[i]}-{path[i + 1]}" for i in range(len(path) - 1)])
            order = formatted_path.count(",") + 1
            group_3_output.append((formatted_path, order))


    group_4_output = []

    for idx, group in enumerate([paths_group_4], start=1):

        for path in group:
            formatted_path = ", ".join([f"{path[i]}-{path[i + 1]}" for i in range(len(path) - 1)])
            order = formatted_path.count(",") + 1
            group_4_output.append((formatted_path, order))

    def create_order_dictionary(group_output):
        order_dictionary = {}
        for path, order in group_output:
            edges = path.split(", ")
            for i, edge in enumerate(edges):
                if i + 1 not in order_dictionary:
                    order_dictionary[i + 1] = set()  #  avoid duplicates
                order_dictionary[i + 1].add(edge)
        return order_dictionary

    order_dictionary_group_1 = create_order_dictionary(group_1_output)
    order_dictionary_group_2 = create_order_dictionary(group_2_output)
    order_dictionary_group_3 = create_order_dictionary(group_3_output)
    order_dictionary_group_4 = create_order_dictionary(group_4_output)

    def join_order_dictionaries(dict1, dict2):
        joined_dict = {}
        for order in set(dict1.keys()).union(dict2.keys()):
            joined_dict[order] = dict1.get(order, set()).union(dict2.get(order, set()))
        return joined_dict

    order_dictionary_group_1_2 = join_order_dictionaries(order_dictionary_group_1, order_dictionary_group_2)
    order_dictionary_group_3_4 = join_order_dictionaries(order_dictionary_group_3, order_dictionary_group_4)



    def replace_inner_edges_with_boots(order_dict, df):
        replaced_dict = {}
        for order, edges in order_dict.items():
            replaced_edges = []
            for edge in edges:
                # checking for inner edge (both parts do not start with 't')
                if not any(part.startswith('t') for part in edge.split('-')):
                    edge_no_star = edge.replace('*', '')
                    edge_reversed = '-'.join(edge_no_star.split('-')[::-1])
                    # finding DataFrame matches
                    match = df[(df['NodeA-NodeB'] == edge_no_star) | (df['NodeA-NodeB'] == edge_reversed)]
                    if not match.empty:
                        # getting UFBoot values
                        boots_value = match.iloc[0]['boots_float']
                        replaced_edges.append(f"{edge}:{boots_value}")
                    else:
                        replaced_edges.append(edge)
                else:
                    replaced_edges.append(edge)
            replaced_dict[order] = replaced_edges
        return replaced_dict

    order_dictionary_group_1_2_replaced = replace_inner_edges_with_boots(order_dictionary_group_1_2, df_for_graph)
    order_dictionary_group_3_4_replaced = replace_inner_edges_with_boots(order_dictionary_group_3_4, df_for_graph)

    print("\nOrder Dictionary Group 1 & 2 with Boots:")
    for order, edges in order_dictionary_group_1_2_replaced.items():
        print(f"{order}: {', '.join(edges)}")

    print("Order Dictionary Group 3 & 4 with Boots:")
    for order, edges in order_dictionary_group_3_4_replaced.items():
        print(f"{order}: {', '.join(edges)}")

    def calculate_order_averages(order_dict, df):
        order_averages = {}
        for order, edges in order_dict.items():
            boot_values = []
            for edge in edges:
                if not any(part.startswith('t') for part in edge.split('-')):
                    edge_no_star = edge.replace('*', '')
                    edge_reversed = '-'.join(edge_no_star.split('-')[::-1])
                    match = df[(df['NodeA-NodeB'] == edge_no_star) | (df['NodeA-NodeB'] == edge_reversed)]
                    if not match.empty:
                        boots_value = match.iloc[0]['boots_float']
                        boot_values.append(boots_value)
            if boot_values:
                average_boot = sum(boot_values) / len(boot_values)
                order_averages[order] = round(average_boot * 100, 1)
        return order_averages

    order_averages_group_1_2 = calculate_order_averages(order_dictionary_group_1_2, df_for_graph)
    order_averages_group_3_4 = calculate_order_averages(order_dictionary_group_3_4, df_for_graph)


    print("\nOrder Dictionary Group 1 & 2 with Averages:")
    for order, avg in order_averages_group_1_2.items():
        print(f"{order}: {avg}")

    print("Order Dictionary Group 3 & 4 with Averages:")
    for order, avg in order_averages_group_3_4.items():
        print(f"{order}: {avg}")

    def analyze_group_paths(paths_group, group_number, orders):
        for path in paths_group:
            order = len(path) - 1  # Each path length - 1 gives the order
            if order in orders[group_number]:
                orders[group_number][order] += 1
            else:
                orders[group_number][order] = 1

    # dictionary for counting branches per order
    orders = {1: {}, 2: {}, 3: {}, 4: {}}
    analyze_group_paths(paths_group_1, 1, orders)
    analyze_group_paths(paths_group_2, 2, orders)
    analyze_group_paths(paths_group_3, 3, orders)
    analyze_group_paths(paths_group_4, 4, orders)

    print("\nTree variations: ")
    for group_number, group_orders in orders.items():
        for order, count in group_orders.items():
            if count == 1:
                print(f"Single external branch in Group {group_number} with Order: {order}")
            else:
                pairwise_count = count // 2
                if pairwise_count > 0:
                    if pairwise_count == 1:
                        print(f"{pairwise_count} pairwise branch in Group {group_number} with Order: {order}")
                    else:
                        print(f"{pairwise_count} pairwise branches in Group {group_number} with Order: {order}")
                if count % 2 != 0:
                    print(f"Single external branch in Group {group_number} with Order: {order}")


if __name__ == "__main__":
    description = (
        "This script processes a CSV file along with specific Newick tree files. \n"
        "It identifies the center branch and calculates the average UFBoot values \n"
        "for branches of the same Order in Groups 1 & 2 and Groups 3 & 4. \n"
        "The results are then used to write out tree variations.\n\n"
        "For example, you can use the file 'four_0.01_1.fa_p4_0.05.satute.csv' "
        "located in the 'four_0.01' test folder.\n\n"
        "Command: python3 tree_variations.py <csv_file>"
    )

    parser = argparse.ArgumentParser(description=description, formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument("csv", type=str, help="Path to the CSV file")
    parser.add_argument("-p", "--print-p-value-output", action="store_true", help="Print p-value output")

    args = parser.parse_args()
    main(args.csv, args.print_p_value_output)