import regex as re

def cutting(t, NodeA, NodeB):
    node1_node2 = t.search_nodes(name=NodeA)[0].search_nodes(name=NodeB)[0]
    subtree = node1_node2.detach()
    return subtree, t

def modify_tree_branch_length(t, NodeA, NodeB):
    # cuts between NodeA and NodeB
    subtree, t = cutting(t, NodeA, NodeB)

    # searches for nodeA
    nodeA = t.search_nodes(name=NodeA)[0]

    all_children = nodeA.get_children()


    neighbors = [child for child in all_children if child != nodeA]

    # gets neighboring nodes
    neighbor_names = [neighbor.name for neighbor in neighbors]

    neighbor_node = t.search_nodes(name=neighbor_names[0])[0]
    branch_length_nodeA = nodeA.dist
    branch_length_neighbor = neighbor_node.dist

    total_branch_length = branch_length_nodeA + branch_length_neighbor

    # deletes NodeA
    t.search_nodes(name=NodeA)[0].delete(prevent_nondicotomic=True, preserve_branch_length=False)

    # new branch length
    neighbor_node.dist = total_branch_length

    return t, subtree


def subtree_to_fasta(subtree):
    taxas = set()
    comma_taxas = re.findall(r',([^:,()]+):', subtree)
    taxas = {taxa for taxa in comma_taxas if not taxa.startswith('(')}

    parenthesis_taxas = re.findall(r'\(([^:(]+):', subtree)
    taxas.update({taxa for taxa in parenthesis_taxas if not taxa.startswith('(')})
    return taxas



def fasta_copy(fasta_file, taxas):
    sequences = []
    pattern = r">(.*?)\n(.*?)\n"
    matches = re.findall(pattern, fasta_file, re.DOTALL)

    for match in matches:
        taxa_name, sequence = match
        if taxa_name in taxas:
            #print(f"Matched taxa: {taxa_name}")
            sequences.append(f">{taxa_name}\n{sequence}")
    return sequences
