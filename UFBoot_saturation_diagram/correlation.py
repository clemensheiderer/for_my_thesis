import os
import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def plot_data(x, y, result_test, title, ylabel):
    x_start = 0
    x_end = 1

    mask = (x >= x_start) & (x <= x_end)
    x_filtered = x[mask]
    y_filtered = y[mask]
    result_test_filtered = result_test[mask]

    x_mean = np.mean(x_filtered)
    y_mean = np.mean(y_filtered)

    x_stddev = np.std(x_filtered, ddof=1)
    y_stddev = np.std(y_filtered, ddof=1)

    cov = np.cov(x_filtered, y_filtered)[0, 1]

    r = cov / (x_stddev * y_stddev)
    slope, intercept = np.polyfit(x_filtered, y_filtered, 1)

    saturated_mask = (result_test_filtered == 'Saturated')

    plt.scatter(x_filtered[~saturated_mask], y_filtered[~saturated_mask], s=30)
    plt.scatter(x_filtered[saturated_mask], y_filtered[saturated_mask], s=30, color='red')

    plt.plot(x_filtered, slope * x_filtered + intercept, linestyle='--', color='gray')

    plt.text(0.1, 0.8, 'r = {:.2f}'.format(r), transform=plt.gca().transAxes,
             bbox=dict(facecolor='white', alpha=0.5), fontsize=40)

    plt.xlabel("UFBoot support values", fontsize=40)
    plt.ylabel(ylabel, fontsize=40)

    plt.xticks(fontsize=35)
    plt.yticks(fontsize=35)

    plt.grid(True, linestyle=':', linewidth=0.5, color='gray')

    plt.title(title, fontsize=40)
    plt.show()




def main(csv_file, print_p_value_output=False):

    directory = os.path.dirname(csv_file)
    b = None  # UFBoot value file
    s = None  #nod file

    for file_name in os.listdir(directory):
        if file_name.endswith(".treefile"):
            treefile_path = os.path.join(directory, file_name)
            with open(treefile_path, 'r') as f:
                tree_content = f.read()
                if "Node" in tree_content:
                    s = tree_content
                else:
                    b = tree_content

    try:
        df2 = pd.read_csv(csv_file)

    except FileNotFoundError:
        print(f"No {csv_file} found.")


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

    df_delta_node_boot = df_b[['c_s', 'delta', 'p_value', 'branch_length', 'result_test', 'NodeA-NodeB', 'boots_float']]
    pd.set_option('display.max_rows', None)
    pd.set_option('display.max_columns', None)


    df_for_graph = df_delta_node_boot.copy()
    df_for_graph['A-B'] = df_delta_node_boot['NodeA-NodeB'].apply(
        lambda x: f"{int(x.split('-')[0][4:])}-{int(x.split('-')[1][4:])}")

    pd.set_option('display.max_colwidth', None)


    pearson_corr = df_for_graph['delta'].corr(df_for_graph['boots_float'], method='pearson')

    first_line = input("Enter the title for the plot: ")
    second_line = input("Enter text for the second line: ")

    if not print_p_value_output:
        ylabel = r"$\hat{E}[C_k^\partial]$"
        y = df_for_graph['delta']
    else:
        ylabel = r"$p_{value}$"
        y = df_for_graph['p_value']

    plot_data(df_for_graph['boots_float'], y, df_for_graph['result_test'], f"{first_line}\n{second_line}", ylabel)





if __name__ == "__main__":
    title = "****************Correlation diagram between saturation and UFBoot values***************"
    description = f"This script reads data from a CSV file and generates plots with correlation analysis. It requires a CSV file containing the data to be plotted. It will search for the .node.treefile and .treefile automatically. It needs all three files for creating the saturation, UFBoot diagram with a regression line."

    parser = argparse.ArgumentParser(description=description)
    parser.add_argument("csv", type=str, help="Path to the CSV file")
    parser.add_argument("-p", "--print-p-value-output", action="store_true", help="Print p-value output")

    # Print the title
    print("\033[1m" + title + "\033[0m\n")

    # Call main function with parsed arguments
    args = parser.parse_args()
    main(args.csv, args.print_p_value_output)
