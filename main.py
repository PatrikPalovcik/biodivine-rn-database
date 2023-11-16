from collections import defaultdict

# methods to go through normalization database


def all_main_names(file):
    file = open(file, 'r')
    original_names = []
    for line in file:
        columns = line.split()
        original_names.append(columns[2])
    return original_names


def read_synonyms(file):
    dictionary = {}
    with open(file, 'r') as file:
        for line in file:
            columns = line.split()
            columns_2 = columns[2]
            columns_4 = set(columns[4].split('|'))
            if columns_4 != {'-'}:
                dictionary[columns_2] = columns_4
    return dictionary


# methods to find tuples from databases


def read_pathway(file, synonyms):
    interactions = []
    file = open(file, 'r')
    for line in file:
        line = line.split()
        interactions.append((line[0], line[2], 'PathwayCommons', '-', '-', '-', '-', '-'))
    return interactions


def read_maayanlab(file, synonyms):
    interactions = []
    file = open(file, 'r')
    for line in file:
        line = line.split()
        first_column = line[0].split('_')
        for i in range(1, len(line)):
            interactions.append((first_column[0], line[i], 'Maayanlab', '-', '-', '-', '-', '-'))
    return interactions


def read_grndb(file, synonyms):
    interactions = []
    file = open(file, 'r')
    for line in file:
        line = line.split()
        if line[4] == 'High' or line[4] == 'Low':
            line.append(line[4])
            line[4] = '-'
        interactions.append((line[0], line[1], 'Grndb', '-', line[2], line[3], line[4], line[5]))
    return interactions


def read_grnpedia(file, synonyms):
    interactions = []
    file = open(file, 'r')
    for line in file:
        line = line.split()
        interactions.append((line[0], line[1], 'Grnpedia', line[2], '-', '-', '-', '-'))
    return interactions


# method to find value in dictionary

def find_value_in_dict(dictionary, value):
    for key, value_list in dictionary.items():
        if value in value_list:
            return key, value_list
    return None


def find_most_known_name(synonyms, name):
    for key, value_list in synonyms.items():
        if key == name:
            return name
        if name in value_list:
            return key
    return name

# methods to find successors and predecessors


def find_successors_interactions(interactions):
    successors = defaultdict(list)
    interactions = open(interactions, 'r')
    for interaction in interactions:
        interaction = interaction.split(';')
        successors[interaction[0]].append((interaction[1], interaction[2]))
    return successors


def find_predecessors_interactions(interactions):
    predecessors = defaultdict(list)
    interactions = open(interactions, 'r')
    for interaction in interactions:
        interaction = interaction.split(';')
        predecessors[interaction[1]].append((interaction[0], interaction[2]))
    return predecessors


# methods to find all target genes for TF

def find_tgs_for_tf(successors, synonyms):
    while True:
        gene = input("Input a transcription factor gene or type exit to exit the application: ")
        main_names = all_main_names('Homo_sapiens.gene_info')
        if gene == "exit":
            break
        if gene not in main_names:
            known_name = find_value_in_dict(synonyms, gene)
            if known_name is not None:
                gene = known_name[0]
        if successors.get(gene) is not None:
            print_info(gene, synonyms, main_names, successors, True)
        else:
            print("This gene is not correct or is not among successors.")


# method to find all transcription factors for target gene

def find_tfs_for_tg(predecessors, synonyms):
    while True:
        gene = input("Input a target gene or type exit to exit the application: ")
        main_names = all_main_names('Homo_sapiens.gene_info')
        if gene == "exit":
            break
        if gene not in main_names:
            known_name = find_value_in_dict(synonyms, gene)
            if known_name is not None:
                gene = known_name[0]
        if predecessors.get(gene) is not None:
            print_info(gene, synonyms, main_names, predecessors, False)
        else:
            print("This gene is not correct or is not among predecessors.")


def print_info(gene, synonyms, main_names, successors_or_predecessors, succ_or_pred_bool):
    if succ_or_pred_bool:
        print("gene " + gene + " regulates these genes")
    else:
        print("gene " + gene + " is regulated by these genes")
    for tf_or_tg in successors_or_predecessors.get(gene):
        print(tf_or_tg, end='')
        synonyms_for_transcription_factor = find_value_in_dict(synonyms, tf_or_tg[0])
        if tf_or_tg[0] in main_names:
            if tf_or_tg[0] in synonyms:
                print(" Synonyms: ", end='')
                for synonym in synonyms[tf_or_tg[0]]:
                    print(synonym, end=' ')
            print()
        elif synonyms_for_transcription_factor is not None:
            print(" Commonly known as: " + synonyms_for_transcription_factor[0], end='')
            print(", Other synonyms: ", end='')
            for synonym in synonyms_for_transcription_factor[1]:
                print(synonym, end=' ')
            print()
        else:
            print()
    print('\n')


def find_info_about_tf_and_tg(interactions, synonyms):
    main_names = all_main_names('Homo_sapiens.gene_info')
    while True:
        found_interaction = False
        transcription_factor = input("Input a transcription factor: ")
        target_gene = input("Input a target gene: ")
        transcription_factor = find_most_known_name(synonyms, transcription_factor)
        target_gene = find_most_known_name(synonyms, target_gene)
        searched_interaction = binary_search_table(transcription_factor, target_gene, interactions)
        if searched_interaction is not None:
            print(searched_interaction)
            found_interaction = True
        # for item in interactions:
        #     item = item.split(';')
        #     if find_most_known_name(synonyms, item[0]) == transcription_factor and \
        #             find_most_known_name(synonyms, item[1]) == target_gene:
        #         print("Interaction between transcription factor and target gene from databases: " + str(item[2]))
        #         found_interaction = True
        #         break
        if not found_interaction:
            print("Interaction between these two genes was not found.")
        continue_inputting = input("If you want to exit, type E, else press enter: ")
        if continue_inputting == "E":
            break


def find_info_about_group(interactions_file, synonyms):
    while True:
        found_interaction = False
        genes = input("Input genes separated by spaces: ")
        genes = genes.split(" ")
        normalized_genes = []
        for gene in genes:
            gene = find_most_known_name(synonyms, gene)
            normalized_genes.append(gene)
        if len(normalized_genes) < 2:
            print("There must be at least two genes.")
        for tf in normalized_genes:
            for tg in normalized_genes:
                searched_interaction = binary_search_table(tf, tg, interactions_file)
                if searched_interaction is not None:
                    print(searched_interaction)
                    found_interaction = True

        # for item in interactions:
        #     item = item.split(';')
        #     first_gene = find_most_known_name(synonyms, item[0])
        #     second_gene = find_most_known_name(synonyms, item[1])
        #     if first_gene in normalized_genes and second_gene in normalized_genes:
        #         print(item[0] + " regulates " + item[1] + " based on data from database " + str(item[2]))
        #     # if item[0] in genes and item[1] in genes:
        #     #     print(item[0] + " regulates " + item[1] + " based on data from database " + str(item[2]))
        #     #     found_interaction = True
        if not found_interaction:
            print("Interaction between these two genes was not found.")
        continue_inputting = input("If you want to exit, type E, else press enter: ")
        if continue_inputting == "E":
            break


def start_application(interactions, successors, predecessors, synonyms):
    while True:
        mode = input("Type TF to choose TF and find TGs for it\nTG to choose TG and find"
                     " all TFs that regulate it \ninfo to find info"
                     " about TF and TG\ngroup to find interactions in group of genes\nexit to end the application: ")
        if mode == "TF":
            find_tgs_for_tf(successors, synonyms)
        elif mode == "TG":
            find_tfs_for_tg(predecessors, synonyms)
        elif mode == "info":
            find_info_about_tf_and_tg(interactions, synonyms)
        elif mode == "group":
            find_info_about_group(interactions, synonyms)
        elif mode == "exit":
            break
        else:
            print("Incorrect mode.")


# method to merge interactions

def merge_lists(lists):
    merged_list = []
    for small_list in lists:
        for item in small_list:
            merged_list.append(item)
    return merged_list


# method to merge interaction from all databases that keeps sources

def merge_triples_lists(list_of_lists, filename):
    final_file = open(filename, 'w')
    lists = merge_lists(list_of_lists)
    triple_dict = defaultdict(list)
    for triple in lists:
        key = (triple[0], triple[1])
        if triple[2] not in triple_dict[key]:
            triple_dict[key].append(triple[2])

    for key, values in triple_dict.items():
        values_str = ','.join(values)
        final_file.write(key[0] + ';' + key[1] + ';' + values_str)
        final_file.write('\n')

    final_file.close()
    return final_file


def merge_eight_columns_lists(list_of_lists, filename):
    final_file = open(filename, 'w')
    lists = merge_lists(list_of_lists)
    eight_columns_dict = defaultdict(lambda: [set() for _ in range(6)])

    for row in lists:
        key = tuple(row[:2])  # Assuming the first two columns are the unique identifiers
        values = row[2:]      # The rest of the columns are treated as values for the tuple

        for i, val in enumerate(values):
            if val != '-':
                eight_columns_dict[key][i].add(val)

    sorted_items = sorted(eight_columns_dict.items(), key=lambda x: (x[0][0], x[0][1]))

    for key, sets_list in sorted_items:
        values_with_dashes = [list(s) if s else ['-'] for s in sets_list]
        values_str = ';'.join(','.join(map(str, s)) for s in values_with_dashes)
        final_file.write(key[0] + ';' + key[1] + ';' + values_str)
        final_file.write('\n')

    final_file.close()
    return final_file


# binary searcher for all lines with given tf

def binary_search_table(tf, tg, table_file):
    table_file = open(table_file, 'r')
    table = []
    for line in table_file:
        line = line.split(';')
        table.append(line)
    low, high = 0, len(table) - 1
    start_index, end_index = -1, -1
    while low <= high:
        mid = (low + high) // 2
        current_row = table[mid]

        if current_row[0] == tf:
            start_index = mid
            end_index = mid

            while start_index > 0 and table[start_index - 1][0] == tf:
                start_index -= 1

            while end_index < len(table) - 1 and table[end_index + 1][0] == tf:
                end_index += 1
            break
        elif current_row[0] < tf:
            low = mid + 1
        else:
            high = mid - 1
    for i in range(start_index, end_index + 1):
        if table[i][0] == tf and table[i][1] == tg:
            return table[i]

    return None


if __name__ == '__main__':
    orig_names = all_main_names('Homo_sapiens.gene_info')
    synonyms_1 = read_synonyms('Homo_sapiens.gene_info')
    interactions_1 = read_pathway('homo-sapiens-9606.sif', synonyms_1)
    interactions_2 = read_grndb('AML_TCGA-regulons.txt', synonyms_1)
    interactions_3 = read_grnpedia('trrust_rawdata.human.tsv', synonyms_1)
    interactions_4 = read_maayanlab('ARCHS4_Coexpression.gmt', synonyms_1)
    final_interactions = merge_eight_columns_lists([interactions_1, interactions_2, interactions_3, interactions_4],
                                           'interactions3.csv')
    # final_interactions = 'interactions2.csv'
    # successors_1 = find_successors_interactions(final_interactions)
    # predecessors_1 = find_predecessors_interactions(final_interactions)
    # start_application(final_interactions, successors_1, predecessors_1, synonyms_1)
