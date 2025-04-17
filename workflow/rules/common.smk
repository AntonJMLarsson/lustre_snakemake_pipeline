def get_final_output(wildcards):
    final_output = []
    for donor in config["custom_references"].keys():
        final_output.extend(["insertion_tables/{project}_KNR_{donor}_insertion_table.csv".format(config["project"], donor), "insertion_tables/{project}_UNK_{donor}_insertion_table.csv".format(config["project"], donor)])
    return final_output