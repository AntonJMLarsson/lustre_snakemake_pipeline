def get_final_output(wildcards):
    final_output = []
    for donor in config["custom_references"].keys():
        final_output.extend(["insertion_tables/{}_KNR_{}_insertion_table.csv".format(config["project"], donor), "insertion_tables/{}_UNK_{}_insertion_table.csv".format(config["project"], donor)])
    return final_output