def get_final_output(wildcards):
    final_output = ["results/KR.done"]
    for donor in config["custom_references"].keys():
        final_output.extend(["results/insertion_tables/{}_KNR_{}_insertion_table.csv".format(config["project"], donor), "results/insertion_tables/{}_UNK_{}_insertion_table.csv".format(config["project"], donor)])
    return final_output