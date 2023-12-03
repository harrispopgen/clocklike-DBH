import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np
import pandas as pd
import string

change_species_name = True
raw_header = ["pid", "date", "group", "group2", "name", "species", "populations", "sid", 
              "reproduction", "type", "unit", "u_mean", "u_lower", "u_upper", "se", 
              "range_low", "range_up", "events", "callable", "method", "sample_size", "generation", 
              "sequencing", "tools", "depth", "comments", "citation", "citation2", "link"]
df_rate = pd.read_excel("data/wang_2023_mutation_rate_estimation.xlsx", sheet_name="raw", names=raw_header)

# drop rows where "type" is not of value "snp"
df_rate = df_rate[df_rate["type"] == "snp"]
df_rate = df_rate[df_rate["unit"] == "bp/generation"]
# remove white space from the beginning and end from the following columns:
df_rate['species'] = df_rate['species'].str.strip()
df_rate['name'] = df_rate['name'].str.strip()
df_rate['group'] = df_rate['group'].str.strip()
df_rate['group2'] = df_rate['group2'].str.strip()

# map species name in the excel file to species name in the tree file
species_excel_to_tree_map = {
        "Ailurus fulgens": "Ailurus",
        "Amphilophus": "Amphilophus labiatus",
        "Apis mellifera": "Apis mellifera mellifera",
        "Bombus terrestris": "Bombus",
        "Brassica rapa": "Brassica juncea",
        "Canis lupus": "Canis lupus lupus",
        "Cavia aperea": "Cavia aperea guianae",
        "Ceratotherium simum simum": "Ceratotherium simum",
        "Cervus elaphus yarkandensis": "Cervus hanglu",
        "Cervus nippon": "Cervus nippon centralis",
        "Chironomus riparius": "Chironomus pallidivittatus", # doesn't have diversity estimate
        "Picochlorum costavermella": "Coccomyxa subellipsoidea",
        "Cyanistes caeruleus": "Cyanistes caeruleus caeruleus",
        "Daphnia galeata": "Daphnia dubia", 
        "Emiliania huxleyi": "Emiliania",
        "Gallus gallus domesticus": "Gallus gallus",
        "Giraffa camelopardalis": "Giraffa reticulata",
        "Gyps fulvus": "Gyps",
        "Heliconius melpomene": "Heliconius ethilla", 
        "Homo sapiens": "Homo sapiens neanderthalensis",
        "Hylobates lar": "Hylobates lar lar",
        "Macaca mulatta": "Macaca mulatta vestita",
        "Rhodotorula toruloides": "Microbotryomycetes", 
        "Mus musculus": "Mus musculus musculus",
        "Oryza sativa": "Oryza longistaminata", 
        "Pan troglodytes": "Pan troglodytes troglodytes",
        "Panthera tigris": "Panthera tigris tigris",
        "Pelecanus crispus": "Pelecanus occidentalis", 
        "Phoenicopterus roseus": "Phoenicopterus",
        "Pristionchus pacificus": "Pristionchus",
        "Sphaerodactylus inigo": "Sphaerodactylus copei",
        "Tupaia chinensis belangeri": "Tupaia belangeri",
        "Turdus merula": "Turdus merula merula",
        "Zea mays": "Zea diploperennis"
    }


if change_species_name:
    # change species in the dataframe based on the excel_to_tree_map
    df_rate = df_rate.replace({"species": species_excel_to_tree_map})
# take the average of values from the "u_mean" column for rows with the same "species" value
df_rate_mean = df_rate.groupby(["species", "name", "group2"], as_index=False)["u_mean"].mean()

# drop house mouse and domestic cat from dataset
name_drop_list = ["house mouse", "domestic cat"]
df_rate_mean = df_rate_mean[~df_rate_mean['name'].isin(name_drop_list)]
print(f"there are {len(df_rate_mean)} species with snp mutation rate estimates (bp/generation) in the original dataset")

# load in the generation time and diversity dataframe
gsgt_header = ["species", "genome_size_mb", "generation_time_year", "log10_popsize", "diversity"]
df_gen_time = pd.read_excel("data/Yiguan_mutation_rate_estimation.xlsx", sheet_name="gs_gt", names=gsgt_header)
df_gen_time = df_gen_time[["species", "genome_size_mb", "generation_time_year", "diversity"]]
df_gen_time['species'] = df_gen_time['species'].str.strip()

if change_species_name:
    df_gen_time = df_gen_time.replace({"species": species_excel_to_tree_map})

# combine df_rate_mean and df_gen_time with matching species names
df_merge = df_rate_mean.merge(df_gen_time, on="species")
# get rid of rows where diversity is 0 or NA
df_merge = df_merge[df_merge["diversity"] > 0]
df_merge["Ne"] = df_merge["diversity"] / (4 * df_merge["u_mean"])
df_merge["Ne*gen"] = df_merge["Ne"] * df_merge["generation_time_year"]
df_merge["u_mean_year"] = df_merge["u_mean"] / df_merge["generation_time_year"]
print(f"there are {len(df_merge)} species with snp mutation rate estimates AND Ne estimates in the dataset")

df_merge.to_csv("output/processed_mutation_rate_estimate.csv", index=False)
