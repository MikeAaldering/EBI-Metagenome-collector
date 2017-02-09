# -*- coding: utf-8 -*-
"""
Author: Mike Aaldering
Company: NIOO Knaw
Email: M.Aaldering2@nioo.knaw.nl
Private Email: mike.aaldering098@gmail.com

Search metadata and species from EBI metagenomics and store it to a .csv file.
1       Download metadata on ENA.
    1.1 Download run, sample and project metadata on EBI metagenomics.
2       Search article id's for every sample on pubmed central and add article id to metadata.
3       Download UTO table for every sample.
4       Save  the metadata and OTU table to a seperate CSV file.
5       Normalize metadata.
6       Join OTU table and metadata and save to .csv file.
"""

from urllib2 import URLError
import urllib, urllib2, json, re, sys, ctypes, xmltodict, requests
import xml.etree.ElementTree as ET
from Bio import Entrez
import pandas as pd
import unicodedata


def main():
    only_normalize = raw_input('Only normalize type: "n", \ndownload data type: "d" \n')
    CSVfile_name = "EBI"
    query = """\"soil\"||\"peat\"||\"sediment\"||\"wetland\""""
    only_species_filename = CSVfile_name + " Only Species.txt"
    only_metadata_filename = CSVfile_name + " Only Metadata.txt"
    normalize = True
    relative = True
    fieldname_order_filename = "Metadata_order.txt"
    fieldname_order = order_fielnames(fieldname_order_filename)

    if only_normalize == "n":
        print "normalizing data"
        join_data(only_species_filename, only_metadata_filename, normalize, fieldname_order, CSVfile_name + " Only Normalized")
    elif only_normalize == "d":
        search_pubmed = False                                   # Check if sample is used in a Article: this slows down the script!
        only_interesting_species = True                         # False = save only species from the list, True = make list with all the species
        interestingspecies_list = get_species_of_interest()     # Make a list with species where you interested in.

        #Check how much runs are available om the EBI metagenomics database
        url = make_run_url(1, 1, query)
        print url
        run_data_dict = get_run_data_ebi(url)
        run_data_dict["run_metadata_url"] = url
        total_runs = run_data_dict["hitCount"]
        print "total runs: ", run_data_dict["hitCount"]


        succes_amount = 0
        download_amount = 100
        all_species = []
        all_metadata = []
        while (succes_amount <= int(total_runs)):
            print "downloading", succes_amount, " - ", succes_amount + download_amount, " from total ", total_runs
            url = make_run_url(download_amount, succes_amount, query)
            run_data_dict = get_run_data_ebi(url)
            for idx, run in enumerate(run_data_dict["entries"]):
                ###############################DOWNLOAD ALL METADATA#######################################
                print run
                sample_data_dict = get_sample_metadata_ebi(                                 # Sample metadata.
                    run["fields"]['METAGENOMICS_SAMPLES'][0])
                metadata_dict = merge_two_dicts(run, sample_data_dict)

                project_data_dict = get_project_metadata(
                    sample_data_dict["entries"][0]['fields']['METAGENOMICS_PROJECTS'][0])   # Project metadata.
                metadata_dict = merge_two_dicts(metadata_dict, project_data_dict)
                metadata_dict = flatten_dict(metadata_dict)                                 # Flatten dict: dict in dict in dict to flatten dict.

                # set all ids to a string
                run_id = metadata_dict['fields_name_0']
                sample_id = metadata_dict['fields_METAGENOMICS_SAMPLES_0']
                project_id = metadata_dict['entries_0_id']
                pipeline_version = metadata_dict['fields_pipeline_version_0']

                ena_metadata = get_sample_metadata_ena(sample_id)                           # ENA project metadata.
                metadata_dict = merge_two_dicts(ena_metadata, metadata_dict)
                metadata_dict = get_experiment_metadata_ena(metadata_dict, run_id)          # ENA experiment metadata.


                metadata_dict = get_pubmed_metadata_ena(project_id, metadata_dict)          # Pubmed id.
                if search_pubmed == True:
                    pmc_ids = search_PMC_article_ids(project_id, sample_id, run_id)         # Check if there are articles available for the project in pubmed central.
                    metadata_dict["PMC_article_ids"] = pmc_ids


                ################################DOWNLOAD OTU TABLE##########################
                OTU_dict = Get_OTU_species(run_id, sample_id, project_id, pipeline_version)    # Download the OTU table.

                total_reads = OTU_dict[1]
                OTU_dict = OTU_dict[0]
                if only_interesting_species == True:
                    OTU_dict = get_only_interesting_species(OTU_dict, interestingspecies_list)

                OTU_dict = clean_otu_dict(OTU_dict)
                if relative == True:
                    OTU_dict = make_relative(OTU_dict, pipeline_version, total_reads)

                #################CLEAN METADATA######################
                metadata_dict = delete_new_lines(metadata_dict)

                #####################WRITING ALL DATA TO FILE##########################
                all_species.append(OTU_dict)
                all_metadata.append(metadata_dict)
                write_csv(all_species, only_species_filename, "0")
                write_csv(all_metadata, only_metadata_filename, "NA")

                join_data(only_species_filename, only_metadata_filename, normalize, fieldname_order, CSVfile_name)


                print succes_amount+idx, " / ", total_runs, " run id: ", run_id
            succes_amount = succes_amount + download_amount
        print "Done!, filename: ", CSVfile_name
    else:
        print only_normalize, "is not validate, try it again"
        main()

def join_data(only_species_filename, only_metadata_filename, normalize, fieldname_order, CSVfile_name):
    df_species = pd.read_csv(only_species_filename, error_bad_lines=False, sep='\t')
    df_metadata = pd.read_csv(only_metadata_filename, error_bad_lines=False, sep='\t')
    MetaSpecies = pd.concat([df_metadata, df_species], axis=1, join_axes=[df_metadata.index])
    MetaSpecies_list_with_dict = MetaSpecies.T.to_dict().values()
    all_data_list_with_dict = make_list_witch_dict_without_nan(MetaSpecies_list_with_dict)

    ######NORMALIZE DATA#######
    if normalize == True:
        final_all_data_list_with_dict = []
        for row in all_data_list_with_dict:
            row_dict = normalize_data(row)
            final_all_data_list_with_dict.append(row_dict)
    else:
        final_all_data_list_with_dict = all_data_list_with_dict

    write_to_csv(CSVfile_name, final_all_data_list_with_dict, fieldname_order)

def write_csv(data, filename, fill):
    DataFrame = pd.DataFrame(data)
    DataFrame = DataFrame.fillna(value=fill)
    DataFrame.to_csv(filename, sep='\t', encoding="utf-8", index=False)


def make_list_witch_dict_without_nan(list_with_dict):
    all_data = []
    for row in list_with_dict:
        row_dict = {}
        for data in row:
            if str(row[data]) != "nan":
                row_dict[data] = row[data]
        all_data.append(row_dict)
    return all_data

def delete_new_lines(meta_OTU):
            for datapoint in meta_OTU:
                try:
                    if isinstance(datapoint, unicode):
                        new_datapoint = unicodedata.normalize('NFKD', datapoint).encode('ascii', 'ignore')
                    else:
                        new_datapoint = datapoint
                    meta_OTU[new_datapoint] = meta_OTU.pop(datapoint)
                    meta_OTU[new_datapoint] = re.sub('\n', '', str(meta_OTU[new_datapoint]))

                except UnicodeEncodeError:
                    pass
            return meta_OTU

def order_data(header_list, fieldname_order):
    new_list = []
    for item in fieldname_order:
        if item in header_list:
            new_list.append(item)
            # fieldname_order = fieldname_order.remove(item)
    for deze in header_list:
        if deze not in new_list:
            new_list.append(deze)
    return new_list

def order_fielnames(filename):
    fieldnames = []
    with open(filename, "r") as f:
        for row in f:
            dataname = row.rstrip('\n')
            # print dataname
            if dataname not in fieldnames:
                fieldnames.append(dataname)
    return fieldnames

def write_to_csv(CSVfile_name, all_data, fieldname_order):
    # print all_data
    DataFrame = pd.DataFrame(all_data)

    header_list = list(DataFrame.columns.values)
    data_order = order_data(header_list, fieldname_order)
    DataFrame = DataFrame[data_order]
    #print DataFrame

    DataFrame = DataFrame.fillna(value="NA")
    all_column = []
    for column_name in DataFrame.columns.values:
        column_names = {}
        column_names["Orig"] = column_name
        #print column_name
        if "k__" in column_name:
            name_split = column_name.split(';')
            for level in name_split:
                level = level.replace(" ", "")
                if "k__" in level:
                    level = level.replace("k__", "")
                    column_names["Kingdom"] = level
                if "p__" in level:
                    level = level.replace("p__", "")
                    if level == "":
                        break
                    else:
                        column_names["Phylum"] = level

                if "c__" in level:
                    level = level.replace("c__", "")
                    if level == "":
                        break
                    else:
                        column_names["Class"] = level

                if "o__" in level:
                    level = level.replace("o__", "")
                    if level == "":
                        break
                    else:
                        column_names["Order"] = level
                if "f__" in level:
                    level = level.replace("f__", "")
                    if level == "":
                        break
                    else:
                        column_names["Family"] = level

                if "g__" in level:
                    level = level.replace("g__", "")
                    if level == "":
                        break
                    else:
                        column_names["Genus"] = level

                if "s__" in level:
                    level = level.replace("s__", "")
                    if level == "":
                        break
                    else:
                        column_names["Species"] = column_names["Genus"] + "_" + level
            #print column_names

            if "Phylum" not in column_names:
                column_names["Phylum"] = column_names["Kingdom"] + "_unclassified"
                column_names["Class"] = column_names["Kingdom"] + "_unclassified"
                column_names["Order"] = column_names["Kingdom"] + "_unclassified"
                column_names["Family"] = column_names["Kingdom"] + "_unclassified"
                column_names["Genus"] = column_names["Kingdom"] + "_unclassified"
                column_names["Species"] = column_names["Kingdom"] + "_unclassified"

            elif "Class" not in column_names:
                column_names["Class"] = column_names["Phylum"] + "_unclassified"
                column_names["Order"] = column_names["Phylum"] + "_unclassified"
                column_names["Family"] = column_names["Phylum"] + "_unclassified"
                column_names["Genus"] = column_names["Phylum"] + "_unclassified"
                column_names["Species"] = column_names["Phylum"] + "_unclassified"

            elif "Order" not in column_names:
                column_names["Order"] = column_names["Class"] + "_unclassified"
                column_names["Family"] = column_names["Class"] + "_unclassified"
                column_names["Genus"] = column_names["Class"] + "_unclassified"
                column_names["Species"] = column_names["Class"] + "_unclassified"

            elif "Family" not in column_names:
                column_names["Family"] = column_names["Order"] + "_unclassified"
                column_names["Genus"] = column_names["Order"] + "_unclassified"
                column_names["Species"] = column_names["Order"] + "_unclassified"

            elif "Genus" not in column_names:
                column_names["Genus"] = column_names["Family"] + "_unclassified"
                column_names["Species"] = column_names["Family"] + "_unclassified"

            elif "Species" not in column_names:
                column_names["Species"] = column_names["Genus"] + "_unclassified"

        all_column.append(column_names)

    DataFrame2 = pd.DataFrame(all_column)

    # print DataFrame2
    new_col_names = ["Orig", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
    final_col_names = []
    for col_name in new_col_names:
        if col_name in DataFrame2.columns.values:
            final_col_names.append(col_name)
    DataFrame2 = DataFrame2[final_col_names]
    DataFrame2 = DataFrame2.set_index(["Orig"])
    DataFrame2 = DataFrame2.T
    DataFrame2 = DataFrame2.append(DataFrame)

    #print "header list voor temperatuur", header_list
    DataFrame2.to_csv(CSVfile_name + ".txt", sep='\t', encoding="utf-8")
    # print DataFrame2


def get_experiment_metadata_ena(metadata_dict, run_id):
    url = "https://www.ebi.ac.uk/ena/data/warehouse/filereport?accession=" + run_id + "&result=read_run&fields=instrument_model"
    model = urllib.urlopen(url)
    for item in model:
        #print item
        metadata_dict["instrument_model"] = item
    #print metadata_dict["instrument_model"]
    return metadata_dict

def normalize_data(meta_OTU):
    #ID
    if "extreme_unusual_properties/heavy metals" in meta_OTU:
        del meta_OTU["extreme_unusual_properties/heavy metals"]
    if "id" in meta_OTU:
        if "," in meta_OTU["id"]:
            value = meta_OTU["id"]
            value = value.replace(" ", "")
            value = value.replace(",", "_")
            meta_OTU["id"] = value
    #print meta_OTU
    if "fields_temperature_0" in meta_OTU:
        meta_OTU["Temp_std"] = meta_OTU["fields_temperature_0"]
    if "temperature (ºC)" in meta_OTU:
        meta_OTU["Temp_std"] = meta_OTU["Temperature (ºC)"]
    if "temp (°C)" in meta_OTU:
        meta_OTU["Temp_std"] = meta_OTU["temp (°C)"]
    if "temp" in meta_OTU:
        meta_OTU["Temp_std"] = meta_OTU["temp"]
    if "temperature (ºC)" in meta_OTU:
        meta_OTU["Temp_std"] = meta_OTU["temperature (ºC)"]
    if "Temperature (\\xbaC)" in meta_OTU:
        meta_OTU["Temp_std"] = meta_OTU["Temperature (\\xbaC)"]
    if "Temperature (ÂºC)" in meta_OTU:
        meta_OTU["Temp_std"] = meta_OTU["Temperature (ÂºC)"]
    if "temp (oC)" in meta_OTU:
        meta_OTU["Temp_std"] = meta_OTU["temp (oC)"]
    if "temperature (oC)" in meta_OTU:
        meta_OTU["Temp_std"] = meta_OTU["temperature (oC)"]
    if "Temperature (oC)" in meta_OTU:
        meta_OTU["Temp_std"] = meta_OTU["Temperature (oC)"]
    if "temp (C)" in meta_OTU:
        meta_OTU["Temp_std"] = meta_OTU["temp (C)"]
    #environment (biome), env_biome
    if "env_biome" in meta_OTU:
        meta_OTU["biome_std"] = meta_OTU["env_biome"]
    if "biome" in meta_OTU:
        meta_OTU["biome_std"] = meta_OTU["biome"]
    if "environment (biome)" in meta_OTU:
        meta_OTU["biome_std"] = meta_OTU["environment (biome)"]

    if "biome_std" in meta_OTU:
        if "ENVO" in meta_OTU["biome_std"]:
            meta_OTU["biome_std"] = meta_OTU["biome_std"].replace("ENVO:", "")
        if "447" in meta_OTU["biome_std"]:
            meta_OTU["biome_std"] = "marine"
        elif "soybean" in meta_OTU["biome_std"]:
            meta_OTU["biome_std"] = "soil"
        elif "Soil Sample" in meta_OTU["biome_std"]:
            meta_OTU["biome_std"] = "soil"
        elif "Zoig wetland-associated" in meta_OTU["biome_std"]:
            meta_OTU["biome_std"] = "wetland"
        elif "39" in meta_OTU["biome_std"]:
            meta_OTU["biome_std"] = "marine"
        elif "39" in meta_OTU["biome_std"]:
            meta_OTU["biome_std"] = "marine"
        elif "126" in meta_OTU["biome_std"]:
            meta_OTU["biome_std"] = "sediment"
        elif "043" in meta_OTU["biome_std"]:
            meta_OTU["biome_std"] = "marine"
        elif "Desert" in meta_OTU["biome_std"]:
            meta_OTU["biome_std"] = "desert"
        elif "desert" in meta_OTU["biome_std"]:
            meta_OTU["biome_std"] = "desert"
        elif "terrestrial" in meta_OTU["biome_std"]:
            meta_OTU["biome_std"] = "terrestrial"
        elif "Terrestrial" in meta_OTU["biome_std"]:
            meta_OTU["biome_std"] = "terrestrial"
        elif "woodland" in meta_OTU["biome_std"]:
            meta_OTU["biome_std"] = "woodland"
        elif "marine" in meta_OTU["biome_std"]:
            meta_OTU["biome_std"] = "marine"
        elif "ocean" in meta_OTU["biome_std"]:
            meta_OTU["biome_std"] = "ocean"
        elif "ocean" in meta_OTU["biome_std"]:
            meta_OTU["biome_std"] = "ocean"
        elif "lake" in meta_OTU["biome_std"]:
            meta_OTU["biome_std"] = "lake"
        elif meta_OTU["biome_std"] == "Temperate":
            meta_OTU["biome_std"] = "soil"
        elif "grasslands" in meta_OTU["biome_std"]:
            meta_OTU["biome_std"] = "grassland"
        elif "Tropical humid" in meta_OTU["biome_std"]:
            meta_OTU["biome_std"] = "forest"
        elif "unknown" in meta_OTU["biome_std"]:
            meta_OTU["biome_std"] = "NA"
        elif "urban" in meta_OTU["biome_std"]:
            meta_OTU["biome_std"] = "urban"
        elif "wastewater treatment plant" in meta_OTU["biome_std"]:
            meta_OTU["biome_std"] = "soil"
        elif "wastewater treatment plant" in meta_OTU["biome_std"]:
            meta_OTU["biome_std"] = "soil"
        elif meta_OTU["biome_std"] == "water":
            meta_OTU["biome_std"] = "marine"
        elif "01000047" in meta_OTU["biome_std"]:
            meta_OTU["biome_std"] = "marine"
        elif "00000446" in meta_OTU["biome_std"]:
            meta_OTU["biome_std"] = "marine"
        if " " in meta_OTU["biome_std"]:
            meta_OTU["biome_std"] = meta_OTU["biome_std"].replace(" ", "_")
        meta_OTU["biome_std"] = meta_OTU["biome_std"].lower()
    #environment (material),env_material,material
    if "env_material" in meta_OTU:
        meta_OTU["material_std"] = meta_OTU["env_material"]
    if "material" in meta_OTU:
        meta_OTU["material_std"] = meta_OTU["material"]
    if "environment (material)" in meta_OTU:
        meta_OTU["material_std"] = meta_OTU["environment (material)"]
    #environment (feature), env_feature,feature
    if "env_feature" in meta_OTU:
        meta_OTU["feature_std"] = meta_OTU["env_feature"]
    if "feature" in meta_OTU:
        meta_OTU["feature_std"] = meta_OTU["feature"]
    if "environment (feature)" in meta_OTU:
        meta_OTU["feature_std"] = meta_OTU["environment (feature)"]
    #geographic location (longitude) (DD), Longitude Start (DD), longitude, geographic location (longitude)
    if "Longitude Start (DD)" in meta_OTU:
        if "Latitude Start (DD)" in meta_OTU:
            meta_OTU["Lat_Lon_std"] = str(meta_OTU["Latitude Start (DD)"] ) + " , " + str(meta_OTU["Longitude Start (DD)"])
    if "longitude" in meta_OTU:
        if "latitude" in meta_OTU:
            meta_OTU["Lat_Lon_std"] = str(meta_OTU["latitude"]) + " , " + str(meta_OTU["longitude"])
    if "geographic location (latitude and longitude)" in meta_OTU:
        meta_OTU["Lat_Lon_std"] = meta_OTU["geographic location (latitude and longitude)"]
    if "lat_lon" in meta_OTU:
        meta_OTU["Lat_Lon_std"] = meta_OTU["lat_lon"]
    if "geographic location (longitude)" in meta_OTU:
        meta_OTU["geographic location (longitude) (DD)"] = meta_OTU["geographic location (longitude)"]
    if "geographic location (latitude)" in meta_OTU:
        if "geographic location (longitude)" in meta_OTU:
            meta_OTU["Lat_Lon_std"] = str(meta_OTU["geographic location (latitude)"]) + " , " + str(meta_OTU["geographic location (longitude)"])
    if "geographic location (latitude) (DD)" in meta_OTU:
        if "geographic location (longitude) (DD)" in meta_OTU:
            meta_OTU["Lat_Lon_std"] = str(meta_OTU["geographic location (latitude) (DD)"]) + " , " + str(meta_OTU["geographic location (longitude) (DD)"])
    #altitude
    if "altitude" in meta_OTU:
        meta_OTU["altitude_std"] = re.sub("[^\d,.]", "", str(meta_OTU["altitude"]))
        #print "altitude_std", meta_OTU["altitude"]
    #geographic location (elevation) (m), elev, elevation
    if "elevation" in meta_OTU:
        if "cm" in str(meta_OTU["elevation"]):
            #print meta_OTU["elevation"]
            meta_OTU["elevation"] = re.sub("[^\d,.]", "", str(meta_OTU["elevation"]))
            #print "elevalion:", meta_OTU["elevation"]
            meta_OTU["elevation_std"] = float(meta_OTU["elevation"]) / 100
        else:
            meta_OTU["elevation_std"] = re.sub("[^\d,.]", "", str(meta_OTU["elevation"]))
    if "geographic location (elevation) (m)" in meta_OTU:
        meta_OTU["elevation_std"] = meta_OTU["geographic location (elevation) (m)"]
    if "elev" in meta_OTU:
        #print meta_OTU
        #print meta_OTU['elev']
        if ("ft" in str(meta_OTU["elev"]) or "FT" in str(meta_OTU["elev"])):
            meta_OTU["elev"] = re.sub("[^\d,.]", "", meta_OTU["elev"])
            foot = meta_OTU["elev"]
            meter = float(foot) * 0.3048
            meta_OTU["elevation_std"] = meter
            #print "normalisation: ", foot, " foot to ", meter, " meter"
            #dprogramPause = raw_input("Press the <ENTER> key to continue...")
        elif "m.a.s.l." in str(meta_OTU["elev"]):
            meta_OTU["elevation_std"] = re.sub("m.a.s.l.", "", meta_OTU["elev"])
        elif "-" in str(meta_OTU["elev"]):
            meta_OTU["elev"] = re.sub("[^\d.,-]", "", meta_OTU["elev"])
            two_digit_list = meta_OTU["elev"].split("-")
            average = (float(two_digit_list[1])-float(two_digit_list[0]))/2+float(two_digit_list[0])
            meta_OTU["elevation_std"] = average
            #print "elev:", meta_OTU["elev"]
        else:
            #print "elev:", meta_OTU["elev"]
            meta_OTU["elevation_std"] = re.sub("[^\d,.]", "", str(meta_OTU["elev"]))
    #geographic location (depth) (m),Depth (m),depth,fields_depth_0
    if "depth" in meta_OTU:
        if "-" in str(meta_OTU["depth"]):
            #print meta_OTU["depth"]
            depth = re.sub("[^\d.,-]", "", meta_OTU["depth"])
            #print meta_OTU["depth"]
            two_digit_list= depth.split("-")
            average = (float(two_digit_list[1]) - float(two_digit_list[0])) / 2 + float(two_digit_list[0])
            meta_OTU["depth_std"] = average
        elif "cm" in str(meta_OTU["depth"]):
            depth = re.sub("[^\d,.]", "", meta_OTU["depth"])
            meta_OTU["depth_std"] = float(depth) / 100
        else:
            meta_OTU["depth_std"] = meta_OTU["depth"]
    if "Depth (m)" in meta_OTU:
        if "-" in str(meta_OTU["Depth (m)"]):
            two_digit_list= meta_OTU["Depth (m)"].split("-")
            average = (float(two_digit_list[1]) - float(two_digit_list[0])) / 2 + float(two_digit_list[0])
            meta_OTU["depth_std"] = average
    if "depth (m)" in meta_OTU:
         meta_OTU["depth_std"] = meta_OTU["depth (m)"]
    if "fields_depth_0" in meta_OTU:
        meta_OTU["depth_std"] = meta_OTU["fields_depth_0"]

    if "depth_std" in meta_OTU:
        value = meta_OTU["depth_std"]
        if "cm" in str(meta_OTU["depth_std"]):
            value = meta_OTU["depth_std"].replace("cm", "")
            value = float(value) / 100
            meta_OTU["depth_std"] = value
        elif "m" in str(meta_OTU["depth_std"]):
            meta_OTU["depth_std"] = value.replace("m", "")
        elif "metagenome" in str(meta_OTU["depth_std"]):
            meta_OTU["depth_std"] = "NA"
        elif "Not" in str(meta_OTU["depth_std"]):
            meta_OTU["depth_std"] = "NA"
    #geographic location, geo_loc_name, geographic location (country and/or sea), geographic location (country and/or sea,region),geographic location (country:region,area)
    if "geo_loc_name" in meta_OTU:
        meta_OTU["geographic_location_std"] = meta_OTU["geo_loc_name"]
    if "geographic location (country and/or sea)" in meta_OTU:
        meta_OTU["geographic_location_std"] = meta_OTU["geographic location (country and/or sea)"]
    if "geographic location (country and/or sea,region)" in meta_OTU:
        meta_OTU["geographic_location_std"] = meta_OTU["geographic location (country and/or sea,region)"]
    if "geographic location (country:region,area)" in meta_OTU:
        meta_OTU["geographic_location_std"] = meta_OTU["geographic location (country:region,area)"]
    if "geographic location" in meta_OTU:
        meta_OTU["geographic_location_std"] = meta_OTU["geographic location"]

    if "geographic_location_std" in meta_OTU:
        value = meta_OTU["geographic_location_std"]
        if "Adelaide, South Australia, Australia" in meta_OTU["geographic_location_std"]:
            meta_OTU["geographic_location_std"] = "Australia"
        if "Argentina" in meta_OTU["geographic_location_std"]:
            meta_OTU["geographic_location_std"] = "Argentina"
        if "Australia: Sydney" in meta_OTU["geographic_location_std"]:
            meta_OTU["geographic_location_std"] = "Australia"
        if "Canada: Axel Heiberg Island" in meta_OTU["geographic_location_std"]:
            meta_OTU["geographic_location_std"] = "Canada"
        if "China" in meta_OTU["geographic_location_std"]:
            meta_OTU["geographic_location_std"] = "China"
        if "France" in meta_OTU["geographic_location_std"]:
            meta_OTU["geographic_location_std"] = "France"
        if "United States of America" in meta_OTU["geographic_location_std"]:
            meta_OTU["geographic_location_std"] = "USA"
        if "India" in meta_OTU["geographic_location_std"]:
            meta_OTU["geographic_location_std"] = "India"
        if "Japan" in meta_OTU["geographic_location_std"]:
            meta_OTU["geographic_location_std"] = "Japan"
        if "China" in meta_OTU["geographic_location_std"]:
            meta_OTU["geographic_location_std"] = "China"
        if "Germany" in meta_OTU["geographic_location_std"]:
            meta_OTU["geographic_location_std"] = "Germany"
        if "Spain" in meta_OTU["geographic_location_std"]:
            meta_OTU["geographic_location_std"] = "Spain"
        if "unknown" in meta_OTU["geographic_location_std"]:
            meta_OTU["geographic_location_std"] = "NA"
        if "USA" in meta_OTU["geographic_location_std"]:
            meta_OTU["geographic_location_std"] = "USA"
        if "not applicable" in meta_OTU["geographic_location_std"]:
            meta_OTU["geographic_location_std"] = "NA"
        if "87826" in meta_OTU["geographic_location_std"]:
            meta_OTU["geographic_location_std"] = "NA"
        if "ddd" in meta_OTU["geographic_location_std"]:
            meta_OTU["geographic_location_std"] = "NA"
        meta_OTU["geographic_location_std"] = meta_OTU["geographic_location_std"].replace(" ", "_")
    #collection date, collection_date
    if "collection_date" in meta_OTU:
        meta_OTU["collection_date_std"] = meta_OTU["collection_date"]
    if "collection date" in meta_OTU:
        meta_OTU["collection_date_std"] = meta_OTU["collection date"]
    #percent_org_carb, total organic carbon (g/kg)-> perc, total organic carbon -> /10
    if "total organic carbon (g/kg)" in meta_OTU:
        meta_OTU["C_percent_std"] = float(meta_OTU["total organic carbon (g/kg)"])/10
    if "total organic carbon" in meta_OTU:
        meta_OTU["C_percent_std"] = float(meta_OTU["total organic carbon"]) / 10
    if "percent_org_carb" in meta_OTU:
        meta_OTU["C_percent_std"] = meta_OTU["percent_org_carb"]
    if "tot_org_carb" in meta_OTU:
        if "%" in str(meta_OTU["tot_org_carb"]):
            new = re.sub("%", "", str(meta_OTU["tot_org_carb"]))
            new = re.sub(",", ".", str(new))
            meta_OTU["C_percent_std"] = new
        else:
            meta_OTU["C_percent_std"] = float(meta_OTU["tot_org_carb"]) / 10
    if "tot_nitro" in meta_OTU:
        meta_OTU["N_percent_std"] = meta_OTU["tot_nitro"]
    if "tot_n" in meta_OTU:
        nitrogen = meta_OTU["tot_n"]
        if "%" in str(nitrogen):
            nitrogen = nitrogen.replace("%", "")
            nitrogen = nitrogen.replace(",", ".")
        if meta_OTU["entries_0_id"] == "SRP041174":
            nitrogen = float(nitrogen) / 10
        meta_OTU["N_percent_std"] = nitrogen
    # ph level (h2o), ph --> ph_std
    if "ph level (h2o)" in meta_OTU:
        meta_OTU["pH_std"] = meta_OTU["ph level (h2o)"]
    if "ph" in meta_OTU:
        if "," in str(meta_OTU["ph"]):
            meta_OTU["ph"] = re.sub(",", ".", meta_OTU["ph"])
        #print meta_OTU["ph"], "phphphphhphphphphpphphphphphphphphphphphphphphphphphphpphph"
        ph = meta_OTU["ph"]
        meta_OTU["pH_std"] = ph

    # Sequencing method
    if "instrument_model" in meta_OTU:
        meth = meta_OTU["instrument_model"]
        if "454" in str(meth):
            meta_OTU["seq_meth_std"] = "pyrosequencing"
        elif "Illumina Genome Analyzer II" in str(meth):
            meta_OTU["seq_meth_std"] = "Illumina_Genome_Analyzer_II"
        elif "HiSeq" in str(meth):
            meta_OTU["seq_meth_std"] = "Illumina_HiSeq"
        elif "MiSeq" in str(meth):
            meta_OTU["seq_meth_std"] = "Illumina_MiSeq"
        elif "instrument_model" in str(meth):
            if "sequencing method" in meta_OTU:
                if "Hi" in meta_OTU["sequencing method"]:
                    meta_OTU["seq_meth_std"] = "Illumina_HiSeq"
                elif "Mi" in meta_OTU["sequencing method"]:
                    meta_OTU["seq_meth_std"] = "Illumina_MiSeq"
            else:
                meta_OTU["seq_meth_std"] = "NA"
        elif "Ion Torrent PGM" in str(meth):
            meta_OTU["seq_meth_std"] = "Ion_Torrent"
        elif "NextSeq 500" in str(meth):
            meta_OTU["seq_meth_std"] = "Illumina_NextSeq_500"
        else:
            meta_OTU["seq_meth_std"] = "NA"
    if "C_percent_std" and "N_percent_std" in meta_OTU:
        meta_OTU["C_N_ratio_std"] = float(meta_OTU["C_percent_std"]) / float(meta_OTU["N_percent_std"])
    return meta_OTU


def search_PMC_article_ids(project_id, sample_id, run_id):
    search_term = project_id + " OR " + sample_id + " OR " + run_id
    Entrez.email = "M.Aaldering2@nioo.knaw.nl"
    handle = Entrez.esearch(db="PMC", retmax=100, term=search_term)
    record = Entrez.read(handle)
    article_ids = ""
    for articleid in record["IdList"]:

        article_ids = article_ids + articleid + " "
    return article_ids



def Get_OTU_species(run_id, sample_id, project_id, version):
    url = 'https://www.ebi.ac.uk/metagenomics//projects/' + project_id + '/samples/' + sample_id + '/runs/' + run_id + '/results/versions/' + version + '/taxonomy/OTU-TSV'
    otudict = {}
    total_reads = 0
    try:
        # Pipeline 1.0 gives uses another method, check the EBI metagenomics documentation.
        if version == "1.0":
            otudict["OTU_URL"] = url
            for item in urllib.urlopen(url):
                x = item.rstrip("\n").split('\t')
                otudict[x[1]] = x[2]
                #print "version 1.0 uses another pipeline, you don't want to compare between different pipelines!"

        # Pipeline 2.0 and 3.0 are the same for OTU prediction.
        else:
            for item in urllib.urlopen(url):
                #print item
                if "#" in item:
                    pass
                    # print "No species info line: ", item
                else:
                    x = item.rstrip("\n").split('\t')
                    reads = x[1].split('.')
                    no_decimal_reads = reads[0]
                    otudict[x[2]] = no_decimal_reads
            for count in otudict:

                total_reads = total_reads + float(otudict[count])
        return otudict, total_reads
    except URLError as url_error:
        print(url_error)
        raise
    except  IOError as io_error:
        print(io_error)
        raise
        # save all species or only species from a list (list = interestingspecies_list)

def make_relative(otudict, version, total_reads):
    if version == "1.0":
        return otudict
    else:
        total_reads_species_of_interest = 0
        for interest in otudict:
            total_reads_species_of_interest = total_reads_species_of_interest + float(otudict[interest])
        if total_reads_species_of_interest != 0:
            otudict["total_reads_species_of_interest"] = total_reads_species_of_interest
            for rel in otudict:
                otudict[rel] = (float(otudict[rel]) / total_reads_species_of_interest)
        otudict["total_reads_species_of_interest"] = total_reads_species_of_interest
        otudict["total_reads"] = total_reads
        return otudict


def clean_otu_dict(otudict):
        for species_name in otudict:
            new_species_name = species_name.replace("Root;", "")
            otudict[new_species_name] = otudict.pop(species_name)
            if new_species_name not in otudict:
                otudict.append(new_species_name)
        to_remove = []
        for name in otudict:
            if "f__Methylobacteriaceae; g__Methylobacterium; s__" in name or "k__Bacteria; p__Proteobacteria; c__Alphaproteobacteria; o__Rhizobiales;" in name:
                to_remove.append(name)
        for remove in to_remove:
            del otudict[remove]
        return otudict


def get_only_interesting_species(species_dict, interestingspecies_list):
    keys = []
    # print species_dict
    for key in species_dict:
        for species in interestingspecies_list:
            keylow = key.lower()
            specieslow = species.lower()
            # print key, species
            if specieslow in keylow:
                # print specieslow, keylow
                # print keylow, specieslow
                keys.append(key)
                # print keys
    species_dict = {x: species_dict[x] for x in keys}
    # print species_dict
    return species_dict


def get_species_of_interest():
    file = "InterestingSpecies.txt"
    with open(file, "r") as f:
        specieslist = []
        for species in f:
            nonewlinespecies = species.rstrip('\n')
            specieslist.append(nonewlinespecies)
            firstword = nonewlinespecies.split(' ')
            specieslist.append(firstword[0])
        specieslist = list(set(specieslist))
        # print specieslist
        return specieslist


def get_project_metadata(project_id):
    # allfields = METAGENOMICS_RUNS,METAGENOMICS_SAMPLES,biome,biome_name,centre_name,creation_date,description,domain_source,id,last_modification_date,name,releaseDate_date
    fields = "biome,biome_name,centre_name,creation_date,description,domain_source,id,last_modification_date,name,releaseDate_date"
    url = "http://www.ebi.ac.uk/ebisearch/ws/rest/metagenomics_projects?format=json&query=" + project_id + "&fields=" + fields
    # print "project_url: ", url
    project_data_dict = get_run_data_ebi(url)
    project_data_dict["project_metadata_url"] = url
    # print "project_data_dict: ", project_data_dict
    return project_data_dict


def get_sample_metadata_ebi(sample_id):
    # allfields =  GO,INTERPRO,METAGENOMICS_PROJECTS,METAGENOMICS_SAMPLES,altitude,biome,completion_date,creation_date,depth,domain_source,experiment_type,go_term,id,interpro_entry,name,organism,pH,pipeline_version,project_name,sample_name,temperature
    fields = "GO,INTERPRO,METAGENOMICS_PROJECTS,METAGENOMICS_SAMPLES,altitude,biome,completion_date,creation_date,depth,domain_source,experiment_type,go_term,id,interpro_entry,name,organism,pH,pipeline_version,project_name,sample_name,temperature"
    url = "http://www.ebi.ac.uk/ebisearch/ws/rest/METAGENOMICS_SAMPLES?format=json&query=" + sample_id + "&fields=" + fields
    #print "sample_url: ", url
    sample_data_dict = get_run_data_ebi(url)
    sample_data_dict["sample_metadata_url"] = url
    # print "sample_data_dict: ", sample_data_dict
    return sample_data_dict


def merge_two_dicts(x, y):
    '''Given two dicts, merge them into a new dict as a shallow copy.'''
    z = x.copy()
    z.update(y)
    return z


def flatten_dict(y):
    out = {}
    def flatten(x, name=''):
        if type(x) is dict:
            for a in x:
                flatten(x[a], name + a + '_')
        elif type(x) is list:
            i = 0
            for a in x:
                flatten(a, name + str(i) + '_')
                i += 1
        else:
            out[name[:-1]] = x
    flatten(y)
    return out


def get_run_data_ebi(url):
    # print "url: ", url
    response = urllib2.urlopen(url)
    data = json.load(response)
    return data


def make_run_url(size, start, query):
    fields = "METAGENOMICS_PROJECTS,METAGENOMICS_SAMPLES,altitude,biome,completion_date,creation_date,depth,domain_source,experiment_type,go_term,id,interpro_entry,name,organism,pH,pipeline_version,project_name,sample_name,temperature"
    # url = "http://www.ebi.ac.uk/ebisearch/ws/rest/metagenomics_runs?format=json&start=" + str(start) + "&size=" + str(size) + "&query="+query+"&fields="+fields
    url = "http://www.ebi.ac.uk/ebisearch/ws/rest/metagenomics_runs?format=json&start=" + str(start) + "&size=" + str(
        size) + "&query=experiment_type:metagenomic+(" + query + ")" + "&fields=" + fields
    # print "run url: ", url
    return url


def get_sample_metadata_ena(identifier):
    url = "https://www.ebi.ac.uk/ena/data/view/" + identifier + "&display=xml"
    requestURL = url
    root = ET.parse(urllib.urlopen(requestURL)).getroot()
    metadata = {}
    metadata["ENA_metadata_url"] = url
    for sample in root:
        for sample_attributes in sample.findall("SAMPLE_ATTRIBUTES"):
            for sampleatribute in sample_attributes.findall('SAMPLE_ATTRIBUTE'):
                tag = sampleatribute.find('TAG')
                value = sampleatribute.find('VALUE')
                unit = sampleatribute.find('UNITS')
                if tag is not None:
                    tag = sampleatribute.find('TAG').text
                    if tag is None:
                        tag = 'unknown'
                if value is not None:
                    value = sampleatribute.find('VALUE').text
                    if value is None:
                        value = 'unknown'
                if unit is not None:
                    unit = sampleatribute.find('UNITS').text
                    if unit is None:
                        unit = 'unknown'
                    metadata[tag + " (" + unit + ")"] = value
                else:
                    if tag in metadata:
                        sys.stderr.write("double metadata: %s\n" % tag)
                        if tag != "BioSampleModel":
                            ctypes.windll.user32.MessageBoxW(0, u"This value is multiple times stored: " + tag,
                                                             u"Error", 0)
                        tag = tag + "2"
                    tag = tag.lower()
                    metadata[tag] = value
    # print "metadict: ", metadata
    metadata.pop("unknown", None)
    return metadata


def get_pubmed_metadata_ena(project, decoded_metadata_complete):
    url = "https://www.ebi.ac.uk/ena/data/view/" + project + "&display=xml"
    #print "prject", url
    requestURL = url
    root = ET.parse(urllib.urlopen(requestURL)).getroot()
    pubmed = ""
    for study in root:
        for study_links in study.findall("STUDY_LINKS"):

            for study_link in study_links.findall("STUDY_LINK"):

                for xref_link in study_link.findall("XREF_LINK"):

                    db = xref_link.find('DB')
                    ids = xref_link.find('ID')
                    if db is not None:
                        db = xref_link.find('DB').text
                        if db is None:
                            db = 'unknown'
                        if db == "pubmed":

                            if ids is not None:
                                ids = xref_link.find('ID').text
                                if ids is None:
                                    ids = 'unknown'
                                if pubmed == "":
                                    pubmed = str(ids)
                                else:
                                    pubmed = pubmed + " , " +str(ids)
    if pubmed != "":
        decoded_metadata_complete["pubmed"] = pubmed
        #print pubmed, "Yes there is a pubmed id!!!!!!!!!"
    return decoded_metadata_complete
main()
