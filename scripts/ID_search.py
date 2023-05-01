import pandas as pd
import copy
import random

safetyreportid_list = []
# makes a list of list of dates between start and end date
def get_dates(strt_dt, end_dt):
    # create a range of all dates between start and end date
    my_range = pd.date_range(start=strt_dt, end=end_dt)
    f_range = []
    for dt in my_range:
        y = str(dt)[0:4]
        m = str(dt)[5:7]
        d = str(dt)[8:10]
        new_dt = y + m + d
        f_range.append(new_dt)
    return f_range

# 2018 European heat wave selected dates
start = '20180726'
end = '20180804'
# 2018 country list
country_list = ["PT", "ES", "FR", "NL", "BE", "IT"]
date_range = get_dates(start, end)


# 2019 European heat wave selected dates
start0 = '20190623'
end0 = '20190628 '
# 2019 country list
country_list0 = ["BE", "FR", "DE", "NL", "GB", "PL", "ES", "CH"]
date_range0 = get_dates(start0, end0)

# dictionary templates
# template for KEGG value
KEGG_template = {"get_URL": None, "D_number": None, "Classes": None, "Target": None, "Pathway": None }  
# template for a drug value
drug_template = {"KEGG": KEGG_template.copy()}

# template for a primaryid value in control group
pid_template_control = { 'group': 0, 'reactions_MedDRA': [], 'drugs': {}}
# template for a primaryid value in heatwave group
pid_template_hw = { 'group': 1, 'reactions_MedDRA': [], 'drugs': {}}

# get heatwave ADE primary IDs
def FAERS_standard_generate_primaryids_hw(fp):
    print ("adding primaryids from heatwave")
    # iterate through each line of the demographics file and store primaryids from heatwaves
    for line in fp:
        split_line = line.split("$")
        if split_line[5] in date_range and split_line[8] in country_list:
            # store primaryid
            primaryid = split_line[1]
            yield primaryid
        elif split_line[5] in date_range0 and split_line[8] in country_list0:
            # store primaryid
            primaryid = split_line[1]
            yield primaryid

# get a list of entries by picking lines with matching countries
def FAERS_standard_generate_primaryids_control(path):
    print("adding primaryids from control")
    # make list of all primaryids in given date range from given countries
    with open("/Users/loaner/Documents/GitHub/Symbolic-Methods-FAERS-Project/text_files/all_countries.txt", "w") as output:
        for line in path:
            spl = line.split("$")
            # all 2018 ADE report primaryids from 2018 hw countries that do NOT fall within the heatwave date range
            if spl[8] in country_list:
                if spl[5] not in date_range:
                    output.write(line)
            # all 2019 ADE report primaryids from 2019 hw countries that do NOT fall within the heatwave date range
            elif spl[8] in country_list0:
                if spl[5] not in date_range0:
                    output.write(line)     
    print("done writing countries")
    # selects a random set of 40000 ADE reports
    filt = open("/Users/loaner/Documents/GitHub/Symbolic-Methods-FAERS-Project/text_files/all_countries.txt").read().splitlines()
    lset = set([])
    while len(lset) < 40000:
        x = random.choice(filt).split("$")
        lset.add(x[1])
    return lset

# generate drug info associated with each primaryid
def FAERS_standard_generate_drugs(id_dict, fp):
    dct = 0
    dset = {}
    # iterate through DRUGS_STANDARDIZED.txt
    for line in fp:
        split_line = line.strip("\n").split("$")
        if split_line[0] in id_dict:
            dct += 1
            # add drug info to corresponding dictionary entry
            dset[split_line[6]] = KEGG_template.copy()
            tmp = copy.deepcopy(id_dict[split_line[0]]['drugs'])
            tmp[split_line[6]] = dset[split_line[6]]
            id_dict[split_line[0]]['drugs'] = tmp
    print("num total drugs:", dct)
    print("num unique drugs:", len(dset))
    return [id_dict, dset]

# generate reactions associated with each primaryid
def FAERS_standard_generate_reactions(id_dict, fp):
    rct = 0
    for line in fp:
        split_line = line.split("$")
        # add reaction info to corresponding dictionary entry
        if split_line[0] in id_dict:
            rlist = id_dict[split_line[0]]['reactions_MedDRA'].copy()
            rlist.append(split_line[2].strip("\n"))
            id_dict[split_line[0]]['reactions_MedDRA'] = rlist
            rct += 1
            #yield [split_line[2].replace("\n", ""), string_list[i]]
    print("num reactions:", rct)
    return id_dict

# adds a drug subdictionary for each drug for each corresponding safetyreportid key in FAERS_standardized
def add_FAERS_standard_data_to_dictionary():
    # generate a dictionary of random primaryids for the control group
    with open("FAERS_standardized/DEMOGRAPHICS.txt", "r") as fp:
        primaryid_gen = FAERS_standard_generate_primaryids_control(fp)
        primaryid_list = []
        for l in primaryid_gen:
            primaryid_list.append(l)
    F_res_dict = dict( (pid, pid_template_control.copy()) for pid in primaryid_list )
    print("found ", len(F_res_dict), " control IDs")

    # find all primaryids from the experimental (heat wave) group and add them to the dictionary
    with open("/Users/loaner/Documents/GitHub/Symbolic-Methods-FAERS-Project/FAERS_standardized/DEMOGRAPHICS.txt", "r") as fp:
        # generate a list of demographic info for each primaryid that fulfills the search requirements
        primaryid_list_hw = FAERS_standard_generate_primaryids_hw(fp)
        hw_ade = 0
        for hw_id in primaryid_list_hw:
            F_res_dict[hw_id] = pid_template_hw.copy()
            hw_ade += 1
    print("found ", hw_ade, " heatwave IDs")
    
    # generate dictionary of all drugs associated with primaryids, update pid dictionary
    print("starting drugs")
    with open("/Users/loaner/Documents/GitHub/Symbolic-Methods-FAERS-Project/FAERS_standardized/DRUGS_STANDARDIZED.txt", "r") as fp:
        F_res_dict_d = FAERS_standard_generate_drugs(F_res_dict, fp)

    print("drugs done, starting reactions")
    # adds reactions for each primary key by searching ADVERSE_REACTIONS.txt
    with open("/Users/loaner/Documents/GitHub/Symbolic-Methods-FAERS-Project/FAERS_standardized/ADVERSE_REACTIONS.txt") as fp:
        rdict = FAERS_standard_generate_reactions(F_res_dict_d[0], fp)
    print("done with reactions")
    F_res_dict_d[0] = rdict
    return(F_res_dict_d)
