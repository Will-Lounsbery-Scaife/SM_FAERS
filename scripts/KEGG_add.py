import requests
import copy
import traceback
import time

# strings for KEGG request URLs
KEGG_find_base = 'https://rest.kegg.jp/find/drug/'
KEGG_get_base = 'https://rest.kegg.jp/get/'


# function to generate the appropriate URLs for KEGG API queries
# uses the find request to search a drug string, then determines the d-number of the drug for the get-url 
def KEGG_find_query(drug_string):
    poss1 = 0
    find_url = KEGG_find_base + drug_string
    find_req = requests.get(find_url)
    a = find_req.text.split("\n")
    b = []
    for l in a:
        b.append(l.split("\t"))
    b.pop()
    # iterate through find entries to find INN names, extract the d number
    # if multiple entries contain INN, store the first one with matching start
    for found_entry in b:
        if "INN" in found_entry[1] and found_entry[1].lower().startswith(drug_string.lower()):
            get_url = KEGG_get_base + found_entry[0]
            x = found_entry[0].split(":")
            d_number = x[1]
            # put d number at end of get URL
            return [find_url, get_url, d_number]
    for found_entry in b:
        if "INN" in found_entry[1]:
            get_url = KEGG_get_base + found_entry[0]
            x = found_entry[0].split(":")
            d_number = x[1]
            poss1 = [find_url, get_url, d_number]
            break
    # if no "INN" entries found, store first that has drug name not preceded by non-whitespace chars
    for found_entry in b:
        if found_entry[1].lower().startswith(drug_string.lower()):
            d_number = found_entry[0].split(":")[1]
            x = [find_url, KEGG_get_base + b[0][0], d_number]
            return x
    # otherwise, just store first entry
    #if poss1: return poss1
    if poss1 != 0: return poss1
    return [find_url, KEGG_get_base + b[0][0], b[0][0].split(":")[1]]

# send get request to Kegg and store drug targets, pathway, and class
def KEGG_get_query(get_URL):
    target = []
    pathway = []
    class_info = []
    get_req = requests.get(get_URL)
    klines = get_req.iter_lines()
    # booleans to keep track of if a line stores info on pathway, class, or target
    is_targ = 0
    is_path = 0
    is_class = 0
    # iterate through get response
    for line in klines:
        dec = line.decode("utf-8")
        # line has class info
        if dec.startswith("CLASS"):
            is_class = 1
            is_path = 0
            is_targ = 0
        if dec[2:9]=="PATHWAY":
            is_path = 1
            is_targ = 0
            is_class = 0
        # line does not have class info
        elif dec.startswith("REMARK"):
            is_class = 0
            is_targ = 0
            is_path = 0
        # line has target info
        elif dec.startswith("TARGET"):
            is_targ = 1
            is_class = 0
            is_path = 0
        # line has pathway info, not target info
        elif dec.startswith("STR_MAP"):
            is_targ = 0
            is_path = 0
            is_class = 0
        elif dec.startswith("BRITE") or dec.startswith("METABOLISM") or dec.startswith("REMARK") or dec.startswith("EFFICACY"):
            is_targ = 0
            is_path = 0
            is_class = 0
        # line no longer has path info
        elif dec.startswith("INTERACTION"):
            is_targ = 0
            is_class =0
            is_path = 0
        # if the line corresponds to a field of interest, store it
        if is_class == 1: 
            # gets the name and dg-number of the KEGG drug classes 
            class_info.append(dec[12:])
        elif is_targ == 1: target.append(dec[12:])
        elif is_path == 1: 
            p = dec[12:].split("(")
            o = p[0].split(" ")
            pathway.append(o[0])
    # store drug groups
    class_list = set()
    for ind in range(len(class_info)):
        if class_info[ind].startswith(" DG"):
            y = class_info[ind].strip(" ").split(" ")[0]
            x = tuple([y, 1])
            class_list.add(x)
        elif class_info[ind].startswith("  DG"):
            y = class_info[ind].strip(" ").split(" ")[0]
            x = tuple([y, 2])
            class_list.add(x)
        elif class_info[ind].startswith("   DG"):
            y = class_info[ind].strip(" ").split(" ")[0]
            x = tuple([y, 3])
            class_list.add(x)
        elif class_info[ind].startswith("    DG"):
            y = class_info[ind].strip(" ").split(" ")[0]
            x = tuple([y, 4])
            class_list.add(x)
    return(target, pathway, class_list)

# update dictionary with info from KEGG API requests
def FAERS_standard_get_KEGG_info(results_dict, F_dlist):

    #print("res dict before:\n", results_dict, "\n\n")
    #print("d dict before:\n", F_dlist, "\n")

    kegg_find_results = 0
    kegg_get_results = 0
    dct = 0
    
    startTime = time.time()
    for a in F_dlist:
        if dct % 10 == 0: print("drugs processed:", dct)
        dct += 1
        # fix for polyethylene glycol
        if a.startswith("polyethylene glycol"):
            kegg_get_results = KEGG_get_query('https://rest.kegg.jp/get/dr:D03370')
            F_dlist[a]['get_URL'] = 'https://rest.kegg.jp/get/dr:D03370'
            F_dlist[a]['D_number'] = 'D03370'
            F_dlist[a]['Target'] = kegg_get_results[0]
            F_dlist[a]['Pathway'] = kegg_get_results[1]
            F_dlist[a]['Classes'] = kegg_get_results[2]
            continue
        # fix for ibrutinib
        if a == 'ibrutinib':
            kegg_get_results = KEGG_get_query('https://rest.kegg.jp/get/dr:D03936')
            F_dlist[a]['get_URL'] = 'https://rest.kegg.jp/get/dr:D03936'
            F_dlist[a]['D_number'] = 'D03936'
            F_dlist[a]['Target'] = kegg_get_results[0]
            F_dlist[a]['Pathway'] = kegg_get_results[1]
            F_dlist[a]['Classes'] = kegg_get_results[2]
            continue
        # fix for "Streptococcus pneumoniae" vaccine
        if a.startswith('Streptococcus pneumoniae'):
            kegg_get_results = KEGG_get_query('https://rest.kegg.jp/get/dr:D10455')
            F_dlist[a]['get_URL'] = 'https://rest.kegg.jp/get/dr:D10455'
            F_dlist[a]['D_number'] = 'D10455'
            F_dlist[a]['Target'] = kegg_get_results[0]
            F_dlist[a]['Pathway'] = kegg_get_results[1]
            F_dlist[a]['Classes'] = kegg_get_results[2]
            continue
        # fix for "insulin, regular, human"
        if a == 'insulin, regular, human':
            kegg_find_results = KEGG_find_query('insulin human')
            kegg_get_results = KEGG_get_query(kegg_find_results[1])
            F_dlist[a]['get_URL'] = kegg_find_results[1]
            F_dlist[a]['D_number'] = kegg_find_results[2]
            F_dlist[a]['Target'] = kegg_get_results[0]
            F_dlist[a]['Pathway'] = kegg_get_results[1]
            F_dlist[a]['Classes'] = kegg_get_results[2]
            continue
        # fix for Ursodiol
        if a == 'ursodeoxycholate':
            kegg_find_results = KEGG_find_query('Ursodiol')
            kegg_get_results = KEGG_get_query(kegg_find_results[1])
            F_dlist[a]['get_URL'] = kegg_find_results[1]
            F_dlist[a]['D_number'] = kegg_find_results[2]
            F_dlist[a]['Target'] = kegg_get_results[0]
            F_dlist[a]['Pathway'] = kegg_get_results[1]
            F_dlist[a]['Classes'] = kegg_get_results[2]
            continue
        # try with original RxNorm drug name
        try:
            kegg_find_results = KEGG_find_query(a)
            kegg_get_results = KEGG_get_query(kegg_find_results[1])
            F_dlist[a]['get_URL'] = kegg_find_results[1]
            F_dlist[a]['D_number'] = kegg_find_results[2]
            F_dlist[a]['Target'] = kegg_get_results[0]
            F_dlist[a]['Pathway'] = kegg_get_results[1]
            F_dlist[a]['Classes'] = kegg_get_results[2]
            continue
        except Exception as e1:
            print("tried original name:", a, "error:", e1)  
            pass        
        # try fix for drugs that end with product
        if str(a).strip(" ").endswith("product"):
            try:
                fixy = str(a).strip(" ")
                fixy = fixy[:-7]
                kegg_find_results = KEGG_find_query(fixy[0])
                kegg_get_results = KEGG_get_query(kegg_find_results[1])
                F_dlist[a]['get_URL'] = kegg_find_results[1]
                F_dlist[a]['D_number'] = kegg_find_results[2]
                F_dlist[a]['Target'] = kegg_get_results[0]
                F_dlist[a]['Pathway'] = kegg_get_results[1]
                F_dlist[a]['Classes'] = kegg_get_results[2]
                print ("fixed product", a)
                continue
            except Exception as e1:
                print("tried removing product", e1, a)  
                pass
        # try fix for drugs with hyphens
        if "-" in str(a):
            try:
                fix2 = str(a).replace("-", " ")
                kegg_find_results = KEGG_find_query(fix2)
                kegg_get_results = KEGG_get_query(kegg_find_results[1])
                F_dlist[a]['get_URL'] = kegg_find_results[1]
                F_dlist[a]['D_number'] = kegg_find_results[2]
                F_dlist[a]['Target'] = kegg_get_results[0]
                F_dlist[a]['Pathway'] = kegg_get_results[1]
                F_dlist[a]['Classes'] = kegg_get_results[2]
                print("fixed hyphens", a)
                continue
            except Exception as e1:
                print("tried removing hyphens", e1, a)  
                pass
        # try fix for drugs that end with ", human"
        if str(a).endswith(", human"):
            try:
                fix3 = str(a).split(",")
                kegg_find_results = KEGG_find_query(fix3[0])
                kegg_get_results = KEGG_get_query(kegg_find_results[1])
                F_dlist[a]['get_URL'] = kegg_find_results[1]
                F_dlist[a]['D_number'] = kegg_find_results[2]
                F_dlist[a]['Target'] = kegg_get_results[0]
                F_dlist[a]['Pathway'] = kegg_get_results[1]
                F_dlist[a]['Classes'] = kegg_get_results[2]
                print("fixed human", a)
                continue
            except Exception as e1:
                print("tried removing ', human", e1, a)  
                pass
            # fix for vaccines with "innactivated"
            if str(a).endswith(", inactivated"):
                try:
                    fix3 = str(a).split(",")
                    kegg_find_results = KEGG_find_query(fix3[0])
                    kegg_get_results = KEGG_get_query(kegg_find_results[1])
                    F_dlist[a]['get_URL'] = kegg_find_results[1]
                    F_dlist[a]['D_number'] = kegg_find_results[2]
                    F_dlist[a]['Target'] = kegg_get_results[0]
                    F_dlist[a]['Pathway'] = kegg_get_results[1]
                    F_dlist[a]['Classes'] = kegg_get_results[2]
                    continue
                except Exception as e1:
                    print("tried removing ' inactivated", e1, a)  
                    pass
        # try fix for drugs that end with "4000"
        if str(a).strip(" ").endswith("4000"):
            try:
                fix3 = str(a).strip(" ")
                fix3 = fix3[:-4]
                kegg_find_results = KEGG_find_query(fix3[0])
                kegg_get_results = KEGG_get_query(kegg_find_results[1])
                F_dlist[a]['get_URL'] = kegg_find_results[1]
                F_dlist[a]['D_number'] = kegg_find_results[2]
                F_dlist[a]['Target'] = kegg_get_results[0]
                F_dlist[a]['Pathway'] = kegg_get_results[1]
                F_dlist[a]['Classes'] = kegg_get_results[2]
                print("fixed 4000", a)
                continue
            except Exception as e1:
                print("tried removing 4000", e1, a)  
                pass 
        # try fix for drugs that end with ", USP"
        if str(a).strip(" ").endswith(", USP") :
            try:
                fixit = str(a).strip(" ")
                fix3 = fixit[:-5]
                kegg_find_results = KEGG_find_query(fix3)
                kegg_get_results = KEGG_get_query(kegg_find_results[1])
                F_dlist[a]['get_URL'] = kegg_find_results[1]
                F_dlist[a]['D_number'] = kegg_find_results[2]
                F_dlist[a]['Target'] = kegg_get_results[0]
                F_dlist[a]['Pathway'] = kegg_get_results[1]
                F_dlist[a]['Classes'] = kegg_get_results[2]
                print("fixed USP:", a, "\nused:", fix3)
                continue
            except Exception as e1:
                print("tried removing ', USP", e1, a)  
                pass
        # try removing unnecessary trailing s
        if str(a)[-1] == 's':
            try:
                fix1 = str(a)
                fix1 = fix1[:-1]
                kegg_find_results = KEGG_find_query(fix1)
                kegg_get_results = KEGG_get_query(kegg_find_results[1])
                F_dlist[a]['get_URL'] = kegg_find_results[1]
                F_dlist[a]['D_number'] = kegg_find_results[2]
                F_dlist[a]['Target'] = kegg_get_results[0]
                F_dlist[a]['Pathway'] = kegg_get_results[1]
                F_dlist[a]['Classes'] = kegg_get_results[2]
                print("fixed trailing s", a)
                continue
            except Exception as e1:
                print("tried removing trailing s", e1, a)  
                pass
        print("No KEGG:", a, "\n")
    executionTime = (time.time() - startTime)
    print('seconds to retrieve all KEGG info: ' + str(executionTime))
    
    startTime = time.time()
    # iterate through the drugs dictionary and add corresponding drugs to the results dictionary
    print("adding KEGG drug info to results dict")
    for pid in results_dict:
        for d in results_dict[pid]['drugs']:
            if d in F_dlist:
                tm = F_dlist[d].copy()
                results_dict[pid]['drugs'][d] = tm
    executionTime = (time.time() - startTime)
    print('seconds to update results dict from dlist: ' + str(executionTime))
    return(results_dict)

