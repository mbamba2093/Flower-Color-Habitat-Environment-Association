import requests
from bs4 import BeautifulSoup
import re
import pandas as pd
import pickle
import time
import random
import sys


def Extract_description(taxon_id):
    base_url = "http://www.efloras.org/florataxon.aspx?flora_id=2&taxon_id="
    positive_words = ["Herbs", "Shrubs", "Trees", "Plants", "Climbers", "Creepers", "Stem", "Subshrubs", "Perennial", "Annuals","Caudex", "Roots", "Main", "Petals", "Sepals", "Rhizome", "Basal", "Rootstock", "Rosette", "Leaf"]
    ng_word = " after publication of the family treatment for the Flora of China."

    response_species = requests.get(f"{base_url}{taxon_id}")
    if response_species.status_code != 200:
        retries = 0
        max_retries = 10
        time.sleep(10)
        while retries < max_retries:
            response_species = requests.get(Species_i_url)
            if response_species.status_code == 200:
                break
            retries += 1
    
    if response_species.status_code != 200:
        print(f"{taxon_id} is Status {response_species.status_code}")
        return "Status Error"
    else:
        soup_species = BeautifulSoup(response_species.text, 'html.parser')
        species_desc = soup_species.find_all('span', id='lblTaxonDesc')
        species_desc = species_desc[0].find_all('p') 
        species_desc = species_desc[0]

        return_key = 0
        herbs_lines = []
        for p1 in species_desc.get_text().split('\n'):
            if p1 != "\r" and p1 != "":
                herbs_lines.append(p1)
                return_key += 1
        result_description = "::".join(herbs_lines)

        if return_key == 0 or ng_word in result_description:
            if any(positive_word in result_description for positive_word in positive_words):
                return result_description
            print(taxon_id)
            return "No_description"
        else:
            return result_description
    
    




def Extract_taxalist(taxon_id, target = "non-volume", maxpages = 5):
    base_url = "http://www.efloras.org/browse.aspx?flora_id=2&start_taxon_id="
    volume_url = f"http://www.efloras.org/"

    if not target == "non-volume":
        base_url = volume_url
    
    # Assuming page listings 1-5.
    # Since FOC are displaying 200 taxa per page
    
    html_list = []
    start_page = 1
    while start_page < maxpages:
        taxa_url = f"{base_url}{taxon_id}&page={start_page}"
        response_taxa = requests.get(taxa_url)

        # Redirect when there is no response from the server
        if response_taxa.status_code != 200:
            retries = 0
            max_retries = 10
            time.sleep(10)
            while retries < max_retries:
                print(f"Retry {taxon_id}")
                response_taxa = requests.get(taxa_url)
                if response_taxa.status_code == 200:
                    break
                retries += 1

        if response_taxa.status_code != 200:
            print(f"{taxon_id} is Status {response_taxa.status_code}")
            return "Empty data"
        
        # Check if the target table exists.
        soup_taxa = BeautifulSoup(response_taxa.text, 'html.parser')
        tables_taxa = soup_taxa.find_all('table')
        rows = tables_taxa[0].find_all('tr', class_='underline') 
        if rows:
            """
            html_list.append(soup_taxa.prettify())#soup_taxa
            """
            rowtemp = []
            for row in rows:
                rowtemp.append(row.prettify())
            html_list.append(" ".join(rowtemp))
            start_page += 1
        else:
            break
    
    if not html_list:
        print("no_html")
        return "Empty data"
    
    # Combine the htmls from multiple pages.
    # Then, analyze and extract taxa tables

    combined_html = '<table>' + "".join(html_list) + '</table>'
    soup_taxa = BeautifulSoup(combined_html, 'html.parser')
    tables_taxa = soup_taxa.find_all('table')
    rows = tables_taxa[0].find_all('tr', class_='underline') 
    taxa_list = []
    for row in rows:
        cols = row.find_all('td')
        cols = [ele.text.strip() for ele in cols]
        taxa_list.append([ele for ele in cols if ele])
    
    if not taxa_list:
        print("no_taxalist")
        return "Empty data"
    
    # Cleaning the table.
    # Since the tables come from multiple htmls, header and fotter were included.
    remove_index = set()
    for i in range(len(taxa_list)):
        if taxa_list[i] == []:
            remove_index.add(i)
            if i + 1 < len(taxa_list):
                remove_index.add(i + 1)
                
    for indexi in sorted(remove_index, reverse=True):
        taxa_list.pop(indexi)
    
    # Transform to dataframe
    taxa_df = pd.DataFrame(taxa_list)
    max_cols = max(len(row) for row in taxa_list)
    column_names = ['Taxon Id', 'Name'] + [f'Info{i}' for i in range(1, max_cols - 1)]
    taxa_df.columns = column_names[:max_cols]
    taxa_df["Name"] = [' '.join(i.replace('\n', ' ').split()) for i in taxa_df["Name"]] # Sometaxa contains 

    temp_dict = {}
    for index, row in taxa_df.iterrows():
        # {Name: {id: Taxon Id}}
        temp_dict[row['Name']] = {'id': row['Taxon Id']}
    
    return temp_dict


###
###
###
###
###
###
###
###
###
###
###
###




# Top page URL
Taxa_data_out = "/media/masarubamba/vol1/Project/WebScraping/FOC_"
url_top = 'http://www.efloras.org/flora_page.aspx?flora_id=2'
response_top = requests.get(url_top)

# Read top page
soup_top = BeautifulSoup(response_top.text, 'html.parser')
#Volume search
url_pattern_volume = re.compile(r'volume_page\.aspx\?volume_id=\d+&flora_id=2')
Volume_url_list = []
for link in soup_top.find_all('a', href=url_pattern_volume):
    href = link.get('href')
    if href:
        Volume_url_list.append(href)

# Volume 1 is an introduction. 
Volume_url_list = Volume_url_list[1:]



# Volume_url_list

#for Volume_i in Volume_url_list:
Volume_url_list = Volume_url_list[5:]
for Volume_i in Volume_url_list:
    taxa_data = {}
    print(f"{Volume_i} volume started")
    
    volume_name = Volume_i.split("volume_id=")[1].split("&")[0]
    taxa_data[volume_name] = {"Family": Extract_taxalist(Volume_i, target = "volume")}

    #tempid = 0
    maxattempt = 5
    retry_interval = 300

    for Family_i in taxa_data[volume_name]["Family"].keys():
        print(f"{Family_i} is started")
        taxon_id_i = taxa_data[volume_name]["Family"][Family_i]["id"]

        for attempt in range(maxattempt):
            try:
                Extract_flags = Extract_taxalist(taxon_id_i)
                break
            except Exception as e:
                print("Error")
                if attempt < maxattempt - 1:
                    time.sleep(retry_interval)
                else:
                    print("Max retry is over")
                    sys.exit()

        if Extract_flags == "Empty data":
            print(f"{Family_i} is empty")
            continue
        else:
            Extract_flags = {key:value for key, value in Extract_flags.items() if len(key.split()) < 2}
            if "Genus" in taxa_data[volume_name]["Family"][Family_i].keys():
                if not set(Extract_flags.keys()) == set(taxa_data[volume_name]["Family"][Family_i]["Genus"].keys()):
                    taxa_data[volume_name]["Family"][Family_i]["Genus"] = Extract_flags
            else:
                taxa_data[volume_name]["Family"][Family_i]["Genus"] = Extract_flags

        for Genus_i in taxa_data[volume_name]["Family"][Family_i]["Genus"].keys():
            print(f"{Genus_i} is started")
            taxon_id_i = taxa_data[volume_name]["Family"][Family_i]["Genus"][Genus_i]["id"]

            for attempt in range(maxattempt):
                try:
                    Extract_flags = Extract_taxalist(taxon_id_i)
                    break
                except Exception as e:
                    print("Error")
                    if attempt < maxattempt - 1:
                        time.sleep(retry_interval)
                    else:
                        print("Max retry is over")
                        sys.exit()

            if Extract_flags == "Empty data":
                print(f"{Genus_i} is empty")
                continue
            else:
                Extract_flags = {key:value for key, value in Extract_flags.items() if '\"' not in key}
                if "Species" in taxa_data[volume_name]["Family"][Family_i]["Genus"][Genus_i].keys():
                    if not set(Extract_flags.keys()) == set(taxa_data[volume_name]["Family"][Family_i]["Genus"][Genus_i]["Species"].keys()):
                        taxa_data[volume_name]["Family"][Family_i]["Genus"][Genus_i]["Species"] = Extract_flags
                else:
                    taxa_data[volume_name]["Family"][Family_i]["Genus"][Genus_i]["Species"] = Extract_flags

            for Species_i in taxa_data[volume_name]["Family"][Family_i]["Genus"][Genus_i]["Species"].keys():
                if "Description" not in taxa_data[volume_name]["Family"][Family_i]["Genus"][Genus_i]["Species"][Species_i].keys() or taxa_data[volume_name]["Family"][Family_i]["Genus"][Genus_i]["Species"][Species_i]["Description"] == "No_description":
                    print(Species_i)
                    taxon_id_i = taxa_data[volume_name]["Family"][Family_i]["Genus"][Genus_i]["Species"][Species_i]["id"]

                    for attempt in range(maxattempt):
                        try:
                            Extract_flags = Extract_description(taxon_id_i)
                            break
                        except Exception as e:
                            print("Error")
                            if attempt < maxattempt - 1:
                                time.sleep(retry_interval)
                            else:
                                print("Max retry is over")
                                sys.exit()

                    if Extract_flags == "Status Error":
                        print(f"{Species_i} is Status Error")
                        taxa_data[volume_name]["Family"][Family_i]["Genus"][Genus_i]["Species"][Species_i]["Description"] = "Status Error"
                    else:
                        taxa_data[volume_name]["Family"][Family_i]["Genus"][Genus_i]["Species"][Species_i]["Description"] = Extract_flags
                    
                    delay = random.randint(3, 10)  # 3秒から10秒のランダムな時間を生成
                    time.sleep(delay)
                    #tempid += 1 
                    #tempflag = False
                else:
                    print(f"{Species_i} is skipped")
                    continue

    Taxa_data_out_temp = f"{Taxa_data_out}{volume_name}.pkl"
    with open(Taxa_data_out_temp, 'wb') as f:
        pickle.dump(taxa_data, f)