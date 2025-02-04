from openai import OpenAI
import pandas as pd
import pickle
import time
import glob
import sys

def convert_to_dataframe(data):
    rows = []
    for volume, family_data in data.items():
        for family, genus_data in family_data['Family'].items():
            for genus, species_data in genus_data['Genus'].items():
                if "Species" not in species_data.keys():
                    print(f"{genus} does not have species data")
                    break
                else:
                    for species, details in species_data['Species'].items():
                        rows.append({
                            'Species': species,
                            'Family': family,
                            'Genus': genus,
                            'Volume': volume,
                            'Description': details['Description']
                        })
    return pd.DataFrame(rows)

def Extract_features(desc1, openai_api_key):
    client = OpenAI(api_key=openai_api_key)
    prompt1 = "Extract the part of the text that describes the flower petal color from [Description] and respond with that exact text."
    prompt2 = "If there is no information about the petal color, respond with 'no_data'."
    chat_completion = client.chat.completions.create(
    messages=[
        {"role": "user",
         "content": f"{prompt1}\n{prompt2}\n [Description]: {desc1}"}],
    model="gpt-4o")
    chat_content = chat_completion.choices[0].message.content     
    modified_string = chat_content.replace("\n", ", ")
    modified_string = modified_string.replace('"', '')
    modified_string = modified_string.replace("'", "")
    return modified_string
  
openai_api_key = 'XXXXXXXXXX'
Data_dir = "/Datadir/"
pkl_list = sorted(glob.glob(f"{Data_dir}FOC*.pkl"))
maxattempt = 5
retry_interval = 300

for i in pkl_list:
    pkl_name = i.split("/")[-1].split(".pkl")[0]
    with open(i, 'rb') as f:
        Desc_data = pickle.load(f)
    Desc_df = convert_to_dataframe(Desc_data)
    
    Return_list = []
    temp_num = 0
    for j in Desc_df["Description"]:
        for attempt in range(maxattempt):
            try:
                features = Extract_features(j, openai_api_key)
                break
            except Exception as e:
                print("Error")
                if attempt < maxattempt - 1:
                    time.sleep(retry_interval)
                else:
                    print("Max retry is over")
                    sys.exit()
        temp_num += 1
        Return_list.append(features)
        print(f"{temp_num} : {features}")

    Desc_df["GPT_return1"] = Return_list
    Desc_df.to_csv(f"{Data_dir}{pkl_name}_GPT4o.csv")
    print(f"{i} is finished")
    
csv_list = sorted(glob.glob(f"{Data_dir}*_GPT4o.csv"))
df_list = []
for i in csv_list:
    temp_df = pd.read_csv(i)
    df_list.append(temp_df)
df_combined = pd.concat(df_list, ignore_index = True)
df_combined.to_csv(f"{Data_dir}GPT4o_all.csv")

df_filtered = df_combined[~df_combined["GPT_return1"].str.lower().str.contains("no_data", na = False)]
df_filtered['Unnamed: 0'] = pd.to_numeric(df_filtered['Unnamed: 0'], errors='coerce')
df_filtered = df_filtered.dropna(subset=['Unnamed: 0'])
df_filtered = df_filtered.dropna(subset=['GPT_return1'])
df_filtered.to_csv(f"{Data_dir}GPT4o_all_filt.csv")

