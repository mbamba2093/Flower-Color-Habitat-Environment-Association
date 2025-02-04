from openai import OpenAI
import pandas as pd
import time
import glob
import sys


def Extract_color(desc1, openai_api_key):
    client = OpenAI(api_key=openai_api_key)
    prompt1 = "[Color] is a description for a specific flower color. Estimate to which category ('White', 'Yellow', 'Red', 'Blue', 'Purple', 'Green', or 'No description') the flower color belongs."
    prompt2 = "Output should be a single word indicating the category. If there are multiple candidates, separate them with ':::'."
    chat_completion = client.chat.completions.create(
    messages=[
        {"role": "user",
         "content": f"{prompt1}\n{prompt2}\n[Color]:{desc1}"}],
    model="gpt-4o")
    chat_content = chat_completion.choices[0].message.content     
    modified_string = chat_content.replace("\n", ", ")
    modified_string = modified_string.replace('"', '')
    modified_string = modified_string.replace("'", "")
    return modified_string




#####
#####
##### For Flora of China datasets
#####
#####



openai_api_key = 'XXXXXXXXXX'
Data_dir = "/media/masarubamba/vol1/Project/WebScraping/"
pkl_list = sorted(glob.glob(f"{Data_dir}FOC*.pkl"))
maxattempt = 5
retry_interval = 300

GPT_df = pd.read_csv(f"{Data_dir}GPT4o_all_filt.csv")

Return_list = []
temp_num = 0
for j in GPT_df["GPT_return1"]:
    for attempt in range(maxattempt):
        try:
            features = Extract_color(j, openai_api_key)
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

GPT_df["GPT_return2"] = Return_list
GPT_df.to_csv(f"{Data_dir}GPT4o_all_filt_categ.csv")
print("finished")

Data_dir = "/media/masarubamba/vol1/Project/WebScraping/"
GPT_df = pd.read_csv(f"{Data_dir}GPT4o_all_filt_categ.csv")

GPT_df["GPT_return2_copy"] = GPT_df["GPT_return2"]
GPT_df["White"] = GPT_df["GPT_return2"].str.lower().str.contains("white").astype(int)
GPT_df["GPT_return2_copy"] = GPT_df["GPT_return2_copy"].str.lower().str.replace("white", "")

GPT_df["Yellow"] = GPT_df["GPT_return2"].str.lower().str.contains("yellow|orange").astype(int)
GPT_df["GPT_return2_copy"] = GPT_df["GPT_return2_copy"].str.lower().str.replace("yellow", "").replace("orange", "")
GPT_df["GPT_return2_copy"] = GPT_df["GPT_return2_copy"].str.lower().str.replace("orange", "")

GPT_df["Red"] = GPT_df["GPT_return2"].str.lower().str.contains("red|pink").astype(int)
GPT_df["GPT_return2_copy"] = GPT_df["GPT_return2_copy"].str.lower().str.replace("red", "").replace("pink", "")
GPT_df["GPT_return2_copy"] = GPT_df["GPT_return2_copy"].str.lower().str.replace("pink", "")

GPT_df["Blue"] = GPT_df["GPT_return2"].str.lower().str.contains("blue").astype(int)
GPT_df["GPT_return2_copy"] = GPT_df["GPT_return2_copy"].str.lower().str.replace("blue", "")

GPT_df["Purple"] = GPT_df["GPT_return2"].str.lower().str.contains("purple|violet").astype(int)
GPT_df["GPT_return2_copy"] = GPT_df["GPT_return2_copy"].str.lower().str.replace("purple", "").replace("violet", "")
GPT_df["GPT_return2_copy"] = GPT_df["GPT_return2_copy"].str.lower().str.replace("violet", "")

GPT_df["Green"] = GPT_df["GPT_return2"].str.lower().str.contains("green").astype(int)
GPT_df["GPT_return2_copy"] = GPT_df["GPT_return2_copy"].str.lower().str.replace("green", "")

GPT_df["NoDescription"] = GPT_df["GPT_return2"].str.lower().str.contains("no description|no descrição").astype(int)
GPT_df["GPT_return2_copy"] = GPT_df["GPT_return2_copy"].str.lower().str.replace("no description", "").replace("no descrição", "")
GPT_df["GPT_return2_copy"] = GPT_df["GPT_return2_copy"].str.lower().str.replace("no descrição", "")

GPT_df["GPT_return2_copy"] = GPT_df["GPT_return2_copy"].str.replace(":", "")

GPT_df.to_csv(f"{Data_dir}GPT4o_all_fin.csv")


#####
#####
##### For TRY dataset
#####
#####



Data_dir = "/media/masarubamba/vol1/Project/WebScraping/"
TRY_df = pd.read_csv(f"{Data_dir}TRY_data.csv", encoding="cp932")

Return_list = []
temp_num = 0
for j in TRY_df["OrigValueStr"]:
    for attempt in range(maxattempt):
        try:
            features = Extract_color(j, openai_api_key)
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

TRY_df["GPT_return2_copy"] = TRY_df["GPT_return2"]
TRY_df["White"] = TRY_df["GPT_return2"].str.lower().str.contains("white").astype(int)
TRY_df["GPT_return2_copy"] = TRY_df["GPT_return2_copy"].str.lower().str.replace("white", "")

TRY_df["Yellow"] = TRY_df["GPT_return2"].str.lower().str.contains("yellow|orange").astype(int)
TRY_df["GPT_return2_copy"] = TRY_df["GPT_return2_copy"].str.lower().str.replace("yellow", "").replace("orange", "")
TRY_df["GPT_return2_copy"] = TRY_df["GPT_return2_copy"].str.lower().str.replace("orange", "")

TRY_df["Red"] = TRY_df["GPT_return2"].str.lower().str.contains("red|pink").astype(int)
TRY_df["GPT_return2_copy"] = TRY_df["GPT_return2_copy"].str.lower().str.replace("red", "").replace("pink", "")
TRY_df["GPT_return2_copy"] = TRY_df["GPT_return2_copy"].str.lower().str.replace("pink", "")

TRY_df["Blue"] = TRY_df["GPT_return2"].str.lower().str.contains("blue").astype(int)
TRY_df["GPT_return2_copy"] = TRY_df["GPT_return2_copy"].str.lower().str.replace("blue", "")

TRY_df["Purple"] = TRY_df["GPT_return2"].str.lower().str.contains("purple|violet").astype(int)
TRY_df["GPT_return2_copy"] = TRY_df["GPT_return2_copy"].str.lower().str.replace("purple", "").replace("violet", "")
TRY_df["GPT_return2_copy"] = TRY_df["GPT_return2_copy"].str.lower().str.replace("violet", "")

TRY_df["Green"] = TRY_df["GPT_return2"].str.lower().str.contains("green").astype(int)
TRY_df["GPT_return2_copy"] = TRY_df["GPT_return2_copy"].str.lower().str.replace("green", "")

TRY_df["NoDescription"] = TRY_df["GPT_return2"].str.lower().str.contains("no description|no descrição").astype(int)
TRY_df["GPT_return2_copy"] = TRY_df["GPT_return2_copy"].str.lower().str.replace("no description", "").replace("no descrição", "")
TRY_df["GPT_return2_copy"] = TRY_df["GPT_return2_copy"].str.lower().str.replace("no descrição", "")

TRY_df["GPT_return2_copy"] = TRY_df["GPT_return2_copy"].str.replace(":", "")

TRY_df.to_csv(f"{Data_dir}TRY_data_categ_fin.csv", encoding="cp932")


