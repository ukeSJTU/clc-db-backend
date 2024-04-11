import sqlite3
import csv
import pandas as pd

# Path to the SQLite database
database_path = "db.sqlite3"

xlsx_file_path = "data.xlsx"

df = pd.read_excel(xlsx_file_path)


# Connect to the SQLite database
conn = sqlite3.connect(database_path)
cursor = conn.cursor()

task_id = [1]

#################
# 1. handle categories
#################
if 1 in task_id:
    delim = "、"
    table_name = "core_category"
    df_column_name = "类别"

    categories = []

    for index, row in df.iterrows():
        for single_category in row[df_column_name].split(delim):
            if single_category != "无":
                categories.append(single_category)

    categories = set(categories)

    for idx, single_category in enumerate(categories):
        print(single_category)
        cursor.execute(
            f"INSERT INTO {table_name} (id, name) VALUES (?, ?)",
            (
                idx + 1,
                single_category,
            ),
        )


#################
# 2. handle smile types
#################
if 2 in task_id:
    delim = "、"
    table_name = "core_smile"
    df_column_name = "手性种类"
    smile_types = []
    for index, row in df.iterrows():
        for single_smile in row[df_column_name].split(delim):
            if single_smile != "无":
                smile_types.append(single_smile)

    smile_types = set(smile_types)

    for idx, single_smile in enumerate(smile_types):
        print(single_smile)
        cursor.execute(
            f"INSERT INTO {table_name} (id, smile) VALUES (?, ?)",
            (
                idx + 1,
                single_smile,
            ),
        )

conn.commit()
conn.close()

# print(set(categories))

print("Data imported successfully.")
