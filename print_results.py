import argparse
import sys
import os
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(add_help=False)
parser.add_argument("--testproblem", type=str, required=True)

args, _ = parser.parse_known_args()
testproblem = args.testproblem

linearisations = ["newton"]

for linearisation in linearisations:
    mypath = f"results/results{linearisation}{testproblem}"
    my_files = [f for f in os.listdir(mypath) if (os.path.isfile(os.path.join(mypath, f)) and f.endswith(".txt"))]
    var1, var2 = my_files[0].split("_")[0:2]
    values_1 = []
    values_2 = []
    for my_file in my_files:
        my_file = my_file.replace(".txt", "")
        _, _, val1, val2 = my_file.split("_")
        values_1.append(val1)
        values_2.append(val2)
    key = lambda x: float(x)
    reverse1 = np.min([float(v) for v in values_1]) < 1.0
    reverse2 = np.min([float(v) for v in values_2]) < 1.0
    values_1 = sorted(np.unique(values_1), key=key, reverse=reverse1)
    values_2 = sorted(np.unique(values_2), key=key, reverse=reverse2)

    df = pd.DataFrame(index = values_1, columns = values_2)

    for my_file in my_files:
        with open(os.path.join(mypath, my_file), 'r') as f:
            it_num = f.read()  
        my_file = my_file.replace(".txt", "")
        _, _, val1, val2 = my_file.split("_")
        df[val2][val1] = it_num
    df = df.T
    df.columns.name = f"{var2}/{var1}"
    df = df.fillna("( 0) -- ")
    print(df)
    print(df.to_latex())
    df.to_csv(os.path.join(mypath, "it_numbers.csv"))
