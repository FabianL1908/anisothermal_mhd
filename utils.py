import itertools
import csv

def get_branches():
    branches = []
    with open('branches.csv', newline='') as csvfile:
        data = csv.reader(csvfile, delimiter=',')
        for row in data:
            int_row = [int(r) for r in row]
            branches.append(int_row)
    return branches

def get_colors():
    return itertools.cycle(('b', 'g', 'r', 'c', 'm', 'y', 'k'))
