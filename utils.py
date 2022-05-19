import itertools
import csv

def get_branches():
    branch_dict = {}
    with open('branches.csv', newline='') as csvfile:
        data = csv.reader(csvfile, delimiter=',')
        for row in data:
            branch_dict[row[0]] = []
    with open('branches.csv', newline='') as csvfile:
        data = csv.reader(csvfile, delimiter=',')
        for row in data:
            int_row = [int(r) for r in row[1:]]
            branch_dict[row[0]].append(int_row)
    return branch_dict

def get_colors():
    return itertools.cycle(('b', 'g', 'r', 'c', 'm', 'y', 'k'))

def get_linestyles():
    return itertools.cycle(('-', '-', '-', '-', '-', '-', '-', '--', '--', '--', '--', '--', '--', '--', '-.', '-.', '-.', '-.', '-.', '-.', '-', ':'))
    
