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

def get_rot_degree_dict():
    rot_degree_dict = {}
    with open('branches.csv', newline='') as csvfile:
        data = csv.reader(csvfile, delimiter=',')
        num_branches = 0
        for row in data:
            num_branches = int(row[0])
        for i in range(0, num_branches+1):
            rot_degree_dict[i] = 0
    with open('rot_degree.csv', newline='') as csvfile:
        data = csv.reader(csvfile, delimiter=',')
        for row in data:
            rot_degree_dict[int(row[0])] = int(row[1])
    return rot_degree_dict


def get_colors():
    return itertools.cycle(('b', 'g', 'r', 'c', 'm', 'y'))

def get_linestyles():
    return itertools.cycle(('-', '-', '-', '-', '-', '-', '-', '--', '--', '--', '--', '--', '--', '--', '-.', '-.', '-.', '-.', '-.', '-.', '-', ':'))
    
