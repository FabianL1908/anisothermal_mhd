import itertools
import csv
import os
from pathlib import Path


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

def get_image_dict():
    image_dict = {}
    Path(os.path.join(os.getcwd(), 'images.csv')).touch(exist_ok=True)
    with open('images.csv', newline='') as csvfile:
        data = csv.reader(csvfile, delimiter=',')
        for row in data:
            image_dict[row[0]] = [r for r in row[1:]]
    return image_dict

def get_num_eigs():
    num_eigs_dict = {}
    Path(os.path.join(os.getcwd(), 'num_eigs.csv')).touch(exist_ok=True)
    with open('num_eigs.csv', newline='') as csvfile:
        data = csv.reader(csvfile, delimiter=',')
        for row in data:
            num_eigs_dict[row[0]] = row[1]
    return num_eigs_dict

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

def get_xybox(xy, xdist, ydist, pos):
    if int(pos) > 9:
        if pos[1] == "0":
            scale = float(f"0.{pos[2:]}")
        else:
            scale = int(pos[1:])
        xdist *= scale
        ydist *= scale
        pos = pos[0]
    if pos == "1": return (xy[0]-xdist, xy[1]+ydist)
    if pos == "2": return (xy[0], xy[1]+ydist)
    if pos == "3": return (xy[0]+xdist, xy[1]+ydist)
    if pos == "4": return (xy[0]-xdist, xy[1])
    if pos == "5": return (xy[0], xy[1])
    if pos == "6": return (xy[0]+xdist, xy[1])
    if pos == "7": return (xy[0]-xdist, xy[1]-ydist)
    if pos == "8": return (xy[0], xy[1]-ydist)
    if pos == "9": return (xy[0]+xdist, xy[1]-ydist)
    raise ValueError("pos has to be string in [1,9]")


def get_colors():
    return ('b', 'g', 'r', 'c', 'm', 'y', 'orange', 'g', 'r', 'c', 'm', 'y')

def get_linestyles():
    return ('-', '-', '-', '-', '-', '--', '--', '--', '--', '--', '--', '--', '--', '--')
    
