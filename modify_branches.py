import os
output_dir = os.path.join(os.getcwd(), "output")

dirs = os.listdir(output_dir)

del_branch = 59

for d in dirs:
#    import ipdb; ipdb.set_trace()
    try:
        Ra = d.split("@")[0].split("=")[1]
        Ra = float(Ra)
        if Ra < 47685.61:
            my_dir = os.path.join(output_dir, d)
            for f in os.listdir(my_dir):
                if f == f"solution-{del_branch}.h5" or f == f"functional-{del_branch}.txt":
                    os.remove(os.path.join(my_dir, f))
                    print(d, f)
    except:
        pass
