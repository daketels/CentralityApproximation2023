import matplotlib.pyplot as plt

# Phi = 0.1
prefix_mc = "phi_harmonic_mc"
prefix_phi = "phi_harmonic_0-1"

avg_ae_0_1 = 0.0
avg_scc_0_1 = 0.0

max_ae_0_1 = 0.0
min_scc_0_1 = 1.0
max_scc_0_1 = 0.0

run_0_1 = 0.0
run_0_2 = 0.0
run_0_3 = 0.0
run_0_4 = 0.0
run_0_5 = 0.0
run_0_6 = 0.0
run_0_7 = 0.0
run_0_8 = 0.0
run_0_9 = 0.0
run_1_0 = 0.0

for i in range(1,21):
    ba_beta_mc_file = open("./" + prefix_mc + "/G_ba_beta_{}.out".format(i))
    ba_beta_mc_vals = ba_beta_mc_file.readlines()[1]


    ba_beta_phi_file = open("./" + prefix_phi + "/G_ba_beta_{}.out".format(i))
    ba_beta_phi_lines = ba_beta_phi_file.readlines()
    ba_beta_phi_vals = ba_beta_phi_lines[0]

    run_0_1 += float(ba_beta_phi_lines[3].split()[1]) / 20.0

    ba_beta_mc_file.close()
    ba_beta_phi_file.close()

    ba_beta_mc_vals = ba_beta_mc_vals.split(",")
    ba_beta_phi_vals = ba_beta_phi_vals.split(",")

    ba_beta_mc_vals = [float(elem) for elem in ba_beta_mc_vals]
    ba_beta_phi_vals = [float(elem) for elem in ba_beta_phi_vals]

    acc_error = 0.0

    for j in range(500):
        error = abs(ba_beta_mc_vals[j]-ba_beta_phi_vals[j])
        avg_ae_0_1 += error / 500.0
        acc_error += error / 500.0

    max_ae_0_1 = max(max_ae_0_1, acc_error)

    L_mc = [(ba_beta_mc_vals[j], j+1) for j in range(500)]
    L_psp = [(ba_beta_phi_vals[j], j+1) for j in range(500)]

    L_mc.sort()
    L_psp.sort()

    Perm_mc = [elem[1] for elem in L_mc]
    Perm_psp = [elem[1] for elem in L_psp]

    scc_sum = sum([(elem1-elem2)**2  for elem1, elem2 in zip(Perm_mc, Perm_psp)])
    scc = 1.0 - (6*scc_sum)/(500*((500**2)-1))

    avg_scc_0_1 += scc
    min_scc_0_1 = min(scc, min_scc_0_1)
    max_scc_0_1 = max(scc, max_scc_0_1)

    #print("for G_ba_beta_{}".format(i) + ", phi = 0.1:\n")
    #print("\tmax ae: " + str(max_ae_0_1))
    #print("\tmin ae: " + str(min_ae_0_1))
    #print("\tavg ae: " + str(avg_ae_0_1))
avg_ae_0_1 /= 20.0
avg_scc_0_1 /= 20.0


# Phi = 0.2
prefix_mc = "phi_harmonic_mc"
prefix_phi = "phi_harmonic_0-2"

avg_ae_0_2 = 0.0
avg_scc_0_2 = 0.0

max_ae_0_2 = 0.0
min_scc_0_2 = 1.0
max_scc_0_2 = 0.0

for i in range(1,21):
    ba_beta_mc_file = open("./" + prefix_mc + "/G_ba_beta_{}.out".format(i))
    ba_beta_mc_vals = ba_beta_mc_file.readlines()[1]

    ba_beta_phi_file = open("./" + prefix_phi + "/G_ba_beta_{}.out".format(i))

    ba_beta_phi_lines = ba_beta_phi_file.readlines()
    ba_beta_phi_vals = ba_beta_phi_lines[0]

    run_0_2 += float(ba_beta_phi_lines[3].split()[1]) / 20.0

    ba_beta_mc_file.close()
    ba_beta_phi_file.close()

    ba_beta_mc_vals = ba_beta_mc_vals.split(",")
    ba_beta_phi_vals = ba_beta_phi_vals.split(",")

    ba_beta_mc_vals = [float(elem) for elem in ba_beta_mc_vals]
    ba_beta_phi_vals = [float(elem) for elem in ba_beta_phi_vals]

    acc_error = 0.0

    for j in range(500):
        error = abs(ba_beta_mc_vals[j]-ba_beta_phi_vals[j])
        avg_ae_0_2 += error / 500.0
        acc_error += error / 500.0
        
    max_ae_0_2 = max(acc_error, max_ae_0_2)

    L_mc = [(ba_beta_mc_vals[j], j+1) for j in range(500)]
    L_psp = [(ba_beta_phi_vals[j], j+1) for j in range(500)]

    L_mc.sort()
    L_psp.sort()

    Perm_mc = [elem[1] for elem in L_mc]
    Perm_psp = [elem[1] for elem in L_psp]

    scc_sum = sum([(elem1-elem2)**2  for elem1, elem2 in zip(Perm_mc, Perm_psp)])
    scc = 1.0 - (6*scc_sum)/(500*((500**2)-1))

    avg_scc_0_2 += scc
    min_scc_0_2 = min(scc, min_scc_0_2)
    max_scc_0_2 = max(scc, max_scc_0_2)

avg_ae_0_2 /= 20.0
avg_scc_0_2 /= 20.0




# Phi = 0.3
prefix_mc = "phi_harmonic_mc"
prefix_phi = "phi_harmonic_0-3"

avg_ae_0_3 = 0.0
avg_scc_0_3 = 0.0

max_ae_0_3 = 0.0
min_scc_0_3 = 1.0
max_scc_0_3 = 0.0

for i in range(1,21):
    ba_beta_mc_file = open("./" + prefix_mc + "/G_ba_beta_{}.out".format(i))
    ba_beta_mc_vals = ba_beta_mc_file.readlines()[1]

    ba_beta_phi_file = open("./" + prefix_phi + "/G_ba_beta_{}.out".format(i))

    ba_beta_phi_lines = ba_beta_phi_file.readlines()
    ba_beta_phi_vals = ba_beta_phi_lines[0]

    run_0_3 += float(ba_beta_phi_lines[3].split()[1]) / 20.0

    ba_beta_mc_file.close()
    ba_beta_phi_file.close()

    ba_beta_mc_vals = ba_beta_mc_vals.split(",")
    ba_beta_phi_vals = ba_beta_phi_vals.split(",")

    ba_beta_mc_vals = [float(elem) for elem in ba_beta_mc_vals]
    ba_beta_phi_vals = [float(elem) for elem in ba_beta_phi_vals]

    acc_error = 0.0

    for j in range(500):
        error = abs(ba_beta_mc_vals[j]-ba_beta_phi_vals[j])
        avg_ae_0_3 += error / 500.0
        acc_error += error / 500.0

    max_ae_0_3 = max(acc_error, max_ae_0_3)

    L_mc = [(ba_beta_mc_vals[j], j+1) for j in range(500)]
    L_psp = [(ba_beta_phi_vals[j], j+1) for j in range(500)]

    L_mc.sort()
    L_psp.sort()

    Perm_mc = [elem[1] for elem in L_mc]
    Perm_psp = [elem[1] for elem in L_psp]

    scc_sum = sum([(elem1-elem2)**2  for elem1, elem2 in zip(Perm_mc, Perm_psp)])
    scc = 1.0 - (6*scc_sum)/(500*((500**2)-1))

    avg_scc_0_3 += scc
    min_scc_0_3 = min(scc, min_scc_0_3)
    max_scc_0_3 = max(scc, max_scc_0_3)

avg_ae_0_3 /= 20.0
avg_scc_0_3 /= 20.0



# Phi = 0.4
prefix_mc = "phi_harmonic_mc"
prefix_phi = "phi_harmonic_0-4"

avg_ae_0_4 = 0.0
avg_scc_0_4 = 0.0

max_ae_0_4 = 0.0
min_scc_0_4 = 1.0
max_scc_0_4 = 0.0

for i in range(1,21):
    ba_beta_mc_file = open("./" + prefix_mc + "/G_ba_beta_{}.out".format(i))
    ba_beta_mc_vals = ba_beta_mc_file.readlines()[1]

    ba_beta_phi_file = open("./" + prefix_phi + "/G_ba_beta_{}.out".format(i))

    ba_beta_phi_lines = ba_beta_phi_file.readlines()
    ba_beta_phi_vals = ba_beta_phi_lines[0]

    run_0_4 += float(ba_beta_phi_lines[3].split()[1]) / 20.0

    ba_beta_mc_file.close()
    ba_beta_phi_file.close()

    ba_beta_mc_vals = ba_beta_mc_vals.split(",")
    ba_beta_phi_vals = ba_beta_phi_vals.split(",")

    ba_beta_mc_vals = [float(elem) for elem in ba_beta_mc_vals]
    ba_beta_phi_vals = [float(elem) for elem in ba_beta_phi_vals]

    acc_error = 0.0

    for j in range(500):
        error = abs(ba_beta_mc_vals[j]-ba_beta_phi_vals[j])
        avg_ae_0_4 += error / 500.0
        acc_error += error / 500.0
    
    max_ae_0_4 = max(acc_error, max_ae_0_4)

    L_mc = [(ba_beta_mc_vals[j], j+1) for j in range(500)]
    L_psp = [(ba_beta_phi_vals[j], j+1) for j in range(500)]

    L_mc.sort()
    L_psp.sort()

    Perm_mc = [elem[1] for elem in L_mc]
    Perm_psp = [elem[1] for elem in L_psp]

    scc_sum = sum([(elem1-elem2)**2  for elem1, elem2 in zip(Perm_mc, Perm_psp)])
    scc = 1.0 - (6*scc_sum)/(500*((500**2)-1))

    avg_scc_0_4 += scc
    min_scc_0_4 = min(scc, min_scc_0_4)
    max_scc_0_4 = max(scc, max_scc_0_4)

avg_ae_0_4 /= 20.0
avg_scc_0_4 /= 20.0



# Phi = 0.5
prefix_mc = "phi_harmonic_mc"
prefix_phi = "phi_harmonic_0-5"

avg_ae_0_5 = 0.0
avg_scc_0_5 = 0.0

max_ae_0_5 = 0.0
min_scc_0_5 = 1.0
max_scc_0_5 = 0.0

for i in range(1,21):
    ba_beta_mc_file = open("./" + prefix_mc + "/G_ba_beta_{}.out".format(i))
    ba_beta_mc_vals = ba_beta_mc_file.readlines()[1]

    ba_beta_phi_file = open("./" + prefix_phi + "/G_ba_beta_{}.out".format(i))

    ba_beta_phi_lines = ba_beta_phi_file.readlines()
    ba_beta_phi_vals = ba_beta_phi_lines[0]

    run_0_5 += float(ba_beta_phi_lines[3].split()[1]) / 20.0

    ba_beta_mc_file.close()
    ba_beta_phi_file.close()

    ba_beta_mc_vals = ba_beta_mc_vals.split(",")
    ba_beta_phi_vals = ba_beta_phi_vals.split(",")

    ba_beta_mc_vals = [float(elem) for elem in ba_beta_mc_vals]
    ba_beta_phi_vals = [float(elem) for elem in ba_beta_phi_vals]

    acc_error = 0.0

    for j in range(500):
        error = abs(ba_beta_mc_vals[j]-ba_beta_phi_vals[j])
        avg_ae_0_5 += error / 500.0
        acc_error += error / 500.0
    
    max_ae_0_5 = max(acc_error, max_ae_0_5)

    L_mc = [(ba_beta_mc_vals[j], j+1) for j in range(500)]
    L_psp = [(ba_beta_phi_vals[j], j+1) for j in range(500)]

    L_mc.sort()
    L_psp.sort()

    Perm_mc = [elem[1] for elem in L_mc]
    Perm_psp = [elem[1] for elem in L_psp]

    scc_sum = sum([(elem1-elem2)**2  for elem1, elem2 in zip(Perm_mc, Perm_psp)])
    scc = 1.0 - (6*scc_sum)/(500*((500**2)-1))

    avg_scc_0_5 += scc
    min_scc_0_5 = min(scc, min_scc_0_5)
    max_scc_0_5 = max(scc, max_scc_0_5)

avg_ae_0_5 /= 20.0
avg_scc_0_5 /= 20.0



# Phi = 0.6
prefix_mc = "phi_harmonic_mc"
prefix_phi = "phi_harmonic_0-6"

avg_ae_0_6 = 0.0
avg_scc_0_6 = 0.0

max_ae_0_6 = 0.0
min_scc_0_6 = 1.0
max_scc_0_6 = 0.0

for i in range(1,21):
    ba_beta_mc_file = open("./" + prefix_mc + "/G_ba_beta_{}.out".format(i))
    ba_beta_mc_vals = ba_beta_mc_file.readlines()[1]

    ba_beta_phi_file = open("./" + prefix_phi + "/G_ba_beta_{}.out".format(i))

    ba_beta_phi_lines = ba_beta_phi_file.readlines()
    ba_beta_phi_vals = ba_beta_phi_lines[0]

    run_0_6 += float(ba_beta_phi_lines[3].split()[1]) / 20.0

    ba_beta_mc_file.close()
    ba_beta_phi_file.close()

    ba_beta_mc_vals = ba_beta_mc_vals.split(",")
    ba_beta_phi_vals = ba_beta_phi_vals.split(",")

    ba_beta_mc_vals = [float(elem) for elem in ba_beta_mc_vals]
    ba_beta_phi_vals = [float(elem) for elem in ba_beta_phi_vals]

    acc_error = 0.0

    for j in range(500):
        error = abs(ba_beta_mc_vals[j]-ba_beta_phi_vals[j])
        avg_ae_0_6 += error / 500.0
        acc_error += error / 500.0
    
    max_ae_0_6 = max(acc_error, max_ae_0_6)

    L_mc = [(ba_beta_mc_vals[j], j+1) for j in range(500)]
    L_psp = [(ba_beta_phi_vals[j], j+1) for j in range(500)]

    L_mc.sort()
    L_psp.sort()

    Perm_mc = [elem[1] for elem in L_mc]
    Perm_psp = [elem[1] for elem in L_psp]

    scc_sum = sum([(elem1-elem2)**2  for elem1, elem2 in zip(Perm_mc, Perm_psp)])
    scc = 1.0 - (6*scc_sum)/(500*((500**2)-1))

    avg_scc_0_6 += scc
    min_scc_0_6 = min(scc, min_scc_0_6)
    max_scc_0_6 = max(scc, max_scc_0_6)

avg_ae_0_6 /= 20.0
avg_scc_0_6 /= 20.0



# Phi = 0.7
prefix_mc = "phi_harmonic_mc"
prefix_phi = "phi_harmonic_0-7"

avg_ae_0_7 = 0.0
avg_scc_0_7 = 0.0

max_ae_0_7 = 0.0
min_scc_0_7 = 1.0
max_scc_0_7 = 0.0

for i in range(1,21):
    ba_beta_mc_file = open("./" + prefix_mc + "/G_ba_beta_{}.out".format(i))
    ba_beta_mc_vals = ba_beta_mc_file.readlines()[1]

    ba_beta_phi_file = open("./" + prefix_phi + "/G_ba_beta_{}.out".format(i))

    ba_beta_phi_lines = ba_beta_phi_file.readlines()
    ba_beta_phi_vals = ba_beta_phi_lines[0]

    run_0_7 += float(ba_beta_phi_lines[3].split()[1]) / 20.0

    ba_beta_mc_file.close()
    ba_beta_phi_file.close()

    ba_beta_mc_vals = ba_beta_mc_vals.split(",")
    ba_beta_phi_vals = ba_beta_phi_vals.split(",")

    ba_beta_mc_vals = [float(elem) for elem in ba_beta_mc_vals]
    ba_beta_phi_vals = [float(elem) for elem in ba_beta_phi_vals]

    acc_error = 0.0

    for j in range(500):
        error = abs(ba_beta_mc_vals[j]-ba_beta_phi_vals[j])
        avg_ae_0_7 += error / 500.0
        acc_error += error / 500.0
        
    max_ae_0_7 = max(acc_error, max_ae_0_7)

    L_mc = [(ba_beta_mc_vals[j], j+1) for j in range(500)]
    L_psp = [(ba_beta_phi_vals[j], j+1) for j in range(500)]

    L_mc.sort()
    L_psp.sort()

    Perm_mc = [elem[1] for elem in L_mc]
    Perm_psp = [elem[1] for elem in L_psp]

    scc_sum = sum([(elem1-elem2)**2  for elem1, elem2 in zip(Perm_mc, Perm_psp)])
    scc = 1.0 - (6*scc_sum)/(500*((500**2)-1))

    avg_scc_0_7 += scc
    min_scc_0_7 = min(scc, min_scc_0_7)
    max_scc_0_7 = max(scc, max_scc_0_7)

avg_ae_0_7 /= 20.0
avg_scc_0_7 /= 20.0



# Phi = 0.8
prefix_mc = "phi_harmonic_mc"
prefix_phi = "phi_harmonic_0-8"

avg_ae_0_8 = 0.0
avg_scc_0_8 = 0.0

max_ae_0_8 = 0.0
min_scc_0_8 = 1.0
max_scc_0_8 = 0.0

for i in range(1,21):
    ba_beta_mc_file = open("./" + prefix_mc + "/G_ba_beta_{}.out".format(i))
    ba_beta_mc_vals = ba_beta_mc_file.readlines()[1]

    ba_beta_phi_file = open("./" + prefix_phi + "/G_ba_beta_{}.out".format(i))

    ba_beta_phi_lines = ba_beta_phi_file.readlines()
    ba_beta_phi_vals = ba_beta_phi_lines[0]

    run_0_8 += float(ba_beta_phi_lines[3].split()[1]) / 20.0

    ba_beta_mc_file.close()
    ba_beta_phi_file.close()

    ba_beta_mc_vals = ba_beta_mc_vals.split(",")
    ba_beta_phi_vals = ba_beta_phi_vals.split(",")

    ba_beta_mc_vals = [float(elem) for elem in ba_beta_mc_vals]
    ba_beta_phi_vals = [float(elem) for elem in ba_beta_phi_vals]

    acc_error = 0.0

    for j in range(500):
        error = abs(ba_beta_mc_vals[j]-ba_beta_phi_vals[j])
        avg_ae_0_8 += error / 500.0
        acc_error += error / 500.0
    
    max_ae_0_8 = max(acc_error, max_ae_0_8)

    L_mc = [(ba_beta_mc_vals[j], j+1) for j in range(500)]
    L_psp = [(ba_beta_phi_vals[j], j+1) for j in range(500)]

    L_mc.sort()
    L_psp.sort()

    Perm_mc = [elem[1] for elem in L_mc]
    Perm_psp = [elem[1] for elem in L_psp]

    scc_sum = sum([(elem1-elem2)**2  for elem1, elem2 in zip(Perm_mc, Perm_psp)])
    scc = 1.0 - (6*scc_sum)/(500*((500**2)-1))

    avg_scc_0_8 += scc
    min_scc_0_8 = min(scc, min_scc_0_8)
    max_scc_0_8 = max(scc, max_scc_0_8)

avg_ae_0_8 /= 20.0
avg_scc_0_8 /= 20.0



# Phi = 0.9
prefix_mc = "phi_harmonic_mc"
prefix_phi = "phi_harmonic_0-9"

avg_ae_0_9 = 0.0
avg_scc_0_9 = 0.0

max_ae_0_9 = 0.0
min_scc_0_9 = 1.0
max_scc_0_9 = 0.0

for i in range(1,21):
    ba_beta_mc_file = open("./" + prefix_mc + "/G_ba_beta_{}.out".format(i))
    ba_beta_mc_vals = ba_beta_mc_file.readlines()[1]

    ba_beta_phi_file = open("./" + prefix_phi + "/G_ba_beta_{}.out".format(i))
    
    ba_beta_phi_lines = ba_beta_phi_file.readlines()
    ba_beta_phi_vals = ba_beta_phi_lines[0]

    run_0_9 += float(ba_beta_phi_lines[3].split()[1]) / 20.0

    ba_beta_mc_file.close()
    ba_beta_phi_file.close()

    ba_beta_mc_vals = ba_beta_mc_vals.split(",")
    ba_beta_phi_vals = ba_beta_phi_vals.split(",")

    ba_beta_mc_vals = [float(elem) for elem in ba_beta_mc_vals]
    ba_beta_phi_vals = [float(elem) for elem in ba_beta_phi_vals]

    acc_error = 0.0

    for j in range(500):
        error = abs(ba_beta_mc_vals[j]-ba_beta_phi_vals[j])
        avg_ae_0_9 += error / 500.0
        acc_error += error / 500.0
        
    max_ae_0_9 = max(acc_error, max_ae_0_9)

    L_mc = [(ba_beta_mc_vals[j], j+1) for j in range(500)]
    L_psp = [(ba_beta_phi_vals[j], j+1) for j in range(500)]

    L_mc.sort()
    L_psp.sort()

    Perm_mc = [elem[1] for elem in L_mc]
    Perm_psp = [elem[1] for elem in L_psp]

    scc_sum = sum([(elem1-elem2)**2  for elem1, elem2 in zip(Perm_mc, Perm_psp)])
    scc = 1.0 - (6*scc_sum)/(500*((500**2)-1))

    avg_scc_0_9 += scc
    min_scc_0_9 = min(scc, min_scc_0_9)
    max_scc_0_9 = max(scc, max_scc_0_9)

avg_ae_0_9 /= 20.0
avg_scc_0_9 /= 20.0



# Phi = 1.0
prefix_mc = "phi_harmonic_mc"
prefix_phi = "phi_harmonic_1-0"

avg_ae_1_0 = 0.0
avg_scc_1_0 = 0.0

max_ae_1_0 = 0.0
min_scc_1_0 = 1.0
max_scc_1_0 = 0.0

for i in range(1,21):
    ba_beta_mc_file = open("./" + prefix_mc + "/G_ba_beta_{}.out".format(i))
    ba_beta_mc_vals = ba_beta_mc_file.readlines()[1]

    ba_beta_phi_file = open("./" + prefix_phi + "/G_ba_beta_{}.out".format(i))
    
    ba_beta_phi_lines = ba_beta_phi_file.readlines()
    ba_beta_phi_vals = ba_beta_phi_lines[0]

    run_1_0 += float(ba_beta_phi_lines[3].split()[1]) / 20.0

    ba_beta_mc_file.close()
    ba_beta_phi_file.close()

    ba_beta_mc_vals = ba_beta_mc_vals.split(",")
    ba_beta_phi_vals = ba_beta_phi_vals.split(",")

    ba_beta_mc_vals = [float(elem) for elem in ba_beta_mc_vals]
    ba_beta_phi_vals = [float(elem) for elem in ba_beta_phi_vals]

    acc_error = 0.0

    for j in range(500):
        error = abs(ba_beta_mc_vals[j]-ba_beta_phi_vals[j])
        avg_ae_1_0 += error / 500.0
        acc_error += error / 500.0
        
    max_ae_1_0 = max(acc_error, max_ae_1_0)

    L_mc = [(ba_beta_mc_vals[j], j+1) for j in range(500)]
    L_psp = [(ba_beta_phi_vals[j], j+1) for j in range(500)]

    L_mc.sort()
    L_psp.sort()

    Perm_mc = [elem[1] for elem in L_mc]
    Perm_psp = [elem[1] for elem in L_psp]

    scc_sum = sum([(elem1-elem2)**2  for elem1, elem2 in zip(Perm_mc, Perm_psp)])
    scc = 1.0 - (6*scc_sum)/(500*((500**2)-1))

    avg_scc_1_0 += scc
    min_scc_1_0 = min(scc, min_scc_1_0)
    max_scc_1_0 = max(scc, max_scc_1_0)

avg_ae_1_0 /= 20.0
avg_scc_1_0 /= 20.0


# Runtime VS PHI
f = plt.figure()
plt.rc("text", usetex=True)
plt.rc("font", family="serif")

x = ["0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"]
y = [run_0_1, run_0_2, run_0_3, run_0_4, run_0_5, run_0_6, run_0_7, run_0_8, run_0_9, run_1_0]

print("max runtime: ")
print([int(elem) for elem in y])

plt.ylabel("Time (ms)", size = 18)
plt.xlabel(r"$\phi$", size = 18)
plt.bar(x,y, width=0.3)
f.savefig("./plots/phi_ba_beta/phi_avg_runtime_ba_beta.pdf")
plt.close()


# AVG MAE VS PHI
f = plt.figure()
plt.rc("text", usetex=True)
plt.rc("font", family="serif")

x = ["0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"]
y = [avg_ae_0_1, avg_ae_0_2, avg_ae_0_3, avg_ae_0_4, avg_ae_0_5, avg_ae_0_6, avg_ae_0_7, avg_ae_0_8, avg_ae_0_9, avg_ae_1_0]
plt.axhline(min(y), color="red")

print("AVG MAE min: " + str(min(y)))

plt.ylabel("MAE", size = 18)
plt.xlabel(r"$\phi$", size = 18)
plt.bar(x,y, width=0.3)
f.savefig("./plots/phi_ba_beta/phi_avg_mae_ba_beta.pdf")
plt.close()

# MAX MAE VS PHI
f = plt.figure()
plt.rc("text", usetex=True)
plt.rc("font", family="serif")

x = ["0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"]
y = [max_ae_0_1, max_ae_0_2, max_ae_0_3, max_ae_0_4, max_ae_0_5, max_ae_0_6, max_ae_0_7, max_ae_0_8, max_ae_0_9, max_ae_1_0]
plt.axhline(min(y), color="red")

print("MAX MAE min: " + str(min(y)))

plt.ylabel("MAE", size = 18)
plt.xlabel(r"$\phi$", size = 18)
plt.bar(x,y, width=0.3)
f.savefig("./plots/phi_ba_beta/phi_max_mae_ba_beta.pdf")
plt.close()


# AVG SCC VS PHI
f = plt.figure()
plt.rc("text", usetex=True)
plt.rc("font", family="serif")
#plt.ylim(0.6, 0.7)
x = ["0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"]
y = [avg_scc_0_1, avg_scc_0_2, avg_scc_0_3, avg_scc_0_4, avg_scc_0_5, avg_scc_0_6, avg_scc_0_7, avg_scc_0_8, avg_scc_0_9, avg_scc_1_0]
#plt.axhline(max(y), color="red")

(print("AVG SCC max: " + str(max(y))))

plt.ylabel("SCC", size = 18)
plt.xlabel(r"$\phi$", size = 18)
plt.bar(x,y, width=0.3)
f.savefig("./plots/phi_ba_beta/phi_avg_scc_ba_beta.pdf")
plt.close()


# Max SCC VS PHI
f = plt.figure()
plt.rc("text", usetex=True)
plt.rc("font", family="serif")
#plt.ylim(0.65, 0.75)
x = ["0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"]
y = [max_scc_0_1, max_scc_0_2, max_scc_0_3, max_scc_0_4, max_scc_0_5, max_scc_0_6, max_scc_0_7, max_scc_0_8, max_scc_0_9, max_scc_1_0]
#plt.axhline(max(y), color="red")

(print("Max SCC max: " + str(max(y))))

plt.ylabel("SCC", size = 18)
plt.xlabel(r"$\phi$", size = 18)
plt.bar(x,y, width=0.3)
f.savefig("./plots/phi_ba_beta/phi_max_scc_ba_beta.pdf")
plt.close()



