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
    rh_uniform_mc_file = open("./" + prefix_mc + "/G_rh_uniform_{}.out".format(i))
    rh_uniform_mc_vals = rh_uniform_mc_file.readlines()[1]


    rh_uniform_phi_file = open("./" + prefix_phi + "/G_rh_uniform_{}.out".format(i))
    rh_uniform_phi_lines = rh_uniform_phi_file.readlines()
    rh_uniform_phi_vals = rh_uniform_phi_lines[0]

    run_0_1 += float(rh_uniform_phi_lines[3].split()[1]) / 20.0

    rh_uniform_mc_file.close()
    rh_uniform_phi_file.close()

    rh_uniform_mc_vals = rh_uniform_mc_vals.split(",")
    rh_uniform_phi_vals = rh_uniform_phi_vals.split(",")

    rh_uniform_mc_vals = [float(elem) for elem in rh_uniform_mc_vals]
    rh_uniform_phi_vals = [float(elem) for elem in rh_uniform_phi_vals]

    acc_rhror = 0.0

    for j in range(500):
        rhror = abs(rh_uniform_mc_vals[j]-rh_uniform_phi_vals[j])
        avg_ae_0_1 += rhror / 500.0
        acc_rhror += rhror / 500.0

    max_ae_0_1 = max(max_ae_0_1, acc_rhror)

    L_mc = [(rh_uniform_mc_vals[j], j+1) for j in range(500)]
    L_psp = [(rh_uniform_phi_vals[j], j+1) for j in range(500)]

    L_mc.sort()
    L_psp.sort()

    Prhm_mc = [elem[1] for elem in L_mc]
    Prhm_psp = [elem[1] for elem in L_psp]

    scc_sum = sum([(elem1-elem2)**2  for elem1, elem2 in zip(Prhm_mc, Prhm_psp)])
    scc = 1.0 - (6*scc_sum)/(500*((500**2)-1))

    avg_scc_0_1 += scc
    min_scc_0_1 = min(scc, min_scc_0_1)
    max_scc_0_1 = max(scc, max_scc_0_1)

    #print("for G_rh_uniform_{}".format(i) + ", phi = 0.1:\n")
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
    rh_uniform_mc_file = open("./" + prefix_mc + "/G_rh_uniform_{}.out".format(i))
    rh_uniform_mc_vals = rh_uniform_mc_file.readlines()[1]

    rh_uniform_phi_file = open("./" + prefix_phi + "/G_rh_uniform_{}.out".format(i))

    rh_uniform_phi_lines = rh_uniform_phi_file.readlines()
    rh_uniform_phi_vals = rh_uniform_phi_lines[0]

    run_0_2 += float(rh_uniform_phi_lines[3].split()[1]) / 20.0

    rh_uniform_mc_file.close()
    rh_uniform_phi_file.close()

    rh_uniform_mc_vals = rh_uniform_mc_vals.split(",")
    rh_uniform_phi_vals = rh_uniform_phi_vals.split(",")

    rh_uniform_mc_vals = [float(elem) for elem in rh_uniform_mc_vals]
    rh_uniform_phi_vals = [float(elem) for elem in rh_uniform_phi_vals]

    acc_rhror = 0.0

    for j in range(500):
        rhror = abs(rh_uniform_mc_vals[j]-rh_uniform_phi_vals[j])
        avg_ae_0_2 += rhror / 500.0
        acc_rhror += rhror / 500.0
        
    max_ae_0_2 = max(acc_rhror, max_ae_0_2)

    L_mc = [(rh_uniform_mc_vals[j], j+1) for j in range(500)]
    L_psp = [(rh_uniform_phi_vals[j], j+1) for j in range(500)]

    L_mc.sort()
    L_psp.sort()

    Prhm_mc = [elem[1] for elem in L_mc]
    Prhm_psp = [elem[1] for elem in L_psp]

    scc_sum = sum([(elem1-elem2)**2  for elem1, elem2 in zip(Prhm_mc, Prhm_psp)])
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
    rh_uniform_mc_file = open("./" + prefix_mc + "/G_rh_uniform_{}.out".format(i))
    rh_uniform_mc_vals = rh_uniform_mc_file.readlines()[1]

    rh_uniform_phi_file = open("./" + prefix_phi + "/G_rh_uniform_{}.out".format(i))

    rh_uniform_phi_lines = rh_uniform_phi_file.readlines()
    rh_uniform_phi_vals = rh_uniform_phi_lines[0]

    run_0_3 += float(rh_uniform_phi_lines[3].split()[1]) / 20.0

    rh_uniform_mc_file.close()
    rh_uniform_phi_file.close()

    rh_uniform_mc_vals = rh_uniform_mc_vals.split(",")
    rh_uniform_phi_vals = rh_uniform_phi_vals.split(",")

    rh_uniform_mc_vals = [float(elem) for elem in rh_uniform_mc_vals]
    rh_uniform_phi_vals = [float(elem) for elem in rh_uniform_phi_vals]

    acc_rhror = 0.0

    for j in range(500):
        rhror = abs(rh_uniform_mc_vals[j]-rh_uniform_phi_vals[j])
        avg_ae_0_3 += rhror / 500.0
        acc_rhror += rhror / 500.0

    max_ae_0_3 = max(acc_rhror, max_ae_0_3)

    L_mc = [(rh_uniform_mc_vals[j], j+1) for j in range(500)]
    L_psp = [(rh_uniform_phi_vals[j], j+1) for j in range(500)]

    L_mc.sort()
    L_psp.sort()

    Prhm_mc = [elem[1] for elem in L_mc]
    Prhm_psp = [elem[1] for elem in L_psp]

    scc_sum = sum([(elem1-elem2)**2  for elem1, elem2 in zip(Prhm_mc, Prhm_psp)])
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
    rh_uniform_mc_file = open("./" + prefix_mc + "/G_rh_uniform_{}.out".format(i))
    rh_uniform_mc_vals = rh_uniform_mc_file.readlines()[1]

    rh_uniform_phi_file = open("./" + prefix_phi + "/G_rh_uniform_{}.out".format(i))

    rh_uniform_phi_lines = rh_uniform_phi_file.readlines()
    rh_uniform_phi_vals = rh_uniform_phi_lines[0]

    run_0_4 += float(rh_uniform_phi_lines[3].split()[1]) / 20.0

    rh_uniform_mc_file.close()
    rh_uniform_phi_file.close()

    rh_uniform_mc_vals = rh_uniform_mc_vals.split(",")
    rh_uniform_phi_vals = rh_uniform_phi_vals.split(",")

    rh_uniform_mc_vals = [float(elem) for elem in rh_uniform_mc_vals]
    rh_uniform_phi_vals = [float(elem) for elem in rh_uniform_phi_vals]

    acc_rhror = 0.0

    for j in range(500):
        rhror = abs(rh_uniform_mc_vals[j]-rh_uniform_phi_vals[j])
        avg_ae_0_4 += rhror / 500.0
        acc_rhror += rhror / 500.0
    
    max_ae_0_4 = max(acc_rhror, max_ae_0_4)

    L_mc = [(rh_uniform_mc_vals[j], j+1) for j in range(500)]
    L_psp = [(rh_uniform_phi_vals[j], j+1) for j in range(500)]

    L_mc.sort()
    L_psp.sort()

    Prhm_mc = [elem[1] for elem in L_mc]
    Prhm_psp = [elem[1] for elem in L_psp]

    scc_sum = sum([(elem1-elem2)**2  for elem1, elem2 in zip(Prhm_mc, Prhm_psp)])
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
    rh_uniform_mc_file = open("./" + prefix_mc + "/G_rh_uniform_{}.out".format(i))
    rh_uniform_mc_vals = rh_uniform_mc_file.readlines()[1]

    rh_uniform_phi_file = open("./" + prefix_phi + "/G_rh_uniform_{}.out".format(i))

    rh_uniform_phi_lines = rh_uniform_phi_file.readlines()
    rh_uniform_phi_vals = rh_uniform_phi_lines[0]

    run_0_5 += float(rh_uniform_phi_lines[3].split()[1]) / 20.0

    rh_uniform_mc_file.close()
    rh_uniform_phi_file.close()

    rh_uniform_mc_vals = rh_uniform_mc_vals.split(",")
    rh_uniform_phi_vals = rh_uniform_phi_vals.split(",")

    rh_uniform_mc_vals = [float(elem) for elem in rh_uniform_mc_vals]
    rh_uniform_phi_vals = [float(elem) for elem in rh_uniform_phi_vals]

    acc_rhror = 0.0

    for j in range(500):
        rhror = abs(rh_uniform_mc_vals[j]-rh_uniform_phi_vals[j])
        avg_ae_0_5 += rhror / 500.0
        acc_rhror += rhror / 500.0
    
    max_ae_0_5 = max(acc_rhror, max_ae_0_5)

    L_mc = [(rh_uniform_mc_vals[j], j+1) for j in range(500)]
    L_psp = [(rh_uniform_phi_vals[j], j+1) for j in range(500)]

    L_mc.sort()
    L_psp.sort()

    Prhm_mc = [elem[1] for elem in L_mc]
    Prhm_psp = [elem[1] for elem in L_psp]

    scc_sum = sum([(elem1-elem2)**2  for elem1, elem2 in zip(Prhm_mc, Prhm_psp)])
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
    rh_uniform_mc_file = open("./" + prefix_mc + "/G_rh_uniform_{}.out".format(i))
    rh_uniform_mc_vals = rh_uniform_mc_file.readlines()[1]

    rh_uniform_phi_file = open("./" + prefix_phi + "/G_rh_uniform_{}.out".format(i))

    rh_uniform_phi_lines = rh_uniform_phi_file.readlines()
    rh_uniform_phi_vals = rh_uniform_phi_lines[0]

    run_0_6 += float(rh_uniform_phi_lines[3].split()[1]) / 20.0

    rh_uniform_mc_file.close()
    rh_uniform_phi_file.close()

    rh_uniform_mc_vals = rh_uniform_mc_vals.split(",")
    rh_uniform_phi_vals = rh_uniform_phi_vals.split(",")

    rh_uniform_mc_vals = [float(elem) for elem in rh_uniform_mc_vals]
    rh_uniform_phi_vals = [float(elem) for elem in rh_uniform_phi_vals]

    acc_rhror = 0.0

    for j in range(500):
        rhror = abs(rh_uniform_mc_vals[j]-rh_uniform_phi_vals[j])
        avg_ae_0_6 += rhror / 500.0
        acc_rhror += rhror / 500.0
    
    max_ae_0_6 = max(acc_rhror, max_ae_0_6)

    L_mc = [(rh_uniform_mc_vals[j], j+1) for j in range(500)]
    L_psp = [(rh_uniform_phi_vals[j], j+1) for j in range(500)]

    L_mc.sort()
    L_psp.sort()

    Prhm_mc = [elem[1] for elem in L_mc]
    Prhm_psp = [elem[1] for elem in L_psp]

    scc_sum = sum([(elem1-elem2)**2  for elem1, elem2 in zip(Prhm_mc, Prhm_psp)])
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
    rh_uniform_mc_file = open("./" + prefix_mc + "/G_rh_uniform_{}.out".format(i))
    rh_uniform_mc_vals = rh_uniform_mc_file.readlines()[1]

    rh_uniform_phi_file = open("./" + prefix_phi + "/G_rh_uniform_{}.out".format(i))

    rh_uniform_phi_lines = rh_uniform_phi_file.readlines()
    rh_uniform_phi_vals = rh_uniform_phi_lines[0]

    run_0_7 += float(rh_uniform_phi_lines[3].split()[1]) / 20.0

    rh_uniform_mc_file.close()
    rh_uniform_phi_file.close()

    rh_uniform_mc_vals = rh_uniform_mc_vals.split(",")
    rh_uniform_phi_vals = rh_uniform_phi_vals.split(",")

    rh_uniform_mc_vals = [float(elem) for elem in rh_uniform_mc_vals]
    rh_uniform_phi_vals = [float(elem) for elem in rh_uniform_phi_vals]

    acc_rhror = 0.0

    for j in range(500):
        rhror = abs(rh_uniform_mc_vals[j]-rh_uniform_phi_vals[j])
        avg_ae_0_7 += rhror / 500.0
        acc_rhror += rhror / 500.0
        
    max_ae_0_7 = max(acc_rhror, max_ae_0_7)

    L_mc = [(rh_uniform_mc_vals[j], j+1) for j in range(500)]
    L_psp = [(rh_uniform_phi_vals[j], j+1) for j in range(500)]

    L_mc.sort()
    L_psp.sort()

    Prhm_mc = [elem[1] for elem in L_mc]
    Prhm_psp = [elem[1] for elem in L_psp]

    scc_sum = sum([(elem1-elem2)**2  for elem1, elem2 in zip(Prhm_mc, Prhm_psp)])
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
    rh_uniform_mc_file = open("./" + prefix_mc + "/G_rh_uniform_{}.out".format(i))
    rh_uniform_mc_vals = rh_uniform_mc_file.readlines()[1]

    rh_uniform_phi_file = open("./" + prefix_phi + "/G_rh_uniform_{}.out".format(i))

    rh_uniform_phi_lines = rh_uniform_phi_file.readlines()
    rh_uniform_phi_vals = rh_uniform_phi_lines[0]

    run_0_8 += float(rh_uniform_phi_lines[3].split()[1]) / 20.0

    rh_uniform_mc_file.close()
    rh_uniform_phi_file.close()

    rh_uniform_mc_vals = rh_uniform_mc_vals.split(",")
    rh_uniform_phi_vals = rh_uniform_phi_vals.split(",")

    rh_uniform_mc_vals = [float(elem) for elem in rh_uniform_mc_vals]
    rh_uniform_phi_vals = [float(elem) for elem in rh_uniform_phi_vals]

    acc_rhror = 0.0

    for j in range(500):
        rhror = abs(rh_uniform_mc_vals[j]-rh_uniform_phi_vals[j])
        avg_ae_0_8 += rhror / 500.0
        acc_rhror += rhror / 500.0
    
    max_ae_0_8 = max(acc_rhror, max_ae_0_8)

    L_mc = [(rh_uniform_mc_vals[j], j+1) for j in range(500)]
    L_psp = [(rh_uniform_phi_vals[j], j+1) for j in range(500)]

    L_mc.sort()
    L_psp.sort()

    Prhm_mc = [elem[1] for elem in L_mc]
    Prhm_psp = [elem[1] for elem in L_psp]

    scc_sum = sum([(elem1-elem2)**2  for elem1, elem2 in zip(Prhm_mc, Prhm_psp)])
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
    rh_uniform_mc_file = open("./" + prefix_mc + "/G_rh_uniform_{}.out".format(i))
    rh_uniform_mc_vals = rh_uniform_mc_file.readlines()[1]

    rh_uniform_phi_file = open("./" + prefix_phi + "/G_rh_uniform_{}.out".format(i))
    
    rh_uniform_phi_lines = rh_uniform_phi_file.readlines()
    rh_uniform_phi_vals = rh_uniform_phi_lines[0]

    run_0_9 += float(rh_uniform_phi_lines[3].split()[1]) / 20.0

    rh_uniform_mc_file.close()
    rh_uniform_phi_file.close()

    rh_uniform_mc_vals = rh_uniform_mc_vals.split(",")
    rh_uniform_phi_vals = rh_uniform_phi_vals.split(",")

    rh_uniform_mc_vals = [float(elem) for elem in rh_uniform_mc_vals]
    rh_uniform_phi_vals = [float(elem) for elem in rh_uniform_phi_vals]

    acc_rhror = 0.0

    for j in range(500):
        rhror = abs(rh_uniform_mc_vals[j]-rh_uniform_phi_vals[j])
        avg_ae_0_9 += rhror / 500.0
        acc_rhror += rhror / 500.0
        
    max_ae_0_9 = max(acc_rhror, max_ae_0_9)

    L_mc = [(rh_uniform_mc_vals[j], j+1) for j in range(500)]
    L_psp = [(rh_uniform_phi_vals[j], j+1) for j in range(500)]

    L_mc.sort()
    L_psp.sort()

    Prhm_mc = [elem[1] for elem in L_mc]
    Prhm_psp = [elem[1] for elem in L_psp]

    scc_sum = sum([(elem1-elem2)**2  for elem1, elem2 in zip(Prhm_mc, Prhm_psp)])
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
    rh_uniform_mc_file = open("./" + prefix_mc + "/G_rh_uniform_{}.out".format(i))
    rh_uniform_mc_vals = rh_uniform_mc_file.readlines()[1]

    rh_uniform_phi_file = open("./" + prefix_phi + "/G_rh_uniform_{}.out".format(i))
    
    rh_uniform_phi_lines = rh_uniform_phi_file.readlines()
    rh_uniform_phi_vals = rh_uniform_phi_lines[0]

    run_1_0 += float(rh_uniform_phi_lines[3].split()[1]) / 20.0

    rh_uniform_mc_file.close()
    rh_uniform_phi_file.close()

    rh_uniform_mc_vals = rh_uniform_mc_vals.split(",")
    rh_uniform_phi_vals = rh_uniform_phi_vals.split(",")

    rh_uniform_mc_vals = [float(elem) for elem in rh_uniform_mc_vals]
    rh_uniform_phi_vals = [float(elem) for elem in rh_uniform_phi_vals]

    acc_rhror = 0.0

    for j in range(500):
        rhror = abs(rh_uniform_mc_vals[j]-rh_uniform_phi_vals[j])
        avg_ae_1_0 += rhror / 500.0
        acc_rhror += rhror / 500.0
        
    max_ae_1_0 = max(acc_rhror, max_ae_1_0)

    L_mc = [(rh_uniform_mc_vals[j], j+1) for j in range(500)]
    L_psp = [(rh_uniform_phi_vals[j], j+1) for j in range(500)]

    L_mc.sort()
    L_psp.sort()

    Prhm_mc = [elem[1] for elem in L_mc]
    Prhm_psp = [elem[1] for elem in L_psp]

    scc_sum = sum([(elem1-elem2)**2  for elem1, elem2 in zip(Prhm_mc, Prhm_psp)])
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

print("runtime: ")
print([int(elem) for elem in y])

plt.ylabel("Time (ms)", size = 18)
plt.xlabel(r"$\phi$", size = 18)
plt.bar(x,y, width=0.3)
f.savefig("./plots/phi_rh_uniform/phi_avg_runtime_rh_uniform.pdf")
plt.close()



# AVG MAE VS PHI
f = plt.figure()
plt.rc("text", usetex=True)
plt.rc("font", family="srhif")

x = ["0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"]
y = [avg_ae_0_1, avg_ae_0_2, avg_ae_0_3, avg_ae_0_4, avg_ae_0_5, avg_ae_0_6, avg_ae_0_7, avg_ae_0_8, avg_ae_0_9, avg_ae_1_0]
plt.axhline(min(y), color="red")

print("AVG MAE min: " + str(min(y)))

plt.ylabel("MAE", size = 18)
plt.xlabel(r"$\phi$", size = 18)
plt.bar(x,y, width=0.3)
f.savefig("./plots/phi_rh_uniform/phi_avg_mae_rh_uniform.pdf")
plt.close()

# MAX MAE VS PHI
f = plt.figure()
plt.rc("text", usetex=True)
plt.rc("font", family="srhif")

x = ["0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"]
y = [max_ae_0_1, max_ae_0_2, max_ae_0_3, max_ae_0_4, max_ae_0_5, max_ae_0_6, max_ae_0_7, max_ae_0_8, max_ae_0_9, max_ae_1_0]
plt.axhline(min(y), color="red")

print("MAX MAE min: " + str(min(y)))

plt.ylabel("MAE", size = 18)
plt.xlabel(r"$\phi$", size = 18)
plt.bar(x,y, width=0.3)
f.savefig("./plots/phi_rh_uniform/phi_max_mae_rh_uniform.pdf")
plt.close()


# AVG SCC VS PHI
f = plt.figure()
plt.rc("text", usetex=True)
plt.rc("font", family="srhif")
#plt.ylim(0.5, 0.6)
x = ["0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"]
y = [avg_scc_0_1, avg_scc_0_2, avg_scc_0_3, avg_scc_0_4, avg_scc_0_5, avg_scc_0_6, avg_scc_0_7, avg_scc_0_8, avg_scc_0_9, avg_scc_1_0]
#plt.axhline(max(y), color="red")

print("AVG SCC max: " + str(max(y)))

plt.ylabel("SCC", size = 18)
plt.xlabel(r"$\phi$", size = 18)
plt.bar(x,y, width=0.3)
f.savefig("./plots/phi_rh_uniform/phi_avg_scc_rh_uniform.pdf")
plt.close()

# Max SCC AE VS PHI
f = plt.figure()
plt.rc("text", usetex=True)
plt.rc("font", family="srhif")
#plt.ylim(0.65, 0.75)
x = ["0.1", "0.2", "0.3", "0.4", "0.5", "0.6", "0.7", "0.8", "0.9", "1.0"]
y = [max_scc_0_1, max_scc_0_2, max_scc_0_3, max_scc_0_4, max_scc_0_5, max_scc_0_6, max_scc_0_7, max_scc_0_8, max_scc_0_9, max_scc_1_0]
#plt.axhline(max(y), color="red")

print("MAX SCC max: " + str(max(y)))

plt.ylabel("SCC", size = 18)
plt.xlabel(r"$\phi$", size = 18)
plt.bar(x,y, width=0.3)
f.savefig("./plots/phi_rh_uniform/phi_max_scc_rh_uniform.pdf")
plt.close()



