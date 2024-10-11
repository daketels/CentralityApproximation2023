import os
import networkit as nk

print("making...")
os.system("cd ../ && make all")

num_of_tests = 30
print_failed_result = True
tests_executed = 0
failed_count = 0

for i in range(1, num_of_tests+1):
    graph_path = "./testcases/testcase_{}/test_graph_{}.txt".format(i, i)
    config_file = open("./testcase_{}/test_config_{}.txt".format(i, i), "r")
    expected_result_file = open("./testcase_{}/test_expected_result_{}.txt".format(i, i), "r")

    config = config_file.read().splitlines()[0]
    expected_result_lines = [line for line in expected_result_file.read().splitlines() if len(line) > 0]

    mc_sample_count = int(config.split(" ")[1])
    graph_file = open("./testcase_{}/test_graph_{}.txt".format(i, i), "r")
    node_count = int(graph_file.read().splitlines()[0])

    print("running test number {}".format(i))
    tests_executed += 1
    os.system("cd ../ && ./UncertainCentrality {} {}".format(config, graph_path))
    current_result_path = "./testcase_{}/test_graph_{}_result.txt".format(i, i)
    current_result_file = open(current_result_path)

    if "distance" in config:
        current_result = [line for line in current_result_file.read().splitlines() if len(line) > 0]
        if len(current_result) != len(expected_result_lines):
            print("Test {} (distance) failed: expected {} lines in result, got {}".format(i,len(expected_result_lines),len(current_result)))
            failed_count += 1
        else:
            lines_failed = 0
            for line_index in range(len(expected_result_lines)):
                if current_result[line_index] != expected_result_lines[line_index]:
                    lines_failed += 1
                    print("Test {} (distance) failed in line {}".format(i, line_index))
                    if print_failed_result:
                        print("\texpected :" + expected_result_lines[line_index])
                        print("\tgot      :" + current_result[line_index])
            if lines_failed > 0:
                failed_count += 1
    else:
        current_result_lines = current_result_file.read().splitlines()
        current_psp_result = current_result_lines[0]
        expected_psp_result = expected_result_lines[0]
        if current_psp_result != expected_psp_result:
            print("test {} failed".format(i))
            if print_failed_result:
                print("\tresult:   {}".format(expected_psp_result))
                print("\texpected: {}".format(expected_psp_result))
            print("skipping monte carlo test for this testcase")
            failed_count += 1
            continue
        if (mc_sample_count == 0):
            continue
        current_mc_result = [float(elem) for elem in current_result_lines[1].split(",") if len(elem) > 0]
        mc_expected_result = [0.0 for i in range(node_count)]
        line_index = 0
        for j in range(mc_sample_count):
            G = nk.graph.Graph(node_count)
            while (line_index < len(current_result_lines)) and (not "sample" in current_result_lines[line_index]):
                line_index +=1
            line_index += 1
            while (line_index < len(current_result_lines)) and ("(" in current_result_lines[line_index]):
                if ("NOT" in current_result_lines[line_index]) or (current_result_lines[line_index].strip()[0] != "(") :
                    line_index += 1
                    continue
                first_node = int(current_result_lines[line_index].split("(")[1].split(",")[0].strip())
                second_node = int(current_result_lines[line_index].split(",")[1].split(")")[0].strip())
                line_index += 1
                G.addEdge(first_node, second_node)
            if "harmonic" in config:
                BC = nk.centrality.HarmonicCloseness(G, True)
            else:
                BC = nk.centrality.Betweenness(G, True)
            BC.run()

            for node in range(node_count):
                mc_expected_result[node] += BC.scores()[node] / float(mc_sample_count)
        mc_expected_result = [ float("{}".format(round(elem,5))) for elem in mc_expected_result]
        for j in range(node_count):
            if abs(current_mc_result[j] - mc_expected_result[j]) > 0.000015:
                failed_count += 1
                print("\tmonte carlo failed for node {} in test {}".format(j,i))
                if print_failed_result:
                    print("\t\tresult:   .{}.".format(current_mc_result[j]))
                    print("\t\texpected: .{}.".format(mc_expected_result[j]))

    config_file.close()
    current_result_file.close()
    expected_result_file.close()

for i in range(1,num_of_tests+1):
    current_result_path = "./testcase_{}/test_graph_{}_result.txt".format(i, i)
    os.system("cp {} ./results/result_{}.txt".format(current_result_path, i))
    os.system("rm {}".format(current_result_path))

if failed_count == 0:
    print("all ({}) tests finished and all tests passed".format(tests_executed))
else:
    print("all ({}) tests finished and {} tests failed".format(tests_executed, failed_count))
