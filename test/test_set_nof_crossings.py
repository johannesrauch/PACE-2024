import os
import subprocess

if __name__ == "__main__":
    test_set = "medium_test_set"
    filenames = os.listdir(test_set)
    filenames.sort()
    for filename in filenames:
        print(filename)
        instance_filepath = test_set + "/" + filename
        solution_filepath = test_set + "/solutions/" + filename[:-3] + ".sol"
        output_filepath = test_set + "/nof_crossings/" + filename[:-3] + ".txt"
        with open(output_filepath, "w") as f:
            subprocess.run(["pace2024verifier", "-c", instance_filepath, solution_filepath], stdout=f)
