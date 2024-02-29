import os
import subprocess

if __name__ == "__main__":
    tiny_test_set = "tiny_test_set"
    filenames = os.listdir(tiny_test_set)
    filenames.sort()
    for filename in filenames:
        print(filename)
        instance_filepath = tiny_test_set + "/" + filename
        solution_filepath = tiny_test_set + "_solutions/" + filename[:-3] + ".sol"
        output_filepath = tiny_test_set + "_nof_crossings/" + filename[:-3] + ".txt"
        with open(output_filepath, "w") as f:
            subprocess.run(["pace2024verifier", "-c", instance_filepath, solution_filepath], stdout=f)
