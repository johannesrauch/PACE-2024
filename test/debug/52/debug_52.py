import re

def read_ordering(filename):
    l = []
    with open(filename) as f:
        for line in f:
            l.append(int(line) - 251)
    # print("read", l[:10])
    return l

def read_fixing(filename):
    d = {}
    with open(filename) as f:
        for line in f:
            split = re.sub("\s\s+", " ", line).split(" ")
            if split[0] == "fixed":
                d[int(split[2])] = int(split[4])
    print("read", list(d.items())[:10])
    return d

def get_variable_index(i, j):
    assert(i < j)
    n1 = 220
    n1_choose_2 = (n1 // 2) * 219
    offset = n1_choose_2 - (n1 - i) * (n1 - i - 1) // 2
    index = offset + j - i
    return index

def inverse(l):
    return [l.index(i) for i in range(len(l))]

if __name__ == "__main__":
    ref = inverse(read_ordering("52.sol"))
    test = read_ordering("52.txt")
    fix = read_fixing("52.log")

    for i in range(len(test)):
        for j in range(i + 1, len(test)):
            x = test[i]
            y = test[j]

            if ref[x] > ref[y]:
                if x < y:
                    j = get_variable_index(x, y)
                    if j in fix:
                        print(f"{j}: {x + 251} ({x}), {y + 251} ({y}), {fix[j]}")
                else:
                    j = get_variable_index(y, x)
                    if j in fix:
                        print(f"{j}: {y + 251} ({y}), {x + 251} ({x}), {fix[j]}")