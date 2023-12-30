import random

if __name__ == "__main__":
    n0 = 10
    n1 = 10

    edges = []
    for i in range(n0):
        threshold = random.random()
        for j in range(n1):
            if random.random() < threshold:
                edges.append((i + 1, j + n0 + 1))
    
    print(f"p ocr {n0} {n1} {len(edges)}")
    for i, j in edges:
        print(f"{i} {j}")
