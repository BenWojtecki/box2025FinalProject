import matplotlib.pyplot as plt

def read_data(filename):
    sizes = []
    times = []
    lengths = []

    with open(filename, "r") as f:
        next(f)
        for line in f:
            parts = line.strip().split()
            if len(parts) != 3:
                continue  

            size, time, length = parts
            try:
                sizes.append(int(size))
                times.append(float(time))
                lengths.append(int(length))
            except ValueError:
                continue

    return sizes, times, lengths


def main():
    input_file = "size.txt"

    sizes, times, lengths = read_data(input_file)

    plt.figure()
    plt.plot(sizes, times, marker='o')
    plt.xlabel("Size (k)")
    plt.ylabel("Time (s)")
    plt.title("Computation Time in function of time")
    plt.grid(True)
    plt.savefig("sizeVTime.png", dpi=300, bbox_inches="tight")
    plt.close()

    plt.figure()
    plt.plot(sizes, lengths, marker='o')
    plt.xlabel("Size (k)")
    plt.ylabel("Number of MAWs")
    plt.title("Number of maws in function of size")
    plt.grid(True)
    plt.savefig("sizeVsLength.png", dpi=300, bbox_inches="tight")
    plt.close()



if __name__ == "__main__":
    main()
