import matplotlib.pyplot as plt

def read_data(filename):
    sizes = []
    times = []

    with open(filename, "r") as f:
        next(f)
        for line in f:
            parts = line.strip().split()
            if len(parts) != 2:
                continue  # ignore incomplete lines
            size, time = parts
            try:
                sizes.append(int(size))
                times.append(float(time))
            except ValueError:
                continue

    return sizes, times

def main():
    input_file = "size.txt"
    output_file = "timeVsSize.png"

    sizes, times = read_data(input_file)

    plt.figure()
    plt.plot(sizes, times, marker='o')
    plt.xlabel("Size (k)")
    plt.ylabel("Time (s)")
    plt.title("Computation Time vs Size")
    plt.grid(True)

    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Figure saved as {output_file}")


if __name__ == "__main__":
    main()