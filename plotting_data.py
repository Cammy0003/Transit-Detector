import sys
import json
import matplotlib.pyplot as plt

def main():
    raw = sys.stdin.read()
    if not raw.strip():
        print("No data received.")
        return

    data = json.loads(raw)

    x = data["x"]
    x_label = data["x_label"]
    y = data["y"]
    y_label = data["y_label"]

    plt.plot(x, y)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.title(f'{y_label} vs. {x_label}')
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    main()
