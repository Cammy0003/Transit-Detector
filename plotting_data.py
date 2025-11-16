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
    y = data["y"]

    plt.plot(x, y)
    plt.xlabel("x")
    plt.ylabel("y")
    plt.title("Rust Vec<f64> vs Vec<f64>")
    plt.grid(True)
    plt.show()

if __name__ == "__main__":
    main()
