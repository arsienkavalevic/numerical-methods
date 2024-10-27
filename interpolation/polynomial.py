import numpy as np
import matplotlib.pyplot as plt

def f1(x):
    return np.exp(np.cos(x))

def f2(x):
    return np.abs(x * np.abs(x) - 1)

def equidistant_nodes(a, b, n):
    return np.linspace(a, b, n)

def chebyshev_nodes(a, b, n):
    nodes = [(a + b) / 2 + (b - a) / 2 * np.cos((2 * i - 1) / (2 * n) * np.pi) for i in range(1, n + 1)]
    return np.array(nodes)

def newton_interpolation(x, nodes, values):
    n = len(nodes)
    coefficients = np.zeros(n)
    coefficients[0] = values[0]

    for j in range(1, n):
        for i in range(n - 1, j - 1, -1):
            values[i] = (values[i] - values[i - 1]) / (nodes[i] - nodes[i - j])
            coefficients[j] += values[i] * np.prod(x - nodes[:i])

    return np.sum(coefficients)

def plot_interpolation(f, nodes_func, nodes_type, a, b, n_values):
    x_vals = np.linspace(a, b, 1000)
    f_vals = f(x_vals)

    plt.figure(figsize=(12, 6))

    for n in n_values:
        nodes = nodes_func(a, b, n)
        values = f(nodes)

        interpolation_vals = [newton_interpolation(x, nodes, values) for x in x_vals]

        plt.plot(x_vals, interpolation_vals, label=f'{nodes_type} nodes (n={n})')

    plt.plot(x_vals, f_vals, label='Original function')
    plt.title(f'Interpolation using {nodes_type} nodes')
    plt.xlabel('x')
    plt.ylabel('y')
    plt.legend()
    plt.grid(True)
    plt.show()

# Определение параметров
a = -3
b = 3
n_values = [3, 10, 20]  # Разные степени интерполяционного многочлена

fig, axs = plt.subplots(5)
axs[0].plot([x for x in np.linspace(-3, 3, 100)], [f1(x) for x in np.linspace(-3, 3, 100)], label='First function')
axs[0].plot([x for x in np.linspace(-3, 3, 100)], [f2(x) for x in np.linspace(-3, 3, 100)], label='Second function')
axs[0].set_title('Functions')
axs[0].legend()
eq_nodes_3 = equidistant_nodes(a, b, 3)
eq_nodes_10 = equidistant_nodes(a, b, 10)
eq_nodes_20 = equidistant_nodes(a, b, 20)
eq_f1_values_3 = f1(eq_nodes_3)
eq_f1_values_10 = f1(eq_nodes_10)
eq_f1_values_20 = f1(eq_nodes_20)
eq_f2_values_3 = f2(eq_nodes_3)
eq_f2_values_10 = f2(eq_nodes_10)
eq_f2_values_20 = f2(eq_nodes_20)
ch_nodes_3 = chebyshev_nodes(a, b, 3)
ch_nodes_10 = chebyshev_nodes(a, b, 10)
ch_nodes_20 = chebyshev_nodes(a, b, 20)
ch_f1_values_3 = f1(ch_nodes_3)
ch_f1_values_10 = f1(ch_nodes_10)
ch_f1_values_20 = f1(ch_nodes_20)
ch_f2_values_3 = f2(ch_nodes_3)
ch_f2_values_10 = f2(ch_nodes_10)
ch_f2_values_20 = f2(ch_nodes_20)
axs[1].plot([x for x in np.linspace(-3, 3, 100)], [newton_interpolation(x, eq_nodes_3, eq_f1_values_3) for x in np.linspace(-3, 3, 100)], label='eq_f1_3')
axs[1].plot([x for x in np.linspace(-3, 3, 100)], [newton_interpolation(x, eq_nodes_10, eq_f1_values_10) for x in np.linspace(-3, 3, 100)], label='eq_f1_10')
axs[1].plot([x for x in np.linspace(-3, 3, 100)], [newton_interpolation(x, eq_nodes_20, eq_f1_values_20) for x in np.linspace(-3, 3, 100)], label='eq_f1_20')
axs[2].plot([x for x in np.linspace(-3, 3, 100)], [newton_interpolation(x, ch_nodes_3, ch_f1_values_3) for x in np.linspace(-3, 3, 100)], label='ch_f1_3')
axs[2].plot([x for x in np.linspace(-3, 3, 100)], [newton_interpolation(x, ch_nodes_10, ch_f1_values_10) for x in np.linspace(-3, 3, 100)], label='ch_f1_10')
axs[2].plot([x for x in np.linspace(-3, 3, 100)], [newton_interpolation(x, ch_nodes_20, ch_f1_values_20) for x in np.linspace(-3, 3, 100)], label='ch_f1_20')
axs[3].plot([x for x in np.linspace(-3, 3, 100)], [newton_interpolation(x, eq_nodes_3, eq_f2_values_3) for x in np.linspace(-3, 3, 100)], label='eq_f2_3')
axs[3].plot([x for x in np.linspace(-3, 3, 100)], [newton_interpolation(x, eq_nodes_10, eq_f2_values_10) for x in np.linspace(-3, 3, 100)], label='eq_f2_10')
axs[3].plot([x for x in np.linspace(-3, 3, 100)], [newton_interpolation(x, eq_nodes_20, eq_f2_values_20) for x in np.linspace(-3, 3, 100)], label='eq_f2_20')
axs[4].plot([x for x in np.linspace(-3, 3, 100)], [newton_interpolation(x, ch_nodes_3, ch_f2_values_3) for x in np.linspace(-3, 3, 100)], label='ch_f2_3')
axs[4].plot([x for x in np.linspace(-3, 3, 100)], [newton_interpolation(x, ch_nodes_10, ch_f2_values_10) for x in np.linspace(-3, 3, 100)], label='ch_f2_10')
axs[4].plot([x for x in np.linspace(-3, 3, 100)], [newton_interpolation(x, ch_nodes_20, ch_f2_values_20) for x in np.linspace(-3, 3, 100)], label='ch_f2_20')
