import numpy as np
import matplotlib.pyplot as plt
from sympy import symbols, Eq, solve

# Define the equations
def line1(x):
    return (2*x - 6) / 3

def line2(x):
    return (42 - 6*x) / 7

def line3(x):
    return x + 3

# Solve the system of equations for intersections
x, y = symbols('x y')

# Equation 1: 2x - 3y = 6
eq1 = Eq(2*x - 3*y, 6)

# Equation 2: 6x + 7y = 42
eq2 = Eq(6*x + 7*y, 42)

# Equation 3: x - y = -3
eq3 = Eq(x - y, -3)

# Solving intersections

# Intersection of line1 and line2
sol1 = solve((eq1, eq2), (x, y))
# Intersection of line1 and line3
sol2 = solve((eq1, eq3), (x, y))
# Intersection of line2 and line3
sol3 = solve((eq2, eq3), (x, y))

# Extracting intersection points
intersection_1 = (sol1[x], sol1[y])
intersection_2 = (sol2[x], sol2[y])
intersection_3 = (sol3[x], sol3[y])

# Define points for intersections (with renamed points)
points = {
    'O (0; 0)': (0, 0),
    'A (0; 3)': (0, 3),
    'B (21/13, 60/13)': (21/13, 60/13),  # Renamed from C to B
    'C (5.25; 1.5)': (5.25, 1.5),  # Renamed from D to C
    'D (3; 0)': (3, 0),  # Renamed from E to D
}

# Define the points for the pentagon (in the new order)
pentagon_points = [
    points['O (0; 0)'],
    points['A (0; 3)'],
    points['B (21/13, 60/13)'],
    points['C (5.25; 1.5)'],
    points['D (3; 0)']
]

# Separate x and y values of pentagon points
pentagon_x = [p[0] for p in pentagon_points]
pentagon_y = [p[1] for p in pentagon_points]

# Close the pentagon loop by adding the first point at the end
pentagon_x.append(pentagon_x[0])
pentagon_y.append(pentagon_y[0])

# Generate x values for plotting
x_values = np.linspace(-1, 10, 400)

# Plot the lines
plt.figure(figsize=(8, 8))

# Plot each line
plt.plot(x_values, line1(x_values), label=r'$2x - 3y = 6$', color='blue')
plt.plot(x_values, line2(x_values), label=r'$6x + 7y = 42$', color='red')
plt.plot(x_values, line3(x_values), label=r'$x - y = -3$', color='green')

# Adding labels for the axes (creating Ox and Oy axes with thicker and bolder black lines)
plt.axhline(0, color='black',linewidth=2)  # Ox axis (horizontal), thicker line
plt.axvline(0, color='black',linewidth=2)  # Oy axis (vertical), thicker line

# Plot the points on the graph
for point, coords in points.items():
    plt.plot(coords[0], coords[1], 'o')
    plt.text(coords[0] + 0.1, coords[1], f'{point}', fontsize=12)

# Plot the pentagon and fill it with color
plt.fill(pentagon_x, pentagon_y, color='lightblue', alpha=0.5, label="Pentagon Area")

# Plot intersection points
plt.scatter(intersection_1[0], intersection_1[1], color='orange', label=f'Intersection 1 {intersection_1}')
plt.scatter(intersection_2[0], intersection_2[1], color='purple', label=f'Intersection 2 {intersection_2}')
plt.scatter(intersection_3[0], intersection_3[1], color='brown', label=f'Intersection 3 {intersection_3}')

# Set graph limits
plt.xlim(-1, 10)
plt.ylim(-5, 10)

# Labeling the axes
plt.xlabel(r'$x$')
plt.ylabel(r'$y$')

# Add title
plt.title('Graph of the lines with intersections on Ox and Oy')

# Add a legend
plt.legend()

# Show grid and plot
plt.grid(True)
plt.show()
