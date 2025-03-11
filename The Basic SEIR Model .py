import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Constants
N = 67000000  # Total population of UK
R0 = 2.4  # Basic reproduction number
t_incubation = 5.1  # Incubation period in days
t_infective = 3.3  # Infectious period in days
alpha = 1/t_incubation  # Rate of progression from exposed to infectious
gamma = 1/t_infective  # Recovery rate
beta = R0*gamma  # Transmission rate

# Initial conditions
S0 = N - 100000  # Susceptible individuals
E0 = 10000  # Exposed individuals
I0 = 10  # Infectious individuals
R0 = 0  # Recovered individuals

# Time points
t_span = (0, 200)  # Time span for simulation
t_eval = np.linspace(t_span[0], t_span[1], 100)  # Evaluation points

# Differential equations
def seir_model(t, y, beta, alpha, gamma):
    S, E, I, R = y
    dSdt = -beta * S * I / N
    dEdt = beta * S * I / N - alpha * E
    dIdt = alpha * E - gamma * I
    dRdt = gamma * I
    return [dSdt, dEdt, dIdt, dRdt]

# Solve ODE
solution = solve_ivp(
    seir_model,
    t_span,
    [S0, E0, I0, R0],
    args=(beta, alpha, gamma),
    t_eval=t_eval,
    method='LSODA',
    max_step=1
)

# Plot results
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 10))
ax1.plot(solution.t, solution.y[0], label='Susceptible')
ax1.plot(solution.t, solution.y[1] + solution.y[2] + solution.y[3], label='Infected', linestyle='--')
ax1.set_title('Population Dynamics')
ax1.legend()

ax2.plot(solution.t, solution.y[1], label='Exposed')
ax2.plot(solution.t, solution.y[2], label='Infectious', linestyle='--')
ax2.plot(solution.t, solution.y[3], label='Recovered', linestyle=':')
ax2.set_title('Disease States')
ax2.legend()

plt.tight_layout()
plt.show()
