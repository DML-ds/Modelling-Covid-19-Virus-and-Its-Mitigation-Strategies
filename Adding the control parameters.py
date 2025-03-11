# Constants
N = 67000000  # Total population of UK
R0 = 2.4  # Basic reproduction number
t_incubation = 5.1  # Incubation period in days
t_infective = 3.3  # Infectious period in days
alpha = 1/t_incubation  # Rate of progression from exposed to infectious
gamma = 1/t_infective  # Recovery rate
beta = R0*gamma  # Transmission rate
delta = 0.01  # Hospitalization rate from asymptomatic cases
mu = 0.001  # Death rate from hospitalized cases
alpha_a = 0.05  # Proportion of infections that become asymptomatic
rho = 0.02  # Average daily infection rate of asymptomatic cases

# Initial conditions
S0 = N - 100000  # Susceptible individuals
E0 = 10000  # Exposed individuals
I0 = 1000  # Infectious individuals
R0 = 0  # Recovered individuals
A0 = 0  # Asymptomatic individuals
H0 = 0  # Hospitalized individuals
D0 = 0  # Dead individuals

# Time points
t_span = (0, 200)  # Time span for simulation
t_eval = np.linspace(t_span[0], t_span[1], 200)  # Evaluation points

# Control parameters
u_contact_tracing = 0.8  # Effectiveness of contact tracing
u_quarantine = 0.7  # Effectiveness of quarantine measures
u_mask_wearing = 0.6  # Effectiveness of mask-wearing
u_vaccination = 0.9  # Effectiveness of vaccination (assuming 90% effective)

def extended_seir_model(t, y, beta, alpha, gamma, delta, mu, alpha_a, rho, u_contact_tracing, u_quarantine, u_mask_wearing, u_vaccination):
    S, E, I, R, A, H, D = y
    dSdt = -(1-u_contact_tracing*u_quarantine*u_mask_wearing*u_vaccination)*beta * S * (I + A) / N
    dEdt = (1-u_contact_tracing*u_quarantine*u_mask_wearing*u_vaccination)*beta * S * (I + A) / N - alpha * E
    dIdt = alpha * E - gamma * I - delta * I
    dRdt = gamma * I
    dAdt = alpha_a * I - gamma * A - rho * A
    dHdt = rho * A - mu * H
    ddDt = mu * H
    return [dSdt, dEdt, dIdt, dRdt, dAdt, dHdt, ddDt]

# Solve ODE
solution = solve_ivp(
    extended_seir_model,
    t_span,
    [S0, E0, I0, R0, A0, H0, D0],
    args=(beta, alpha, gamma, delta, mu, alpha_a, rho, u_contact_tracing, u_quarantine, u_mask_wearing, u_vaccination),
    t_eval=t_eval,
    method='LSODA',
    max_step=1
)

# Plot results
fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(10, 15))
ax1.plot(solution.t, solution.y[0], label='Susceptible')
ax1.plot(solution.t, solution.y[1] + solution.y[2] + solution.y[3] + solution.y[4], label='Infected', linestyle='--')
ax1.set_title('Population Dynamics')
ax1.legend()

ax2.plot(solution.t, solution.y[1], label='Exposed')
ax2.plot(solution.t, solution.y[2], label='Infectious', linestyle='--')
ax2.plot(solution.t, solution.y[3], label='Recovered', linestyle=':')
ax2.plot(solution.t, solution.y[4], label='Asymptomatic', linestyle=':')
ax2.set_title('Disease States')
ax2.legend()

ax3.plot(solution.t, solution.y[5], label='Hospitalized')
ax3.plot(solution.t, solution.y[6], label='Dead', linestyle='--')
ax3.set_title('Severe Cases')
ax3.legend()

plt.tight_layout()
plt.show()
