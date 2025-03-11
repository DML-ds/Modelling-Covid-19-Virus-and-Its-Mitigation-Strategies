# Modelling-Covid-19-Virus-and-Its-Mitigation-Strategies


This project simulates the **COVID-19 spread and mitigation strategies** using the **SEIR (Susceptible-Exposed-Infectious-Recovered) model** in Python. The SEIR model is an epidemiological tool used to understand and predict disease dynamics.

## Overview

The model helps analyze the spread of COVID-19 in a population the size of the UK (\~67 million people) and explores potential mitigation strategies such as isolation, social distancing, and vaccination.

## SEIR Model Explanation

- **S (Susceptible):** Individuals who can contract the virus.
- **E (Exposed):** Individuals who have been infected but are not yet infectious.
- **I (Infectious):** Individuals who are actively spreading the virus.
- **R (Recovered):** Individuals who have recovered and gained immunity.

## Features

- Implements **basic SEIR differential equations** to simulate COVID-19 transmission.
- Uses **SciPy's **`` for solving differential equations numerically.
- Visualizes the **spread of the virus and the effects of interventions**.
- Can be modified to test different **mitigation strategies** (e.g., lockdowns, vaccination, mask-wearing).


## Parameters

The model is initialized with:

- **Total Population (N):** 67,000,000 (UK approximation)
- **Basic Reproduction Number (R0):** 2.4
- **Incubation Period:** 5.1 days
- **Infectious Period:** 3.3 days

## Plots Generated

- **Population Dynamics:** Shows susceptible and total infected individuals over time.
- **Disease States:** Displays exposed, infectious, and recovered individuals over time.

## Output

After running the script, you should see two graphs representing the disease spread and population response over 200 days.

 

## Contributions

Contributions are welcome! Feel free to fork this repository and submit pull requests for improvements.

## Contact

For questions or suggestions, reach out via email or open an issue in the repository.

