# Sensitivity Analysis

Sensitivity analysis on some sets of statistical tests of randomness

For each test, we will measure its ability to detect vulnerabilities. We will use three biased generators: they depend on a parameter that indicates the probability of giving bad results. We will carry out different experiments, varying the size and bias of each generator, to obtain the following results:

- Sensitivity frontier, defined as the bias and size values from which each test begins to find errors.
- Percentage of sensitivity, defined as the area (integral) of the region delimited by the previous boundary.

We will also show the percentage of failed tests of each test with these generators (significance α = 0.01)

Finally, we include the code used (in Python) to obtain these calculations so that anyone can use it in future analyzes.
This work was done by Brian Leiva, as part of his final year thesis in Mathematics at Universidad Complutense de Madrid. Some of these results are also shown in the article "Sensitivity and uniformity in statistical randomness tests" written by Marcos Brian Leiva Cerna, Elena Salomé Almaraz Luengo, Luis Javier García Villalba, and Julio Hernandez-Castro.
