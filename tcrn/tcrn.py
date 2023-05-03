"""PySB version of the general three-component repressive network.

Adapted from the generic three-component repressive network model described 
Figure 1 of Mogilner et al. Developmental Cell 11, 279-287, 2006 (https://doi.org/10.1016/j.devcel.2006.08.004).

"""

# PySB components
from pysb import (
    Model,
    Monomer,
    Parameter,
    Initial,
    Rule,
    Observable,
    Expression,
    Annotation,
    ANY,
    WILD,
)

from pysb.macros import catalyze_state, equilibrate

## Initialize the model.
Model()

## Define the monomers.

# state: U = unphosphorylated, P = phosphorylated
Monomer("A", ["state"], {"state": ["U", "P"]})
Monomer("B", ["state"], {"state": ["U", "P"]})
Monomer("C", ["state"], {"state": ["U", "P"]})

## Define the initial states and amounts:

Parameter("Atot", 0.5) # uM
Parameter("Btot", 0.5) # uM
Parameter("Ctot", 0.5) # uM

Initial(A(state="U"), Atot)
Initial(B(state="U"), Btot)
Initial(C(state="U"), Ctot)

## Define the rules:

# Phosphorylation of A catalyzed by C:
# Michaelis-Menten Eq. for process:
#     A_U + C_U | A_U:C_U >> C_U + A_P
Parameter('k_d1', 1.33e-2) # uM/s
Parameter('K_d1', 1e-2) # uM
Observable("A_u", A(state="U"))
Observable("C_u", C(state="U"))
Expression('k_Au_Ap', k_d1*C_u / (K_d1 + A_u))
Rule("phospho_A", A(state="U") >> A(state="P"), k_Au_Ap)
# Dephosphorylation of A:
# Michaelis-Menten Eq. for process:
#     A_P >> A_U
Parameter('k_c1', 8.33e-3) # uM/s
Parameter('K_1', 5e-3) # uM
Observable("A_p", A(state="P"))
Expression('k_Ap_Au', k_c1 / (K_1 + A_p))
Rule("dephospho_A", A(state='P') >> A(state='U'), k_Ap_Au)


# Dephosphorylation of B catalyzed by A:
# Michaelis-Menten Eq. for process:
#     B_P + A_U | B_P:A_U >> B_U + A_U
Parameter('k_c2', 9.33e-3) # uM/s
Parameter('K_2', 1e-2) # uM
Observable("B_p", B(state="P"))
Expression('k_Bp_Bu', k_c2*A_u / (K_2 + B_p))
Rule("dephospho_B", B(state="P") >> B(state="U"), k_Bp_Bu)
# Phosphorylation of B:
# Michaelis-Menten Eq. for process:
#     B_U >> B_P
Parameter('k_d2', 8.33e-3) # uM/s
Parameter('K_d2', 5e-3) # uM
Observable("B_u", B(state="U"))
Expression('k_Bu_Bp', k_d2 / (K_d2 + B_u))
Rule("phospho_B", B(state='U') >> B(state='P'), k_Bu_Bp)

# Dephosphorylation of C catalyzed by B:
# Michaelis-Menten Eq. for process:
#     C_P + B_U | C_P:B_U >> C_U + B_U
Parameter('k_c3', 1.33e-3) # uM/s
Parameter('K_3', 5e-3) # uM
Observable("C_p", C(state="P"))
Expression('k_Cp_Cu', k_c3*B_u / (K_3 + C_p))
Rule("dephospho_C", C(state="P") >> C(state="U"), k_Cp_Cu)
# Phosphorylation of C:
# Michaelis-Menten Eq. for process:
#     C_U >> C_U
Parameter('k_d3', 8.33e-3) # uM/s
Parameter('K_d3', 5e-3) # uM
Expression('k_Cu_Cp', k_d3 / (K_d3 + C_u))
Rule("phospho_C", C(state='U') >> C(state='P'), k_Cp_Cu)