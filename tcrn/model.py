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


Model()

# state: U = unphosphorylated, P = phosphorylated
Mononer('A', ['state'], {'state':['U','P']})
Mononer('B', ['state'], {'state':['U','P']})
Mononer('C', ['state'], {'state':['U','P']})

Parameter('Atot', 10.)
Parameter('Btot', 10.)
Parameter('Ctot', 10.)

Parameter('A_U_0', 1.)
Parameter('A_U_0', 1.)
Parameter('A_U_0', 1.)
