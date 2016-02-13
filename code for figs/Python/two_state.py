from pysb import *
import numpy as np

Model()

Monomer('Cell', ['d'], {'d' : ['0', '1']})

Parameter('kdiv0', 0.065*np.log(2))     # /hr
Parameter('kdth0', 0.005*np.log(2))     # /hr
Parameter('kon',   1000.)                # /uM-hr
Parameter('koff',  10.)                  # /hr
Parameter('kdiv1', 0.01*np.log(2))      # /hr
Parameter('kdth1', 0.04*np.log(2))      # /hr

Parameter('Drug0', 0) # uM

Initial(Cell(d='0'), Parameter('Cell0_init', 1000)) # number
Initial(Cell(d='1'), Parameter('Cell1_init', 0)) # number

Observable('Cell_total', Cell())
Observable('Cell_free',  Cell(d='0'))
Observable('Cell_drug',  Cell(d='1'))

def full_model():
    Rule('Division', Cell(d='0') >> Cell(d='0') + Cell(d='0'), kdiv0)
    Rule('Death', Cell(d='0') >> None, kdth0)

    Rule('DrugActivation', Cell(d='0') >> Cell(d='1'), Expression('kon_Drug', kon*Drug0))
    Rule('DrugInactivation', Cell(d='1') >> Cell(d='0'), koff)
    
    Rule('DivisionInDrug', Cell(d='1') >> Cell(d='1') + Cell(d='1'), kdiv1)
    Rule('DeathInDrug', Cell(d='1') >> None, kdth1)

def reduced_model():
    Expression('kdiv_eff', (koff*kdiv0 + kon*Drug0*kdiv1) / (koff + kon*Drug0))
    Expression('kdth_eff', (koff*kdth0 + kon*Drug0*kdth1) / (koff + kon*Drug0))
    Rule('Division', Cell(d='0') >> Cell(d='0') + Cell(d='0'), kdiv_eff)
    Rule('Death', Cell(d='0') >> None, kdth_eff)

full_model()
# reduced_model()