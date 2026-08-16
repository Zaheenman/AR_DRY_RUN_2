MMHG_TO_PA = 133.33
UM_TO_M = 1e-6
CM_TO_M = 1e-2
CP_TO_PA_S = 1e-3
M3_S_TO_UL_MIN = 60e9


def mmhg_to_pa(x):
    return x * MMHG_TO_PA


def pa_to_mmhg(x):
    return x / MMHG_TO_PA


def um_to_m(x):
    return x * UM_TO_M


def cm_to_m(x):
    return x * CM_TO_M


def cp_to_pa_s(x):
    return x * CP_TO_PA_S


def m3_s_to_ul_min(x):
    return x * M3_S_TO_UL_MIN


def ul_min_to_m3_s(x):
    return x / M3_S_TO_UL_MIN

def m_s_to_cm_s(x):
    return x * 100.0
    
def resistance_table_to_si(x):
    return x * 1e13
    

def oxygen_demand_to_si(x):
    return x / 100 / 60
    
