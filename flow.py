import math

def resistance(length, viscosity, diameter, number):
    return 128 * viscosity * length / (
        math.pi * number * diameter**4
    )

def total_resistance(resistances):
    return sum(resistances)

def total_flow(p_in, p_out, resistance_total):
    return (p_in - p_out) / resistance_total

def vessel_flow(flow_total, number):
    return flow_total / number

def pressure_drop(flow_total, resistance):
    return flow_total * resistance

def mean_pressure(p_in, p_out):
    return (p_in + p_out) / 2

def velocity(flow, diameter):
    return flow / (math.pi * diameter**2 / 4)

def wall_shear_stress(viscosity, flow, diameter):
    return 32 * viscosity * flow / (
        math.pi * diameter**3
    )
