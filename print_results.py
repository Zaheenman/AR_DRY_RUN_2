from units import m3_s_to_ul_min, pa_to_mmhg


def print_results(results, summary):

    print()
    print(results.to_string(index=False))

    print()
    print(f"Total resistance = {summary['R_total']:.6e} Pa.s/m3")
    print(f"Total flow = {m3_s_to_ul_min(summary['Q_total']):.6f} uL/min")
    print(f"Inlet pressure = {pa_to_mmhg(summary['p_in']):.3f} mmHg")
    print(f"Outlet pressure = {pa_to_mmhg(summary['p_out']):.3f} mmHg")
