TAU_D = 1.0
TAU_A = 60.0

C_SHEAR = 0.258
C_META = 3.0e6


AR_PARAMS = {

    "Large Arteriole": {

        "D_control": 65.2e-6,
        "A_control": 0.5,

        "D0": 156.49e-6,

        "C_pass": 1042.99e-3,
        "C_pass_prime": 8.293,

        "C_act": 1596.3e-3,
        "C_act_prime": 0.6804,
        "C_act_double_prime": 0.2905,

        "C_myo": 10.1,

        # Carlson et al. (2008), Table 1:
        # Eq. 6 constant for the non-metabolic (myogenic/shear) model.
        "C_tone_myo": -2.22,

        # Eq. 7 constant when the metabolic term is included.
        "C_tone_meta": 10.11,
    },

    "Small Arteriole": {

        "D_control": 14.8e-6,
        "A_control": 0.5,

        "D0": 38.99e-6,

        "C_pass": 259.90e-3,
        "C_pass_prime": 11.467,

        "C_act": 274.193e-3,
        "C_act_prime": 0.750,
        "C_act_double_prime": 0.384,

        "C_myo": 35.9,

        # Carlson et al. (2008), Table 1:
        # Eq. 6 constant for the non-metabolic (myogenic/shear) model.
        "C_tone_myo": -0.53,

        # Eq. 7 constant when the metabolic term is included.
        "C_tone_meta": 10.66,
    },
}
