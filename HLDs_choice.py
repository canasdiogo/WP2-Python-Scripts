import numpy as np

#Function to calculate airfoil cl_max
def Clmax(clmax_airfoil):
    return ratio * clmax_airfoil

#Function to calculate the delta Cl obtained by the combination of LE and TE HLDs
def Delta_CL_Max_obtained(Swf_S_te, Swf_S_le, deltaCl_te, deltaCl_le, hinge_sweep_angle_TE, hinge_sweep_angle_LE):
    delta_CL_Max_TE = 0.9 * deltaCl_te * Swf_S_te*np.cos(hinge_sweep_angle_TE * (np.pi / 180))
    delta_CL_Max_LE = 0.9 * deltaCl_le * Swf_S_le*np.cos(hinge_sweep_angle_LE * (np.pi / 180))
    delta_CL_Max = delta_CL_Max_LE + delta_CL_Max_TE
    return delta_CL_Max, delta_CL_Max_TE, delta_CL_Max_LE

#Constants
leading_edge_sweep = 31 #[deg] #EXPLAIN THE CALCULATION OF THIS NUMBER IN OVERLEAF
hinge_sweep_TE = 24.6 #[deg]
hinge_sweep_LE = 29.5 #[deg]
take_off_mass = 200535 #[kg] #EXPLAIN THE CALCULATION OF THIS VALUE IN OVERLEAF
g = 9.80665 #m/s^2
rho = 1.225 #[kg/m^3]
S = 333.43 #m^2
lift_at_take_off = take_off_mass * g
half_span = 28 #[m]

#Parameters to calculate CL_max of clean wing
ratio = 0.8 #Cl/CL ratio taken from graph in ADSEE slides (don't forget to explain the process to get this ratio)
cl_airfoil_take_off = 1.732 #Cl of airfoil at 14 aoa (aoa assumed for take-off)
cl_airfoil_landing = 1.588 # Cl of airfoil at 12 aoa (aoa assumed for landing)

#Parameters to calculate total increment in Cl by the HLDs of design option 1 (TE: fowler flaps + LE: kruger flap)
deltaCl_airfoil_TE_design_option1_landing = 1.6 #value taken from type of TE HLDs
deltaCl_airfoil_TE_design_option1_take_off = deltaCl_airfoil_TE_design_option1_landing * 0.60 #only 60% of landing values, as TE HLD's are not fully deployed
deltaCl_airfoil_LE_design_option1_landing = 0.3 #value taken from type of LE HLDs
deltaCl_airfoil_LE_design_option1_take_off = deltaCl_airfoil_LE_design_option1_landing * 0.60 #only 60% of landing values, as LE HLD's are not fully deployed
Swf_S_LE_design_option1 = 0.75
Swf_S_TE_design_option1 = 0.5

#Parameters to calculate total increment in Cl by the HLDs of design option 2 (TE: single-slotted flaps + LE: kruger flap)
deltaCl_airfoil_TE_design_option2_landing = 1.3 #value taken from type of TE HLDs
deltaCl_airfoil_TE_design_option2_take_off = deltaCl_airfoil_TE_design_option2_landing * 0.60 #only 60% of landing values, as TE HLD's are not fully deployed
deltaCl_airfoil_LE_design_option2_landing = 0.3 #value taken from type of LE HLDs
deltaCl_airfoil_LE_design_option2_take_off = 0.3 * 0.60 #only 60% of landing values, as LE HLD's are not fully deployed
Swf_S_LE_design_option2 = 0.8
Swf_S_TE_design_option2 = 0.57

#Parameters to calculate size of TE HLDs
c = 6.79 #[m] (MAC)
ratio_c_option1 = deltaCl_airfoil_TE_design_option1_landing / 1.3
c_prime_option1 = ratio_c_option1 * c #[m]
delta_c_option1 = c_prime_option1 - c
cf_TE = 0.38 * c #[m]
ratio_deltac_cf_option1 = delta_c_option1 / cf_TE #value which is needed to insert into the plot in the ADSEE slides to then calculate deflection angle
ratio_deltac_cf_option2 = 0.31
delta_c_option2 = ratio_deltac_cf_option2 * cf_TE
c_prime_option2 = delta_c_option2 + c
deflection_angle_option1 = 35 #[deg] value obtained visually from plot in ADSEE slides
deflection_angle_option2 = 40 #[deg]

#Calculate area of TE HLDs
ratio_S_prime_S_option1 = 1 + Swf_S_TE_design_option1 * (c_prime_option1 / c - 1)
print(f'H{ratio_S_prime_S_option1}')
ratio_S_prime_S_option2 = 1 + Swf_S_TE_design_option2 * (c_prime_option2 / c - 1)
area_TE_HLD_option1 = S * (ratio_S_prime_S_option1 - 1) / 2 #Area of each Fowler flap
area_TE_HLD_option2 = S * (ratio_S_prime_S_option2 - 1) / 2 #Area of each Single-Sloted flap

#Calculate CL of the clean wing at take-off
CLmax_clean_wing_take_off = Clmax(cl_airfoil_take_off)
print(CLmax_clean_wing_take_off)
#Calculate CL of the clean wing at approach
CLmax_clean_wing_approach = Clmax(cl_airfoil_landing)

#Calculate take-off velocity
v_take_off = 77 #[m/s]
v_take_off_eff = v_take_off * np.cos((leading_edge_sweep*np.pi)/180) #[m/s]

#Approach velocity
v_approach = 75 #[m/s]
v_approach_eff = v_approach * np.cos((leading_edge_sweep*np.pi)/180)
v_stall = v_approach/1.23
v_stall_eff = v_stall * np.cos((leading_edge_sweep*np.pi)/180)

#Required lift coefficient to take_off
CL_required_to_take_off = 1.7
#Required lift coefficient to land
CL_required_to_land = 2.04

#Calculate Delta_CL that the HLDs will have to provide at take-off
delta_CL_take_off_required = CL_required_to_take_off - CLmax_clean_wing_take_off
#Calculate Delta_CL that the HLDs will have to provide at approach
delta_CL_approach_required = CL_required_to_land - CLmax_clean_wing_approach

_,_,delta_CL_max_obtained_option1_LE_landing = Delta_CL_Max_obtained(Swf_S_te=Swf_S_TE_design_option1, Swf_S_le=Swf_S_LE_design_option1, deltaCl_te=deltaCl_airfoil_TE_design_option1_landing, deltaCl_le=deltaCl_airfoil_LE_design_option1_landing, hinge_sweep_angle_TE=hinge_sweep_TE, hinge_sweep_angle_LE=hinge_sweep_LE)
_,_,delta_CL_max_obtained_option1_LE_take_off = Delta_CL_Max_obtained(Swf_S_te=Swf_S_TE_design_option1, Swf_S_le=Swf_S_LE_design_option1, deltaCl_te=deltaCl_airfoil_TE_design_option1_landing, deltaCl_le=deltaCl_airfoil_LE_design_option1_take_off, hinge_sweep_angle_TE=hinge_sweep_TE, hinge_sweep_angle_LE=hinge_sweep_LE)
_,delta_CL_max_obtained_option1_TE_landing,_ = Delta_CL_Max_obtained(Swf_S_te=Swf_S_TE_design_option1, Swf_S_le=Swf_S_LE_design_option1, deltaCl_te=deltaCl_airfoil_TE_design_option1_landing, deltaCl_le=deltaCl_airfoil_LE_design_option1_landing, hinge_sweep_angle_TE=hinge_sweep_TE, hinge_sweep_angle_LE=hinge_sweep_LE)
_,delta_CL_max_obtained_option1_TE_take_off,_ = Delta_CL_Max_obtained(Swf_S_te=Swf_S_TE_design_option1, Swf_S_le=Swf_S_LE_design_option1, deltaCl_te=deltaCl_airfoil_TE_design_option1_take_off, deltaCl_le=deltaCl_airfoil_LE_design_option1_take_off, hinge_sweep_angle_TE=hinge_sweep_TE, hinge_sweep_angle_LE=hinge_sweep_LE)
print(delta_CL_max_obtained_option1_LE_landing)
print(delta_CL_max_obtained_option1_TE_landing)

delta_CL_max_obtained_option1_landing, _, _ = Delta_CL_Max_obtained(Swf_S_te=Swf_S_TE_design_option1, Swf_S_le=Swf_S_LE_design_option1, deltaCl_te=deltaCl_airfoil_TE_design_option1_landing, deltaCl_le=deltaCl_airfoil_LE_design_option1_landing, hinge_sweep_angle_TE=hinge_sweep_TE, hinge_sweep_angle_LE=hinge_sweep_LE)
delta_CL_max_obtained_option1_take_off, _, _ = Delta_CL_Max_obtained(Swf_S_te=Swf_S_TE_design_option1, Swf_S_le=Swf_S_LE_design_option1, deltaCl_te=deltaCl_airfoil_TE_design_option1_take_off, deltaCl_le=deltaCl_airfoil_LE_design_option1_take_off, hinge_sweep_angle_TE=hinge_sweep_TE, hinge_sweep_angle_LE=hinge_sweep_LE)

_,_,delta_CL_max_obtained_option2_LE_landing = Delta_CL_Max_obtained(Swf_S_te=Swf_S_TE_design_option2, Swf_S_le=Swf_S_LE_design_option2, deltaCl_te=deltaCl_airfoil_TE_design_option2_landing, deltaCl_le=deltaCl_airfoil_LE_design_option2_landing, hinge_sweep_angle_TE=hinge_sweep_TE, hinge_sweep_angle_LE=hinge_sweep_LE)
_,_,delta_CL_max_obtained_option2_LE_take_off = Delta_CL_Max_obtained(Swf_S_te=Swf_S_TE_design_option2, Swf_S_le=Swf_S_LE_design_option2, deltaCl_te=deltaCl_airfoil_TE_design_option2_landing, deltaCl_le=deltaCl_airfoil_LE_design_option2_take_off, hinge_sweep_angle_TE=hinge_sweep_TE, hinge_sweep_angle_LE=hinge_sweep_LE)
_,delta_CL_max_obtained_option2_TE_landing,_ = Delta_CL_Max_obtained(Swf_S_te=Swf_S_TE_design_option2, Swf_S_le=Swf_S_LE_design_option2, deltaCl_te=deltaCl_airfoil_TE_design_option2_landing, deltaCl_le=deltaCl_airfoil_LE_design_option2_landing, hinge_sweep_angle_TE=hinge_sweep_TE, hinge_sweep_angle_LE=hinge_sweep_LE)
_,delta_CL_max_obtained_option2_TE_take_off,_ = Delta_CL_Max_obtained(Swf_S_te=Swf_S_TE_design_option2, Swf_S_le=Swf_S_LE_design_option2, deltaCl_te=deltaCl_airfoil_TE_design_option2_take_off, deltaCl_le=deltaCl_airfoil_LE_design_option2_take_off, hinge_sweep_angle_TE=hinge_sweep_TE, hinge_sweep_angle_LE=hinge_sweep_LE)
print(delta_CL_max_obtained_option2_LE_landing)
print(delta_CL_max_obtained_option2_TE_landing)

delta_CL_max_obtained_option2_landing, _, _ = Delta_CL_Max_obtained(Swf_S_te=Swf_S_TE_design_option2, Swf_S_le=Swf_S_LE_design_option2, deltaCl_te=deltaCl_airfoil_TE_design_option2_landing, deltaCl_le=deltaCl_airfoil_LE_design_option2_landing, hinge_sweep_angle_TE=hinge_sweep_TE, hinge_sweep_angle_LE=leading_edge_sweep)
delta_CL_max_obtained_option2_take_off, _, _ = Delta_CL_Max_obtained(Swf_S_te=Swf_S_TE_design_option2, Swf_S_le=Swf_S_LE_design_option2, deltaCl_te=deltaCl_airfoil_TE_design_option2_take_off, deltaCl_le=deltaCl_airfoil_LE_design_option2_take_off, hinge_sweep_angle_TE=hinge_sweep_TE, hinge_sweep_angle_LE=leading_edge_sweep)
print(f'HERE {delta_CL_max_obtained_option1_take_off + CLmax_clean_wing_take_off}')
print(f'TAKE-OFF: required delta CL (in relation to clean wing cenario) of {round(delta_CL_take_off_required, 3)}')
print(f'LANDING: required delta CL (in relation to clean wing cenario) of {round(delta_CL_approach_required,3)}')

print(f'The delta CL provided by combination 1 of HLDs on take-off is {round(delta_CL_max_obtained_option1_take_off,3)}')
print(f'The delta CL provided by combination 1 of HLDs on landing is {round(delta_CL_max_obtained_option1_landing,3)}')
print(f'The delta CL provided by combination 2 of HLDs on take-off is {round(delta_CL_max_obtained_option2_take_off,3)}')
print(f'The delta CL provided by combination 2 of HLDs on landing is {round(delta_CL_max_obtained_option2_landing,3)}')


if delta_CL_max_obtained_option1_landing > delta_CL_approach_required:
    print('Combination 1 of HLDs works for landing')
if delta_CL_max_obtained_option1_take_off > delta_CL_take_off_required:
    print('Combination 1 of HLDs works for take_off')
if delta_CL_max_obtained_option1_landing < delta_CL_approach_required:
        print('Combination 1 of HLDs does not work for landing')
if delta_CL_max_obtained_option1_take_off < delta_CL_take_off_required:
        print('Combination 1 of HLDs does not work for take_off')

if delta_CL_max_obtained_option2_landing > delta_CL_approach_required:
    print('Combination 2 of HLDs works for landing')
if delta_CL_max_obtained_option2_take_off > delta_CL_take_off_required:
    print('Combination 2 of HLDs works for take_off')
if delta_CL_max_obtained_option2_landing < delta_CL_approach_required:
        print('Combination 2 of HLDs does not work for landing')
if delta_CL_max_obtained_option2_take_off < delta_CL_take_off_required:
        print('Combination 2 of HLDs does not work for take_off')

print(f'The area of each TE HLD corresponding to option 1 is: {round(area_TE_HLD_option1,2)}')
print(f'The area of each TE HLD corresponding to option 2 is: {round(area_TE_HLD_option2,2)}')