#!/usr/bin/env python
import binary_tools as bt
periods = ["02.00","02.50"]
for y_label in ["mass","age"]:
   for z_label in ["total_mass_lost","omega_crit","lg_mstar_dot_1","lg_mstar_dot_2","n14"]:
       for mass_accretor in ["15.00","15.50","16.00","16.50","17.00","17.50","18.00","18.50","19.00","19.50"]:
           bt.plot_fixed_masses("20",mass_accretor,periods,z_label=z_label,y_label=y_label)
   
