verbose=True

# migration methods should be "matrix" or "functional"
proton_energy_migration_method="matrix"
electron_energy_migration_method="matrix"
gamma_energy_migration_method="matrix"

# smooth the resulting Aeff and background_rate curves before
# calculating the sensitivity:
enable_smoothing=False
smooth_parameter=1.0

# cut off data below this fraction of the effictive area:
effective_area_fraction_min = 0.1 
