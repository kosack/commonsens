# migration methods should be "matrix" or "functional"
proton_energy_migration_method="matrix"
electron_energy_migration_method="matrix"
gamma_energy_migration_method="matrix"

# smooth the resulting Aeff and background_rate curves before
# calculating the sensitivity:
enable_smoothing=False
smooth_window_size=5 # must be odd
