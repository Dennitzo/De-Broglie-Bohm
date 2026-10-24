"""Historical regularisation sweep is invalid for the corrected model."""
if __name__ == '__main__':
    raise SystemExit('Density-floor sweeps are disabled: they change the continuity equation. Use slits_density.py with eps=0 and refine dx, dy and domain independently.')
