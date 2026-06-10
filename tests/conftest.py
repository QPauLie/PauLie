"""
Shared test configuration.
"""
import matplotlib

# Use a non-interactive backend so the animation tests run headless.
matplotlib.use("Agg")
