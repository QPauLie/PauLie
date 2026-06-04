"""Shared pytest configuration."""
import matplotlib

# Use a non-interactive backend so the plotting and animation tests run headless.
matplotlib.use("Agg")
