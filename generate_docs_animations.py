import os
import math
from paulie import get_pauli_string as p
from paulie.classifier.recording_canonicalizer import (
    RecordingCanonicalizer,
)
from paulie.helpers._renderer import (
    render_gradient_animation,
)

os.makedirs(
    "docs/source/_static/images",
    exist_ok=True
)

def make_stack(strings):
    return [p(s) for s in strings]

# --- TYP A ---
rc_a = RecordingCanonicalizer()
gens_a = [
    "IYZI", "IIXX", "IIYZ",
    "IXXI", "XXII", "YZII"
]
rc_a._set_vertex_stack(make_stack(gens_a))
if hasattr(rc_a, "classify"):
    rc_a.classify()

print("Generuje petle kamery dla Typu A...")
# Pętla obrotu: prosta do góry i kółko 360 stopni
for i in range(24):
    kat_gora = 30 + (i * 1.5)
    kat_obrot = (i / 24) * 360
    
    # Podmieniamy pozycje kamery w rejestratorze
    if hasattr(rc_a.recorder, "camera"):
        rc_a.recorder.camera.elevation = kat_gora
        rc_a.recorder.camera.azimuth = kat_obrot

    sciezka_a = (
        f"docs/source/_static/images/"
        f"animation_type_a_{i:02d}.png"
    )
    render_gradient_animation(
        rc_a.recorder,
        sciezka_a,
    )

print("Typ A wygenerowany kompletnie!")

# --- TYP B ---
rc_b = RecordingCanonicalizer()
gens_b = ["XY", "XZ"]
rc_b._set_vertex_stack(make_stack(gens_b))
if hasattr(rc_b, "classify"):
    rc_b.classify()

print("Generuje petle kamery dla Typu B...")
for i in range(24):
    kat_gora = 30 + (i * 1.5)
    kat_obrot = (i / 24) * 360
    
    if hasattr(rc_b.recorder, "camera"):
        rc_b.recorder.camera.elevation = kat_gora
        rc_b.recorder.camera.azimuth = kat_obrot

    sciezka_b = (
        f"docs/source/_static/images/"
        f"animation_type_b_{i:02d}.png"
    )
    render_gradient_animation(
        rc_b.recorder,
        sciezka_b,
    )

print("====================")
print("KURWA, SUKCES ANIMACJI!")
