import math
from zadanie_1 import calculate_a

def test_calculate_a():
    delta_r = math.radians(30)
    fi_r = math.radians(45)
    t_r = math.radians(60)
    expected_az = math.radians(150)
    assert math.isclose(calculate_a(delta_r, fi_r, t_r), expected_az, rel_tol=1e-9)