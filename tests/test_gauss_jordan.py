import importlib.util
from pathlib import Path
from unittest.mock import patch

import pytest
import sympy as sp


SCRIPT = (Path(__file__).resolve().parents[1] / "01_Numerical_Methods" /
          "03_Linear_Systems" / "Gauss-Jordan" / "Gauss-Jordan" / "Gauss_Jordan.py")
spec = importlib.util.spec_from_file_location("gauss_jordan", SCRIPT)
solver = importlib.util.module_from_spec(spec)
spec.loader.exec_module(solver)


@pytest.mark.parametrize("matrix", [
    [[0, 2, 4], [1, 1, 3]],
    [[1, 1, 1, 6], [2, 2, 3, 13], [1, 2, 3, 14]],
    [[4, 8]],
    solver.A,
])
def test_both_methods_match_exact_solution(matrix):
    original = [row[:] for row in matrix]
    expected = list(sp.Matrix([row[:-1] for row in matrix]).inv() * sp.Matrix([row[-1] for row in matrix]))
    with patch("builtins.input", return_value=""):
        reduced = solver.gauss_jordan_method([row[:] for row in matrix])
        gaussian = solver.gauss_method([row[:] for row in matrix])
    assert [row[-1] for row in reduced] == pytest.approx([float(x) for x in expected])
    variables = sp.symbols(f"x1:{len(matrix) + 1}")
    assert [gaussian[var] for var in variables] == pytest.approx([float(x) for x in expected])
    for i, row in enumerate(reduced):
        assert row[:-1] == pytest.approx([int(i == j) for j in range(len(matrix))])
    assert matrix == original


@pytest.mark.parametrize("matrix", [[], [[1, 2, 3]], [[1, 2, 3], [2, 4, 6]], [[1, 2, 3], [2, 4, 7]]])
def test_invalid_or_nonunique_system_is_rejected(matrix):
    with patch("builtins.input", return_value=""), pytest.raises(ValueError):
        solver.gauss_jordan_method(matrix)


@pytest.mark.parametrize("choice", ["1", "2"])
def test_menu_options_preserve_original_system(choice):
    original = [row[:] for row in solver.A]
    with patch("builtins.input", side_effect=lambda *args: choice), patch.object(solver, "clear_screen"):
        solver.main()
    assert solver.A == original
