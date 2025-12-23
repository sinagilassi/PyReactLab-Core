# chem_balance.py
from __future__ import annotations

import math
import re
from dataclasses import dataclass
from fractions import Fraction
from typing import Dict, List, Tuple, Iterable, Optional

import numpy as np


# ---------------------------
# Integer helpers
# ---------------------------

def _gcd(a: int, b: int) -> int:
    return math.gcd(abs(a), abs(b))


def _lcm(a: int, b: int) -> int:
    a, b = abs(a), abs(b)
    if a == 0 or b == 0:
        return 0
    return a // _gcd(a, b) * b


def _lcmm(nums: Iterable[int]) -> int:
    out = 1
    for n in nums:
        out = _lcm(out, n)
    return out


# ---------------------------
# Parsing
# ---------------------------

_CHARGE_CURLY_RE = re.compile(r"^(.*)\{([0-9]*)([+-])\}$")  # Fe{3+}, SO4{2-}
_CHARGE_CARET_RE = re.compile(r"^(.*)\^([0-9]*)([+-])$")   # SO4^2-
_CHARGE_TRAIL_RE = re.compile(r"^(.*?)([0-9]*)([+-])$")    # Fe2+, NH4+


def _normalize_token(tok: str) -> str:
    tok = tok.strip()
    tok = tok.replace("·", "*").replace(".", "*")  # hydrate dot to '*'
    if tok in ("e", "e-", "e⁻"):
        return "e-"
    return tok


def _parse_charge(token: str) -> Tuple[str, int]:
    """
    Return (formula_without_charge, charge).
    Supports:
      - Fe{3+}, Cr2O7{2-}
      - SO4^2-
      - Fe2+, NH4+
    """
    t = token.strip()

    m = _CHARGE_CURLY_RE.match(t)
    if m:
        base, num, sign = m.group(1), m.group(2), m.group(3)
        n = int(num) if num else 1
        return base, (n if sign == "+" else -n)

    m = _CHARGE_CARET_RE.match(t)
    if m:
        base, num, sign = m.group(1), m.group(2), m.group(3)
        n = int(num) if num else 1
        return base, (n if sign == "+" else -n)

    # trailing charge (ambiguous with separator '+', so prefer brace style in your UI)
    if t.endswith(("+", "-")):
        m = _CHARGE_TRAIL_RE.match(t)
        if m:
            base, num, sign = m.group(1), m.group(2), m.group(3)
            if base:
                n = int(num) if num else 1
                return base, (n if sign == "+" else -n)

    return t, 0


def _smart_split_plus(side: str) -> List[str]:
    """
    Split by '+' only when NOT inside {...} charge blocks.
    Fixes: H{+}, Cr{3+}, etc.
    """
    out: List[str] = []
    buf: List[str] = []
    brace_depth = 0

    for ch in side:
        if ch == "{":
            brace_depth += 1
            buf.append(ch)
        elif ch == "}":
            brace_depth = max(0, brace_depth - 1)
            buf.append(ch)
        elif ch == "+" and brace_depth == 0:
            token = "".join(buf).strip()
            if token:
                out.append(token)
            buf = []
        else:
            buf.append(ch)

    token = "".join(buf).strip()
    if token:
        out.append(token)
    return out


def _split_equation(eq: str) -> Tuple[List[str], List[str]]:
    s = eq.strip()
    if "->" in s:
        left, right = s.split("->", 1)
    elif "=" in s:
        left, right = s.split("=", 1)
    else:
        raise ValueError("Equation must contain '->' or '='.")

    left_parts = _smart_split_plus(left)
    right_parts = _smart_split_plus(right)

    if not left_parts or not right_parts:
        raise ValueError("Invalid equation: empty side.")
    return left_parts, right_parts


def _tokenize_formula(formula: str) -> List[str]:
    f = formula.replace("[", "(").replace("]", ")")
    return re.findall(r"[A-Z][a-z]?|\d+|\(|\)", f)


def _parse_atoms_single(formula: str) -> Dict[str, int]:
    toks = _tokenize_formula(formula)
    stack: List[Dict[str, int]] = [dict()]
    i = 0

    while i < len(toks):
        tok = toks[i]

        if tok == "(":
            stack.append({})
            i += 1
            continue

        if tok == ")":
            if len(stack) == 1:
                raise ValueError(f"Unmatched ')' in formula: {formula}")
            group = stack.pop()
            i += 1
            mult = 1
            if i < len(toks) and toks[i].isdigit():
                mult = int(toks[i])
                i += 1
            top = stack[-1]
            for el, cnt in group.items():
                top[el] = top.get(el, 0) + cnt * mult
            continue

        if re.match(r"^[A-Z][a-z]?$", tok):
            el = tok
            i += 1
            mult = 1
            if i < len(toks) and toks[i].isdigit():
                mult = int(toks[i])
                i += 1
            top = stack[-1]
            top[el] = top.get(el, 0) + mult
            continue

        if tok.isdigit():
            raise ValueError(
                f"Unexpected number '{tok}' in formula: {formula}")

        raise ValueError(f"Unexpected token '{tok}' in formula: {formula}")

    if len(stack) != 1:
        raise ValueError(f"Unmatched '(' in formula: {formula}")

    return stack[0]


@dataclass(frozen=True)
class Species:
    raw: str
    formula: str
    charge: int
    atoms: Dict[str, int]

    def is_electron(self) -> bool:
        return self.raw == "e-" and self.charge == -1 and self.atoms == {}


def _is_electron_token(tok: str) -> bool:
    """
    Detect electron tokens:
    e, e-, e{−}, e{1-}, e^- , etc.
    """
    t = tok.strip().lower()

    # Normalize unicode minus
    t = t.replace("−", "-")

    # Direct forms
    if t in ("e", "e-", "e^-"):
        return True

    # Curly-brace charge forms: e{−}, e{1-}
    m = re.fullmatch(r"e\{([0-9]*)(-)\}", t)
    if m:
        return True

    return False


def parse_species(token: str) -> Species:
    tok = _normalize_token(token)

    if _is_electron_token(tok):
        return Species(
            raw="e-",
            formula="e",
            charge=-1,
            atoms={}
        )

    base, charge = _parse_charge(tok)
    base = base.strip()

    # hydrates: CuSO4*5H2O
    parts = [p.strip() for p in base.split("*") if p.strip()]
    atoms_total: Dict[str, int] = {}

    for part in parts:
        m = re.match(r"^(\d+)(.*)$", part)
        mult = 1
        frag = part
        if m:
            mult = int(m.group(1))
            frag = m.group(2).strip()
        frag_atoms = _parse_atoms_single(frag)
        for el, cnt in frag_atoms.items():
            atoms_total[el] = atoms_total.get(el, 0) + cnt * mult

    return Species(raw=tok, formula=base, charge=charge, atoms=atoms_total)


# ---------------------------
# Exact rational nullspace (NO floating/SciPy)
# ---------------------------

def _rref_fraction(mat: List[List[Fraction]]) -> Tuple[List[List[Fraction]], List[int]]:
    """Return RREF and pivot column indices."""
    A = [row[:] for row in mat]
    m = len(A)
    n = len(A[0]) if m else 0
    pivots: List[int] = []
    r = 0

    for c in range(n):
        pivot_row = None
        for rr in range(r, m):
            if A[rr][c] != 0:
                pivot_row = rr
                break
        if pivot_row is None:
            continue

        A[r], A[pivot_row] = A[pivot_row], A[r]
        piv = A[r][c]
        A[r] = [x / piv for x in A[r]]

        for rr in range(m):
            if rr == r:
                continue
            factor = A[rr][c]
            if factor != 0:
                A[rr] = [A[rr][j] - factor * A[r][j] for j in range(n)]

        pivots.append(c)
        r += 1
        if r == m:
            break

    return A, pivots


def _nullspace_rational(A_int: np.ndarray) -> List[List[Fraction]]:
    """
    Exact rational nullspace basis for A x = 0.
    Returns list of basis vectors (each is list[Fraction], length n).
    """
    m, n = A_int.shape
    mat = [[Fraction(int(A_int[i, j]), 1) for j in range(n)] for i in range(m)]
    rref, pivots = _rref_fraction(mat)
    pivot_set = set(pivots)
    free_cols = [c for c in range(n) if c not in pivot_set]
    if not free_cols:
        return []

    pivot_row_for_col = {pivots[i]: i for i in range(len(pivots))}
    basis: List[List[Fraction]] = []

    for free in free_cols:
        x = [Fraction(0, 1) for _ in range(n)]
        x[free] = Fraction(1, 1)
        for pc in pivots:
            row = pivot_row_for_col[pc]
            x[pc] = -rref[row][free]
        basis.append(x)

    return basis


def _fraction_vec_to_small_int(v: List[Fraction]) -> List[int]:
    den_lcm = _lcmm([f.denominator for f in v])
    ints = [int(f * den_lcm) for f in v]

    # make sign consistent
    for a in ints:
        if a != 0:
            if a < 0:
                ints = [-k for k in ints]
            break

    g = 0
    for a in ints:
        g = _gcd(g, a)
    if g > 1:
        ints = [a // g for a in ints]
    return ints


# ---------------------------
# Matrix + formatting
# ---------------------------

def _build_balance_matrix(
    reactants: List[Species],
    products: List[Species],
    include_charge: bool
) -> np.ndarray:
    all_species = reactants + products
    elements = sorted({el for sp in all_species for el in sp.atoms.keys()})

    rows: List[List[int]] = []
    for el in elements:
        row = [sp.atoms.get(el, 0) for sp in reactants] + \
            [-sp.atoms.get(el, 0) for sp in products]
        rows.append(row)

    if include_charge:
        crow = [sp.charge for sp in reactants] + \
            [-sp.charge for sp in products]
        rows.append(crow)

    return np.array(rows, dtype=int)


def _format_side(items: List[Tuple[int, Species]]) -> str:
    parts = []
    for c, sp in items:
        if c == 1:
            parts.append(sp.raw)
        else:
            parts.append(f"{c} {sp.raw}")
    return " + ".join(parts) if parts else "0"


def _format_equation(lhs: List[Tuple[int, Species]], rhs: List[Tuple[int, Species]]) -> str:
    return f"{_format_side(lhs)} -> {_format_side(rhs)}"


def _cancel_both_sides(
    lhs: List[Tuple[int, Species]],
    rhs: List[Tuple[int, Species]],
) -> Tuple[List[Tuple[int, Species]], List[Tuple[int, Species]]]:
    lm = {sp.raw: [c, sp] for c, sp in lhs}
    rm = {sp.raw: [c, sp] for c, sp in rhs}
    for k in list(lm.keys()):
        if k in rm:
            m = min(lm[k][0], rm[k][0])
            lm[k][0] -= m
            rm[k][0] -= m
    lhs2 = [(c, sp) for c, sp in (v for v in lm.values()) if c != 0]
    rhs2 = [(c, sp) for c, sp in (v for v in rm.values()) if c != 0]
    return lhs2, rhs2


# ---------------------------
# Method 1: Algebraic (core)
# ---------------------------

def balance_algebraic(equation: str, include_charge: str = "auto") -> str:
    left_tokens, right_tokens = _split_equation(equation)
    reactants = [parse_species(t) for t in left_tokens]
    products = [parse_species(t) for t in right_tokens]

    if include_charge == "auto":
        use_charge = any(sp.charge != 0 for sp in reactants +
                         products) or any(sp.is_electron() for sp in reactants + products)
    elif include_charge == "yes":
        use_charge = True
    elif include_charge == "no":
        use_charge = False
    else:
        raise ValueError("include_charge must be 'auto', 'yes', or 'no'.")

    A = _build_balance_matrix(reactants, products, include_charge=use_charge)
    basis = _nullspace_rational(A)
    if not basis:
        raise ValueError(
            "No non-trivial solution found. Equation may be impossible to balance.")

    # choose simplest integer basis vector
    candidates = [_fraction_vec_to_small_int(v) for v in basis]
    best = min(candidates, key=lambda ints: sum(abs(x) for x in ints))

    nL = len(reactants)
    lhs = [(best[i], reactants[i]) for i in range(nL) if best[i] != 0]
    rhs = [(best[nL + j], products[j])
           for j in range(len(products)) if best[nL + j] != 0]

    # overall positive
    if any(c < 0 for c, _ in lhs + rhs):
        lhs = [(-c, sp) for c, sp in lhs]
        rhs = [(-c, sp) for c, sp in rhs]

    # gcd normalize
    allc = [c for c, _ in lhs + rhs]
    g = 0
    for c in allc:
        g = _gcd(g, c)
    if g > 1:
        lhs = [(c // g, sp) for c, sp in lhs]
        rhs = [(c // g, sp) for c, sp in rhs]

    return _format_equation(lhs, rhs)


# ---------------------------
# Method 2: Oxidation-number (kept as wrapper; robust solver is algebraic)
# ---------------------------

_COMMON_FIXED_OX = {
    "F": -1,
    "Li": +1, "Na": +1, "K": +1, "Rb": +1, "Cs": +1,
    "Be": +2, "Mg": +2, "Ca": +2, "Sr": +2, "Ba": +2,
}


def _infer_ox_state_single_unknown(sp: Species, target_el: str) -> Optional[Fraction]:
    """Simple inference if only target is unknown and others follow common rules."""
    if target_el not in sp.atoms:
        return None
    total_charge = Fraction(sp.charge, 1)
    known_sum = Fraction(0, 1)
    unknown_count = 0

    for el, cnt in sp.atoms.items():
        cntF = Fraction(cnt, 1)
        if el == target_el:
            unknown_count += cnt
            continue

        if el in _COMMON_FIXED_OX:
            known_sum += Fraction(_COMMON_FIXED_OX[el], 1) * cntF
        elif el == "O":
            known_sum += Fraction(-2, 1) * cntF
        elif el == "H":
            known_sum += Fraction(+1, 1) * cntF
        else:
            return None

    if unknown_count == 0:
        return None
    return (total_charge - known_sum) / Fraction(unknown_count, 1)


def balance_oxidation_number(equation: str) -> str:
    # You can extend this to output steps; final robust answer comes from algebraic solve.
    return balance_algebraic(equation, include_charge="auto")


# ---------------------------
# Method 3: Half-reaction (computational redox mode)
# ---------------------------

def _detect_medium(tokens: List[str]) -> str:
    toks = set(tokens)
    if "H{+}" in toks:
        return "acid"
    if "OH{-}" in toks:
        return "base"
    return "unknown"


def _looks_ionic_or_redox(species: List[Species]) -> bool:
    return any(sp.charge != 0 for sp in species) or any(sp.is_electron() for sp in species)


def balance_half_reaction(equation: str, medium: str = "auto") -> str:
    """
    Redox/ionic-friendly mode:
      - decides acid/base for medium='auto' without allowing BOTH H+ and OH-
      - adds only needed helpers on LHS, then solves with charge conservation
      - moves negative coefficients across, cancels both sides, normalizes
      - falls back to algebraic for non-ionic reactions
    """
    left_tokens, right_tokens = _split_equation(equation)
    base_left = [_normalize_token(t) for t in left_tokens]
    base_right = [_normalize_token(t) for t in right_tokens]

    all_base = base_left + base_right
    parsed_base = [parse_species(t) for t in all_base]

    if medium not in ("auto", "acid", "base"):
        raise ValueError("medium must be 'auto', 'acid', or 'base'.")

    if medium == "auto":
        md = _detect_medium(all_base)
        if md == "unknown":
            # If not ionic/redox, do NOT invent H+/OH-/e-
            if not _looks_ionic_or_redox(parsed_base):
                return balance_algebraic(equation, include_charge="auto")
            md = "acid"  # default for ionic redox if no hint
        medium = md

    if medium == "acid":
        helpers = ["H2O", "H{+}", "e-"]
    else:
        helpers = ["H2O", "OH{-}", "e-"]

    present = set(all_base)
    add_left = [h for h in helpers if h not in present]

    reactants = [parse_species(t) for t in (base_left + add_left)]
    products = [parse_species(t) for t in base_right]

    A = _build_balance_matrix(reactants, products, include_charge=True)
    basis = _nullspace_rational(A)
    if not basis:
        return balance_algebraic(equation, include_charge="auto")

    # choose simplest, but also prefer smaller helper usage
    helper_idxs = list(range(len(base_left), len(base_left) + len(add_left)))

    best_ints = None
    best_score = None
    for v in basis:
        ints = _fraction_vec_to_small_int(v)
        score = sum(abs(x) for x in ints) + 3 * \
            sum(abs(ints[i]) for i in helper_idxs)
        if best_score is None or score < best_score:
            best_score = score
            best_ints = ints

    if best_ints is None:
        return balance_algebraic(equation, include_charge="auto")

    coeffs = best_ints
    nL = len(reactants)
    lhs_pairs = [(coeffs[i], reactants[i]) for i in range(nL)]
    rhs_pairs = [(coeffs[nL + j], products[j]) for j in range(len(products))]

    # Move negative coefficients across arrow
    lhs: List[Tuple[int, Species]] = []
    rhs: List[Tuple[int, Species]] = []

    for c, sp in lhs_pairs:
        if c > 0:
            lhs.append((c, sp))
        elif c < 0:
            rhs.append((-c, sp))

    for c, sp in rhs_pairs:
        if c > 0:
            rhs.append((c, sp))
        elif c < 0:
            lhs.append((-c, sp))

    # Cancel same species on both sides
    lhs, rhs = _cancel_both_sides(lhs, rhs)

    # GCD normalize
    allc = [c for c, _ in lhs + rhs]
    g = 0
    for c in allc:
        g = _gcd(g, c)
    if g > 1:
        lhs = [(c // g, sp) for c, sp in lhs]
        rhs = [(c // g, sp) for c, sp in rhs]

    return _format_equation(lhs, rhs)


# ---------------------------
# Unified API
# ---------------------------

def balance(equation: str, method: str = "algebraic", medium: str = "auto") -> str:
    m = method.lower().strip()
    if m in ("algebraic", "matrix"):
        return balance_algebraic(equation, include_charge="auto")
    if m in ("half", "half-reaction", "ion", "ion-electron"):
        return balance_half_reaction(equation, medium=medium)
    if m in ("oxidation", "oxidation-number", "ox"):
        return balance_oxidation_number(equation)
    raise ValueError("method must be: algebraic | half | oxidation")
