# MAUVE Fixed Galaxy Color Scheme

Date: 2026-05-30

Files:

- `MAUVE_galaxy_colors.dat`
- `20260530_MAUVE_galaxy_colors.md`

Copies are saved in both `ICRAR/further` and the parent `ICRAR` folder so the active notebook workspace and the project-level folder share the same convention.

## Purpose

From this point onward, MAUVE plots that encode galaxy identity should use one fixed color per galaxy. The old `tab20` convention is not sufficient because the active sample contains 40 component galaxies, while `tab20` only has 20 categorical colors.

The goal is consistency across all current and future notebooks:

- the same galaxy should always have the same color;
- different galaxies should be visually distinguishable where possible;
- black should remain reserved for the total assembled or pooled galaxy sample;
- `NGC4567_8` should be treated as a combined-system alias and use `NGC4567`'s color.

## Official Palette

Official palette name:

```python
MAUVE_GALAXY_GLASBEY_PACKAGE_40
```

The component-galaxy colors were generated with:

```python
glasbey.create_palette(
    palette_size=40,
    as_hex=True,
    lightness_bounds=(25, 88),
    chroma_bounds=(20, 90),
    optimize_palette=True,
)
```

The environment used when fixing the scheme had:

- `colorcet 3.2.1`
- `glasbey 0.3.0`
- `distinctipy 1.3.4`

The mapping is now persisted in the `.dat` file and should be treated as fixed. For reproducibility, future plotting code should load the saved mapping or copy its fixed dictionary rather than regenerating a new palette unless the whole convention is intentionally revised.

## Ordering Rule

The 40 component galaxies are read from `new_redshifts.txt` and ordered as:

1. IC galaxies first.
2. NGC galaxies second.
3. Increasing numeric ID within each prefix.

This gives one deterministic list and one deterministic color assignment.

## Special Cases

### `NGC4567_8`

`NGC4567_8` is not assigned a new independent color. It is an alias for the combined system and falls back to the `NGC4567` color:

```python
FIXED_GALAXY_COLORS["NGC4567_8"] = FIXED_GALAXY_COLORS["NGC4567"]
```

Current values:

- `NGC4567`: `#04aa75`
- `NGC4568`: `#aa592d`
- `NGC4567_8`: `#04aa75`

The alias redshift in the `.dat` file is the mean of the two component redshifts:

```text
(0.007495 + 0.007446) / 2 = 0.0074705
```

### Total Assembled Sample

Black is reserved for total assembled, pooled, or all-galaxy reference curves and points:

```python
TOTAL_ASSEMBLED_COLOR = "black"
```

No individual component galaxy should use black.

## Data File Columns

The `.dat` file is whitespace-delimited. Comment lines begin with `#`.

Columns:

1. `order`: fixed palette order for component galaxies; alias/reserved rows use text labels.
2. `galaxy`: galaxy identifier.
3. `redshift`: redshift used for the current sample; `NGC4567_8` uses the component mean.
4. `color_hex`: fixed plotting color.
5. `r_float`, `g_float`, `b_float`: RGB channels in the 0-1 Matplotlib convention.
6. `r_255`, `g_255`, `b_255`: RGB channels in 0-255 integer convention.
7. `role`: `component`, `alias_of_NGC4567`, or `reserved_total_sample`.

## Recommended Loader

```python
from pathlib import Path
from matplotlib.colors import ListedColormap

COLOR_TABLE = Path("MAUVE_galaxy_colors.dat")

def load_mauve_galaxy_colors(path=COLOR_TABLE):
    colors = {}
    order = []
    total_color = "black"
    with Path(path).open() as handle:
        for line in handle:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            row_id, galaxy, redshift, color_hex = parts[:4]
            role = parts[-1]
            if role == "component":
                order.append(galaxy)
                colors[galaxy] = color_hex
            elif role == "alias_of_NGC4567":
                colors[galaxy] = color_hex
            elif role == "reserved_total_sample":
                total_color = color_hex
    cmap = ListedColormap([colors[g] for g in order], name="MAUVE_GALAXY_GLASBEY_PACKAGE_40")
    return order, colors, cmap, total_color
```

For a galaxy-colored scatter plot:

```python
order, galaxy_colors, galaxy_cmap, total_color = load_mauve_galaxy_colors()

for galaxy in order:
    ax.scatter(x[galaxy], y[galaxy], color=galaxy_colors[galaxy], label=galaxy)

ax.plot(x_total, y_total, color=total_color, lw=2.0, label="total assembled")
```

## Fixed Mapping

| Order | Galaxy | Redshift | Color |
|---:|---|---:|---|
| 01 | IC3392 | 0.0055660 | `#d21820` |
| 02 | NGC4064 | 0.0030250 | `#6175ff` |
| 03 | NGC4189 | 0.0069970 | `#008a00` |
| 04 | NGC4192 | -0.0004740 | `#ff79fb` |
| 05 | NGC4216 | 0.0004370 | `#aafb00` |
| 06 | NGC4222 | 0.0007670 | `#713971` |
| 07 | NGC4254 | 0.0080260 | `#00dbc2` |
| 08 | NGC4293 | 0.0030030 | `#ff9e04` |
| 09 | NGC4294 | 0.0011840 | `#35595d` |
| 10 | NGC4298 | 0.0037490 | `#6d4904` |
| 11 | NGC4302 | 0.0036630 | `#aa657d` |
| 12 | NGC4321 | 0.0052400 | `#a2aaff` |
| 13 | NGC4330 | 0.0052140 | `#92aa65` |
| 14 | NGC4351 | 0.0077550 | `#a241ff` |
| 15 | NGC4380 | 0.0032120 | `#00ff96` |
| 16 | NGC4383 | 0.0057040 | `#c60079` |
| 17 | NGC4388 | 0.0084190 | `#ffaeca` |
| 18 | NGC4394 | 0.0030750 | `#697996` |
| 19 | NGC4396 | -0.0004170 | `#a67d00` |
| 20 | NGC4402 | 0.0007910 | `#697955` |
| 21 | NGC4405 | 0.0057880 | `#ffe392` |
| 22 | NGC4419 | -0.0008710 | `#ff6182` |
| 23 | NGC4424 | 0.0014580 | `#863135` |
| 24 | NGC4450 | 0.0065250 | `#a6e7ff` |
| 25 | NGC4457 | 0.0029420 | `#1cc600` |
| 26 | NGC4501 | 0.0076190 | `#35518e` |
| 27 | NGC4522 | 0.0077690 | `#ff5d00` |
| 28 | NGC4535 | 0.0065510 | `#a682be` |
| 29 | NGC4548 | 0.0016210 | `#04a6ef` |
| 30 | NGC4567 | 0.0074950 | `#04aa75` |
| 31 | NGC4568 | 0.0074460 | `#aa592d` |
| 32 | NGC4569 | -0.0007840 | `#cac200` |
| 33 | NGC4579 | 0.0050600 | `#08612d` |
| 34 | NGC4580 | 0.0034530 | `#dfd2fb` |
| 35 | NGC4606 | 0.0055000 | `#ca9679` |
| 36 | NGC4607 | 0.0075350 | `#d200df` |
| 37 | NGC4654 | 0.0034560 | `#8eaaba` |
| 38 | NGC4689 | 0.0053500 | `#bae3c2` |
| 39 | NGC4694 | 0.0038690 | `#6500ff` |
| 40 | NGC4698 | 0.0033660 | `#008a8a` |
| alias | NGC4567_8 | 0.0074705 | `#04aa75` |
| reserved | TOTAL_ASSEMBLED | n/a | `black` |

## Policy For Future MAUVE Plots

Use this fixed color scheme whenever galaxy identity is encoded by color in MAUVE notebooks, scripts, figures, or reports. Do not use `tab20` for the 40-galaxy MAUVE identity mapping. Use black only for total assembled or pooled all-galaxy references.
