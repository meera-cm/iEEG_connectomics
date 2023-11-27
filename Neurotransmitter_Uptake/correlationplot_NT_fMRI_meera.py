import os
import pathlib
import random
import re
from typing import Union

import cycler
import matplotlib as mpl
from matplotlib import figure
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.stats
import seaborn as sns
import gaussianize

import sys

sys.path.append(r'C:\\Users\\ICN_guest\\Downloads\\data_meera')

os.chdir(r'C:\Users\ICN_guest\Downloads\data_meera')
print(os.getcwd())
print(sys.argv[0])
print(os.path.dirname(os.path.realpath('__file__')))

#IN_PATH = pathlib.Path(__file__).parent
IN_PATH = r"C:\Users\ICN_guest\Downloads\data_meera"
FNAME_FMAP = r"sbeta_str_network_n982_compound_atlas_HCPex_SUIT_ABGT.csv"
#FNAMES_FDOPA = {
 #   "fdopa_hc12": "FDOPA_fluorodopa_hc12_gomez_compound_atlas_HCPex_SUIT_ABGT.csv",
  #  "d1_hc13": "D1_SCH23390_hc13_kaller_compound_atlas_HCPex_SUIT_ABGT.csv",
   # "d2_alakurtti": "D2_raclopride_hc7_alakurtti_compound_atlas_HCPex_SUIT_ABGT.csv",
    #"DAT_hc6": "DAT_fepe2i_hc6_sasaki_compound_atlas_HCPex_SUIT_ABGT.csv",
    #"DAT_hc174": "DAT_fpcit_hc174_dukart_spect_compound_atlas_HCPex_SUIT_ABGT.csv",
    #"d2_smith": "D2_flb457_hc37_smith_compound_atlas_HCPex_SUIT_ABGT.csv",
    #"DA_aggregate_ALL": "Dopamine_aggregate_map_compound_atlas_HCPex_SUIT_ABGT.csv",
    #"D2_aggregate": "",
    #"DAT_aggregate": ""
#}
FNAMES_FDOPA = {"Dopamine_aggregate":"Dopamine_aggregate_compound_atlas_HCPex_SUIT_ABGT.csv"}
X_LIMS = (-3.5, 3.5)
Y_LIMS = (-3.5, 3.5)
FIGSIZE = 3  # height in inches, width will be adapted accordingly


def activate_plotting_settings() -> None:
    PALETTE_WJN_2022 = np.array(
        (
            (55 / 255, 110 / 255, 180 / 255),  # blue
            (113 / 255, 188 / 255, 173 / 255),  # green
            (223 / 255, 74 / 255, 74 / 255),  # red
            (125 / 255, 125 / 255, 125 / 255),  # medium grey
            (150 / 255, 150 / 255, 150 / 255),  # medium grey
            (215 / 255, 178 / 255, 66 / 255),  # orange
            (100 / 255, 100 / 255, 100 / 255),  # dark grey
            (215 / 255, 213 / 255, 203 / 255),  # medium/light grey
            (234 / 255, 234 / 255, 234 / 255),  # light grey
        )
    )
    mpl.rcParams["font.family"] = "Arial"
    mpl.rcParams["font.size"] = 8

    mpl.rcParams["pdf.fonttype"] = "TrueType"
    mpl.rcParams["svg.fonttype"] = "none"

    mpl.rcParams["savefig.dpi"] = 300
    mpl.rcParams["savefig.bbox"] = "tight"
    mpl.rcParams["savefig.pad_inches"] = 0
    mpl.rcParams["savefig.transparent"] = False

    mpl.rcParams["axes.prop_cycle"] = cycler.cycler(color=PALETTE_WJN_2022)
    mpl.rcParams["axes.xmargin"] = 0

    mpl.rcParams["legend.title_fontsize"] = mpl.rcParams["legend.fontsize"]

    mpl.rcParams["lines.linewidth"] = 1


def gaussianize_data(data: pd.DataFrame, column_name: str) -> pd.DataFrame:
    out = gaussianize.Gaussianize(strategy="brute", max_iter=200, verbose=True)
    out.fit(data[column_name])
    data[column_name] = out.transform(data[column_name])
    return data


def patch_affinity_svg(svg_text: str) -> str:
    """Patch Matplotlib SVG so that it can be read by Affinity Designer."""
    matches = [
        x
        for x in re.finditer(
            'font:( [0-9][0-9]?[0-9]?[0-9]?)? ([0-9.]+)px ([^;"]+)[";]',
            svg_text,
        )
    ]
    if not matches:
        return svg_text
    svg_pieces = [svg_text[: matches[0].start()]]
    for i, match in enumerate(matches):
        # Change "font" style property to separate "font-size" and
        # "font-family" properties because Affinity ignores "font".
        group = match.groups()
        if len(group) == 2:
            font_weight, font_size_px, font_family = match.groups()
            new_font_style = (
                f"font-size: {float(font_size_px):.1f}px; "
                f"font-family: {font_family}"
            )
        else:
            font_weight, font_size_px, font_family = match.groups()
            if font_weight is not None:
                new_font_style = (
                    f"font-weight: {font_weight}; "
                    f"font-size: {float(font_size_px):.1f}px; "
                    f"font-family: {font_family}"
                )
            else:
                new_font_style = (
                    f"font-size: {float(font_size_px):.1f}px; "
                    f"font-family: {font_family}"
                )
        svg_pieces.append(new_font_style)
        if i < len(matches) - 1:
            svg_pieces.append(
                svg_text[match.end() - 1 : matches[i + 1].start()]
            )
        else:
            svg_pieces.append(svg_text[match.end() - 1 :])
    return "".join(svg_pieces)


def save_fig(fig: figure.Figure, outpath: Union[str, os.PathLike]) -> None:
    outpath = str(outpath)
    fig.savefig(outpath)  # , bbox_inches="tight")
    if outpath.endswith(".svg"):
        with open(outpath, "r", encoding="utf-8") as f:
            svg_text = f.read()
        patched_svg = patch_affinity_svg(svg_text)
        with open(outpath, "w", encoding="utf-8") as f:
            f.write(patched_svg)


def spearmans_rho_permutation(
    x: Union[np.ndarray, pd.Series],
    y: Union[np.ndarray, pd.Series],
    n_perm: int = 10000,
) -> tuple[float, float]:
    """Calculate permutation test for multiple repetitions of Spearmans Rho

    https://towardsdatascience.com/how-to-assess-statistical-significance-in-your-data-with-permutation-tests-8bb925b2113d

    Parameters
    ----------
    x (np array) : first distibution
    y (np array) : second distribution
    n_permp (int): number of permutations

    Returns
    -------
    gT (float) : estimated ground truth, here spearman's rho
    p_val (float) : p value of permutation test
    """

    # compute ground truth difference
    gT = scipy.stats.spearmanr(x, y)[0]
    #
    pV = np.array((x, y))
    # Initialize permutation:
    pD = []
    # Permutation loop:
    args_order = np.arange(0, pV.shape[1], 1)
    args_order_2 = np.arange(0, pV.shape[1], 1)
    for _ in range(n_perm):
        # Shuffle the data:
        random.shuffle(args_order)
        random.shuffle(args_order_2)
        # Compute permuted absolute difference of your two sampled
        # distributions and store it in pD:
        pD.append(
            scipy.stats.spearmanr(pV[0, args_order], pV[1, args_order_2])[0]
        )

    # calculate p value
    if gT < 0:
        p_val = (len(np.where(pD <= gT)[0]) + 1) / (n_perm + 1)
    else:
        p_val = (len(np.where(pD >= gT)[0]) + 1) / (n_perm + 1)

    return gT, p_val


def main() -> None:
    activate_plotting_settings()
    x = "F-Map [AU]"
    hue = "Structure"

    # read fmap values
    IN_PATH = pathlib.Path(r'C:\Users\ICN_guest\Downloads\data_meera')
    #fpath_fmap = IN_PATH / FNAME_FMAP
    # = pathlib.Path(r'C:\Users\ICN_guest\Downloads\data_meera\beta_spmT_0001_compound_atlas_HCPex_SUIT_ABGT.csv')
    fpath_fmap = pathlib.Path(r'C:\Users\ICN_guest\Downloads\data_meera\sbeta_str_network_n982_compound_atlas_HCPex_SUIT_ABGT.csv')
    fmap = pd.read_csv(fpath_fmap, index_col="Name").rename(columns={"Value": x})

    for descr, fname_fdopa in FNAMES_FDOPA.items():
        if descr == "Dopamine_aggregate":
            y = "Dopamine Aggregate [a.u.]"
       # elif descr == "fdopa_hc12":
        #    y = "Fluorodopa [a.u.]"
        #elif descr == "d1_hc13":
         #   y = "D1 Receptor [a.u.]"
        #elif descr == "d2_alakurtti":
        #    y = "D2 Alakurtti [a.u.]"
        #elif descr == "d2_smith":
        #    y = "D2 Smith [a.u.]"
        #elif descr == "DAT_hc6":
         #   y = "DAT HC 6 [a.u.]"
        #elif descr == "DAT_hc174":
         #   y = "DAT HC 174 [a.u.]"

        else:
            raise ValueError(
                f"Unknown description: {descr}. Please provide a value for y."
            )

        fpath_fdopa = IN_PATH / fname_fdopa
        # read fdopa values
        fdopa = pd.read_csv(fpath_fdopa, index_col="Name").rename(
            columns={"Value": y}
        )
        # preprocess data
        data = pd.concat((fmap, fdopa), axis="columns", join="outer")
        data[hue] = data[hue].str.capitalize()
        data = data[data[x] != 0]
        data = data[data[y] != 0]
        data = gaussianize_data(data=data, column_name=x)
        data = gaussianize_data(data=data, column_name=y)

        # do statistics
        new_labels = {}
        for item in data[hue].unique():
            data_xy = data[data[hue] == item]
            data_x = data_xy[x].to_numpy()
            data_y = data_xy[y].to_numpy()
            rho, p = spearmans_rho_permutation(data_x, data_y, n_perm=10000)
            res_lin = scipy.stats.linregress(data_x, data_y)
            p_lin = res_lin.pvalue
            r_lin = res_lin.rvalue
            print(
                f"\n{item}:\n",
                f"Correlation {x} - {y}:\n"
                f" - Rho={'{:.2f}'.format(rho)},  P={'{:.4f}'.format(p)}\n"
                f" - r = {r_lin:.2f}, P={p_lin:.4f}",
            )
            new_label = (
                f"{item}\n"
                f" \u03C1={rho:.2f}, P={p:.4f}\n"
                f" r={r_lin:.2f}, P={p_lin:.4f}"
            )
            new_labels[item] = new_label
        data.replace(new_labels, inplace=True)

        # visualise with seaborn
        grid = sns.lmplot(
            data=data, x=x, y=y, hue=hue, palette="mako", height=3
        )
        grid.legend.set_title(None)

        # make plot look nicer
        ax = grid.ax
        ax.set_xlim([X_LIMS[0], X_LIMS[1]])
        ax.set_xticks([X_LIMS[0], 0, X_LIMS[1]])
        ax.set_ylim([Y_LIMS[0], Y_LIMS[1]])
        ax.set_yticks([Y_LIMS[0], 0, Y_LIMS[1]])
        # ax.set_ybound(lower=Y_LIMS[0], upper=Y_LIMS[1])
        ax.spines["left"].set_position(("outward", 3))
        ax.spines["bottom"].set_position(("outward", 3))

        basename = f"correlation_fmap_{descr}"
        save_fig(grid, IN_PATH / f"{basename}.png")
        save_fig(grid, IN_PATH / f"{basename}.svg")
    plt.show(block=True)


if __name__ == "__main__":
    main()
