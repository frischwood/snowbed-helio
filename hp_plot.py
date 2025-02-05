from pathlib import Path
from HP_plotter import HP_plotter
import argparse

def get_parser():
    parser = argparse.ArgumentParser(
                    prog='ProgramName',
                    description='What the program does',
                    epilog='Text at the bottom of help')
    parser.add_argument("-sf","--simufolder", type=str, 
                        required=True,
                        help="Provide the path to the snowBedFoam simulations folder")
    parser.add_argument("-v","--variable", type=str, 
                        choices=["mean","median"],
                        required=False, default="mean",
                        help="Provide variable to plot")
    return parser

if __name__ == "__main__":
    # mycmd: python hp_plot.py -sf "/Users/frischho/Documents/sunwell/gondo/simulations/snowbedfoam/experiments" -v "mean"
    parser = get_parser()
    args = parser.parse_args()
    simu_path = args.simufolder
    var = args.variable
    # set path to stats and to output location for plots
    stat_path = Path(simu_path) / "analysis" / "stats"
    out_path = Path(simu_path) / "analysis" / "plots"/ "pub"
    out_path.mkdir(exist_ok=True, parents=True)

    # get plotter instance
    hp_plotter = HP_plotter(stat_path, out_path, var)

    # generic topology
    hp_plotter.plot_generic_topo()

    # all simulations
    hp_plotter.plot_exp_group()

    # differences
    hp_plotter.plot_sectors_diff()

    # correlations
    hp_plotter.plot_corr()

    # plot sectorial SMD of reference simulations
    hp_plotter.plot_sector_singleRef()
    hp_plotter.plot_sector_groupRef()