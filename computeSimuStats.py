import os
import sys
import pickle
import pandas as pd
from ParserSnowBedFoam import ParserSnowBedFoam
from pathlib import Path
import argparse

def nestedToMultiIndexDict(nested_dict):
    """outerkey must be tuple of len==3"""
    reformed_dict = {} 
    for outerKey, innerDict in nested_dict.items(): 
        for innerKey, values in innerDict.items(): 
            reformed_dict[(outerKey[0],outerKey[1],outerKey[2],
                        innerKey)] = values 
    
    return reformed_dict

def get_parser():
    parser = argparse.ArgumentParser(
                    prog='ProgramName',
                    description='What the program does',
                    epilog='Text at the bottom of help')
    parser.add_argument("-sf","--simufolder", type=str, 
                        required=True,
                        help="Provide the path to the snowBedFoam simulations folder")
    return parser

if __name__=="__main__":
    # mycmd: python computeSimuStats.py -sf /Users/frischho/Documents/sunwell/gondo/simulations/snowbedfoam/experiments
    parser = get_parser()
    args = parser.parse_args()
    simu_path = Path(args.simufolder)

    coordsCode = 0 # snowBedFoam 
    timestepCode = 300 # snowBedFoam
    samplingDistance = 0.3 #[m] size of sampling cells 
    inOutRatio = 2 # ratio between inner and outer analysis area radius length
    analysisName = f"mesh_{samplingDistance}"
    analysisDir = simu_path / "analysis" / "stats" / analysisName
    circleDir = analysisDir / "polar"
    squareDir = analysisDir / "cart"
    analysisDir.mkdir(parents=True, exist_ok=True)
    circleDir.mkdir(parents=True, exist_ok=True)
    squareDir.mkdir(parents=True, exist_ok=True)

    # dimensions of a helioplant unit
    helio_dims = {"W":1.8915, "L":5.4, "t":0.2, "p":0.381, "h":0.7, "H":6.203} 
    # helioplant unit number to plot by exp and param value. An exaustive list of pv numbers can also be given. 
    helio_number = {
                    "SINGLE":{1:1},
                    "AZIMUTH":{0:9,30:9,45:9},
                    "FOREST":{9:9,16:16,25:25},
                    "HEIGHT":{0.15:9, 0.2:9, 0.3:9, 0.4:9, 0.5:9, 0.6:9, 0.7:9},
                    "INTERSPACE":{8:9,10:9,12:9},
                    "QUINCUX":{8:8,16:[1,3,4,5,6,7,8,9,10,11,12,14,15,16],23:[1,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20,21,22,23,24,25]}
                    }
    # experiment paths
    simuFolders = {
                "SINGLE":{
                            1:simu_path / "SINGLE" / "helioplant_n1_h07_sp5_ARL45_BL0015_windDir5"
                        },
                "AZIMUTH":{
                            0:simu_path / "REFERENCE" / "helioplant_n9_h07_sp5_ARL45_BL0015_windDir5",
                            30:simu_path / "AZIMUTH" / "helioplant_n9_h07_sp5_ARL45_BL0015_az30_windDir5",
                            45:simu_path / "AZIMUTH" / "helioplant_n9_h07_sp5_ARL45_BL0015_az45_windDir5"
                        },
                "FOREST":{
                            9:simu_path / "REFERENCE" / "helioplant_n9_h07_sp5_ARL45_BL0015_windDir5",
                            16:simu_path / "FOREST" / "helioplant_n16_h07_sp5_ARL45_BL0015_windDir5",
                            25:simu_path / "FOREST" / "helioplant_n25_h07_sp5_ARL45_BL0015_windDir5"
                            },
                "HEIGHT":{ 
                            0.15:simu_path / "HEIGHT" / "helioplant_n9_h015_sp5_ARL45_BL0015_windDir5",
                            0.2:simu_path / "HEIGHT" / "helioplant_n9_h02_sp5_ARL45_BL0015_windDir5",
                            0.3:simu_path / "HEIGHT" / "helioplant_n9_h03_sp5_ARL45_BL0015_windDir5",
                            0.4:simu_path / "HEIGHT" / "helioplant_n9_h04_sp5_ARL45_BL0015_windDir5",
                            0.5:simu_path / "HEIGHT" / "helioplant_n9_h05_sp5_ARL45_BL0015_windDir5",
                            0.6:simu_path / "HEIGHT" / "helioplant_n9_h06_sp5_ARL45_BL0015_windDir5",
                            0.7:simu_path / "REFERENCE" / "helioplant_n9_h07_sp5_ARL45_BL0015_windDir5",
                        },   
                "INTERSPACE":{
                            8:simu_path / "INTERSPACE" / "helioplant_n9_h07_sp8_ARL45_BL0015_windDir5",
                            10:simu_path / "REFERENCE" / "helioplant_n9_h07_sp5_ARL45_BL0015_windDir5",
                            12:simu_path / "INTERSPACE" / "helioplant_n9_h07_sp12_ARL45_BL0015_windDir5"
                        },
                "QUINCUX":{
                            8:simu_path / "QUINCUNX" / "qc_helioplant_n9_h07_sp5_ARL45_BL0015_windDir5",
                            16:simu_path / "QUINCUNX" / "qc_helioplant_n16_h07_sp5_ARL45_BL0015_windDir5",
                            23:simu_path / "QUINCUNX" / "qc_helioplant_n25_h07_sp5_ARL45_BL0015_windDir5",
                        }
                    }
    # loop on all exp
    simus= {}
    for key, simuFolder in simuFolders.items():
        simus_setups = {}
        for key_setup, simuSetup in simuFolder.items():
            print(f"Parsing {key} ---> {key_setup}")
            simus_setups[key_setup] = ParserSnowBedFoam(simu_root=simuSetup, coordsCode=coordsCode, timestepCode=timestepCode,sd=samplingDistance, nb_PV=helio_number[key][key_setup], dim_PV=helio_dims, inOutRatio=inOutRatio)
            # build df from stats dict and dump it
            df_setup_square = pd.DataFrame.from_dict(nestedToMultiIndexDict(simus_setups[key_setup].squareDriftSectors),orient="index")
            df_setup_circle = pd.DataFrame.from_dict(nestedToMultiIndexDict(simus_setups[key_setup].circDriftSectors),orient="index")
            df_setup_square.to_csv(squareDir / f"{key}_{key_setup}.csv")
            df_setup_circle.to_csv(circleDir / f"{key}_{key_setup}.csv")
        # store in all stats dict
        simus[key] = simus_setups
    # dump all stats dict 
    with open(analysisDir / "all_simus_stats.pkl","wb") as file:
        pickle.dump(simus, file)
        
        
            