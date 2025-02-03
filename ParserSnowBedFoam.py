import os
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

class ParserSnowBedFoam():

    def __init__(self, simu_root, coordsCode, timestepCode, sd, nb_PV, dim_PV, inOutRatio):

        """Parses some ascii openfoam output specific files at res=res and builds a ndarray of it.
        res: resolution of grid in meters"""
        
        self.simu_root = simu_root
        self.coordsCode = coordsCode
        self.timestepCode = timestepCode
        self.sd = sd
        self.nb_PV = nb_PV
        self.dim_PV = dim_PV 
        self.PVcenters = self.get_PVcenters()
        self.snowDriftGrid_x, self.snowDriftGrid_y, self.snowDriftGrid = self.get_snowDriftGrid()
        self.dump_grid(self.snowDriftGrid)
        # compute stats
        helio_radius = self.dim_PV["W"] + self.dim_PV["p"]
        self.squareDriftSectors = self.computeSquareDriftSectors(inSquareSide=2*helio_radius,inOutRatio=inOutRatio)
        self.circDriftSectors = self.computeCircDriftSectors(inCircleRad=helio_radius, inOutRatio=inOutRatio)
        
    def get_PVcenters(self):
        """computes center of helioplant units (mean of Cx and Cy). If self.nb_pv is an int: then takes all int from 1 to self.nb_pv included. IF self.nb_pv is a list it takes only the elements of the list.
        Returns:
            dict(tuple): (Cx,Cy) of all helioplants units
        """
        if type(self.nb_PV)==list:
            pv_nb_iter = self.nb_PV # in case a list of pv number is expressively given
        else:
            pv_nb_iter = range(1,self.nb_PV+1) 
        
        pvCenters = {}
        for i in pv_nb_iter:
            grid_dict = {"Cx":None, "Cy":None}
            for key in grid_dict.keys():
                grid_dict[key] = self.extract_openfoamascii_values(self.simu_root / str(self.coordsCode) / key, f"solarPV_{i}")  
            pvCenters[i] = (round(grid_dict["Cx"].mean(),3),
                            round(grid_dict["Cy"].mean(),3))
        return pvCenters
    
    def get_snowDriftGrid(self):
        """_summary_

        Args:
            sd (float, optional): sampling distance (m). Defaults to 0.1.

        Returns:
            _type_: _description_
        """
        gridMass_dict = {Path(str(self.coordsCode)) / "Cx":None, Path(str(self.coordsCode)) / "Cy":None, Path(str(self.timestepCode)) / "kinematicCloud-massCheckPatterns":None}
        for key in gridMass_dict.keys():
            gridMass_dict[key] = self.extract_openfoamascii_values(self.simu_root /  key, "snowBed")

        # Bin the data onto a res x res grid
        # Have to reverse x & y due to row-first indexing
        (_,x),(_,y),(_,z) = gridMass_dict.items()
        # y,x = gridMass_dict[f"{self.coordsCode}/Cy"], gridMass_dict[f"{self.coordsCode}/Cx"]
        # z = gridMass_dict[f"{self.timestepCode}/kinematicCloud-massCheckPatterns"]

        y_extent = y.max()-y.min()
        x_extent = x.max()-x.min()
        domainFormatHW = (y_extent)/(x_extent) # to compute square cells
        bins_x =  x_extent / self.sd # number of bins to get 
        zi, yi, xi = np.histogram2d(y,x, bins=(round(domainFormatHW*bins_x),round(bins_x)), weights=z, density=False)
        counts, _, _ = np.histogram2d(y,x, bins=(round(domainFormatHW*bins_x),round(bins_x)), density=False)
        zi = zi / counts # mean value
        # zi = np.nan_to_num(zi)

        # plot for debug
        self.plot_grids(x,y,z,xi,yi,zi)
        return xi, yi, zi 
    @staticmethod
    def extract_openfoamascii_values(filepath, matcherSection):
        """ reads openfoam ascii folder and return value in matchersection"""
        
        values = []
        with open(filepath,'r') as file:
            lines = file.readlines()
            eof = False
            linenum = 0
            inMatcherSection = False
            inMatcherValuesSection = False
            while not eof:
                line = lines[linenum]

                if inMatcherSection & inMatcherValuesSection:
                    if ")" in line:
                        eof = True
                    else:
                        values.append(np.fromstring(line, dtype=float, sep=','))
                    linenum += 1
                    continue # pass following unnecessary checks

                if (matcherSection in line) or inMatcherSection: # check if got to the right section
                    inMatcherSection = True
                    if inMatcherSection & ("(" in line): # Once in section, check if got to the values part
                        inMatcherValuesSection = True 
                
                # go to next line
                linenum += 1
            
            values_vec = np.concatenate(values)
            return values_vec

    def computeSquareDriftSectors(self, inSquareSide, inOutRatio):
        """Gets stats for quadratic sector"""
        PVdriftSectors = {}
        inOutSquareSideRatio = inOutRatio # how large is largesquare
         
        for i, (key, item) in enumerate(self.PVcenters.items()):
            s1_mask_in = self.get_cartGridMasks(item[1],(item[1] + inSquareSide/2), item[0], (item[0] + inSquareSide/2))
            s1_mask_out = self.get_cartGridMasks(item[1], (item[1] + inOutSquareSideRatio*inSquareSide/2), item[0], (item[0] + inOutSquareSideRatio*inSquareSide/2))# trigo sector 1
            s2_mask_in = self.get_cartGridMasks(item[1], (item[1] + inSquareSide/2), (item[0] - inSquareSide/2), item[0])
            s2_mask_out = self.get_cartGridMasks(item[1],(item[1] + inOutSquareSideRatio*inSquareSide/2), (item[0] - inOutSquareSideRatio*inSquareSide/2), item[0])# trigo sector 2
            s3_mask_in = self.get_cartGridMasks((item[1] - inSquareSide/2), item[1] ,(item[0] - inSquareSide/2), item[0])
            s3_mask_out = self.get_cartGridMasks((item[1] - inOutSquareSideRatio*inSquareSide/2), item[1] ,(item[0] - inOutSquareSideRatio*inSquareSide/2), item[0])# trigo sector 3
            s4_mask_in = self.get_cartGridMasks((item[1] - inSquareSide/2), item[1] , item[0], (item[0] + inSquareSide/2))
            s4_mask_out = self.get_cartGridMasks((item[1] - inOutSquareSideRatio*inSquareSide/2), item[1], item[0], (item[0] + inOutSquareSideRatio*inSquareSide/2))# trigo sector 4 
           
            PVdriftSectors[(key,item[0],item[1])] = self.computeSectorStats([s1_mask_in,s2_mask_in,s3_mask_in,s4_mask_in],
                                                         [s1_mask_out,s2_mask_out,s3_mask_out,s4_mask_out],sectorType="cart")
        return PVdriftSectors

    def computeCircDriftSectors(self, inCircleRad, inOutRatio):
        """Gets stats for polar sectors"""
        PVdriftSectors = {}
        inOutCircRadRatio = inOutRatio # how large is largesquare
         
        for i, (key, item) in enumerate(self.PVcenters.items()):
            s1_mask_in = self.get_polarGridMasks(item[0],item[1], phiStart=0, phiEnd=90,r=inCircleRad)
            s1_mask_out = self.get_polarGridMasks(item[0],item[1], phiStart=0, phiEnd=90,r=inOutCircRadRatio*inCircleRad)# trigo sector 1
            s2_mask_in= self.get_polarGridMasks(item[0],item[1], phiStart=90, phiEnd=180,r=inCircleRad)
            s2_mask_out= self.get_polarGridMasks(item[0],item[1], phiStart=90, phiEnd=180,r=inOutCircRadRatio*inCircleRad)# trigo sector 2
            s3_mask_in = self.get_polarGridMasks(item[0],item[1], phiStart=-180, phiEnd=-90,r=inCircleRad)
            s3_mask_out = self.get_polarGridMasks(item[0],item[1], phiStart=-180, phiEnd=-90,r=inOutCircRadRatio*inCircleRad)# trigo sector 3
            s4_mask_in = self.get_polarGridMasks(item[0],item[1], phiStart=-90, phiEnd=0,r=inCircleRad)
            s4_mask_out = self.get_polarGridMasks(item[0],item[1], phiStart=-90, phiEnd=0,r=inOutCircRadRatio*inCircleRad)# trigo sector 4
            
            PVdriftSectors[(key,item[0],item[1])]= self.computeSectorStats([s1_mask_in,s2_mask_in,s3_mask_in,s4_mask_in],[s1_mask_out,s2_mask_out,s3_mask_out,s4_mask_out],sectorType="polar")
        return PVdriftSectors
        
    def computeSectorStats(self,masks_in, masks_out, sectorType="polar"):
        """compute stats for each sector"""
        s1_mask_in, s2_mask_in, s3_mask_in, s4_mask_in = masks_in
        s1_mask_out, s2_mask_out, s3_mask_out, s4_mask_out = masks_out
        # get values in circle
        s1=self.snowDriftGrid[s1_mask_in]
        s2=self.snowDriftGrid[s2_mask_in]
        s3=self.snowDriftGrid[s3_mask_in]
        s4=self.snowDriftGrid[s4_mask_in]
        inCircleMask = s1_mask_in + s2_mask_in + s3_mask_in + s4_mask_in
        inCircle = self.snowDriftGrid[inCircleMask]# in circle
        inWwMask = s2_mask_in + s3_mask_in
        inWw = self.snowDriftGrid[inWwMask]# in windward
        inLeeMask = s1_mask_in + s4_mask_in
        inLee = self.snowDriftGrid[inLeeMask]# in lee
        # get values out circle
        largeCircleMask = s1_mask_out + s2_mask_out + s3_mask_out + s4_mask_out
        outCircle = self.snowDriftGrid[largeCircleMask * (~inCircleMask)]
        largeWwMask = s2_mask_out + s3_mask_out
        outWw = self.snowDriftGrid[largeWwMask * (~inWwMask)]
        largeLeeMask = s1_mask_out + s4_mask_out
        outLee = self.snowDriftGrid[largeLeeMask * (~inLeeMask)]
        # get stats 
        sectorStats_dict = {f"s1_{sectorType}":self.getArrayStats(s1),f"s2_{sectorType}":self.getArrayStats(s2),
                                f"s3_{sectorType}":self.getArrayStats(s3),f"s4_{sectorType}":self.getArrayStats(s4),
                            f"inWw_{sectorType}":self.getArrayStats(inWw), f"inLee_{sectorType}":self.getArrayStats(inLee), 
                            f"in_{sectorType}":self.getArrayStats(inCircle),
                            f"outWw_{sectorType}":self.getArrayStats(outWw), f"outLee_{sectorType}":self.getArrayStats(outLee), 
                            f"out_{sectorType}": self.getArrayStats(outCircle) } 
        return sectorStats_dict

    def getArrayStats(self, array):
        """compute stats"""
        stats = {}
        stats["median"] = np.nanmedian(array)
        stats["mean"] = np.nanmean(array)
        stats["std"] = np.nanstd(array)
        stats["min"] = np.nanmin(array)
        stats["max"] = np.nanmax(array)
        stats["p25"] = np.nanpercentile(array,25)
        stats["p75"] = np.nanpercentile(array,75)
        stats["eros_median"] = np.nan_to_num(np.nanmedian(array[array<0]))
        stats["eros_mean"] = np.nan_to_num(np.nanmean(array[array<0]))
        stats["eros_sum"] = np.nan_to_num(np.nansum(array[array<0]))
        stats["depo_median"] = np.nan_to_num(np.nanmedian(array[array>0]))
        stats["depo_mean"] = np.nan_to_num(np.nanmean(array[array>0]))
        stats["depo_sum"] = np.nan_to_num(np.nansum(array[array>0]))
        stats["sector_area"] = array.flatten().shape[0] * self.sd**2 # [m^2]
        return stats

    def get_cartGridSliceIndices(self,  yCoordStart, yCoordEnd, xCoordStart, xCoordEnd):
        """computes cartesian grid indices slices from cartesian grid coordinates extents""" 
        x_indices = np.argwhere((self.snowDriftGrid_x >= xCoordStart) * (self.snowDriftGrid_x <= xCoordEnd)).squeeze(-1)
        y_indices = np.argwhere((self.snowDriftGrid_y >= yCoordStart) * (self.snowDriftGrid_y <= yCoordEnd)).squeeze(-1)
        xSlice =  slice(x_indices[0],x_indices[-1]+1) 
        ySlice = slice(y_indices[0], y_indices[-1]+1)
        return ySlice, xSlice
    
    def get_cartGridMasks(self,  yCoordStart, yCoordEnd, xCoordStart, xCoordEnd):
        """computes cartesian grid masks from coords masks""" 
        x_mask= (self.snowDriftGrid_x[np.newaxis,1:] >= xCoordStart) * (self.snowDriftGrid_x[np.newaxis,1:] <= xCoordEnd)
        y_mask = (self.snowDriftGrid_y[1:,np.newaxis] >= yCoordStart) * (self.snowDriftGrid_y[1:,np.newaxis] <= yCoordEnd)
        mask_sector = x_mask & y_mask
        # # plot for debug:
        # self.snowDriftGrid[mask_sector] = 123.
        # plt.figure(figsize=(6, 6))
        # plt.pcolormesh(self.snowDriftGrid_x, self.snowDriftGrid_y, self.snowDriftGrid, alpha=0.5)
        # plt.colorbar()
        # plt.close()
 
        return mask_sector

    def get_polarGridMasks(self,xCenter, yCenter, phiStart, phiEnd, r):
        """computes cart grid masks from polar coords masks

        Args:
            xCenter (_type_): Helioplant's centroid x coord
            yCenter (_type_): Helioplant's centroid y coord
            phiStart (_type_): start polar sector angle in trigo angle[-180,180]
            phiEnd (_type_): end polar sector angle in trigo angle[-180,180]
            r (_type_): mask radius

        Returns:
            _type_: cart mask
        """
        # mask for r (must skip first index because coords are grid cells boundaries (there're one more than cells))
        mask_r = (self.snowDriftGrid_x[np.newaxis,1:]-xCenter)**2 + (self.snowDriftGrid_y[1:,np.newaxis]-yCenter)**2 < r**2
        # mask for phi
        atan2_phi = np.rad2deg(np.arctan2(self.snowDriftGrid_y[1:,np.newaxis]-yCenter, self.snowDriftGrid_x[np.newaxis,1:]-xCenter))
        mask_phi = (atan2_phi >= phiStart) & (atan2_phi < phiEnd)
        # combine
        mask_sector = mask_r & mask_phi

        # # plot for debug:
        # self.snowDriftGrid[mask_sector] = 223.
        # plt.figure(figsize=(6, 6))
        # plt.pcolormesh(self.snowDriftGrid_x, self.snowDriftGrid_y, self.snowDriftGrid,alpha=0.5)
        # plt.colorbar()
        # plt.close()

        return mask_sector

    def plot_grids(self, x,y,z, xi,yi,zi):
        """plots parsed data for debug"""
        # make output folder
        out_folder = self.simu_root / "plots"
        out_folder.mkdir(exist_ok=True, parents=True)
        # global scatter
        fig, ax = plt.subplots(figsize=(12,12))
        scat = ax.scatter(x, y, c=z, s=10, cmap="coolwarm", vmin=-1.5, vmax=1.5)
        ax.set_xlim(-20,20)
        ax.set_ylim(-20,20)
        cb = fig.colorbar(scat)
        cb.set_label("snow mass distribution [kg/$m^2$]")
        ax.set_title(f"Scatter: {x.shape[0]} particles\n {os.path.basename(self.simu_root)}")
        plt.savefig(out_folder / f"scatter_global.png")
        plt.close()
        # scatter vs grid
        figmesh, axmesh = plt.subplots(1,2,sharey=True,figsize=(12,6))
        scat = axmesh[0].scatter(x, y, c=z, s=10, cmap="coolwarm", vmin=-1.5, vmax=1.5)
        axmesh[0].set_xlim(-15,15)
        axmesh[0].set_ylim(-15,15)
        cmesh = axmesh[1].pcolormesh(xi, yi, zi, edgecolors='white', lw=0.01, cmap="coolwarm",
                                  vmin=-1.5,vmax=1.5) 
        axmesh[1].set_xlim(-7,7)
        axmesh[1].set_ylim(-15,7)        
        cbmesh = figmesh.colorbar(cmesh)
        cbmesh.set_label("snow mass distribution [kg/$m^2$]")
        figmesh.suptitle(f"Mesh {self.sd} [m] vs. Scatter\n {os.path.basename(self.simu_root)}")
        plt.savefig(out_folder / f"scatter_vs_grid_sd{self.sd}.png")
        plt.close()

    def dump_grid(self, np_array):
        """writes out grids"""
        # make output folder
        out_folder = self.simu_root / "grids"
        out_folder.mkdir(exist_ok=True, parents=True)
        np.savetxt(out_folder / "snowDriftGrid.txt", X=np_array)

if __name__=="__main__":
    """for debug
    """
    simu_root = "simulations/snowbedfoam/experiments/REFERENCE/helioplant_n9_h07_sp5_ARL45_BL0015_windDir5"
    sd = 0.25 # sampling_distance in meter
    inOutRatio=2
    nb_PV = 9 # number of helioplant units in simu
    helio_dims = {"W":1.8915, "L":5.4, "t":0.2, "p":0.381, "h":0.7, "H":6.203} # dims of helioplant
    ParserSnowBedFoam(simu_root, coordsCode=0, timestepCode=300, sd=sd, nb_PV=nb_PV, dim_PV=helio_dims, inOutRatio=inOutRatio)