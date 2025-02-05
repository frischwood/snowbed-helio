import os
import pandas as pd
import numpy as np
import shapely
from shapely.plotting import plot_polygon
from sklearn import preprocessing
import seaborn as sns
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D


class HP_plotter():
    def __init__(self, stat_path, out_path, var):
        
        # I/O
        self.out_path = out_path
        self.stat_path = stat_path
        # plot cosmetics
        self.fs_suplabels = 18
        self.fs_labels = 14
        self.fs_titles = 14 
        self.fs_legend = 14
        self.fs_cb = 16

        # select experiments type
        self.var = var
        self.sd = 0.3
        self.sector_type = "polar"
        # path of all exps
        self.exp_df_dict = {
            "SINGLE":{
                1:self.stat_path / f"mesh_{self.sd}" / self.sector_type / "SINGLE_1.csv"
                },
            "AZIMUTH":{
                0:self.stat_path / f"mesh_{self.sd}" / self.sector_type / "AZIMUTH_0.csv",
                30:self.stat_path / f"mesh_{self.sd}" / self.sector_type / "AZIMUTH_30.csv",
                45:self.stat_path / f"mesh_{self.sd}" / self.sector_type / "AZIMUTH_45.csv"
                },
            "HEIGHT":{
                0.2:self.stat_path / f"mesh_{self.sd}" / self.sector_type / "HEIGHT_0.2.csv",
                0.3:self.stat_path / f"mesh_{self.sd}" / self.sector_type / "HEIGHT_0.3.csv",
                0.4:self.stat_path / f"mesh_{self.sd}" / self.sector_type / "HEIGHT_0.4.csv",
                0.5:self.stat_path / f"mesh_{self.sd}" / self.sector_type / "HEIGHT_0.5.csv",
                0.6:self.stat_path / f"mesh_{self.sd}" / self.sector_type / "HEIGHT_0.6.csv",
                0.7:self.stat_path / f"mesh_{self.sd}" / self.sector_type / "HEIGHT_0.7.csv",
                },
            "INTERSPACE":{
                8:self.stat_path / f"mesh_{self.sd}" / self.sector_type / "INTERSPACE_8.csv",
                10:self.stat_path / f"mesh_{self.sd}" / self.sector_type / "INTERSPACE_10.csv",
                12:self.stat_path / f"mesh_{self.sd}" / self.sector_type / "INTERSPACE_12.csv"
                },
            "GROUP SIZE":{
                3:self.stat_path / f"mesh_{self.sd}" / self.sector_type / "GROUPSIZE_9.csv",
                4:self.stat_path / f"mesh_{self.sd}" / self.sector_type / "GROUPSIZE_16.csv",
                5:self.stat_path / f"mesh_{self.sd}" / self.sector_type / "GROUPSIZE_25.csv"},
            "ALIGNMENT":{
                "in-line":self.stat_path / f"mesh_{self.sd}" / self.sector_type / "GROUPSIZE_9.csv",
                "staggered":self.stat_path / f"mesh_{self.sd}" / self.sector_type / "ALIGNMENT_8.csv",
            #   15:f"{self.stat_path} / f"mesh_{self.sd}" / self.sector_type / QUINCUX_16.csv",
            #   23:f"{self.stat_path} / f"mesh_{self.sd}" / self.sector_type / QUINCUX_23.csv"       
                }
            }
        # path for generic exps
        self.exp_df_dict_generic = {
            "GROUP SIZE":{
                3:self.stat_path / f"mesh_{self.sd}" / self.sector_type / "GROUPSIZE_9.csv",
                4:self.stat_path / f"mesh_{self.sd}" / self.sector_type / "GROUPSIZE_16.csv",
                5:self.stat_path / f"mesh_{self.sd}" / self.sector_type / "GROUPSIZE_25.csv"},
            "ALIGNMENT":{
                8:self.stat_path / f"mesh_{self.sd}" / self.sector_type / "ALIGNMENT_8.csv",
                15:self.stat_path / f"mesh_{self.sd}" / self.sector_type / "ALIGNMENT_14.csv",
                23:self.stat_path / f"mesh_{self.sd}" / self.sector_type / "ALIGNMENT_23.csv"
                }
            }
        
        # some containers helpers
        # sector naming dict
        self.sector_name_dict = {"in_polar":"IN", "inLee_polar": "IN_LEE", "inWw_polar": "IN_WINDWARD", "out_polar": "OUT", "outLee_polar": "OUT_LEE", "outWw_polar": "OUT_WINDWARD"}  
        # parameter's unit
        self.param_units = {"HEIGHT": "m","AZIMUTH": "°",  "ALIGNMENT": "-", "INTERSPACE": "m",
               "GROUP SIZE": "row",} 
        # selected units number in snowBedFoam topology
        self.helio_topology_dict = {"AZIMUTH": [4,5,6],
                  "HEIGHT": [4,5,6],
                   "INTERSPACE": [4,5,6],
                    "GROUP SIZE": [7,8,9,15,23],
                    "ALIGNMENT": [4,5,6]}
        # reference value for each parameter
        self.ref_values = {"AZIMUTH":0,"HEIGHT":0.7,"ALIGNMENT":9, "INTERSPACE":10,"GROUP SIZE":9}
        # quantitative parameters
        self.quant_params = ["AZIMUTH", "HEIGHT", "INTERSPACE"]
        
        # for descr-respvar table
        self.idx = pd.IndexSlice
        self.df = self.get_df()
        # data table of group reference simulation.
        self.data_ref = None # 
        
    def get_df(self):
        """reads all stat csv into a single df (var x [param, sector, row])
        """
        dfs = []
        for f_num, (key, item_param) in enumerate(self.exp_df_dict.items()):
            for i, (df_key, df_path) in enumerate(item_param.items()):
                    df_temp = pd.read_csv(df_path, index_col=[0,3]).loc[self.idx[:,["in_polar","out_polar","inLee_polar","outLee_polar","inWw_polar","outWw_polar"]],["eros_mean","depo_mean","mean"]]
                    df_temp["param"] = key
                    df_temp["param_value"] = float(os.path.basename(df_path).split("_")[-1].split(".csv")[0])
                    df_temp = df_temp.set_index("param", append=True).swaplevel(i=0,j=-1)
                    dfs.append(df_temp.T)

        df = pd.concat(dfs,axis=1)
        return df

    def plot_generic_topo(self,):
        # analysed units for generic setups
        selected_units = {"GROUP SIZE":[[4,5,6],[7,8,9,15],[7,8,9,15,23]], "ALIGNMENT":[[4,5,6],[7,8,9,15],[7,11,9,16,23]]}

        fig, axes = plt.subplots(2,3,figsize=(12,8))#, sharex=True, sharey=True)
        alignments = ["IN-LINE", "STAGGERED"]
        rows_number = [3,4,5]
        W=1.89
        p=0.3
        rho=2
        r = W + p/2 
        for i, (key, item) in enumerate(self.exp_df_dict_generic.items()):
            for j, (df_key,df_path) in enumerate(item.items()):
                df_topo = pd.read_csv(df_path)
                
                axes[i,j].set_ylim([-30,30])
                axes[i,j].set_xlim([-20,40])
                axes[i,j].set_frame_on(True)
                axes[i,j].yaxis.set_ticks_position('none')
                axes[i,j].set_yticklabels("")
                axes[i,j].xaxis.set_ticks_position('none')
                axes[i,j].set_xticklabels("")
                if j == 0:
                    axes[i,j].set_ylabel(alignments[i], fontsize=self.fs_labels)
                if i > 0:
                    axes[i,j].set_xlabel(rows_number[j], fontsize=self.fs_labels)
                for x,y,unit_number in zip(df_topo.iloc[:,1].values,df_topo.iloc[:,2].values, df_topo.iloc[:,0].values):
                    axes[i,j].scatter(x, y, s=220, marker="+", color="k")
                    if unit_number in selected_units[key][j]:
                        inner = shapely.buffer(shapely.Point(x,y),r)
                        outer = shapely.buffer(shapely.Point(x,y),rho*r)
                        plot_polygon(inner,ax=axes[i,j], fc=(0,0,0,0), ec="k", linewidth=1, add_points=False)
                        plot_polygon(outer,ax=axes[i,j], fc=(0,0,0,0), ec="k", linewidth=1, add_points=False)

        # manual hack
        axes[1,1].set_xlim([-25,35])
        fig.legend(handles=[
            Line2D([0], [0], marker='+', lw=0, markerfacecolor=(0,0,0,0), markeredgecolor="k", markersize=15, label="Helioplant"),
            Line2D([0], [0], marker='o', lw=0, markerfacecolor=(0,0,0,0), markeredgecolor="k", markersize=15, label="Analysed sector")], loc="lower right", ncol=2, fontsize=self.fs_legend)
        fig.supylabel("Alignment", fontweight="bold", fontsize=self.fs_suplabels)
        fig.supxlabel("Number of rows", fontweight="bold", fontsize=self.fs_suplabels)


        plt.tight_layout()
        fig.savefig(self.out_path / "0_generic_topo.jpeg",dpi=150,bbox_inches='tight')

    def plot_exp_group(self):
        vmin_cmap = -1.3
        vmax_cmap = 1.3
        regions = ["in_polar","inWw_polar","inLee_polar","out_polar","outWw_polar","outLee_polar"]
        markers = ["o","s","x","*","^"]
        linestyles = ["-","--",":","-.",(0,(3,1,1,1))]
        custom_leg = [Line2D([0], [0],linestyle=linestyle, color="k", lw=1, label=i+1) for i, (linestyle, marker) in enumerate(zip(linestyles,markers))]# marker=marker
        fig, axes = plt.subplots(6,5, sharey=False, sharex=False, figsize=(len(regions) * 11/6,len(self.helio_topology_dict.keys()) * 12/5))
        for i, region in enumerate(regions):
            vmin, vmax = [], []
            for j, param in enumerate(self.param_units.keys()):
                df_sel = self.df.T.loc[self.idx[param, region, self.helio_topology_dict[param]]].reset_index()
                sns.lineplot(data=df_sel, x="param_value", y=self.var,  style="level_2", ax=axes[i,j], color="k",markers=False, legend=False) # hue="param"
                # sns.move_legend(axes[i,j], "upper left") #, bbox_to_anchor=(1, 1)
                axes[i,j].hlines(y=0, xmin=-1, xmax=10, color="w", lw=2, alpha=0.5)
                axes[i,j].grid(axis="y",alpha=0.5)
                axes[i,j].grid(axis="x",alpha=0.5)

                #store min max values
                vmin.append(df_sel[self.var].min())
                vmax.append(df_sel[self.var].max())

                # cosmetics
                if i == 0:
                    axes[i,j].set_title(param + f" [{self.param_units[param]}]", fontsize=self.fs_titles)

                if j == 0:
                    axes[i,j].set_ylabel(self.sector_name_dict[region], fontsize=self.fs_labels)
                else:
                    axes[i,j].set_ylabel("")
                
                # background
                x = np.arange(vmin_cmap,vmax_cmap,0.01)
                y = np.arange(df_sel["param_value"].min(), df_sel["param_value"].max()+0.1,0.1)
                axes[i,j].set_xticks(df_sel["param_value"].unique(),df_sel["param_value"].unique()) # y+1
                axes[i,j].set_xlim([y[0],y[-1]])  
                grad_grid = np.repeat(np.expand_dims(x,axis=1),repeats=len(y), axis=1)
                axes[i,j].pcolormesh(y,x,grad_grid,  cmap="coolwarm", shading="nearest")

            # set y_extent, background and xticks
            vlim_max = np.stack(vmax).max()
            vlim_min = np.stack(vmin).min()
            for j, param in enumerate(self.param_units.keys()):
                # for j, param in enumerate(self.param_units.keys()):
                axes[i,j].set_ylim([vlim_min-0.05,vlim_max+0.05])  
                # x and y ticks and labels
                if i < (len(regions)-1):
                                axes[i,j].xaxis.set_ticks_position('none')
                                axes[i,j].set_xticklabels("")
                                axes[i,j].set_xlabel("")
                else:
                    axes[i,j].set_xlabel(f"[{self.param_units[param]}]", fontsize=self.fs_labels)
                if j > 0:
                        axes[i,j].yaxis.set_ticks_position('none')
                        axes[i,j].set_yticklabels("")

        # some manual hacks
        axes[-1,2].set_xticklabels([ "in-line", "staggered"],rotation=0)

        # colorbar
        cbar_ax = fig.add_axes([1, 0.061, 0.03, 0.91])
        cmap = mpl.cm.coolwarm
        norm = mpl.colors.Normalize(vmin=x[0], vmax=x[-1])
        cb1 = mpl.colorbar.ColorbarBase(cbar_ax, cmap=cmap,
                                        norm=norm,
                                        orientation='vertical',
                                        )
        cb1.ax.text(.5, x[0]+0.05, "EROSION", ha='center', va='bottom',rotation=90, fontsize=16, fontweight="normal")
        cb1.ax.text(.5, x[-1]-0.05, "DEPOSITION", ha='center', va='top',rotation=90, fontsize=16, fontweight="normal")
        cb1.set_ticks([])

        fig.supxlabel("Parameter value", fontweight="bold", fontsize=self.fs_suplabels)
        fig.supylabel("Mean snow mass distribution $[kg/m^2]$", fontweight="bold", fontsize=self.fs_suplabels)

        # legend
        cbar_ax.legend(handles=custom_leg, title="Row", bbox_to_anchor=(1.5,0.98), ncols=1, fontsize=self.fs_legend)

        plt.tight_layout()
        plt.savefig(self.out_path / f"1_all_exp_{self.var}.jpeg",dpi=150,bbox_inches='tight')

    def plot_sectors_diff(self):
        cat_params = ["HEIGHT","AZIMUTH","ALIGNMENT","INTERSPACE", "GROUP SIZE",]
        diff_values = {"ALIGNMENT":[8], "GROUP SIZE":[16,25], "AZIMUTH":[30,45],"HEIGHT":[0.4,0.6],"INTERSPACE":[8,12] } 
        locs = {"ALIGNMENT":[4,5,6],"GROUP SIZE":[4,5,6], "AZIMUTH":[4,5,6],"HEIGHT":[4,5,6],"INTERSPACE":[4,5,6]}
        self.idx = pd.IndexSlice

        fig_circ, axes_circ = plt.subplots(len(cat_params),2, sharey=True, sharex=True,figsize=(11,10))

        for i,cat_param in enumerate(cat_params):
            for j,diff_value in enumerate(diff_values[cat_param]):
                df_cat = self.df.T.set_index("param_value",append=True).loc[cat_param,self.var]
                df_cat_diff = (df_cat.loc[self.idx[:,:,diff_value]] - df_cat.loc[self.idx[:,:,self.ref_values[cat_param]]])
                # keep ref data
                df_ref = df_cat.loc[self.idx[:,:,self.ref_values[cat_param]]]
                self.data_ref = df_ref.loc[self.idx[:,locs[cat_param]]].reset_index().pivot(index="level_0", columns="level_1")
                # plots
                data = df_cat_diff.loc[self.idx[:,locs[cat_param]]].reset_index().pivot(index="level_0", columns="level_1")
                # circles
                if j==0:
                    ylabel=f"{cat_param}"
                else:
                    ylabel=""
                sectors_dict = self.get_sectors(n=3)
                title=f"{diff_value} [{self.param_units[cat_param]}]"
                if cat_param == "ALIGNMENT": title = "Staggered" # manual hack
                if cat_param == "GROUP SIZE": title = f"{int(np.sqrt(int(diff_value)))} [{self.param_units[cat_param]}]"
                cbar_norm, cbar_cmap = self.plot_sectors(data=data, sectors_dict=sectors_dict, scale_f=1, cmap_dyn_range=1, ax=axes_circ[i,j], ylabel=ylabel, cmap="coolwarm", x_ticks=[0,10,20], x_ticklabels=[1,2,3], title=title)
                # cosmetics
                axes_circ[i,j].yaxis.set_ticks_position('none')
                axes_circ[i,j].set_yticklabels("")
        axes_circ[2,-1].set_axis_off()
        cbar_ax = fig_circ.add_axes([1, 0.061, 0.03, 0.91])
        cb = fig_circ.colorbar(mpl.cm.ScalarMappable(norm=cbar_norm, cmap=cbar_cmap),
                    cax=cbar_ax, orientation='vertical',extend='both')
        cb.set_label(label='$SMD_{mean} - SMD_{mean,ref}$ [$kg/m^2$]',fontsize=self.fs_cb)
        # plt.subplots_adjust(top=1.3)
        plt.tight_layout()
        plt.savefig(self.out_path / f"3_diff_{self.var}.jpeg", dpi=150,bbox_inches='tight')

    def plot_sector_singleRef(self):
        param = "SINGLE"
        fig, ax = plt.subplots(figsize=(1.9,2))
        df_single = self.df.T.set_index("param_value",append=True).loc[param,self.var]
        # plots
        data = df_single.reset_index().pivot(index="level_0", columns="level_1")[self.var]
        sectors_dict = self.get_sectors(n=1)
        cbar_norm, cbar_cmap = self.plot_sectors(data=data, sectors_dict=sectors_dict, scale_f=1, cmap_dyn_range=1, ax=ax, cmap="coolwarm", cb_vmax=1.5)
        # cosmetics
        cbar_ax = fig.add_axes([1, 0.061, 0.03, 0.91])
        cb = fig.colorbar(mpl.cm.ScalarMappable(norm=cbar_norm, cmap=cbar_cmap),
                    cax=cbar_ax, orientation='vertical',extend='both')
        cb.set_label(label='$SMD_{mean}$ [$kg/m^2$]',fontsize=self.fs_cb)
        ax.set_axis_off()
        fig.supxlabel("")
        plt.tight_layout()
        fig.savefig(self.out_path / f"2_sector_singleRef_{self.var}.png", dpi=150, transparent=True, bbox_inches="tight")

    def plot_sector_groupRef(self):
        fig_circ, axes_circ = plt.subplots(sharey=True, sharex=True,figsize=(11.4,4))
        # ref data sector plot
        sectors_dict = self.get_sectors(n=3)
        cbar_norm, cbar_cmap = self.plot_sectors(data=self.data_ref, sectors_dict=sectors_dict, scale_f=1, cmap_dyn_range=1, ax=axes_circ, cmap="coolwarm")
        # cosmetics
        axes_circ.yaxis.set_ticks_position('none')
        axes_circ.set_yticklabels("")
        axes_circ.set_axis_off()
        cbar_ax = fig_circ.add_axes([1, 0.061, 0.03, 0.91])
        cb = fig_circ.colorbar(mpl.cm.ScalarMappable(norm=cbar_norm, cmap=cbar_cmap),
                    cax=cbar_ax, orientation='vertical',extend='both')
        cb.set_label(label='$SMD_{mean}$ [$kg/m^2$]',fontsize=self.fs_cb)
        # plt.subplots_adjust(top=1.3)
        plt.tight_layout()
        fig_circ.savefig(self.out_path / f"2bis_sector_groupRef_{self.var}.png", dpi=150, transparent=True, bbox_inches="tight")

    def plot_corr(self):
        regions = ["in_polar","inWw_polar","inLee_polar","out_polar","outWw_polar","outLee_polar"]
        nlocs_list = [[4],[5],[6]] # ,[4,5,6],[5,6]
        param_values = {"HEIGHT":[0.7,0.7,0.7,0.2,0.3,0.4,0.5,0.6,0.7,0.7,0.7,0.7],
                        "AZIMUTH":[0,30,45,0,0,0,0,0,0,0,0,0],
                        "INTERSPACE":[10,10,10,10,10,10,10,10,10,8,10,12],}
        # corr
        method= "pearson"
        # prepare containers 
        smd_corr_locs_dict, smd_corr_param_dict = {}, {}
        for key in param_values.keys(): 
            smd_corr_param_dict[key] = {}
            for reg in regions:
                smd_corr_param_dict[key][reg] = [] 
        # main loop: rebuild df for pca analysis
        for r, region in enumerate(regions):
            for l, nlocs in enumerate(nlocs_list):
                df_pca_cont = []
                for i, nloc in enumerate(nlocs):
                    param_values["Row"] = len(param_values["HEIGHT"])*[i+1]
                    df_pca_temp = pd.DataFrame(param_values)
                    df_quant = self.df.T.set_index("param_value", append=True).loc[self.quant_params,self.var]
                    df_ = df_quant.loc[self.idx[:,region,nloc,:]].reset_index().rename(columns={"level_0":"region", "level_1": "location"})
                    df_pca_temp[self.var] = df_[self.var]
                    df_pca_cont.append(df_pca_temp)

                df_pca = pd.concat(df_pca_cont,axis=0)
                # normalize data
                data_scaled = pd.DataFrame(preprocessing.scale(df_pca),columns = df_pca.columns) 
                # correlation 
                m_corr = data_scaled.corr(method=method)
                
                # fill containers: store smd corrs
                m_corr_proc = m_corr.loc[self.var].dropna().drop(self.var).rename(self.sector_name_dict[region])
                if str(nlocs) not in smd_corr_locs_dict.keys(): smd_corr_locs_dict[str(nlocs)] = []
                smd_corr_locs_dict[str(nlocs)].append(m_corr_proc)
                for key in smd_corr_param_dict.keys():
                    smd_corr_param_dict[key][region].append(m_corr_proc[key])


        # plot locs smd corr
        fig_smd, axes_smd = plt.subplots(1, len(nlocs_list), sharex=False, sharey=True, figsize=(10,4))
        for l, (locs, smd_corrs) in enumerate(smd_corr_locs_dict.items()):
            df_smd_corrs = pd.concat(smd_corrs, axis=1).T
            # df_smd_corrs.reset_index()["sectors"] = regions
            axes_smd[l].set_title(f"Row {l+1}",fontsize=self.fs_titles)
            sns.heatmap(data=df_smd_corrs, cmap="coolwarm", vmin=-1, vmax=1,annot=True, fmt=".2f", cbar_kws={"label":f"{method} correlation"}, ax = axes_smd[l],  cbar=False) 

        fig_smd.supxlabel("Parameter", fontweight="bold",fontsize=self.fs_suplabels)
        fig_smd.supylabel("Sector", fontweight="bold",fontsize=self.fs_suplabels)
        plt.tight_layout()
        fig_smd.savefig(self.out_path / f"4_corr_param_{self.var}.jpeg", dpi=150, bbox_inches="tight")

        # plot params smd corr
        fig_smd, axes_smd = plt.subplots(1, len(param_values.keys())-1, sharex=False, sharey=True, figsize=(10,4))
        for l, (param, smd_corrs) in enumerate(smd_corr_param_dict.items()):
            df_smd_corrs = pd.DataFrame(smd_corrs).T
            # update index
            df_smd_corrs.index = [self.sector_name_dict[sect] for sect in df_smd_corrs.index]
            # update columns
            df_smd_corrs.columns = df_smd_corrs.columns + 1
            # title
            axes_smd[l].set_title(param,fontsize=self.fs_titles)
            # plot
            sns.heatmap(data=df_smd_corrs, cmap="coolwarm", vmin=-1, vmax=1,annot=True, fmt=".2f", cbar_kws={"label":f"{method} correlation"}, ax = axes_smd[l],  cbar=False) 
            
        fig_smd.supxlabel("Row", fontweight="bold", fontsize=self.fs_suplabels)
        fig_smd.supylabel("Sector", fontweight="bold", fontsize=self.fs_suplabels)
        plt.tight_layout()
        fig_smd.savefig(self.out_path / f"4_corr_locs_{self.var}.jpeg", dpi=150, bbox_inches="tight")

    def plot_sectors(self, data, sectors_dict, scale_f=1, cmap_dyn_range = 1, **kwargs):
        # text placement hacks
        text_loc_inc = {"inLee_polar":-0.1,"inWw_polar":0.1,"outLee_polar":1.1,"outWw_polar":-1.1}
        # rescale
        data = scale_f * data
        # cmap
        if "cmap" not in kwargs.keys():
            cmap = "coolwarm"
        else:
            cmap = kwargs["cmap"]
        cmap = mpl.cm.get_cmap(cmap)
        # norm = mpl.colors.Normalize(vmin=-max(abs(data).max()) * cmap_dyn_range, vmax=max(abs(data).max()) * cmap_dyn_range)
        if "cb_vmax" not in kwargs.keys():
            cb_vmax = 0.4
        else:
            cb_vmax = kwargs["cb_vmax"]
        
        norm = mpl.colors.Normalize(vmin=-cb_vmax, vmax=cb_vmax)

        if "ax" not in kwargs.keys():
            fig, ax = plt.subplots(figsize=(20,5))
        else:
            ax = kwargs["ax"]
            fig = ax.get_figure()
        for i, (key_loc, item_loc) in enumerate(sectors_dict.items()):
            for key_sector, item_sector in item_loc.items():
                # erosion (blue) or deposition (red)
                val = data.loc[key_sector].iloc[i]
                color = cmap(norm(val))
                plot_polygon(item_sector,ax=ax, alpha=1, facecolor="w", color="k", linewidth=1, add_points=False)
                plot_polygon(item_sector,ax=ax, alpha=1, color=color, linewidth=0, add_points=False)
                # etiquettes
                ax.text(item_sector.centroid.x + text_loc_inc[key_sector], item_sector.centroid.y, s=f"{round(val,2)}", rotation=90, fontweight="bold", horizontalalignment="center", verticalalignment="center")
                if "title" in kwargs.keys():
                    ax.set_title(kwargs["title"], fontsize=self.fs_titles)
                if "ylabel" in kwargs.keys():
                    ax.set_ylabel(kwargs["ylabel"], fontsize=self.fs_labels)
                if "x_ticks" in kwargs.keys():
                    ax.set_xticks(kwargs["x_ticks"])
                if "x_ticklabels" in kwargs.keys():
                    ax.set_xticklabels(kwargs["x_ticklabels"])
        fig.supxlabel("Row", fontsize=self.fs_suplabels, fontweight="bold")
        
        return norm, cmap

    ###
    # utils 
    ###
    @staticmethod
    def get_sectors(W=1.89,p=0.3,rho=2,d=10,n=5):
        """ returns shapely geometries for plotting
        W = 1.89 # wing size
        p = 0.3 # post size
        rho = 2 # ratio outer / inner 
        d = 10 # interspace
        n = 5 # number of units
        """
        r = W + p/2 # radius
        sectors_dict = {}
        centers = np.arange(0, n*d, d)
        for row, c in enumerate(centers):
            right_mask = shapely.Polygon(((c,rho*r),(c,-rho*r),(c+rho*r,-rho*r),(c+rho*r,rho*r)))
            left_mask = shapely.Polygon(((c,rho*r),(c,-rho*r),(c-rho*r,-rho*r),(c-rho*r,rho*r)))
            inner = shapely.buffer(shapely.Point(c,0),r)
            inner_lee = shapely.difference(inner, left_mask)
            inner_ww = shapely.difference(inner, right_mask) 
            outer_circle = shapely.buffer(shapely.Point(c,0),rho*r)
            outer = shapely.difference(outer_circle,inner)
            outer_lee = shapely.difference(outer, left_mask)
            outer_ww = shapely.difference(outer, right_mask) 
            sectors_dict[row] = {"inLee_polar":inner_lee,"inWw_polar":inner_ww,"outLee_polar":outer_lee, "outWw_polar":outer_ww}
            # sectors_dict[row] = {"in_polar":inner,"inLee_polar":inner_lee,"inWw_polar":inner_ww,"out_polar":outer,"outLee_polar":outer_lee, "outWw_polar":outer_ww}
        return sectors_dict


        