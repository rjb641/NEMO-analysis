"""Useful functions for creating maps of the Labrador Sea.

        Functions:
                xrLSminmax
                LSmap
"""

import matplotlib.pyplot as plt
import matplotlib.path as mpath
import cartopy.crs as ccrs
import cartopy.feature as feature

def xrLSminmax(xrData,lats,lons):
        """ Finds min and max values in xr data containing the Labrador Sea.
        These are principally for later use in defining a colourbar range.
        
                Parameters:
                        xrData: DataArray containing data to be plotted
                        lats:   2D array of cell latitudes from .nc grid file
                        lons:   2D array of cell longitudes from .nc grid file
                        
                Returns:
                        Tuple containing min and max values in the DataArray
        """

        #define Labrador Sea region
        cond = ((lats > 40) & (lats < 85) & (lons < -30) & (lons > -80))

        #find min and max values within Labrador Sea region
        max = xrData.where(cond).max(skipna=True).to_numpy()
        min = xrData.where(cond).min(skipna=True).to_numpy()

        return min, max

def LSmap(xrData,lons,lats,minmax,CBlabel,title,fileName):
        """ Saves PNG map of Labrador Sea with data from an xarray DataArray. 
        If there are multiple days (i.e. len(time_counter)>1), then multiple 
        PNGs are created and named according to the date.
        
                Parameters:
                        xrData:         DataArray containing data to be plotted
                                        (Should contain time_counter dim of len>=1)
                        lats:           2D array of cell latitudes from .nc grid file
                        lons:           2D array of cell longitudes from .nc grid file
                        minmax:         Tuple with min and max values in the DataArray
                                        (From xrLSminmax function)
                        CBlabel:        String for the colourbar label
                        title:          String for the title (to display in the PNG)
                                        (Will have date appended to it)
                                        E.g., 'Sea surface height (EPM151)' + ' ' + date
                        fileName:       String for the name of the PNG(s)
                                        E.g., 'LS_convective_energy_EPM151_' + date + '.png'
                        
                Returns:
                        None
        """

        #shapefile of land with 1:50,000,000 scale
        land_50m = feature.NaturalEarthFeature('physical', 'land', '50m',edgecolor='black', facecolor='gray')

        #defining the projection, note that standard parallels are the parallels of correct scale
        projection = ccrs.AlbersEqualArea(central_longitude=-55, central_latitude=50,standard_parallels=(40,85))

        #create figure (using the specified projection)
        ax = plt.subplot(1, 1, 1, projection=projection)

        #define map dimensions (using Plate Carree coordinate system)
        ax.set_extent([-80, -30, 40, 85], crs=ccrs.PlateCarree())

        #add land to map
        ax.add_feature(land_50m, color=[0.8, 0.8, 0.8])

        #add coast lines 
        ax.coastlines(resolution='50m')

        #updating boundary shape
        vertices = [(lon, 40) for lon in range(-80, -30, 1)] + [(lon, 85) for lon in range(-30, -80, -1)]
        boundary = mpath.Path(vertices)
        ax.set_boundary(boundary, transform=ccrs.PlateCarree())

        #unpacking tuple
        min, max = minmax

        #ticks
        gls = ax.gridlines(draw_labels=True, dms=True)
        gls.top_labels=False #suppress top labels
        gls.xlabel_style = {'rotation': 0}
        gls.ylabel_style = {'rotation': 0}

        #get list of times
        times = xrData["time_counter"].dt.strftime("%Y-%m-%d").to_numpy()

        if times.size == 1:
                #plotting data
                p1 = ax.pcolormesh(lons, lats, xrData, transform=ccrs.PlateCarree(), cmap='gist_ncar', vmin=min, vmax=max)

                #colourbar 
                ax_cb = plt.axes([0.8, 0.25, 0.015, 0.5])
                cb = plt.colorbar(p1,cax=ax_cb, orientation='vertical')
                cb.ax.set_ylabel(CBlabel)

                #title
                date = times
                ax.set_title(title + ' ' + date)#,fontdict={'fontsize': 12})

                #save figure
                plt.savefig(fileName + date + '.png',dpi=900, bbox_inches="tight")

        else:
                for i,t in enumerate(times): 
                        #plotting data
                        p1 = ax.pcolormesh(lons, lats, xrData[i,:,:], transform=ccrs.PlateCarree(), cmap='gist_ncar', vmin=min, vmax=max)

                        #colourbar 
                        ax_cb = plt.axes([0.8, 0.25, 0.015, 0.5])
                        cb = plt.colorbar(p1,cax=ax_cb, orientation='vertical')
                        cb.ax.set_ylabel(CBlabel) 

                        #title
                        date = times[i]
                        ax.set_title(title + ' ' + date)#,fontdict={'fontsize': 12})

                        #save figure
                        plt.savefig(fileName + date + '.png',dpi=900, bbox_inches="tight")


        

