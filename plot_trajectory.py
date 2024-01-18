def plot_geo(lat1,lon1):
     import cartopy.feature as cfeature
     import cartopy.crs as ccrs
     import matplotlib.pyplot as plt
     ax = plt.axes(projection=ccrs.PlateCarree())
     # ax.set_extent([125, 150, 35, 63])
#     ax.stock_img()
#
#     ax.add_feature(cfeature.LAND)  # If I comment this => all ok, but I need
#     ax.add_feature(cfeature.LAKES)
#     ax.add_feature(cfeature.RIVERS)
     ax.coastlines()
#     ax.coastlines()
     ax.scatter(lon1, lat1, transform=ccrs.PlateCarree())
     return