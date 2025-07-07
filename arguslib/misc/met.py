from csat2.ECMWF import download

def download_era5_winds(dtime):
    levels = [
        150,175,200, 225,250, 300, 350
    ]
    
    for level in levels:
        print(f"downloading {level}hPa...")
        try:
            download(dtime.year,dtime.month, 'U-wind-component',f"{level}hPa", "0.25grid", [dtime.day])
            download(dtime.year,dtime.month, 'V-wind-component',f"{level}hPa", "0.25grid", [dtime.day])
        except PermissionError as e:
            print(f"Failed at {level}: {e}")
            raise e