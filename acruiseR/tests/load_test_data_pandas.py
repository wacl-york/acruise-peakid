import datetime
import pandas as pd
from acruisepy import peakid

df = pd.read_csv(
    "../../../data/C258_FGGA_FAAM.txt",
    names=["conc", "a", "b"],
    parse_dates=True,
    header=1,
    date_parser=lambda x: datetime.datetime.strptime(x, "%d/%m/%Y %H:%M:%S.%f"),
)
start_time = datetime.datetime(2021, 10, 4, 9, 40, 0)
end_time = datetime.datetime(2021, 10, 4, 13, 30, 0)
df = df.loc[start_time:end_time]
