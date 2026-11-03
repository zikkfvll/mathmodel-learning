import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def convert(src, out):
    with np.load(src, allow_pickle=False) as d:
        path = pd.DataFrame({
            'latitude': d['lat'],
            'longitude': d['lon'],
            'elevation_m': d['elev'],
            'author_chainage_m': d['dist_m'],
            'segment_azimuth_deg': d['seg_azimuth_deg'],
        })
        index = int(d['i_confluence'])
        assert np.all(np.diff(path.author_chainage_m) >= 0)
        path['author_confluence_node'] = False
        path.loc[index, 'author_confluence_node'] = True
    path.to_csv(out, index=False)


def sol():
    p = argparse.ArgumentParser()
    p.add_argument('--input', type=Path, default=Path(__file__).resolve().parent / 'route_data/path.npz')
    p.add_argument('--output', type=Path, default=Path('route_new.csv'))
    a = p.parse_args()
    if a.output.exists():
        raise FileExistsError(a.output)
    a.output.parent.mkdir(parents=True, exist_ok=True)
    convert(a.input, a.output)


if __name__ == '__main__':
    sol()
