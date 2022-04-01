# coding=utf-8
import pandas as pd
import numpy as np
from itertools import product
import os

os.chdir('./original')

names = {
    'region': 'region',
    'grödgrupp': 'crop',
    'gödselslag': 'strategy',
    'Areal, hektar 2012/2013': 'a',
    'Tillförd mängd kväve från handelsgödsel, kg per hektar 2012/2013': 'f'
}

vals = {
    'crop': {
        'samtliga åkergrödor': 'all_crops',
        'spannmål': 'cereals',
        'slåttervall': 'ley_silage',
        'betesvall': 'ley_pasture',
        'samtliga åkergrödor utom spannmål och slåtter- och betesvall': 'other'
    },
    'strategy': {
        'enbart handelsgödsel': 'syn',
        'både handels- och stallgödsel': 'syn_and_man',
        'enbart stallgödsel': 'man',
        'varken handels- eller stallgödsel': 'none'
    },
    'region': {
        "PO1 Götalands södra slättbygder": "PO1",
        "PO2 Götalands mellanbygder": "PO2",
        "PO3 Götalands norra slättbygder": "PO3",
        "PO4 Svealands slättbygder": "PO4",
        "PO5 Götalands skogsbygder": "PO5",
        "PO6 Mellersta Sveriges skogsbygder": "PO6",
        "PO7 Nedre Norrland": "PO7",
        "PO8 Övre Norrland": "PO8",
    }
}

mf_nan_replacement = {'m': {'..': float('nan')}, 'f': {'..': float('nan')}}


datasets = {}

datasets['mRC'] = (
    pd.read_csv('mRC_orig.csv', skiprows=[0,1], encoding='cp1252')
    .rename(columns=names)
    .rename(columns={'2012/2013': 'm'})
    .replace(vals)
    .replace(mf_nan_replacement)
    .set_index(['region', 'crop'])
    .drop('all_crops', level='crop')
    ['m']
    )

datasets['mR'] = (
    pd.read_csv('mRC_orig.csv', skiprows=[0,1], encoding='cp1252')
    .rename(columns=names)
    .rename(columns={'2012/2013': 'm'})
    .replace(vals)
    .replace(mf_nan_replacement)
    .set_index(['region', 'crop'])
    .xs('all_crops', level='crop')
    ['m']
    )

datasets['mC'] = (
    pd.read_csv('mC_orig.csv', skiprows=[0,1], encoding='cp1252')
    .rename(columns=names)
    .rename(columns={'2012/2013': 'm'})
    .replace(vals)
    .set_index(['crop'])
    .drop('all_crops')
    ['m']
    )

afCS = (
    pd.read_csv('afCS_orig.csv', skiprows=[0,1], encoding='cp1252')
    .rename(columns=names)
    .replace(vals)
    .set_index(['crop', 'strategy'])
    )
afCS['f'] = afCS['f'].astype(float)

mCS = ((afCS['f'] * afCS['a']).drop('all_crops', level='crop') * 1e-3).round()
mCS.name = 'm'
datasets['mCS'] = mCS


afRCS = (
    pd.read_csv('afRCS_orig.csv', skiprows=[0,1], encoding='cp1252')
    .rename(columns=names)
    .replace(vals)
    .replace(mf_nan_replacement)
    .set_index(['region', 'crop', 'strategy'])
    )
afRCS['f'] = afRCS['f'].astype(float)

datasets['aRCS'] = (
    afRCS['a']
    .drop('all_crops', level='crop')
    )

mRCS = (
    (afRCS['f'] * afRCS['a'])
    .drop('all_crops', level='crop')
    * 1e-3
    ).round()
mRCS.name = 'm'
datasets['mRCS'] = mRCS

os.chdir('./..')

for name, ds in datasets.items():
    print(name)
    ds.astype(float).to_csv('{}.csv'.format(name), header=True)
