# -*- coding: utf-8 -*-
"""
Created on Wed Mar 22 10:58:30 2023

@author: Rasmus
"""
G10_flux = np.array([evt_data['GALEX_FUV'][1:10],evt_data['GALEX_NUV'][1:10],evt_data['CFHT_u_FLUX'][1:10],
                     evt_data['CFHT_ustar_FLUX'][1:10],evt_data['IB427'][1:10],evt_data['IB464'][1:10],
                     evt_data['HSC_g'][1:10],evt_data['IA484'][1:10],evt_data['IB505'][1:10],
                     evt_data['IA527'][1:10],evt_data['IB574'][1:10],evt_data['HSC_r'][1:10],
                     evt_data['IA624'][1:10],evt_data['IA679'][1:10],evt_data['IB709'][1:10],
                     evt_data['NB711'][1:10],evt_data['IA738'][1:10],evt_data['IA767'][1:10],
                     evt_data['HSC_i'][1:10],evt_data['NB816'][1:10],evt_data['IB827'][1:10],
                     evt_data['F814W'][1:10],evt_data['HSC_z'][1:10],evt_data['HSC_y'][1:10],
                     evt_data['UVISTA_Y'][1:10],evt_data['UVISTA_J'][1:10],evt_data['UVISTA_H'][1:10],
                     evt_data['UVISTA_Ks'][1:10],evt_data['IRAC_CH1'][1:10],evt_data['IRAC_CH2'][1:10],
                     evt_data['SPLASH_CH2'][1:10],evt_data['SPLASH_CH3']])