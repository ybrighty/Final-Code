import os
import numpy as np
import calc
import structure as st


class Ag(st.Layer):
    '''
    Creates a Layer of silver. Minimum user input is the thickness.
    '''
    def __init__(self, thickness, unit='nm'):
        wl = _Ag_wl / 1e3 if unit == 'micron' else _Ag_wl
        n = _Ag_nk
        label = 'Ag'
        st.Layer.__init__(self, wl, n, thickness, label=label, unit=unit)


class Al(st.Layer):
    '''
    Creates a Layer of aluminum. Minimum user input is the thickness.
    '''
    def __init__(self, thickness, unit='nm'):
        wl = _Al_wl / 1e3 if unit == 'micron' else _Al_wl
        n = _Al_nk
        label = 'Al'
        st.Layer.__init__(self, wl, n, thickness, label=label, unit=unit)


class Cr(st.Layer):
    '''
    Creates a Layer of aluminum. Minimum user input is the thickness.
    '''
    def __init__(self, thickness, unit='nm'):
        wl = _Cr_wl / 1e3 if unit == 'micron' else _Cr_wl
        n = _Cr_nk
        label = 'Cr'
        st.Layer.__init__(self, wl, n, thickness, label=label, unit=unit)


class Mo(st.Layer):
    '''
    Creates a Layer of aluminum. Minimum user input is the thickness.
    '''
    def __init__(self, thickness, unit='nm'):
        wl = _Mo_wl / 1e3 if unit == 'micron' else _Mo_wl
        n = _Mo_nk
        label = 'Mo'
        st.Layer.__init__(self, wl, n, thickness, label=label, unit=unit)     


      
class AlN(st.Layer):
    '''
    Creates a Layer of aluminum nitride. Minimum user input is the thickness.
    '''
    def __init__(self, thickness, unit='nm'):
        wl = _AlN_wl / 1e3 if unit == 'micron' else _AlN_wl
        n = _AlN_nk
        label = 'AlN'
        st.Layer.__init__(self, wl, n, thickness, label=label, unit=unit)

class AlNH(st.Layer):
    '''
    Creates a Layer of hydrogenated aluminum nitride.
    '''
    def __init__(self, thickness, unit='nm'):
        wl = _AlNH_wl / 1e3 if unit == 'micron' else _AlNH_wl
        n = _AlNH_nk
        label = 'AlNH'
        st.Layer.__init__(self, wl, n, thickness, label=label, unit=unit)
 

class BK7(st.Layer):
    '''
    Creates a Layer of BK7 glass. Minimum user input is the thickness.
    '''
    def __init__(self, thickness, unit='nm'):
        wl = _BK7_wl / 1e3 if unit == 'micron' else _BK7_wl
        n = _BK7_nk
        label = 'BK7'
        st.Layer.__init__(self, wl, n, thickness, label=label, unit=unit)


class DLC3W(st.Layer):
    '''
    Creates a Layer of DLC3W. Minimum user input is the thickness.
    extended=False allows user to only use valid experimental data in
    the visible range. nk value from 800nm to 2600nm is exponentially
    extrapolated. See associated image files in the nklib folder.
    '''
    def __init__(self, thickness, unit='nm', extended=True):
        wl = _DLC3W_ext_wl if extended else _DLC3W_wl
        wl = wl / 1e3 if unit == 'micron' else wl
        n = _DLC3W_ext_nk if extended else _DLC3W_nk
        label = 'DLC3W'
        st.Layer.__init__(self, wl, n, thickness, label=label, unit=unit)


class DLC5W(st.Layer):
    '''
    Creates a Layer of DLC5W. Minimum user input is the thickness.
    extended=False allows user to only use valid experimental data in
    the visible range. nk value from 800nm to 2600nm is exponentially
    extrapolated. See associated image files in the nklib folder.
    '''
    def __init__(self, thickness, unit='nm', extended=True):
        wl = _DLC5W_ext_wl if extended else _DLC5W_wl
        wl = wl / 1e3 if unit == 'micron' else wl
        n = _DLC5W_ext_nk if extended else _DLC5W_nk
        label = 'DLC5W'
        st.Layer.__init__(self, wl, n, thickness, label=label, unit=unit)


class DLC10W(st.Layer):
    '''
    Creates a Layer of DLC10W. Minimum user input is the thickness.
    extended=False allows user to only use valid experimental data in
    the visible range. nk value from 800nm to 2600nm is exponentially
    extrapolated. See associated image files in the nklib folder.
    '''
    def __init__(self, thickness, unit='nm', extended=True):
        wl = _DLC10W_ext_wl if extended else _DLC10W_wl
        wl = wl / 1e3 if unit == 'micron' else wl
        n = _DLC10W_ext_nk if extended else _DLC10W_nk
        label = 'DLC10W'
        st.Layer.__init__(self, wl, n, thickness, label=label, unit=unit)


class DLC15W(st.Layer):
    '''
    Creates a Layer of DLC15W. Minimum user input is the thickness.
    extended=False allows user to only use valid experimental data in
    the visible range. nk value from 800nm to 2600nm is exponentially
    extrapolated. See associated image files in the nklib folder.
    '''
    def __init__(self, thickness, unit='nm', extended=True):
        wl = _DLC15W_ext_wl if extended else _DLC15W_wl
        wl = wl / 1e3 if unit == 'micron' else wl
        n = _DLC15W_ext_nk if extended else _DLC15W_nk
        label = 'DLC15W'
        st.Layer.__init__(self, wl, n, thickness, label=label, unit=unit)


class DLC20W(st.Layer):
    '''
    Creates a Layer of DLC20W. Minimum user input is the thickness.
    extended=False allows user to only use valid experimental data in
    the visible range. nk value from 800nm to 2600nm is exponentially
    extrapolated. See associated image files in the nklib folder.
    '''
    def __init__(self, thickness, unit='nm', extended=True):
        wl = _DLC20W_ext_wl if extended else _DLC20W_wl
        wl = wl / 1e3 if unit == 'micron' else wl
        n = _DLC20W_ext_nk if extended else _DLC20W_nk
        label = 'DLC20W'
        st.Layer.__init__(self, wl, n, thickness, label=label, unit=unit)


class DLC40W(st.Layer):
    '''
    Creates a Layer of DLC40W. Minimum user input is the thickness.
    extended=False allows user to only use valid experimental data in
    the visible range. nk value from 800nm to 2600nm is exponentially
    extrapolated. See associated image files in the nklib folder.
    '''
    def __init__(self, thickness, unit='nm', extended=True):
        wl = _DLC40W_ext_wl if extended else _DLC40W_wl
        wl = wl / 1e3 if unit == 'micron' else wl
        n = _DLC40W_ext_nk if extended else _DLC40W_nk
        label = 'DLC40W'
        st.Layer.__init__(self, wl, n, thickness, label=label, unit=unit)


class DLC60W(st.Layer):
    '''
    Creates a Layer of DLC60W. Minimum user input is the thickness.
    extended=False allows user to only use valid experimental data in
    the visible range. nk value from 800nm to 2600nm is exponentially
    extrapolated. See associated image files in the nklib folder.
    '''
    def __init__(self, thickness, unit='nm', extended=True):
        wl = _DLC40W_ext_wl if extended else _DLC40W_wl
        wl = wl / 1e3 if unit == 'micron' else wl
        n = _DLC40W_ext_nk if extended else _DLC40W_nk
        label = 'DLC60W'
        st.Layer.__init__(self, wl, n, thickness, label=label, unit=unit)


class DLC80WA(st.Layer):
    '''
    Creates a Layer of DLC80WA. Minimum user input is the thickness.
    '''
    def __init__(self, thickness, unit='nm'):
        wl = _DLC80WA_wl / 1e3 if unit == 'micron' else _DLC80WA_wl
        n = _DLC80WA_nk
        label = 'DLC80WA'
        st.Layer.__init__(self, wl, n, thickness, label=label, unit=unit)
        

class ITO(st.Layer):
    '''
    Creates a Layer of ITO. Minimum user input is the thickness.
    '''
    def __init__(self, thickness, unit='nm'):
        wl = _ITO_wl / 1e3 if unit == 'micron' else _ITO_wl
        n = _ITO_nk
        label = 'ITO'
        st.Layer.__init__(self, wl, n, thickness, label=label, unit=unit)


class PDMS(st.Layer):
    '''
    Creates a Layer of PDMS. Minimum user input is the thickness.
    '''
    def __init__(self, thickness, unit='nm'):
        wl = _PDMS_wl / 1e3 if unit == 'micron' else _PDMS_wl
        n = _PDMS_nk
        label = 'PDMS'
        st.Layer.__init__(self, wl, n, thickness, label=label, unit=unit)


class TiO2(st.Layer):
    '''
    Creates a Layer of TiO2. Minimum user input is the thickness.
    '''
    def __init__(self, thickness, unit='nm'):
        wl = _TiO2_wl / 1e3 if unit == 'micron' else _TiO2_wl
        n = _TiO2_nk
        label = 'TiO2'
        st.Layer.__init__(self, wl, n, thickness, label=label, unit=unit)

class SiO2(st.Layer):
    '''
    Creates a Layer of SiO2. Minimum user input is the thickness.
    '''
    def __init__(self, thickness, unit='nm'):
        wl = _SiO2_wl / 1e3 if unit == 'micron' else _SiO2_wl
        n = _SiO2_nk
        label = 'SiO2'
        st.Layer.__init__(self, wl, n, thickness, label=label, unit=unit)        

class Al2O3(st.Layer):
    '''
    Creates a Layer of Al2O3.Minimum user input is the thickness.
    '''
    def __init__(self, thickness, unit='nm'):
        wl = _Al2O3_wl / 1e3 if unit == 'micron' else _Al2O3_wl
        n = _Al2O3_nk
        label = 'Al2O3'
        st.Layer.__init__(self, wl, n, thickness, label=label, unit=unit)


class ZnO(st.Layer):
    '''
    Creates a Layer of ZnO. Minimum user input is the thickness.
    '''
    def __init__(self, thickness, unit='nm'):
        wl = _ZnO_wl / 1e3 if unit == 'micron' else _ZnO_wl
        n = _ZnO_nk
        label = 'ZnO'
        st.Layer.__init__(self, wl, n, thickness, label=label, unit=unit)


# Import data
os.chdir(os.getcwd() + "\\opt_sim\\nklib_data") # possibly generalize

'''bulk silver
_Ag_ndata = np.loadtxt("Ag_n.txt", skiprows=1)
_Ag_kdata = np.loadtxt("Ag_k.txt", skiprows=1)
_Ag_wl = _Ag_ndata[:,0] * 1000
_Ag_nk = _Ag_ndata[:,1] - 1j * _Ag_kdata[:,1]
'''

# 102A silver
'''
_Ag_data = np.loadtxt("Ag 102A nk.txt", skiprows=1)
_Ag_wl = _Ag_data[:,0] * 1000
_Ag_nk = _Ag_data[:,1] - 1j * _Ag_data[:,2]
'''

# 87A silver (up to 800 nm, then used 102A silver's 800 - 2480)

_Ag_data = np.loadtxt("Ag 87 + Ag102A to full (with adjustments).txt", skiprows=1)
_Ag_wl = _Ag_data[:,0] * 1000
_Ag_nk = _Ag_data[:,1] - 1j * _Ag_data[:,2]

_Al_ndata = np.loadtxt("Al_n_Rakic1998.txt", skiprows=1)
_Al_kdata = np.loadtxt("Al_k_Rakic1998.txt", skiprows=1)
_Al_wl = _Al_ndata[:,0] * 1000
_Al_nk = _Al_ndata[:,1] - 1j * _Al_kdata[:,1]

_Cr_ndata = np.loadtxt("Chromium_n.txt", skiprows=1)
_Cr_kdata = np.loadtxt("Chromium_k.txt", skiprows=1)
_Cr_wl = _Cr_ndata[:,0] * 1000
_Cr_nk = _Cr_ndata[:,1] - 1j * _Cr_kdata[:,1]

_Mo_nkdata = np.loadtxt("Mo_nk.txt", skiprows=1)
_Mo_wl = _Mo_nkdata[:,0]
_Mo_nk = _Mo_nkdata[:,1] - 1j * _Mo_nkdata[:,2]

_AlN_nkdata = np.loadtxt("aln_ellipsometry_sai_nk.txt", skiprows=1)
_AlN_wl = _AlN_nkdata[:,0] * 1000
_AlN_nk = _AlN_nkdata[:,1] - 1j * _AlN_nkdata[:,2]

_AlNH_nkdata = np.loadtxt("alN-H-nk.txt")
_AlNH_wl = _AlNH_nkdata[:,0] * 1000
_AlNH_nk = _AlNH_nkdata[:,1] - 1j * _AlNH_nkdata[:,2]

_BK7_ndata = np.loadtxt("BK7_n.txt", skiprows=1)
_BK7_kdata = np.loadtxt("BK7_k.txt", skiprows=1)
_BK7_wl = _BK7_ndata[:,0] * 1000
# BK7 doesn't have the same number of data points for n and k
_wl_temp = _BK7_kdata[:,0] * 1000
_k_temp = _BK7_kdata[:,1]
_BK7_k = calc.interpolate(_wl_temp, _k_temp, _BK7_wl)
_BK7_nk = _BK7_ndata[:,1] - 1j * _BK7_k

_DLC_nkdata = [np.loadtxt("DLC{}W_nk.txt".format(n), skiprows=1) for \
               n in [3, 5, 10, 15, 20, 40, 60]]
_DLC_wl = [_DLC_nkdata[i][:,0] for i in range(7)]
_DLC_nk = [_DLC_nkdata[i][:,1] - 1j * _DLC_nkdata[i][:,2] for i in range(7)]
_DLC3W_wl, _DLC5W_wl, _DLC10W_wl, _DLC15W_wl, _DLC20W_wl, _DLC40W_wl, _DLC60W_wl = _DLC_wl
_DLC3W_nk, _DLC5W_nk, _DLC10W_nk, _DLC15W_nk, _DLC20W_nk, _DLC40W_nk, _DLC60W_nk = _DLC_nk

_DLC_ext_nkdata = [np.loadtxt("DLC{}W_extended_nk.txt".format(n), skiprows=1) \
                   for n in [3, 5, 10, 15, 20, 40, 60]]
_DLC_ext_wl = [_DLC_ext_nkdata[i][:,0] for i in range(7)]
_DLC_ext_nk = [_DLC_ext_nkdata[i][:,1] - 1j * _DLC_ext_nkdata[i][:,2] \
               for i in range(7)]
_DLC3W_ext_wl, _DLC5W_ext_wl, _DLC10W_ext_wl, _DLC15W_ext_wl, _DLC20W_ext_wl, _DLC40W_ext_wl, _DLC60W_ext_wl = _DLC_ext_wl
_DLC3W_ext_nk, _DLC5W_ext_nk, _DLC10W_ext_nk, _DLC15W_ext_nk, _DLC20W_ext_nk, _DLC40W_ext_nk, _DLC60W_ext_nk = _DLC_ext_nk
DLC_list = [DLC3W, DLC5W, DLC10W, DLC15W, DLC20W, DLC40W, DLC60W]

_DLC80WA_nkdata = np.loadtxt("DLC80WAnode_nk.txt", skiprows=1)
_DLC80WA_wl = _DLC80WA_nkdata[:,0] * 1000
_DLC80WA_nk = _DLC80WA_nkdata[:,1] - 1j * _DLC80WA_nkdata[:,2]

_ITO_ndata = np.loadtxt("ITO_n_Konig.txt", skiprows=1)
_ITO_kdata = np.loadtxt("ITO_k_Konig.txt", skiprows=1)
_ITO_wl = _ITO_ndata[:,0] * 1000
_ITO_nk = _ITO_ndata[:,1] - 1j * _ITO_kdata[:,1]

_PDMS_nkdata = np.loadtxt("PDMS Refractive Index.csv",
                          skiprows=1, delimiter=',')
_PDMS_wl = _PDMS_nkdata[:,0] * 1000
_PDMS_nk = _PDMS_nkdata[:,1] - 1j * _PDMS_nkdata[:,2]

_TiO2_nkdata = np.loadtxt("TiO2_nk.txt", skiprows=1)
_TiO2_wl = _TiO2_nkdata[:,0] * 1000
_TiO2_nk = _TiO2_nkdata[:,1] - 1j * _TiO2_nkdata[:,2]

_Al2O3_ndata = np.loadtxt("2015Al2O3_n.txt", skiprows=1)
_Al2O3_kdata = np.loadtxt("2015Al2O3_k.txt", skiprows=1)
_Al2O3_wl = _Al2O3_ndata[:,0] * 1000
_Al2O3_nk = _Al2O3_ndata[:,1] - 1j * _Al2O3_kdata[:,1]

_SiO2_ndata = np.loadtxt("SiO2_n.txt", skiprows=1)
_SiO2_kdata = np.loadtxt("SiO2_k.txt", skiprows=1)
_SiO2_wl = _SiO2_ndata[:,0] * 1000
_SiO2_nk = _SiO2_ndata[:,1] - 1j * _SiO2_kdata[:,1]

_ZnO_nkdata = np.loadtxt("ZnO_nk.txt", skiprows=1)
_ZnO_wl = _ZnO_nkdata[:,0] * 1000
_ZnO_nk = _ZnO_nkdata[:,1] - 1j * _ZnO_nkdata[:,2]

os.chdir("..")
