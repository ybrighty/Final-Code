import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import numpy as np
import structure as st
import calc
import os

def TR(structure_list, min_wl=None, max_wl=None, curves='TR', unit='nm',
       show_solar=True, show_photopic=True, legend=True, new_figure=True):
    '''
    Plot the reflectance and transmittance curve of given structures in the
    window. If given, min_wl and max_wl sets the wavelength range in the unit
    the LayerStructure instances are initialized in. Default is nm.
    Inputting 'R' or 'T' for the curves argument plots only the reflectance
    or transmittance, respectively. Default is both.
    '''
    # Alternate call signature for convenience
    if isinstance(structure_list, st.MultiLayer):
        structure_list = [structure_list]
        
    if not unit in ['nm', 'micron']:
            raise DataFormatException(
                "unit should be given as 'nm' or 'micron'")

    for structure in structure_list:
        if not structure._TR_calculated:
            structure.calculate_TR()

    if new_figure:
        plt.figure()

    for structure in structure_list:
        wl = structure.wl * 1e3 if unit == 'nm' and structure.unit == 'micron' else\
             structure.wl / 1e3 if unit == 'micron' and structure.unit == 'nm' else\
             structure.wl
        T = structure.T
        R = structure.R
        if 'T' in curves:
            plt.plot(wl, T, label=structure.label + ' T')
        if 'R' in curves:
            plt.plot(wl, R, label=structure.label + ' R')
        #if 'A' in curves:
         #   A = 1 - R - T
          #  plt.plot(wl, A, label=structure.label + ' A')
        A = 1 - R - T
        plt.plot(wl, A, label=structure.label + ' A') 
    if show_solar:
        AM1p5_data = np.loadtxt("plot support files\\ASTMG173.txt", skiprows=2)
        solar_wl = AM1p5_data[:,0]
        solar_intens = AM1p5_data[:,3]
        solar_intens /= max(solar_intens)
        plt.plot(solar_wl, solar_intens, label="AM1.5", color=(0.55,0.55,0.55))
    if show_photopic:
        photopic_data = np.loadtxt("plot support files\\Photopic_luminosity_function.txt", skiprows=1)
        photopic_wl = photopic_data[:,0]
        photopic_intens = photopic_data[:,1]
        plt.plot(photopic_wl, photopic_intens, label="photopic", color=(1,0,0))
    plt.title('Reflectance and Transmittance plots')
    plt.xlabel('Wavelength ({})'.format(unit))
    if min_wl:
        plt.xlim(xmin=min_wl)
    if max_wl:
        plt.xlim(xmax=max_wl)
    plt.ylim(0, 1)
    if legend:
        plt.legend(loc=0)

def nk(layer_list, min_wl=None, max_wl=None, curves='nk',sep=False, unit='nm'):
    '''
    Plot the complex refractive index of the given layers in the same window.
    If given, min_wl and max_wl sets the wavelength range in the unit
    the Layer instances are initialized in. Default is nm.
    '''
    # Alternate call signature for convenience
    if isinstance(layer_list, st.Layer):
        layer_list = [layer_list]
    if isinstance(layer_list, st.MultiLayer):
        layer_list = layer_list._layers_list
        
    if not unit in ['nm', 'micron']:
            raise DataFormatException(
                "unit should be given as 'nm' or 'micron'")
    unit = layer_list[0].unit
            
    # Reduce redundant code merging everything below
    if curves == "nk" and sep:
        material_list = []
        plt.figure()
        plt.title("Real refractive index plot")
        plt.xlabel("Wavelength ({})".format(unit))
        if min_wl:
            plt.xlim(xmin=min_wl)
        if max_wl:
            plt.xlim(xmax=max_wl)
        for layer in layer_list:
            if layer.label in material_list:
                continue
            material_list.append(layer.label)
            wl = layer.wl * 1e3 if unit == 'nm' and layer.unit == 'micron' else\
                 layer.wl / 1e3 if unit == 'micron' and layer.unit == 'nm' else\
                 layer.wl
            n = layer.n.real
            plt.plot(wl, n, label=layer.label + ' n')
        plt.legend(loc=0)

        material_list = []
        plt.figure()
        plt.title("Imaginary refractive index plot")
        plt.xlabel("Wavelength ({})".format(unit))
        if min_wl:
            plt.xlim(xmin=min_wl)
        if max_wl:
            plt.xlim(xmax=max_wl)
        for layer in layer_list:
            if layer.label in material_list:
                continue
            material_list.append(layer.label)
            wl = layer.wl * 1e3 if unit == 'nm' and layer.unit == 'micron' else\
                 layer.wl / 1e3 if unit == 'micron' and layer.unit == 'nm' else\
                 layer.wl
            k = -layer.n.imag
            plt.plot(wl, k, label=layer.label + ' k')
        plt.legend(loc=0)
        return None

    material_list = []
    plt.figure()
    for layer in layer_list:
        if layer.label in material_list:
            continue
        material_list.append(layer.label)
        wl = layer.wl * 1e3 if unit == 'nm' and layer.unit == 'micron' else\
             layer.wl / 1e3 if unit == 'micron' and layer.unit == 'nm' else\
             layer.wl
        n = layer.n.real
        k = -layer.n.imag
        if 'n' in curves:
            plt.plot(wl, n, label=layer.label + ' n')
        if 'k' in curves:
            plt.plot(wl, k, label=layer.label + ' k')
    plt.title('Refractive index plots')
    plt.xlabel('Wavelength ({})'.format(unit))
    if min_wl:
        plt.xlim(xmin=min_wl)
    if max_wl:
        plt.xlim(xmax=max_wl)
    plt.legend(loc=0)
    return None

def view(structure, original=True):
    '''
    Shows 4 pictures.
    '''
    if not structure._color_calculated:
        structure.calculate_color()
        
    os.chdir("plot support files")
    T_image = mpimg.imread("test_outdoor.png")
    R_image = mpimg.imread("test_indoor.png")
    os.chdir("..")
    
    T_filter = np.array(structure.T_color, float) / 255
    R_filter = np.array(structure.R_color, float) / 255

    T_image_after = T_image * T_filter
    R_image_after = R_image * R_filter
    overlay = T_image_after + R_image_after

    plt.figure()
    plt.axis("off")
    plt.imshow(overlay)

    if original:
        plt.figure()
        plt.axis("off")
        plt.imshow(T_image)
        plt.figure()
        plt.axis("off")
        plt.imshow(R_image)

##    import Image
##    T = T_image_after * 255
##    T = T.astype(int)
##    R = R_image_after * 255
##    R = R.astype(int)
##    I = overlay * 255
##    I = I.astype(int)
##    im_T = Image.fromarray(T_image_after, 'RGB')
##    im_R = Image.fromarray(R_image_after, 'RGB')
##    im = Image.fromarray(overlay, 'RGB')
##    os.chdir("..")
##    im_T.save("outdoors.png")
##    im_R.save("indoors.png")
##    im.save("window_view.png")
##    mpimg.imsave("outdoors.png", T_image)
##    mpimg.imsave("indoors.png", R_image)
##    mpimg.imsave("window_view.png", overlay)
##    os.chdir("opt_sim")

##    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col', sharey='row')
##    ax1.axis("off")
##    ax1.imshow(T_image)
##    ax2.axis("off")
##    ax2.imshow(T_image * T_filter)
##    ax3.axis("off")
##    ax3.imshow(R_image)
##    ax4.axis("off")
##    ax4.imshow(R_image * R_filter)
##    plt.tight_layout()
    
##    plt.subplot(2, 2, 1)
##    plt.axis("off")
##    plt.imshow(T_image)
##
##    plt.subplot(2, 2, 2)
##    plt.axis("off")
##    plt.imshow(T_image * T_filter)
##
##    plt.subplot(2, 2, 3)
##    plt.axis("off")
##    plt.imshow(R_image)
##
##    plt.subplot(2, 2, 4)
##    plt.axis("off")
##    plt.imshow(R_image * R_filter)
    
    
def show():
    plt.show()


class DataFormatException(Exception):
    pass
