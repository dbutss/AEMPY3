import time
import numpy as np

def mt1dfwd(frequency, resistivity, thickness, rtype="rhophas", out=True):
    
    """
    1D magnetotellutic forward modelling'
    based on A. Pethik's script at www.digitalearthlab.com
    Last chanre vr Aug 2020
    """
     
    mu0 = 4*np.pi*1E-7   		# Magnetic Permeability (H/m)
    
    n = np.size(resistivity);
    
    print('freq\tares\t\t\tphase');
    for freq in frequency:   
        w =  2*np.pi*freq       
        impedances = np.array(range(n));
        #compute basement impedance
        impedances[n-1] = np.sqrt(w*mu0*resistivity[n-1]*1j);
       
        for layer in range(n-2,-1,-1):
            res     = resistivity[layer];
            thick   = thickness[layer];
      
            # 3. Compute apparent resistivity from top layer impedance
            #Step 2. Iterate from bottom layer to top(not the basement) 
            # Step 2.1 Calculate the intrinsic impedance of current layer
            dj = np.sqrt(w * mu0 * 1.0/res*1j);
            wj = dj * res;
            # Step 2.2 Calculate Exponential factor from intrinsic impedance
            ej = np.exp(-2*thick*dj);                     
        
            # Step 2.3 Calculate reflection coeficient using current layer
            #          intrinsic impedance and the below layer impedance 
            belowImpedance = impedances[layer + 1];
            rj = (wj - belowImpedance)/(wj + belowImpedance);
            re = rj*ej; 
            Zj = wj * ((1 - re)/(1 + re));
            impedances[layer] = Zj;    
    
        # Step 3. Compute apparent resistivity from top layer impedance
        Z = impedances[0];
    
    
        if 	 rtype.lower() == "imped":
            if out: print(frequency, '\t', Z.real, '\t', Z.imag)
            return Z
    
        elif rtype.lower() == "rhophas":
            rhoa 	= np.abs(Z)**2 / (mu0*w)
            phase   = np.rad2deg(np.arctan(Z.imag / Z.real))
            if out: print(frequency, '\t', rhoa, '\t', phase)
            return rhoa, phase
    
        else:
            raise Exception("rtype must be 'impedance' or 'rhophas', not {}".format(rtype.lower()))
  