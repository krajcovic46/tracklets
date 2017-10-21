/*
 * Copyright (C) 2006 Matthew Funk
 * Licensed under the Academic Free License version 1.2
 * 
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
 * Academic Free License version 1.2 for more details.
 */
//NOTES: package org.jscience.astronomy.solarsystem.artificialsatellites;

package Compute.Silha;

/**
 *
 */
 class C1 {

    final static double E6A = 1.E-6;
    final static double TOTHRD = .66666667;
    final static double XJ2 = 1.082616E-3;
    final static double XJ4 = -1.65597E-6;
    final static double XJ3 = -.253881E-5;
    final static double XKE = .743669161E-1;
    final static double XKMPER = 6378.135;
    final static double XMNPDA = 1440.;
    final static double AE = 1.;
    final static double QO = 120.0;
    final static double SO = 78.0;
       
    final static double CK2 =.5*XJ2*Math.pow(AE,2);
    final static double CK4 = -.375*XJ4*Math.pow(AE,4);
    final static double QOMS2T = Math.pow(((QO-SO)*AE/XKMPER),4);
    final static double S = AE*(1.+SO/XKMPER);

    /**
     * Deny the ability to construct an instance of this class.
     */
    private C1() {
    }

}
