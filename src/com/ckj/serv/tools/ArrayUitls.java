package com.ckj.serv.tools;

public class ArrayUitls {

    public static double[][] to2dArray(double[] array, int rows) {
        double[][] arrar2d = new double[rows][];
        for (int i = 0; i < rows; i++) {
            double[] temp = new double[array.length / rows];
            for (int j = 0; j < temp.length; j++) {
                temp[j] = array[i * temp.length + j];
            }
            arrar2d[i] = temp;
        }
        return arrar2d;
    }
}
