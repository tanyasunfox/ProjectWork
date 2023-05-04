package in.sunfox.healthcare.java.core;

import java.text.DecimalFormat;
import java.text.DecimalFormatSymbols;
import java.util.ArrayList;
import java.util.Locale;



    public class Filters {

        static double[] lowPassFilter(double[] signal, double cutOffFreq) {

            double dt = 0.002;
            double RC = 1 / (2 * Math.PI * cutOffFreq);
            double alpha = dt / (dt + RC);
            double[] filteredSignal = new double[signal.length];

            filteredSignal[0] = alpha * signal[0];
            for (int i = 1; i < signal.length; i++) {
                filteredSignal[i] = filteredSignal[i - 1] + alpha * (signal[i] - filteredSignal[i - 1]);
            }
            return filteredSignal;
        }


        static double[] highPassFilter(double[] signal, double cutOffFreq) {

            double dt = 0.002;
            double RC = 1 / (2 * Math.PI * cutOffFreq);
            double alpha = RC / (dt + RC);
            double[] filteredSignal = new double[signal.length];

            filteredSignal[0] = signal[0];
            for (int i = 1; i < signal.length; i++) {
                filteredSignal[i] = alpha * (filteredSignal[i - 1] + signal[i] - signal[i - 1]);
            }
            return filteredSignal;
        }

        public static double[] movingAverage(double[] signal) {

            double[] filtered = new double[signal.length - 8];

            for (int i = 4, j = 0; i < signal.length - 4; i++, j++) {
                filtered[j] = (signal[i - 4] + signal[i - 3] + signal[i - 2] + signal[i - 1] + signal[i] + signal[i + 1] + signal[i + 2] + signal[i + 3] + signal[i + 4]) / 9.0;
            }
            return filtered;
        }

        public static ArrayList<Double> movingAverage(ArrayList<Double> signal) {

            ArrayList<Double> filtered = new ArrayList<>();

            for (int i = 4, j = 0; i < signal.size() - 4; i++, j++) {
                filtered.add((signal.get(i + 4) + signal.get(i + 3) + signal.get(i + 2) + signal.get(i + 1) + signal.get(i) + signal.get(i - 1) + signal.get(i - 2) + signal.get(i - 3) + signal.get(i - 4)) / 9.0);
            }
            return filtered;
        }

        static double[] derivative(double[] signal) {
            double[] processed = new double[signal.length - 1];

            for (int i = 0; i < signal.length - 1; i++) {
                processed[i] = signal[i + 1] - signal[i];
            }
            return processed;
        }

        static double[] integral(double[] signal) {
            final double N = 16.0, T = 100.0;
            ArrayList<Double> processed = new ArrayList<>();
            for (int i = 9; i < signal.length / 2 + 1; i++) {
                double val = (signal[2 * i - 1 - ((int) N)] * T) / N;

                processed.add(val);
            }
            double[] arr = new double[processed.size()];
            for (int i = 0; i < processed.size(); i++)
                arr[i] = processed.get(i);
            return arr;
        }

        static double[] integral2(double[] signal) {
            final double N = 16.0, T = 100.0;
            ArrayList<Double> processed = new ArrayList<>();

            for (int i = 0; i < 15; i++) {
                processed.add(signal[i]);
            }

            for (int i = 15; i < signal.length; i++) {
                double sum = 0.0;
                for (int j = 0; j < 16; j++) {
                    sum += signal[i - j];
                }
                double val = sum / N;
                processed.add(val);
            }

            double[] arr = new double[processed.size()];
            for (int i = 0; i < processed.size(); i++)
                arr[i] = processed.get(i);
            return arr;
        }

        public static double mean(ArrayList<Double> list) {

            double total = 0.0;

            for (Double aDouble : list) {
                total += aDouble;
            }

            if (list.size() == 0)
                return 0.0;
            return total / (double) list.size();

        }

        public static double mean(double[] list) {

            double total = 0.0;

            for (double v : list) {
                total += v;
            }
            if (list.length == 0)
                return 0.0;

            return total / (double) list.length;

        }

        public static double[] arraySubtract(double[] from, double[] array, int length) {
            double[] result = new double[length];
            for (int i = 0; i < length; i++) {
                result[i] = from[i] - array[i];
            }
            return result;
        }

        public static double[] arraySubtract(ArrayList<Double> from, ArrayList<Double> array, int length) {
            double[] result = new double[length];
            for (int i = 0; i < length; i++) {
                result[i] = from.get(i) - array.get(i);
            }
            return result;
        }

//    public static ArrayList<Double> DFT( ArrayList<Double> array)
//    {
//        ArrayList<Double> X = new ArrayList<>();
//        int k, l;
//        double c_real = 0;
//        double c_img = 0;
//        for (k = 0; k < array.size(); k++)
//        {
//            for (l = 0; l < array.size(); l++)
//            {
//                double real = Math.cos(-2 * Math.PI * l * k / (array.size()));
//                double img = Math.sin(-2 * Math.PI * l * k / (array.size()));
//                c_real += array.get(l) * real;
//                c_img += array.get(l) * img;
//            }
//            X.add(Math.sqrt((c_real * c_real) + (c_img * c_img)));
//            c_real = 0;
//            c_img = 0;
//        }
//        return X;
//    }

//        public static ArrayList<Double> FFT(double[] signal) {
//            double valueOfPower = Math.ceil(Math.log(signal.length) / Math.log(2.0));
//            int MaxOrEqualToSignalLength = (int) Math.pow(2, valueOfPower);
//            ArrayList<Double> fftMagnitude = new ArrayList<>(signal.length);
//
//            if (MaxOrEqualToSignalLength == signal.length) {
////            System.out.println("radix2 Algorithm");
//                RadixTwoFFTUsingDIF p = new RadixTwoFFTUsingDIF(signal, false);
//                p.process();
//                double[] fftRealSection = p.getFftRealSection();
//                double[] fftImaginarySection = p.getFftImaginarySection();
//
//                for (int i = 0; i < signal.length; i++) {
//                    fftMagnitude.add(Math.sqrt(fftRealSection[i] * fftRealSection[i] + fftImaginarySection[i] * fftImaginarySection[i]));
//                }
//
//            } else {
////            System.out.println("Blustein Algorithm");
//                BlueSteinMethodForFFT rs = new BlueSteinMethodForFFT(signal);
//                rs.process();
//                double[] fftRS = rs.getFftOfSignalReal();
//                double[] fftIS = rs.getFftOfSignalImaginary();
//
//                for (int i = 0; i < signal.length; i++) {
//                    fftMagnitude.add(Math.sqrt(fftRS[i] * fftRS[i] + fftIS[i] * fftIS[i]));
//                }
//            }
//            return fftMagnitude;
//        }

        public static double max(ArrayList<Double> array) {
            double max = 0.0;
            for (Double aDouble : array) {
                if (aDouble > max) {
                    max = aDouble;
                }
            }

            return max;
        }

        public static double min(ArrayList<Double> array) {
            double min = max(array);
            for (Double aDouble : array) {
                if (aDouble < min) {
                    min = aDouble;
                }
            }

            return min;
        }

        public static String localDecimalFormat(double value, String pattern) {
            DecimalFormatSymbols symbols = new DecimalFormatSymbols(Locale.US);
            DecimalFormat decimalFormat = new DecimalFormat(pattern, symbols);
            decimalFormat.setMinimumFractionDigits(2);
            return decimalFormat.format(value);
        }
}

