package `in`.sunfox.healthcare.java.core

import java.util.ArrayList

class EcgCharacteristicsCalculator (

    private val signal: ArrayList<Double>,
    private val augmentedLead: Int ,       //7
    private val applyFilter: Boolean,
    private val adjustRPeaks:Boolean
    ) {
        /*
        Responsible for calculating ECG in.sunfox.healthcare.java.commons.ecg_processor.conclusions.characteristics i.e. PR, QRS, QT, QTc, RR, ST-Elev, HeartRate
         */

        fun getCharacteristics(msDifference: Double): EcgCharacteristics? {
            val ecgProcessing = ECGProcessing(signal, augmentedLead, applyFilter,adjustRPeaks)
            ecgProcessing.msDifference = msDifference
            if (!ecgProcessing.isEcgSignalCompatibleForProcessing)
//            println("Not an ECG")
                return null
            ecgProcessing.process()
            return EcgCharacteristics(
                ecgProcessing.prInterval,
                ecgProcessing.qrsInterval,
                ecgProcessing.qtInterval,
                ecgProcessing.qTcInterval,
                ecgProcessing.rrInterval,
                ecgProcessing.heartRate,
                ecgProcessing.stSegmentInMv,
                ecgProcessing.qRSIntervals,
                ecgProcessing.rrIntervals,
                ecgProcessing.prStopIndices,
                ecgProcessing.prStartIndices,
                ecgProcessing.pPoints,
                ecgProcessing.qPoints,
                ecgProcessing.sPoints,
                ecgProcessing.tPoints,
                ecgProcessing.rPeakPoints,
                ecgProcessing.tWaveEndPoints,
                ecgProcessing.averagePAmplitudeInLead,// added For New Algorithm
                ecgProcessing.averageQAmplitudeInLead,// added For New Algorithm
                ecgProcessing.averageSAmplitudeInLead,// added For New Algorithm
                ecgProcessing.averageTAmplitudeInLead,// added For New Algorithm
                ecgProcessing.averageRAmplitudeInLead,// added For New Algorithm
                ecgProcessing.pWidth,// added For New Algorithm
                ecgProcessing.tWidth,// added For New Algorithm
                ecgProcessing.qrsDirectionUpward,// added For New Algorithm
                ecgProcessing.valueRS,// added For New Algorithm
                ecgProcessing.ventricularActivationLOR,// added For New Algorithm
                ecgProcessing.ventricularActivationROR,// added For New Algorithm
                ecgProcessing.concavity,// added For New Algorithm
                ecgProcessing.frequencyOfPatternInQRS,// added For New Algorithm
                ecgProcessing.frequencyOfPatternInRR,// added For New Algorithm
                ecgProcessing.pAmplitudeArrayInMv,
                ecgProcessing.TRRatioSatisfy,
                ecgProcessing.TSRatioSatisfy,
                ecgProcessing.multiPWave,
                ecgProcessing.pPeaksIndex,
                ecgProcessing.prIntervalInMilliSec,
                ecgProcessing.arrayRateOfMultiPWave,
                ecgProcessing.pWavePresent,
                ecgProcessing.tentedTWave
            )
        }

}