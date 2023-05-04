package `in`.sunfox.healthcare.java.core

import com.google.gson.annotations.SerializedName
import java.io.Serializable


data class EcgCharacteristics(
    @SerializedName("pr") var pr: Int,
    @SerializedName("qrs") var qrs: Int,
    @SerializedName("qt") var qt: Int,
    @SerializedName("qtc") var qtc: Double, // added For New Algorithm: change is in return type Integer to Double
    @SerializedName("rr") var rr: Int,
    @SerializedName("bpm") var heartRate: Int,
    @SerializedName("st_elevation") var stElevation: Double,
    @SerializedName("qrs_intervals") var qrsIntervals: ArrayList<Double>,
    @SerializedName("rr_intervals") var rrIntervals: ArrayList<Double>,
    @SerializedName("pr_stop_indices") var prStopIndices: ArrayList<Double>,
    @SerializedName("pr_start_indices") var prStartIndices: ArrayList<Double>,
    @SerializedName("p_wave_points") var pWavePoints: ArrayList<Double>,
    @SerializedName("q_wave_points") var qWavePoints: ArrayList<Double>,
    @SerializedName("s_wave_points") var sWavePoints: ArrayList<Double>,
    @SerializedName("t_wave_points") var tWavePoints: ArrayList<Double>,
    @SerializedName("r_peak_points") var rPeakPoints: ArrayList<Double>,
    @SerializedName("t_wave_end_points") var tWaveEndPoints: ArrayList<Double>,

    //==================================== added For New Algorithm ====================================//

    @SerializedName("averagePAmplitude") var averagePAmplitudeInLead: Double,
    @SerializedName("averageQAmplitude") var averageQAmplitudeInLead: Double,
    @SerializedName("averageSAmplitude") var averageSAmplitudeInLead: Double,
    @SerializedName("averageTAmplitude") var averageTAmplitudeInLead: Double,
    @SerializedName("averageRAmplitude") var averageRAmplitudeInLead: Double,
    @SerializedName("pWidth") var pWidth: Double,
    @SerializedName("tWidth") var tWidth: Double,
    @SerializedName("qrsDirectionUpward") var qrsDirectionUpward:Boolean,
    @SerializedName("ratioRS") var ratioRS: Double,
    @SerializedName("ventricularActivationLOR") var ventricularActivationLOR: Double,
    @SerializedName("ventricularActivationROR") var ventricularActivationROR: Double,
    @SerializedName("indication") var concavity: Boolean,
    @SerializedName("frequencyOfPatterninQrsArray") var frequencyOfPatternInQRS: Int,
    @SerializedName("frequencyOfPatterninRRArray") var frequencyOfPatternInRR: Int,
    @SerializedName("pAmplitudeArrayInMv") var pAmplitudeArrayInMv: ArrayList<Double>,
    @SerializedName("TRRatioSatisfy") var TRRatioSatisfy:Boolean,
    @SerializedName("TSRatioSatisfy") var TSRatioSatisfy:Boolean,

    ///////// For Hyperklemia and Arrhythmia //////////////////////////////////////////////
    @SerializedName("isMultiPWave") var multiPWave:Boolean,
    @SerializedName("pPeaksIndex") var pPeaksIndex: ArrayList<Double>,
    @SerializedName("prIntervalInMilliSec") var prIntervalInMilliSec: ArrayList<Double>,
    @SerializedName("arrayRateOfMultiPWave") var arrayRateOfMultiPWave: ArrayList<Double>,
    @SerializedName("isPwavePresent") var  pWavePresent:Boolean,
    @SerializedName("tentedTWave") var tentedTWave:Boolean,
    ///////// For Hyperklemia and Arrhythmia //////////////////////////////////////////////


    //==================================== added For New Algorithm ====================================//

) : Serializable {
    override fun toString(): String
    {
        return "EcgCharacteristics(pr=$pr, qrs=$qrs, qt=$qt, qtc=$qtc, rr=$rr, heartRate=$heartRate, stElevation=$stElevation, qrsIntervals=$qrsIntervals, rrIntervals=$rrIntervals, prStopIndices=$prStopIndices, pWavePoints=$pWavePoints, qWavePoints=$qWavePoints, sWavePoints=$sWavePoints, tWavePoints=$tWavePoints, rPeakPoints=$rPeakPoints, tWaveEndPoints=$tWaveEndPoints, averagePAmplitudeInLead=$averagePAmplitudeInLead, averageQAmplitudeInLead=$averageQAmplitudeInLead, averageSAmplitudeInLead=$averageSAmplitudeInLead, averageTAmplitudeInLead=$averageTAmplitudeInLead, averageRAmplitudeInLead=$averageRAmplitudeInLead,pWdth=$pWidth,tWidth=$tWidth,qrsDirectionUpward=$qrsDirectionUpward,ratioRS=$ratioRS,ventricularActivationLOR=$ventricularActivationLOR,ventricularActivationROR=$ventricularActivationROR,concavity=$concavity,frequencyOfPatternInQRS=$frequencyOfPatternInQRS,frequencyOfPatternInRR=$frequencyOfPatternInRR,pAmplitudeArrayInMv=$pAmplitudeArrayInMv,TRRatioSatisfy=$TRRatioSatisfy,TSRatioSatisfy=$TSRatioSatisfy,multiPWave,=$multiPWave,, pPeaksIndex=$pPeaksIndex,prIntervalInMilliSec=$prIntervalInMilliSec,arrayRateOfMultiPWave=$arrayRateOfMultiPWave,pWavePresent=$pWavePresent,tentedTWave=$tentedTWave)"
    }

}
